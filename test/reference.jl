# Independent reference implementation of the Posit system, used as an "oracle"
# by the tests. It does NOT use any function from PositPR: it works directly on
# the bit pattern (as an unsigned integer) and uses BigFloat to obtain the exact
# value of each operation before rounding.
#
# Conventions (Posit Standard 2022, Gustafson & Yonemoto 2017):
#   * 0b000…0          -> zero
#   * 0b100…0          -> NaR (Not a Real)
#   * negatives         -> two's complement of the positive pattern
#   * rounding          -> to nearest, ties to the even pattern
#                          (done on the bit string, as in the standard)
#   * no underflow/overflow: nonzero values saturate at minpos / maxpos

using PositPR
using PositPR: Posit, PositEncoder, PositDecoder, x2Posit, Posit2x

setprecision(BigFloat, 1024)

include("watchdog.jl")

"Thrown by `gcall` when the function does not finish within the time limit."
struct HangError <: Exception
    call::String
end

"""
    gcall(f, args...)

Call `f(args...)` under the watchdog (see watchdog.jl). Use it for library
functions that contain `while` loops (eps, nextPosit, exp, log, sin, ...).
"""
function gcall(f, args...; timeout = WATCHDOG_TIMEOUT)
    r = guarded(f, args...; timeout = timeout)
    r isa TimedOut && throw(HangError("$(f)($(join(map(_short, args), ", ")))"))
    return r
end

_short(p::Posit) = string(Posit2x(p))
_short(x) = string(x)

const NAR_VALUE = NaN

mask(ps) = (UInt64(1) << ps) - 1
nar_pattern(ps) = UInt64(1) << (ps - 1)

"Convert the pattern `u` (ps bits) to a signed integer (two's complement)."
function signed_pattern(u::Integer, ps)
    u = UInt64(u)
    u &= mask(ps)
    return (u >> (ps - 1)) == 1 ? Int64(u) - (Int64(1) << ps) : Int64(u)
end

"Distance in ULPs (positions on the posit line) between two patterns."
ulp_distance(u1::Integer, u2::Integer, ps) = abs(signed_pattern(u1, ps) - signed_pattern(u2, ps))

bitstr(u::Integer, ps) = string(UInt64(u) & mask(ps); base = 2, pad = ps)

"Bit pattern of a library Posit, via PositEncoder."
pattern(p::Posit{ps, es}) where {ps, es} = parse(UInt64, PositEncoder(p); base = 2)

"All patterns of ps bits."
all_patterns(ps) = UInt64(0):mask(ps)

"All patterns that represent real numbers (NaR excluded)."
real_patterns(ps) = (u for u in all_patterns(ps) if u != nar_pattern(ps))

"""
    ref_value(u, ps, es) -> Float64

Exact value of the pattern `u` in Posit{ps,es} (NaN for NaR). For ps ≤ 32 and
es ≤ 3 the value is always exactly representable as a Float64.
"""
function ref_value(u::Integer, ps, es)
    u = UInt64(u)
    u &= mask(ps)
    u == 0 && return 0.0
    u == nar_pattern(ps) && return NAR_VALUE
    neg = (u >> (ps - 1)) == 1
    neg && (u = (~u + 1) & mask(ps))

    nb = ps - 1                                   # bits after the sign
    bit(i) = i <= nb ? Int((u >> (nb - i)) & 1) : 0
    r0 = bit(1)
    m = 1
    while m < nb && bit(m + 1) == r0
        m += 1
    end
    k = r0 == 1 ? m - 1 : -m
    pos = m + 2                                   # skip the regime terminating bit
    e = 0
    for _ in 1:es                                 # missing exponent bits are 0
        e = 2e + bit(pos)
        pos += 1
    end
    f = 0.0
    w = 0.5
    for i in pos:nb
        f += bit(i) * w
        w /= 2
    end
    v = ldexp(1.0 + f, k * 2^es + e)
    return neg ? -v : v
end

ref_value(p::Posit{ps, es}) where {ps, es} = ref_value(pattern(p), ps, es)

"""
    ref_encode(x, ps, es) -> UInt64

Bit pattern of the Posit{ps,es} nearest to `x`, with RNE rounding applied to the
infinite bit string (regime, exponent and fraction), as in the standard.
"""
function ref_encode(x::Real, ps, es)
    (isnan(x) || isinf(x)) && return nar_pattern(ps)
    iszero(x) && return UInt64(0)
    y = abs(BigFloat(x))
    E = exponent(y)                               # y = m·2^E, 1 ≤ m < 2
    useed_pow = 2^es
    k = fld(E, useed_pow)
    e = E - k * useed_pow
    body = Int[]
    if k >= 0
        append!(body, fill(1, k + 1)); push!(body, 0)
    else
        append!(body, fill(0, -k)); push!(body, 1)
    end
    for j in es-1:-1:0
        push!(body, (e >> j) & 1)
    end
    frac = y / BigFloat(2)^E - 1
    while length(body) < ps                       # ps-1 kept bits + 1 guard bit
        frac *= 2
        b = frac >= 1 ? 1 : 0
        b == 1 && (frac -= 1)
        push!(body, b)
    end
    kept = body[1:ps-1]
    guard = body[ps]
    sticky = any(==(1), body[ps+1:end]) || !iszero(frac)
    v = UInt64(0)
    for b in kept
        v = (v << 1) | UInt64(b)
    end
    if guard == 1 && (sticky || isodd(v))
        v += 1
    end
    v == 0 && (v = UInt64(1))                     # never rounds to zero
    v >= nar_pattern(ps) && (v = nar_pattern(ps) - 1)  # nor to NaR
    x < 0 && (v = (~v + 1) & mask(ps))
    return v
end

ref_encode(x::Real, ::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = ref_encode(x, ps, es)

"Build a library Posit from a bit pattern."
from_pattern(u::Integer, ps, es) = PositDecoder(bitstr(u, ps), ps, es)
from_pattern(u::Integer, ::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = from_pattern(u, ps, es)

"Exact value (BigFloat) used as the reference for operations."
big_value(u::Integer, ps, es) = BigFloat(ref_value(u, ps, es))

"""
    check_all(f, itr; maxshow=5)

Run `f(x)` for every `x` in `itr`. `f` must return `nothing` when the case
passes or a `String` describing the failure. Exceptions count as failures.
Returns `(nfail, ntotal, examples)` so that a single `@test` can summarize
thousands of cases without flooding the output.
"""
function check_all(f, itr; maxshow = 5, maxhangs = 2)
    nfail = 0
    ntotal = 0
    nhangs = 0
    examples = String[]
    items = collect(itr)
    for (i, x) in enumerate(items)
        ntotal += 1
        msg = try
            f(x)
        catch err
            if err isa HangError
                nhangs += 1
                "infinite loop (> $(WATCHDOG_TIMEOUT) s) in $(err.call)"
            else
                "exception $(typeof(err)) in $(x): " * first(sprint(showerror, err), 120)
            end
        end
        if msg !== nothing
            nfail += 1
            length(examples) < maxshow && push!(examples, msg)
        end
        if nhangs >= maxhangs && i < length(items)
            skipped = length(items) - i
            nfail += skipped
            ntotal += skipped
            push!(examples, "$skipped remaining cases skipped after $nhangs infinite loops")
            break
        end
    end
    return nfail, ntotal, examples
end

"Print a readable summary when there are failures; return `true` if everything passed."
function report(name, (nfail, ntotal, examples))
    if nfail > 0
        println("  ✗ $name: $nfail of $ntotal cases failed. Examples:")
        foreach(ex -> println("      • ", ex), examples)
    end
    return nfail == 0
end

"Patterns used in the 16/32-bit tests: extremes + an evenly spaced sample."
function sample_patterns(ps; step = 7)
    n = mask(ps)
    edge = UInt64[0, 1, 2, 3, n >> 1, (n >> 1) - 1, (n >> 1) - 2, nar_pattern(ps) + 1,
                  nar_pattern(ps) + 2, n, n - 1, UInt64(1) << (ps - 2)]
    return unique(vcat(edge, collect(UInt64(0):UInt64(step):n)))
end

"{ps, es} configurations used in the exhaustive 8-bit tests."
const CONFIGS_8 = [(8, 0), (8, 1), (8, 2), (8, 3)]
