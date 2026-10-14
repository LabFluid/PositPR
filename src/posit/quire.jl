# The quire: an exact accumulator for sums of posits and of products of posits,
# as defined by the Posit Standard (2022). Sums and dot products are accumulated
# with no rounding at all and rounded only once, when converted back to a posit.
#
# Every product of two Posit{ps,es} is an integer multiple of minpos² =
# 2^(-2^(es+1)·(ps-2)), so the quire is a fixed-point number whose least
# significant bit is worth minpos². It is stored as an integer `acc` with
#
#     value = acc · 2^(-quire_scale(ps, es)),   quire_scale = 2^(es+1)·(ps-2).
#
# The standard fixes the quire at 16·ps bits (and turns overflow into NaR); here
# `acc` is a BigInt, so the quire never overflows. Whenever the standard quire
# does not overflow, both give exactly the same results.

"""
    Quire{ps, es}

Exact accumulator for posits of type `Posit{ps, es}`. Create it with
`Quire(Posit{ps,es})` (zero) or `Quire(p)` (value of a posit), accumulate with
`qadd!`, `qsub!`, `qmuladd!`, `qmulsub!` and round it back with `Posit(q)`.

```julia
q = Quire(Posit{16,2})
for (x, y) in zip(xs, ys)
    qmuladd!(q, x, y)        # q += x·y, exactly
end
Posit(q)                     # Σ xᵢ·yᵢ, rounded once
```
"""
mutable struct Quire{ps, es}
    acc::BigInt     # value = acc · 2^(-quire_scale(ps, es))
    nar::Bool       # NaR is contagious: once true, the result is NaR
end

"Exponent of the least significant quire bit: its weight is 2^(-quire_scale) = minpos²."
quire_scale(ps::Integer, es::Integer) = (1 << (es + 1)) * (ps - 2)

"Number of bits the Posit Standard reserves for the quire of a ps-bit posit."
standard_quire_bits(ps::Integer) = 16 * ps

Quire(::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = Quire{ps, es}(big(0), false)
Quire{ps, es}() where {ps, es} = Quire{ps, es}(big(0), false)
Quire(p::Posit{ps, es}) where {ps, es} = qadd!(Quire(typeof(p)), p)

Base.copy(q::Quire{ps, es}) where {ps, es} = Quire{ps, es}(copy(q.acc), q.nar)
Base.zero(::Type{Quire{ps, es}}) where {ps, es} = Quire{ps, es}()

isNaR(q::Quire) = q.nar
Base.iszero(q::Quire) = !q.nar && iszero(q.acc)

"Reset the quire to zero."
function qclear!(q::Quire)
    q.acc = big(0)
    q.nar = false
    return q
end

"Add the exact value (-1)^s · M · 2^P to the quire."
function quire_accumulate!(q::Quire{ps, es}, s::Integer, M::Integer, P::Integer) where {ps, es}
    shift = P + quire_scale(ps, es)      # ≥ 0 for any posit or product of two posits
    term = big(M) << shift
    q.acc = s == 1 ? q.acc - term : q.acc + term
    return q
end

"`q += a`, exactly."
function qadd!(q::Quire{ps, es}, a::Posit{ps, es}) where {ps, es}
    if isNaR(a)
        q.nar = true
    elseif !iszero(a)
        quire_accumulate!(q, exact_parts(a)...)
    end
    return q
end

"`q -= a`, exactly."
qsub!(q::Quire{ps, es}, a::Posit{ps, es}) where {ps, es} = qadd!(q, -a)

"`q += a·b`, exactly (the product is not rounded)."
function qmuladd!(q::Quire{ps, es}, a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    if isNaR(a) || isNaR(b)
        q.nar = true
    elseif !(iszero(a) || iszero(b))
        sa, Ma, Pa = exact_parts(a)
        sb, Mb, Pb = exact_parts(b)
        quire_accumulate!(q, sa ⊻ sb, big(Ma) * Mb, Pa + Pb)
    end
    return q
end

"`q -= a·b`, exactly."
qmulsub!(q::Quire{ps, es}, a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es} = qmuladd!(q, -a, b)

"`q += r` for two quires of the same type."
function qadd!(q::Quire{ps, es}, r::Quire{ps, es}) where {ps, es}
    q.nar |= r.nar
    q.acc += r.acc
    return q
end

"`q -= r` for two quires of the same type."
function qsub!(q::Quire{ps, es}, r::Quire{ps, es}) where {ps, es}
    q.nar |= r.nar
    q.acc -= r.acc
    return q
end

"`q = -q`."
function qneg!(q::Quire)
    q.acc = -q.acc
    return q
end

"`q = |q|`."
function qabs!(q::Quire)
    q.acc = abs(q.acc)
    return q
end

Base.:-(q::Quire) = qneg!(copy(q))
Base.:+(q::Quire{ps, es}, r::Quire{ps, es}) where {ps, es} = qadd!(copy(q), r)
Base.:-(q::Quire{ps, es}, r::Quire{ps, es}) where {ps, es} = qsub!(copy(q), r)

"""
    Posit(q::Quire)

Round the quire to the nearest posit (RNE, a single rounding).
"""
function Posit(q::Quire{ps, es}) where {ps, es}
    T = PositType(ps, es)
    q.nar && return NaR(T)
    iszero(q.acc) && return zero(T)
    return posit_round(T, q.acc < 0 ? 1 : 0, abs(q.acc), -quire_scale(ps, es))
end
Posit{ps, es}(q::Quire{ps, es}) where {ps, es} = Posit(q)
Posit{ps, es, U}(q::Quire{ps, es}) where {ps, es, U} = Posit(q)

"Exact value of the quire as a BigFloat."
function Base.BigFloat(q::Quire{ps, es}) where {ps, es}
    q.nar && return BigFloat(NaN)
    prec = max(256, iszero(q.acc) ? 1 : ndigits(q.acc; base = 2))
    return ldexp(setprecision(() -> BigFloat(q.acc), BigFloat, prec), -quire_scale(ps, es))
end

function Base.show(io::IO, q::Quire{ps, es}) where {ps, es}
    print(io, "Quire{", ps, ",", es, "}(")
    q.nar ? print(io, "NaR") : print(io, Float64(BigFloat(q)))
    print(io, ")")
end

# ------------------------------------------------------------------ fused operations

"""
    fdp(xs, ys)

Fused dot product: Σ xᵢ·yᵢ accumulated exactly in a quire and rounded once.
"""
function fdp(xs::AbstractVector{<:Posit{ps, es}}, ys::AbstractVector{<:Posit{ps, es}}) where {ps, es}
    length(xs) == length(ys) || throw(DimensionMismatch("vectors must have the same length"))
    q = Quire{ps, es}()
    for (x, y) in zip(xs, ys)
        qmuladd!(q, x, y)
    end
    return Posit(q)
end

"""
    fsum(xs)

Fused sum: Σ xᵢ accumulated exactly in a quire and rounded once.
"""
function fsum(xs::AbstractArray{<:Posit{ps, es}}) where {ps, es}
    q = Quire{ps, es}()
    for x in xs
        qadd!(q, x)
    end
    return Posit(q)
end

"""
    fma(a, b, c)

Fused multiply-add: a·b + c with a single rounding (computed in a quire).
"""
function Base.fma(a::Posit{ps, es}, b::Posit{ps, es}, c::Posit{ps, es}) where {ps, es}
    q = Quire(c)
    qmuladd!(q, a, b)
    return Posit(q)
end
