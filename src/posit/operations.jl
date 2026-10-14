# Posit arithmetic.
#
# Every operation computes the EXACT result from the integer mantissas (with the
# hidden bit) and the expanded exponents E = 2^es·k + e of the operands, and only
# then rounds once with `posit_round` (rounding.jl). Hence all basic operations
# are correctly rounded (RNE), as required by the standard. Intermediate integers
# are UInt128 when they fit in 127 bits (always for ps ≤ 32), BigInt otherwise.

import Base: >>, <<
# Multiplication/division by 2^n: only the exponent changes, but the regime may
# change size, so the result also goes through posit_round.
function >>(a::Posit{ps, es}, n::Int64) where {ps, es}
    is_special(a) && return a
    s, M, P = exact_parts(a)
    return posit_round(Posit{ps, es}, s, M, P - n)
end

function <<(a::Posit{ps, es}, n::Int64) where {ps, es}
    is_special(a) && return a
    s, M, P = exact_parts(a)
    return posit_round(Posit{ps, es}, s, M, P + n)
end

"Exact sum of two values `(s, M, P)`, rounded to `T`."
function add_exact(::Type{T}, (sa, Ma, Pa), (sb, Mb, Pb)) where {ps, es, T <: Posit{ps, es}}
    iszero(Ma) && return posit_round(T, sb, Mb, Pb)
    iszero(Mb) && return posit_round(T, sa, Ma, Pa)

    # If one operand is far below the rounding position of the other, it only
    # matters as a sticky bit: with g extra bits, x ± tiny lies strictly between
    # (x·2^g - 1) and x·2^g in units of 2^(P-g). g = ps + 4 makes sure that, for
    # x - tiny, the run of ones of (x·2^g - 1) reaches past the rounding position.
    Ta, Tb = Pa + bitlength(Ma), Pb + bitlength(Mb)
    g = ps + 4
    if Ta - Tb > ps + 8
        A = fit_int(Ma, g) << g
        return posit_round(T, sa, sa == sb ? A : A - one(A), Pa - g, true)
    elseif Tb - Ta > ps + 8
        B = fit_int(Mb, g) << g
        return posit_round(T, sb, sa == sb ? B : B - one(B), Pb - g, true)
    end

    # Alignment: bring both mantissas to the smaller exponent (no bits are lost)
    P = min(Pa, Pb)
    need = max(bitlength(Ma) + Pa - P, bitlength(Mb) + Pb - P) + 1
    A = need <= 127 ? UInt128(Ma) << (Pa - P) : big(Ma) << (Pa - P)
    B = need <= 127 ? UInt128(Mb) << (Pb - P) : big(Mb) << (Pb - P)
    if sa == sb
        S, s = A + B, sa
    elseif A >= B
        S, s = A - B, sa
    else
        S, s = B - A, sb
    end
    iszero(S) && return zero(T)
    return posit_round(T, s, S, P)
end

import Base: +
function +(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    (isNaR(a) || isNaR(b)) && return NaR(Posit{ps, es})
    iszero(a) && return b
    iszero(b) && return a
    return add_exact(Posit{ps, es}, exact_parts(a), exact_parts(b))
end

import Base: -
-(a::Posit{ps, es}) where {ps, es} = from_bits(Posit{ps, es}, negate_bits(bits(a), ps))

function -(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    (isNaR(a) || isNaR(b)) && return NaR(Posit{ps, es})
    iszero(b) && return a
    iszero(a) && return -b
    sb, Mb, Pb = exact_parts(b)
    return add_exact(Posit{ps, es}, exact_parts(a), (1 - sb, Mb, Pb))
end

"Exact product of two integers, as UInt128 when it fits."
mul_int(A::Integer, B::Integer) =
    bitlength(A) + bitlength(B) <= 127 ? UInt128(A) * UInt128(B) : big(A) * big(B)

import Base: *
function *(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    (isNaR(a) || isNaR(b)) && return NaR(Posit{ps, es})
    (iszero(a) || iszero(b)) && return zero(Posit{ps, es})
    sa, Ma, Pa = exact_parts(a)
    sb, Mb, Pb = exact_parts(b)
    # (Ma·2^Pa)(Mb·2^Pb) = (Ma·Mb)·2^(Pa+Pb), exact
    return posit_round(Posit{ps, es}, sa ⊻ sb, mul_int(Ma, Mb), Pa + Pb)
end

import Base: /
function /(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    (isNaR(a) || isNaR(b)) && return NaR(Posit{ps, es})
    iszero(b) && throw(DivideError())
    iszero(a) && return a
    sa, Ma, Pa = exact_parts(a)
    sb, Mb, Pb = exact_parts(b)
    return posit_round_div(Posit{ps, es}, sa ⊻ sb, Ma, Pa, Mb, Pb)
end

# ------------------------------------------------------------- operations with integers
# The integer enters exactly (it is not rounded to a Posit first).

function +(a::Posit{ps, es}, n::Integer) where {ps, es}
    isNaR(a) && return a
    return add_exact(Posit{ps, es}, exact_parts(a), exact_parts(n))
end
+(n::Integer, a::Posit{ps, es}) where {ps, es} = a + n

function -(a::Posit{ps, es}, n::Integer) where {ps, es}
    isNaR(a) && return a
    sn, Mn, Pn = exact_parts(n)
    return add_exact(Posit{ps, es}, exact_parts(a), (1 - sn, Mn, Pn))
end
-(n::Integer, a::Posit{ps, es}) where {ps, es} = (-a) + n

function *(a::Posit{ps, es}, n::Integer) where {ps, es}
    isNaR(a) && return a
    (iszero(a) || iszero(n)) && return zero(Posit{ps, es})
    sa, Ma, Pa = exact_parts(a)
    sn, Mn, _ = exact_parts(n)
    return posit_round(Posit{ps, es}, sa ⊻ sn, mul_int(Ma, Mn), Pa)
end
*(n::Integer, a::Posit{ps, es}) where {ps, es} = a * n

function /(a::Posit{ps, es}, n::Integer) where {ps, es}
    iszero(n) && throw(DivideError())
    isNaR(a) && return a
    iszero(a) && return a
    sa, Ma, Pa = exact_parts(a)
    sn, Mn, _ = exact_parts(n)
    return posit_round_div(Posit{ps, es}, sa ⊻ sn, Ma, Pa, Mn, 0)
end

function /(n::Integer, a::Posit{ps, es}) where {ps, es}
    isNaR(a) && return a
    iszero(a) && throw(DivideError())
    iszero(n) && return zero(Posit{ps, es})
    sa, Ma, Pa = exact_parts(a)
    sn, Mn, _ = exact_parts(n)
    return posit_round_div(Posit{ps, es}, sa ⊻ sn, Mn, 0, Ma, Pa)
end

# ------------------------------------------------------------------- square root
import Base: sqrt
function sqrt(a::Posit{ps, es}) where {ps, es}
    is_special(a) && return a
    sign_bit(a) == 1 && throw(DomainError(a, "square root not defined for negative numbers."))
    _, M, P = exact_parts(a)
    # sqrt(M·2^P) with even P: sqrt(M)·2^(P/2)
    shift = isodd(P) ? 1 : 0
    P -= shift
    K = ps + 8                             # extra result bits
    X = fit_int(M, shift + 2K) << (shift + 2K)
    r = isqrt(X)                           # integer square root (digit by digit)
    return posit_round(Posit{ps, es}, 0, r, (P - 2K) ÷ 2, r * r != X)
end

# ---------------------------------------------------------------------- powers
import Base: ^
const EXACT_POW_LIMIT = 4096

function ^(a::Posit{ps, es}, n::Integer) where {ps, es}
    T = Posit{ps, es}
    isNaR(a) && return a
    n == 0 && return one(T)
    if iszero(a)
        n > 0 && return zero(T)
        throw(DivideError())
    end
    s, M, P = exact_parts(a)
    sr = isodd(n) ? s : 0
    m = abs(n)
    if m <= EXACT_POW_LIMIT
        # a^m = M^m · 2^(P·m), exact; rounded once
        Mm = bitlength(M) * m <= 127 ? UInt128(M)^m : big(M)^m
        n > 0 && return posit_round(T, sr, Mm, P * m)
        return posit_round_div(T, sr, UInt128(1), 0, Mm, P * m)
    end
    # very large exponents: exponentiation by squaring in extended precision
    result = one(wide_type(T))
    base = widen_posit(abs(a))
    while m > 0
        isodd(m) && (result *= base)
        base *= base
        m ÷= 2
    end
    n < 0 && (result = one(result) / result)
    r = narrow_posit(result, T)
    return sr == 1 ? -r : r
end

function ^(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    T = Posit{ps, es}
    (isNaR(a) || isNaR(b)) && return NaR(T)
    iszero(b) && return one(T)
    iszero(a) && return zero(T)

    # integer exponent: exact path
    b_int = round(b, RoundToZero)
    if b_int == b && abs(Posit2x(b)) <= EXACT_POW_LIMIT
        return a ^ Int(b)
    end

    # a^(±1/2) = sqrt
    if abs(b) == x2Posit(0.5, T)
        return sign_bit(b) == 0 ? sqrt(a) : one(T) / sqrt(a)
    end

    if sign_bit(a) == 1
        throw(DomainError((a, b), "Exponentiation a^b is undefined for a < 0 when b is not an integer."))
    end
    w = exp_series(widen_posit(b) * log_series(widen_posit(a)))
    return narrow_posit(w, T)
end

# ------------------------------------------------------------ exponential and logarithm
import Base: exp
function exp(a::Posit{ps, es}) where {ps, es}
    isNaR(a) && return a
    return narrow_posit(exp_series(widen_posit(a)), Posit{ps, es})
end

"Taylor series of e^a, evaluated in the type of `a` itself."
function exp_series(a::Posit{ps, es}) where {ps, es}
    T = Posit{ps, es}
    iszero(a) && return one(T)

    # Negative exponent: use exp(-x) = 1 / exp(x)
    sign_bit(a) == 1 && return one(T) / exp_series(-a)

    # Above ln(maxpos) the result saturates
    Posit2x(a) > (1 << es) * (ps - 2) * log(2.0) && return typemax(T)

    # Range reduction: e^a = (e^(a/2^j))^(2^j), with a/2^j < 1
    j = max(0, expand_exponent(a) + 1)
    r = a >> j

    exp_r = one(T)
    term = one(T)
    n = 1
    while true
        term = term * r / n                # term = r^n / n!
        new = exp_r + term
        new == exp_r && break
        exp_r = new
        n += 1
    end
    for _ in 1:j
        exp_r = exp_r * exp_r
    end
    return exp_r
end

"Series for e^a - 1 (no cancellation for small |a|)."
function expm1_series(a::Posit{ps, es}) where {ps, es}
    T = Posit{ps, es}
    iszero(a) && return a
    sum = zero(T)
    term = one(T)
    n = 1
    while true
        term = term * a / n
        new = sum + term
        new == sum && break
        sum = new
        n += 1
    end
    return sum
end

# Use log(2) in high precision
const LN2 = big"0.693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687542001481020570"
const SQRT2 = big"1.414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534327641572735013846230912"

import Base: log
function log(a::Posit{ps, es}) where {ps, es}
    isNaR(a) && return a
    if sign_bit(a) == 1 || iszero(a)
        throw(DomainError(a, "log is only defined for positive real inputs."))
    end
    return narrow_posit(log_series(widen_posit(a)), Posit{ps, es})
end

"ln(a) = E·ln 2 + 2·atanh((y-1)/(y+1)), with a = 2^E·y and y ∈ [√2/2, √2)."
function log_series(a::Posit{ps, es}) where {ps, es}
    T = Posit{ps, es}
    _, M, P = exact_parts(a)
    E = P + bitlength(M) - 1
    y = posit_round(T, 0, M, P - E)        # mantissa in [1, 2)
    if y > x2Posit(SQRT2, T)               # center y around 1 to avoid cancellation
        y = y >> 1
        E += 1
    end
    z = (y - one(T)) / (y + one(T))
    z2 = z * z
    sum = zero(T)
    term = z
    k = 0
    while k < 500
        new = sum + term / (2k + 1)
        new == sum && break
        sum = new
        term *= z2
        k += 1
    end
    return x2Posit(E * LN2, T) + (sum << 1)
end

# --------------------------------------------------------------- rounding to integers
import Base: round, floor
"""
    round(p, r = RoundNearest)

Round the posit to an integer (returned as a Posit), like `Base.round` does
for Float64. Also provides `floor`, `ceil` and `trunc` (via RoundDown, RoundUp,
RoundToZero).
"""
function round(p::Posit{ps, es}, r::RoundingMode = RoundNearest) where {ps, es}
    is_special(p) && return p
    s, M, P = exact_parts(p)
    P >= 0 && return p                     # already an integer
    v = s == 1 ? -big(M) : big(M)
    q = div(v, big(1) << -P, r)
    iszero(q) && return zero(Posit{ps, es})
    return posit_round(Posit{ps, es}, q < 0 ? 1 : 0, abs(q), 0)
end

floor(p::Posit{ps, es}) where {ps, es} = round(p, RoundDown)

import Base: mod
function mod(a::Posit{ps, es}, b::Posit{ps, es})::Posit{ps, es} where {ps, es}
    if iszero(b)
        error("Division by zero in mod function")
    end
    # Compute a % b as: a - b * floor(a / b)
    return a - (b * floor(a / b))
end

import Base: max
function max(a::Posit{ps, es}, b::Posit{ps2, es2}) where {ps, es, ps2, es2}
    return a >= b ? a : b
end

function max(a::Posit{ps, es}, b::T) where {ps, es, T <: Real}
    b_posit = x2Posit(b, Posit{ps, es})
    return max(a, b_posit)
end

function max(a::T, b::Posit{ps, es}) where {ps, es, T <: Real}
    return max(b, a)
end
