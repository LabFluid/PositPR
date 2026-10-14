"""
    x2Posit(x, T)

Convert a real number to `T = Posit{ps,es}` with correct rounding (RNE).
The value of `x` is decomposed EXACTLY as `(-1)^s · M · 2^P` (integer M) and
rounded once by `posit_round`. Out-of-range values saturate at
±maxpos / ±minpos; `NaN` and `±Inf` become NaR.
"""
function x2Posit(x::Real, ::Val{ps}, ::Val{es}) where {ps, es}
    T = Posit{ps, es}

    # infinity / NaN
    if x isa AbstractFloat && (isnan(x) || isinf(x))
        return NaR(T)
    end

    # zero
    iszero(x) && return zero(T)

    s = x < 0 ? 1 : 0

    if x isa Integer
        return posit_round(T, s, fit_int(abs(big(x))), 0)
    elseif x isa Rational
        return posit_round_div(T, s, abs(numerator(x)), 0, abs(denominator(x)), 0)
    elseif x isa Union{Float16, Float32, Float64}
        # IEEE floats: y = 1.f · 2^E with a 53-bit integer significand, exact
        y = abs(Float64(x))
        E = exponent(y)
        M = UInt128(Int64(ldexp(significand(y), 52)))
        return posit_round(T, s, M, E - 52)
    end

    # Other reals (BigFloat, π, ℯ, ...): decomposition in BigFloat. Irrationals
    # are approximated with many more bits than the posit has.
    prec = x isa BigFloat ? precision(x) : max(256, 4ps)
    y = setprecision(() -> abs(BigFloat(x)), BigFloat, prec)
    E = exponent(y)                              # y = 1.f · 2^E
    M = BigInt(ldexp(y, prec - 1 - E))           # integer mantissa, exact
    return posit_round(T, s, fit_int(M), E - (prec - 1))
end

x2Posit(x::Real, ::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = x2Posit(x, Val(ps), Val(es))
x2Posit(x::Posit, ::Val{ps}, ::Val{es}) where {ps, es} = Posit2Posit(x, Val(ps), Val(es))

function Posit2x(p::Posit) # float representation of a Posit
    if is_special(p)
        return sign_bit(p) == 0 ? 0.0 : NaN     # zero or NaR
    end
    s, M, P = exact_parts(p)
    v = ldexp(Float64(M), P)
    return s == 1 ? -v : v
end

"Exact value of the posit as a BigFloat."
function Base.BigFloat(p::Posit)
    isNaR(p) && return BigFloat(NaN)
    s, M, P = exact_parts(p)
    v = ldexp(setprecision(() -> BigFloat(M), BigFloat, max(256, bitlength(M))), P)
    return s == 1 ? -v : v
end

import Base: convert
import Base: promote_rule
# Promotion always goes to the concrete type Posit{ps, es, U}.
promote_rule(::Type{PT}, ::Type{Bool}) where {ps, es, PT <: Posit{ps, es}} = PositType(ps, es)
promote_rule(::Type{PT}, ::Type{T}) where {ps, es, PT <: Posit{ps, es}, T <: AbstractFloat} = PositType(ps, es)
promote_rule(::Type{PT}, ::Type{T}) where {ps, es, PT <: Posit{ps, es}, T <: Integer} = PositType(ps, es)
promote_rule(::Type{PT}, ::Type{T}) where {ps, es, PT <: Posit{ps, es}, T <: AbstractIrrational} = PositType(ps, es)
promote_rule(::Type{PT}, ::Type{BigFloat}) where {ps, es, PT <: Posit{ps, es}} = PositType(ps, es)

convert(::Type{PT}, x::Real) where {ps, es, PT <: Posit{ps, es}} = x2Posit(x, Val(ps), Val(es))
convert(::Type{PT}, x::Bool) where {ps, es, PT <: Posit{ps, es}} = x ? one(PT) : zero(PT)
convert(::Type{PT}, x::PT) where {ps, es, PT <: Posit{ps, es}} = x

# Constructors from real numbers, for the family (Posit{16,2}(1.5)) and for the
# concrete type (Posit{16,2,UInt16}(1.5) or Posit16{2}(1.5)).
Posit{ps,es}(x::Real) where {ps,es} = x2Posit(x, Val(ps), Val(es))
Posit{ps,es,U}(x::Real) where {ps,es,U} = x2Posit(x, Val(ps), Val(es))
Posit{ps,es}(x::Posit{ps,es}) where {ps,es} = x
Posit{ps,es,U}(x::Posit{ps,es,U}) where {ps,es,U} = x

"""
    Posit2Posit(p, ps_new, es_new)

Convert between configurations. The exact value of `p` is rounded (RNE) to the
new format: exact when widening, correctly rounded when narrowing.
"""
function Posit2Posit(p0::Posit{ps0, es0}, ::Val{ps_new}, ::Val{es_new}) where {ps0, es0, ps_new, es_new}
    T = Posit{ps_new, es_new}
    isNaR(p0) && return NaR(T)
    iszero(p0) && return zero(T)
    s, M, P = exact_parts(p0)
    return posit_round(T, s, M, P)
end

Posit2Posit(p::Posit{ps0, es0}, ps_new::Integer, es_new::Integer) where {ps0, es0} = Posit2Posit(p, Val(ps_new), Val(es_new))

function PositConverter(p1::Posit{ps1, es1}, p2::Posit{ps2, es2}) where {ps1, es1, ps2, es2}
    if ps1 > ps2
        # Convert p2 to the type of p1
        p2 = Posit2Posit(p2, Val(ps1), Val(es1))
    elseif ps1 < ps2 || es1 != es2
        # Convert p1 to the type of p2
        p1 = Posit2Posit(p1, Val(ps2), Val(es2))
    end

    return p1, p2
end

import Base: Int

"Nearest integer (ties away from zero)."
function Int(a::Posit{ps, es}) where {ps, es}
    iszero(a) && return 0
    isNaR(a) && throw(InexactError(:Int, Int, a))
    s, M, P = exact_parts(a)
    M = big(M)
    if P >= 0
        v = M << P
    else
        d = big(1) << -P
        q, r = divrem(M, d)
        v = 2r >= d ? q + 1 : q
    end
    return Int(s == 1 ? -v : v)
end

import Base: Integer
Base.Integer(a::Posit{ps,es}) where {ps,es} = Int(a)

convert(::Type{Int}, a::Posit{ps, es}) where {ps, es} = Int(a)

"Largest integer ≤ a, as an Int."
function int_floor(a::Posit{ps, es}) where {ps, es}
    iszero(a) && return 0
    s, M, P = exact_parts(a)
    v = s == 1 ? -big(M) : big(M)
    return Int(P >= 0 ? v << P : fld(v, big(1) << -P))
end

import Base: real
Base.real(p::Posit{ps,es}) where {ps,es} = p

import Base: Float64
Float64(p::Posit) = Posit2x(p)

import Base: float
float(p::Posit) = Posit2x(p)
