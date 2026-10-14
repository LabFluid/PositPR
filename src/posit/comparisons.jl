# Comparisons. Bit patterns read as signed (two's complement) integers are in the
# same order as the values they represent, and NaR (100…0) is the most negative
# pattern, so comparing posits is comparing integers. NaR == NaR and NaR is
# smaller than every number, as in the Posit standard.

import Base: ==, >, >=, <, <=, isless

"Bit pattern as a signed integer (two's complement over ps bits)."
@inline function signed_bits(p::Posit{ps, es}) where {ps, es}
    u = bits(p)
    ps == 128 && return reinterpret(Int128, u)
    return u >= nar_bits(ps) ? Int128(u) - (Int128(1) << ps) : Int128(u)
end

Base.:(==)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es} = bits(a) == bits(b)
Base.:(==)(::Posit{ps1, es1}, ::Posit{ps2, es2}) where {ps1, es1, ps2, es2} = false

Base.:(<)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es} = signed_bits(a) < signed_bits(b)
Base.:(>)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es} = signed_bits(a) > signed_bits(b)
Base.:(<=)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es} = signed_bits(a) <= signed_bits(b)
Base.:(>=)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es} = signed_bits(a) >= signed_bits(b)
Base.isless(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es} = a < b

function Base.:(>=)(a::Posit{ps,es}, b::Posit{ps2,es2}) where {ps,es,ps2,es2}
    if ps > ps2 || (ps == ps2 && es >= es2)
        return a >= x2Posit(b, Posit{ps,es})
    else
        return x2Posit(a, Posit{ps2,es2}) >= b
    end
end

# Posits of different types: compare in the wider format (PositConverter).
for op in (:<, :>, :<=, :isless)
    @eval Base.$op(a::Posit{ps,es}, b::Posit{ps2,es2}) where {ps,es,ps2,es2} =
        Base.$op(PositConverter(a, b)...)
end

<(a::Real, b::Posit{ps,es}) where {ps,es} = Posit{ps,es}(a) < b
>(a::Real, b::Posit{ps,es}) where {ps,es} = Posit{ps,es}(a) > b
<=(a::Real, b::Posit{ps,es}) where {ps,es} = Posit{ps,es}(a) <= b
Base.:(>=)(a::Real, b::Posit{ps,es}) where {ps,es} = Posit{ps,es}(a) >= b
Base.isless(a::Real, b::Posit{ps, es}) where {ps, es} = Posit{ps,es}(a) < b

<(a::Posit{ps,es}, b::Real) where {ps,es} = a < Posit{ps, es}(b)
>(a::Posit{ps,es}, b::Real) where {ps,es} = a > Posit{ps, es}(b)
<=(a::Posit{ps,es}, b::Real) where {ps,es} = a <= Posit{ps, es}(b)
Base.:(>=)(a::Posit{ps,es}, b::Real) where {ps,es} = a >= Posit{ps, es}(b)
Base.isless(a::Posit{ps, es}, b::Real) where {ps, es} = a < Posit{ps, es}(b)

# Base defines isless(::Real, ::AbstractFloat) and isless(::AbstractFloat, ::Real);
# these methods resolve the ambiguity with Posit.
Base.isless(a::Posit{ps, es}, b::AbstractFloat) where {ps, es} = a < Posit{ps, es}(b)
Base.isless(a::AbstractFloat, b::Posit{ps, es}) where {ps, es} = Posit{ps, es}(a) < b

function useed(es)
	return 2^(2^es)
end

# ---------------------------------------------------------------- constants

import Base: one, zero, isone, iszero, typemin, typemax
Base.one(::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = from_bits(PT, UInt128(1) << (ps - 2))
Base.one(p::Posit{ps, es}) where {ps, es} = one(Posit{ps, es})
Base.isone(a::Posit{ps, es}) where {ps, es} = bits(a) == UInt128(1) << (ps - 2)

Base.zero(::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = from_bits(PT, UInt128(0))
Base.zero(p::Posit{ps, es}) where {ps, es} = zero(Posit{ps, es})
Base.iszero(a::Posit) = iszero(bits(a))

"Smallest positive posit (minpos)."
smallest(::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = from_bits(PT, UInt128(1))

"Most negative finite posit (-maxpos)."
typemin(::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = from_bits(PT, nar_bits(ps) + UInt128(1))

"Largest posit (maxpos)."
typemax(::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = from_bits(PT, maxpos_bits(ps))

"Not a Real: the single exception value of the posit system."
NaR(::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = from_bits(PT, nar_bits(ps))
isNaR(a::Posit{ps,es}) where {ps,es} = bits(a) == nar_bits(ps)

import Base: zeros, ones
function Base.zeros(::Type{T}, dims::Integer...) where {ps, es, T <: Posit{ps, es}}
    fill(zero(Posit{ps, es}), dims...)
end

function Base.ones(::Type{T}, dims::Integer...) where {ps, es, T <: Posit{ps, es}}
    fill(one(Posit{ps, es}), dims...)
end

import Base.:abs
abs(a::Posit{ps, es}) where {ps, es} = sign_bit(a) == 1 ? -a : a   # -NaR == NaR

# ------------------------------------------------------------- neighbours

"""
    nextPosit(a)

Next posit greater than `a`: adds 1 to the bit pattern (patterns are in the same
order as the values). Saturates at `maxpos`; `nextPosit(NaR) == NaR`.
"""
function nextPosit(a::Posit{ps, es}) where {ps, es}
    (isNaR(a) || bits(a) == maxpos_bits(ps)) && return a
    return from_bits(Posit{ps, es}, (bits(a) + UInt128(1)) & posit_mask(ps))
end

"Posit preceding `a` (subtracts 1 from the pattern). Saturates at `-maxpos`."
function prevPosit(a::Posit{ps, es}) where {ps, es}
    (isNaR(a) || bits(a) == nar_bits(ps) + UInt128(1)) && return a
    return from_bits(Posit{ps, es}, (bits(a) - UInt128(1)) & posit_mask(ps))
end

import Base.eps

eps(::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = eps(one(PT))

"""
    eps(a)

Distance from `a` to the next posit of larger magnitude (with the sign of `a`),
rounded to the format. At `maxpos` the distance to the previous posit is used.
"""
function eps(a::Posit{ps, es}) where {ps, es}
    iszero(a) && return smallest(typeof(a))
    isNaR(a) && return a
    m = abs(a)
    nxt = nextPosit(m)
    d = nxt == m ? m - prevPosit(m) : nxt - m
    return sign_bit(a) == 1 ? -d : d
end

import Base: decompose
"`x == num · 2^pow / den`, used by Julia to compare and hash reals of different types."
function decompose(x::Posit{ps, es}) where {ps, es}
    isNaR(x) && return Int128(0), 0, 0
    iszero(x) && return Int128(0), 0, 1
    s, M, P = exact_parts(x)
    return (s == 1 ? -Int128(M) : Int128(M)), P, 1
end
