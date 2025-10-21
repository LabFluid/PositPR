import Base: ==, >, >=, <, <=, isless
function Base.:(==)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    if a.sn == 1 && b.sn == 1
        return a.s == b.s
    end
    for i in 1:nfields(a)
        if getfield(a, i) != getfield(b, i)
            return false
        end
    end
    return true
end

Base.:(==)(::Posit{ps1, es1}, ::Posit{ps2, es2}) where {ps1, es1, ps2, es2} = false

function parse_UInt(s::String)
	return parse(UInt,s, base=2)
end

function Base.:(>)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    if a.sn == 1 && a.s == 1
        return false
    end

    if b.sn == 1 && b.s == 1
        return true
    end

    if a.s != b.s
        return a.s == 0
    end

    if a.sn == 1 && a.s == 0
        return false
    end

    if b.sn == 1 && b.s == 0
        return true
    end

    positives = (a.s == 0)

    if a.k > b.k
        return positives
    elseif a.k == b.k
        if a.e > b.e
            return positives
        elseif a.e == b.e
            parse_UInt(bits::AbstractString) = isempty(bits) ? UInt(0) : parse(UInt, bits; base=2)
            if parse_UInt(a.fi) > parse_UInt(b.fi)
                return positives
            end
        end
    end

    return !positives
end

function Base.:(>=)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    a == b || a > b
end

function Base.:(<)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    !(a >= b)
end

Base.isless(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es} = a < b

function Base.:(<=)(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    a == b || a < b
end

function Base.:(>=)(a::Posit{ps,es}, b::Posit{ps2,es2}) where {ps,es,ps2,es2} 
    if ps > ps2 || (ps == ps2 && es >= es2)
        return a >= x2Posit(b, Posit{ps,es})
    else 
        return x2Posit(a, Posit{ps2,es2}) >= b
    end
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

function useed(es)
	return 2^(2^es)
end

import Base: one
function Base.one(::Type{Posit{ps, es}}) where {ps, es}
    fs = max(ps - 3 - es, 0)
    fi = fs == 0 ? "" : "0"^fs
    ers = min(ps - 3, es)
    return Posit{ps, es}(0, 0, 0, 2, ers, 0, fi, fs, 0.0, useed(es))
end

Base.one(p::Posit{ps, es}) where {ps, es} = one(Posit{ps, es})

import Base: isone
function Base.isone(a::Posit{ps, es}) where {ps, es}
    a == one(Posit{ps, es})
end

import Base: zero
function Base.zero(::Type{Posit{ps, es}}) where {ps, es}
    return Posit{ps, es}(0, 1, 0, 0, 0, 0, "0", 0, 0.0, 2.0^(2^es))
end

Base.zero(p::Posit{ps, es}) where {ps, es} = zero(Posit{ps, es})

import Base: iszero
function Base.iszero(a::Posit{ps, es}) where {ps, es}
    a == zero(Posit{ps, es})
end

import Base: zeros, ones
function Base.zeros(::Type{T}, dims::Integer...) where {ps, es, T <: Posit{ps, es}}
    fill(zero(Posit{ps, es}), dims...)
end

function Base.ones(::Type{T}, dims::Integer...) where {ps, es, T <: Posit{ps, es}}
    fill(one(Posit{ps, es}), dims...)
end

function ispower_oftwo(n::Integer)
    return n > 0 && (n & (n - 1)) == 0
end

function ispower_oftwo(str::String)
	return all(x -> x == '0', str)
end

import Base.:abs
function abs(a::Posit{ps, es}) where {ps, es}
	return (a.s != 0) ? twos_complement(a) : a
end

import Base.eps

function eps(::Type{Posit{ps, es}}) where {ps, es}
    one_posit = one(Posit{ps, es})
    epsilon = one_posit

    while one_posit + epsilon/2 > one_posit
        epsilon /= 2
    end

    return epsilon
end

function eps(a::Posit{ps, es}) where {ps, es}
    if iszero(a)
        return abs(smallest(typeof(a)))
    end

    if a.s == 1  # número negativo
        return -eps(abs(a))
    end

    ε = one(Posit{ps, es})

    while a + ε/2 > a
        ε /= 2
    end

    return ε
end

function nextPosit(a::Posit{ps, es}) where {ps, es}
    if iszero(a) 
        return eps(a)
    end
    if a.s == 1 && ispower_oftwo(a.fi)
        return a + eps(a) / 2
    end
    return a + eps(a)
end

function prevPosit(a::Posit{ps, es}) where {ps, es}
    if iszero(a) 
        return -eps(a)
    end
    if a.s == 0 && ispower_oftwo(a.fi)
        return a - eps(a) / 2
    end
    
    return a - eps(a)
end

isNaR(a::Posit{ps,es}) where {ps,es} = a.sn == 1 && a.s == 1

NaR(::Type{Posit{ps,es}}) where {ps,es} = Posit{ps, es}(1, 1, 0, 0, 0, 0, "0", 0, 0.0, 2.0^(2^es))