import Base: rand, randn
import Random: AbstractRNG, Xoshiro

default_rng = Xoshiro(1234)

"Uniform in [0, 1): draws a Float64 and rounds it; draws that round to 1 are discarded."
function posit_rand(rng::AbstractRNG, ::Type{T}) where {ps, es, T <: Posit{ps, es}}
    while true
        p = x2Posit(Base.rand(rng, Float64), T)
        isone(p) || return p
    end
end

function Base.rand(rng::AbstractRNG, ::Type{T}, dims::Integer...) where {ps, es, T<:Posit{ps, es}}
    result = Array{PositType(ps, es)}(undef, dims...)
    for i in eachindex(result)
        result[i] = posit_rand(rng, T)
    end
    return result
end

function Base.rand(::Type{T}, dims::Integer...) where {ps, es, T <: Posit{ps, es}}
    return rand(default_rng, T, dims...)
end

function Base.rand(::Type{T}) where {ps, es, T <: Posit{ps, es}}
    return posit_rand(default_rng, T)
end

function Base.rand(rng::AbstractRNG, ::Type{T}) where {ps, es, T <: Posit{ps, es}}
    return posit_rand(rng, T)
end

function posit_rand_nonzero(rng::AbstractRNG, ::Type{T}) where {ps, es, T <: Posit{ps, es}}
    p = posit_rand(rng, T)
    while iszero(p)
        p = posit_rand(rng, T)
    end
    return p
end

"Standard normal via the Box-Muller transform, in Posit arithmetic."
function Base.randn(rng::AbstractRNG, ::Type{T}, dims::Integer...) where {ps, es, T <: Posit{ps, es}}
    N = prod(dims)
    result = Vector{PositType(ps, es)}(undef, N)

    two    = x2Posit(2.0, T)
    two_pi = two * pi(T)

    i = 1
    while i ≤ N
        u1 = posit_rand_nonzero(rng, T)
        u2 = posit_rand_nonzero(rng, T)

        r = sqrt(-two * log(u1))
        θ = two_pi * u2

        z0 = r * cos(θ)
        z1 = r * sin(θ)

        result[i] = z0
        if i + 1 ≤ N
            result[i + 1] = z1
        end

        i += 2
    end

    return reshape(result, dims...)
end

function Base.randn(::Type{T}, dims::Integer...) where {ps, es, T <: Posit{ps, es}}
    return randn(default_rng, T, dims...)
end

function Base.randn(::Type{T}) where {ps, es, T <: Posit{ps, es}}
    return randn(default_rng, T, 1)[1]
end

function Base.randn(rng::AbstractRNG, ::Type{T}) where {ps, es, T <: Posit{ps, es}}
    return randn(rng, T, 1)[1]
end
