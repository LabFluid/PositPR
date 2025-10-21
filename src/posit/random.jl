import Base: rand, randn
import Random: AbstractRNG, Xoshiro

default_rng = Xoshiro(1234)

function Base.rand(rng::AbstractRNG, ::Type{T}, dims::Integer...) where {ps, es, T<:Posit{ps, es}}
    float_vals = Base.rand(rng, Float64, dims...)
    return map(x -> x2Posit(x, T), float_vals)
end

function Base.rand(::Type{T}, dims::Integer...) where {ps, es, T <: Posit{ps, es}}
    return rand(default_rng, T, dims...)
end

function Base.rand(::Type{T}) where {ps, es, T <: Posit{ps, es}}
    return x2Posit(Base.rand(default_rng, Float64), T)
end

function posit_rand_nonzero(rng::AbstractRNG, ::Type{T}) where {ps, es, T <: Posit{ps, es}}
    p = rand(rng, T)
    while iszero(p)
        p = rand(rng, T)
    end
    return p
end

function Base.randn(rng::AbstractRNG, ::Type{T}, dims::Integer...) where {ps, es, T <: Posit{ps, es}}
    N = prod(dims)
    result = Vector{T}(undef, N)

    two    = x2Posit(2.0, T)
    pi     = pi(T)
    two_pi = two * pi

    i = 1
    while i ≤ N
        u1 = posit_rand_nonzero(rng, T)
        u2 = posit_rand_nonzero(rng, T)

        r = sqrt(-two * ln(u1))
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
