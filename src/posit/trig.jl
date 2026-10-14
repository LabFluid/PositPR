# Trigonometric functions. As in exp/log, the series are evaluated in a posit with
# GUARD_BITS extra bits (rounding.jl) and the result is rounded back.

const _PI_CACHE = Dict{DataType, Any}()

"π from a BBP-type series, computed (once per type) in the type `T` itself."
function pi_series(::Type{PT}) where {ps, es, PT <: Posit{ps, es}}
    T = PositType(ps, es)
    return get!(_PI_CACHE, T) do
        p = zero(T)
        n = 0
        m_one_fourth = x2Posit(-1 // 4, T)
        tol = eps(T)
        while true
            four_n = 4 * n
            term = ((2 / x2Posit(four_n + 1, T)) +
                    (2 / x2Posit(four_n + 2, T)) +
                    (1 / x2Posit(four_n + 3, T))) * (m_one_fourth ^ n)
            p += term
            abs(term) <= tol && break
            n += 1
        end
        p
    end::T
end

"π rounded to `Posit{ps,es}`."
function pi(::Type{PT}) where {ps, es, PT <: Posit{ps, es}}
    return narrow_posit(pi_series(wide_type(Posit{ps, es})), Posit{ps, es})
end

"Reduce `a` to `[-π/2, π/2]`; returns the reduced argument and whether the sign flips."
function reduce_half_pi(a::Posit{ps, es}) where {ps, es}
    T = Posit{ps, es}
    p = pi_series(T)
    half_pi = p >> 1
    two_pi = p << 1

    # [-π, π]: a - q·2π with q the integer nearest to a/2π
    q = floor(a / two_pi + x2Posit(0.5, T))
    x = iszero(q) ? a : a - q * two_pi

    # [-π/2, π/2]: sin(x ± π) = -sin(x), cos(x ± π) = -cos(x)
    flip = false
    if x > half_pi
        x -= p
        flip = true
    elseif x < -half_pi
        x += p
        flip = true
    end
    return x, flip
end

"Taylor series of sin on [-π/2, π/2]: x - x³/3! + x⁵/5! - …"
function sin_taylor(x::Posit{ps, es}) where {ps, es}
    T = Posit{ps, es}
    iszero(x) && return x
    x2 = x * x
    term = x
    s = x
    n = 1
    while n < 60
        term = -(term * x2 / ((2n) * (2n + 1)))
        new = s + term
        new == s && break
        s = new
        n += 1
    end
    return s
end

"Taylor series of cos on [-π/2, π/2]: 1 - x²/2! + x⁴/4! - …"
function cos_taylor(x::Posit{ps, es}) where {ps, es}
    T = Posit{ps, es}
    x2 = x * x
    term = one(T)
    s = one(T)
    n = 1
    while n < 60
        term = -(term * x2 / ((2n - 1) * (2n)))
        new = s + term
        new == s && break
        s = new
        n += 1
    end
    return s
end

function sin_series(a::Posit)
    x, flip = reduce_half_pi(a)
    s = sin_taylor(x)
    return flip ? -s : s
end

function cos_series(a::Posit)
    x, flip = reduce_half_pi(a)
    c = cos_taylor(x)
    return flip ? -c : c
end

import Base: sin
function sin(a::Posit{ps, es}) where {ps, es}
    (isNaR(a) || iszero(a)) && return a
    return narrow_posit(sin_series(widen_posit(a)), Posit{ps, es})
end

import Base: cos
function cos(a::Posit{ps, es}) where {ps, es}
    isNaR(a) && return a
    return narrow_posit(cos_series(widen_posit(a)), Posit{ps, es})
end

import Base: tan
function tan(a::Posit{ps, es}) where {ps, es}
    (isNaR(a) || iszero(a)) && return a
    w = widen_posit(a)
    x, _ = reduce_half_pi(w)
    s, c = sin_taylor(x), cos_taylor(x)     # tan(x ± π) = tan(x)

    # points where tan is undefined: cos(a) ≈ 0
    if iszero(c)
        error("tan(a) is undefined for a ≈ (2n+1)·π/2")
    end

    return narrow_posit(s / c, Posit{ps, es})
end

import Base: sincos
sincos(a::Posit{ps, es}) where {ps, es} = (sin(a), cos(a))

import Base: tanh
"tanh(x) = t/(t+2) with t = e^(2x) - 1 (no cancellation near 0)."
function Base.tanh(a::Posit{ps, es}) where {ps, es}
    (isNaR(a) || iszero(a)) && return a
    T = Posit{ps, es}
    w = widen_posit(abs(a))
    W = typeof(w)
    if w < x2Posit(0.5, W)
        t = expm1_series(w << 1)
        r = t / (t + 2)
    else
        r = one(W) - 2 / (exp_series(w << 1) + 1)
    end
    r = narrow_posit(r, T)
    return a.s == 1 ? -r : r
end
