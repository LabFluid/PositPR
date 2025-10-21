function pi(::Type{Posit{ps, es}}) where {ps, es}
    π = zero(Posit{ps, es})
    term = one(Posit{ps, es})

    n = 0

    n_one = x2Posit(1, Posit{ps, es})
    n_two = x2Posit(2, Posit{ps, es})
    m_one_fourth = x2Posit(-1/4, Posit{ps, es})

    while abs(term) > eps(Posit{ps, es})
        four_n = x2Posit(4 * n, Posit{ps, es})
        denom1 = four_n + one(Posit{ps, es})
        denom2 = four_n + x2Posit(2, Posit{ps, es})
        denom3 = four_n + x2Posit(3, Posit{ps, es})

        term = ((n_two / denom1) +
                (n_two / denom2) +
                (n_one / denom3)) * (m_one_fourth ^ n)

        π += term
        n += 1
    end

    return π
end

function reduce_to_pi(a::Posit{ps, es}) where {ps, es}
    pi = pi(Posit{ps, es})
    two_pi = pi * x2Posit(2, Posit{ps, es})
    r = mod(a + pi, two_pi)
    return r - pi
end

function reduce_to_pi_sym(a::Posit{ps, es}) where {ps, es}
    T = Posit{ps, es}
    pi = pi(T)
    two_pi = pi * x2Posit(2, T)

    q = floor(a / two_pi + x2Posit(0.5, T))  # arredonda para o inteiro mais próximo
    return a - q * two_pi
end

import Base: sin
function sin(a::Posit{ps, es}) where {ps, es}
    pi = pi(Posit{ps, es})
    half_pi = pi / x2Posit(2, Posit{ps, es})
    return cos(a - half_pi)
end

import Base: cos
function cos(a::Posit{ps, es}) where {ps, es}
    pi = pi(Posit{ps, es})
    half_pi = pi / x2Posit(2, Posit{ps, es})

    # [-π, π] reduction
    temp_a = reduce_to_pi_sym(a)

    # [-π/2, π/2] symmetry
    flip_sign = false
    if temp_a > half_pi
        temp_a -= pi
        flip_sign = true
    elseif temp_a < -half_pi
        temp_a += pi
        flip_sign = true
    end

    # Taylor series: cos(x) = 1 - x^2/2! + x^4/4! - ...
    term = one(Posit{ps, es})
    cos_a = term
    a_squared = temp_a ^ 2

    ε = eps(cos_a)
    n = 1

    while abs(term) > abs(ε) && n < 20
        denom = x2Posit((2n - 1) * (2n), Posit{ps, es})
        term *= -a_squared / denom
        cos_a += term
        n += 1
    end

    return flip_sign ? -cos_a : cos_a
end
import Base: tan

function tan(a::Posit{ps, es}) where {ps, es}
    sin_a = sin(a)
    cos_a = cos(a)

    # points where tan is undefined: cos(a) ≈ 0
    if abs(cos_a) < eps(cos_a)
        error("tan(a) is undefined for a ≈ (2n+1)·π/2")
    end

    return sin_a / cos_a
end

import Base: sincos
sincos(a::Posit{ps, es}) where {ps, es} = (sin(a), cos(a))

import Base: tanh
Base.tanh(x::Posit{ps, es}) where {ps, es} = 
    (exp(x) - exp(-x)) / (exp(x) + exp(-x))