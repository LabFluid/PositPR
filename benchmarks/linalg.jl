# Álgebra linear densa genérica, escrita com laços para que TODA operação seja
# feita no tipo T dos dados (sem BLAS, sem promoção para Float64).
#
# `fused = true` troca cada expressão c - Σ aₖbₖ por um acúmulo exato no quire
# com um único arredondamento final (só para posits).

"Valor inválido do formato: Inf/NaN em ponto flutuante, NaR em posit."
isbad(x::Posit) = isNaR(x)
isbad(x::AbstractFloat) = !isfinite(x)

"""
    subdot(c, a, b, fused)

`c - Σₖ a[k]·b[k]` no tipo de `c`. Sem fusão: `acc = acc - a[k]*b[k]` (dois
arredondamentos por termo, na ordem natural). Com fusão (posits): quire.
"""
function subdot(c::T, a, b, fused::Bool) where {T}
    if fused && T <: Posit
        q = Quire(c)
        @inbounds for k in eachindex(a, b)
            qmulsub!(q, a[k], b[k])
        end
        return Posit(q)
    end
    acc = c
    @inbounds for k in eachindex(a, b)
        acc = acc - a[k] * b[k]
    end
    return acc
end

"""
    LUFact{T}

Fatoração `P·A = L·U` (Doolittle: `L` com diagonal unitária, guardada abaixo da
diagonal de `LU`). `perm` é a permutação de linhas; `status` é `:ok`,
`:pivo_zero` ou `:nao_finito`.
"""
struct LUFact{T}
    LU::Matrix{T}
    perm::Vector{Int}
    status::Symbol
end

"""
    lu_generic(A; pivot = false, fused = false, tracker = nothing)

Eliminação de Gauss em forma de Doolittle/Crout: cada entrada de `L` e `U` é
`aᵢⱼ - Σ lᵢₖuₖⱼ` (produto interno, que o quire calcula com um arredondamento).
Com `pivot = true`, pivoteamento parcial; empates mantêm a ordem atual.
"""
function lu_generic(A::AbstractMatrix{T}; pivot::Bool = false, fused::Bool = false,
                    tracker = nothing) where {T}
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("matriz não quadrada"))
    M = Matrix{T}(A)
    perm = collect(1:n)
    s = Vector{T}(undef, n)
    for k in 1:n
        # coluna k de U (candidatos a pivô): sᵢ = aᵢₖ - Σ_{m<k} lᵢₘ uₘₖ
        for i in k:n
            s[i] = subdot(M[i, k], view(M, i, 1:k-1), view(M, 1:k-1, k), fused)
            track!(tracker, s[i])
        end
        p = k
        if pivot
            for i in k+1:n
                abs(s[i]) > abs(s[p]) && (p = i)
            end
            if p != k
                M[k, :], M[p, :] = M[p, :], M[k, :]
                s[k], s[p] = s[p], s[k]
                perm[k], perm[p] = perm[p], perm[k]
            end
        end
        M[k, k] = s[k]
        isbad(s[k]) && return LUFact(M, perm, :nao_finito)
        iszero(s[k]) && return LUFact(M, perm, :pivo_zero)
        # linha k de U: uₖⱼ = aₖⱼ - Σ_{m<k} lₖₘ uₘⱼ
        for j in k+1:n
            M[k, j] = subdot(M[k, j], view(M, k, 1:k-1), view(M, 1:k-1, j), fused)
            track!(tracker, M[k, j])
        end
        # multiplicadores lᵢₖ = sᵢ / uₖₖ
        for i in k+1:n
            M[i, k] = s[i] / s[k]
            track!(tracker, M[i, k])
            isbad(M[i, k]) && return LUFact(M, perm, :nao_finito)
        end
    end
    return LUFact(M, perm, :ok)
end

"""
    lu_solve(F, b; fused = false, tracker = nothing)

Resolve `A x = b` com a fatoração `F` (substituições progressiva e regressiva no
tipo de `b`). Devolve `(x, status)`.
"""
function lu_solve(F::LUFact{T}, b::AbstractVector{T}; fused::Bool = false,
                  tracker = nothing) where {T}
    F.status === :ok || return (fill(zero(T), length(b)), F.status)
    LU, n = F.LU, length(b)
    y = b[F.perm]
    for i in 2:n
        y[i] = subdot(y[i], view(LU, i, 1:i-1), view(y, 1:i-1), fused)
        track!(tracker, y[i])
    end
    x = similar(y)
    for i in n:-1:1
        t = subdot(y[i], view(LU, i, i+1:n), view(x, i+1:n), fused)
        x[i] = t / LU[i, i]
        track!(tracker, x[i])
    end
    any(isbad, x) && return (x, :nao_finito)
    return (x, :ok)
end

"Fatora e resolve; devolve `(x, status, F)`."
function gauss_solve(A::AbstractMatrix{T}, b::AbstractVector{T}; pivot = false, fused = false,
                     tracker = nothing) where {T}
    F = lu_generic(A; pivot, fused, tracker)
    x, st = lu_solve(F, b; fused, tracker)
    return x, st, F
end
