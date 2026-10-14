# Experimento I — sistemas lineares 2×2 quase singulares
#
#     A_δ = [1  1 ; 1  1+δ],   b = A_δ x*,   det = δ,   κ₂(A_δ) ≈ 4/δ.
#
# Famílias de δ:
#   :diadica   δ = 2^-j               (dados exatos enquanto 1+δ for representável)
#   :generica  δ = 2^-j·(√5-1)/2      (dados exigem arredondamento em qualquer formato)
# Estruturas:
#   :canonica     A_δ como acima: a eliminação é exata (l₂₁ = 1, u₂₂ = (1+δ) - 1 por
#                 Sterbenz), então todo erro vem da representação dos dados;
#   :rotacionada  G(θ)·A_δ e G(θ)·b para NS_NROT ângulos θ ∈ (0, π/2) sorteados com
#                 semente fixa: mesmo κ, mesma solução e mesmo det = δ, mas dados
#                 genéricos — a eliminação sofre cancelamento. Os resultados são
#                 resumidos pela mediana e pelos quartis sobre os ângulos.
# Soluções exatas: :uns = [1, 1] e :alternado = [1, -1].
# Modos de dados:
#   :original  cada formato arredonda o problema construído em BigFloat;
#   :comum     A e b são trocados por valores representáveis em TODOS os formatos
#              comparados; a solução de referência passa a ser a exata desses dados.
#
# Algoritmos (ver linalg.jl):
#   :lu        Gauss sem pivoteamento, tudo em T
#   :lu_pp     Gauss com pivoteamento parcial, tudo em T
#   :ri_*      refinamento iterativo em precisão mista: fatoração, correção e
#              atualização em T; resíduo em BigFloat, contra os dados originais
#              (:ri_orig) ou contra os dados armazenados (:ri_arm); `_esc` indica
#              resíduo escalado por potência de 2 antes de ser arredondado para T.

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "linalg.jl"))

const NS_FAMILIES = (:diadica, :generica)
const NS_STRUCTS = (:canonica, :rotacionada)
const NS_NROT = 32
"Ângulos das rotações (o 1º é π/7; os demais, sorteados com semente fixa)."
const NS_THETAS = let rng = MersenneTwister(314159)
    [BigFloat(pi) / 7; [BigFloat(pi) / 2 * BigFloat(rand(rng)) for _ in 2:NS_NROT]]
end
const NS_XEXACT = (:uns, :alternado)
const NS_MODES = (:original, :comum)
const NS_JS = 1:20
const RI_MAXIT = 30

delta_value(fam::Symbol, j::Integer) =
    fam === :diadica ? BigFloat(2)^(-j) :
    fam === :generica ? BigFloat(2)^(-j) * (sqrt(big(5)) - 1) / 2 :
    error("família desconhecida: $fam")

xexact_value(k::Symbol) = k === :uns ? BigFloat[1, 1] : k === :alternado ? BigFloat[1, -1] :
    error("solução desconhecida: $k")

"Problema em BigFloat: dados, solução de referência (ou `nothing`) e δ efetivo."
struct NSProblem
    A::Matrix{BigFloat}
    b::Vector{BigFloat}
    x::Union{Nothing, Vector{BigFloat}}
    delta::BigFloat
    ok::Bool              # false se o modo :comum não achou dados comuns
end

"Rotação de Givens G(θ) em BigFloat."
givens(θ) = BigFloat[cos(θ) -sin(θ); sin(θ) cos(θ)]

function ns_problem(fam, j, xk; estrutura = :canonica, angulo = 1)
    δ = delta_value(fam, j)
    A = BigFloat[1 1; 1 1+δ]
    estrutura === :rotacionada && (A = givens(NS_THETAS[angulo]) * A)
    estrutura in NS_STRUCTS || error("estrutura desconhecida: $estrutura")
    x = xexact_value(xk)
    return NSProblem(A, matvec(A, x), x, δ, true)
end

"""
Versão do problema com dados representáveis em todos os tipos `Ts`. `ok = false`
se não há ponto fixo comum ou se o arredondamento comum descaracteriza o problema
(|det(A_c) - δ| > δ/2): aí o "mesmo problema" já não é quase singular do mesmo
modo e a comparação perde sentido — essas linhas ficam no CSV mas saem dos gráficos.
"""
function ns_common(P::NSProblem, Ts)
    A, okA = common_value(P.A, Ts)
    b, okb = common_value(P.b, Ts)
    d = det2(A)
    return NSProblem(A, b, solve2_exact(A, b), d, okA && okb && abs(d - P.delta) <= P.delta / 2)
end

# ----------------------------------------------------------------------------
# Refinamento iterativo em precisão mista

"""
    refine(F, bT, Aref, bref, T; fused, scaled, maxit)

Refinamento iterativo: x₀ = solução de `F` com `bT`; a cada passo
`r = bref - Aref·x` em BigFloat, `r` é arredondado para `T` (com `scaled`,
antes divide-se `r` por `2^e`, `‖r‖∞ ∈ [2^e, 2^(e+1))`), resolve-se `F d = r_T`
em `T`, a correção volta à escala com um arredondamento e `x ← x + d` em `T`.

Devolve `(x, status, iters, hist)`, com `hist` o vetor dos iterados. Status:
`:estacionario` (x não mudou), `:ciclo` (x repetiu um iterado anterior),
`:residuo_zero`, `:residuo_underflow` (r ≠ 0 mas r_T = 0), `:max_iter`,
`:nao_finito`, `:pivo_zero`.
"""
function refine(F::LUFact{T}, bT, Aref, bref; fused = false, scaled = false,
                maxit = RI_MAXIT) where {T}
    x, st = lu_solve(F, bT; fused)
    hist = [copy(x)]
    st === :ok || return (x, st, 0, hist)
    status, it = :max_iter, 0
    for k in 1:maxit
        r = bref .- matvec(Aref, exact(x))
        if all(iszero, r)
            status = :residuo_zero; break
        end
        s = scaled ? BigFloat(2)^exponent(maximum(abs, r)) : big(1)
        rT = tofmt(T, r ./ s)
        if all(iszero, rT)
            status = :residuo_underflow; break
        end
        d, st = lu_solve(F, rT; fused)
        if st !== :ok
            status = st; break
        end
        c = scaled ? tofmt(T, exact(d) .* s) : d      # ldexp com um arredondamento
        xn = x .+ c
        it = k
        if any(isbad, xn)
            push!(hist, xn); status = :nao_finito; break
        end
        if xn == x
            status = :estacionario; break
        end
        if any(h -> h == xn, hist)
            push!(hist, xn); status = :ciclo; break
        end
        x = xn
        push!(hist, copy(x))
    end
    return (x, status, it, hist)
end

# ----------------------------------------------------------------------------
# Execução de uma configuração

const NS_ALGS = (:lu, :lu_pp, :ri_orig, :ri_orig_esc, :ri_arm, :ri_arm_esc)

"Métricas de uma solução `x` (no formato T) contra o problema `P` e o armazenado."
function ns_metrics(x, ok::Bool, P::NSProblem, As, bs, xs)
    nanrow = (erro_orig = NaN, erro_arm = NaN, residuo = NaN, eta_orig = NaN, eta_arm = NaN,
              x1 = NaN, x2 = NaN)
    ok || return nanrow
    xe = exact(x)
    return (erro_orig = P.x === nothing ? NaN : f64(forward_error(xe, P.x)),
            erro_arm = xs === nothing ? NaN : f64(forward_error(xe, xs)),
            residuo = f64(residual(P.A, xe, P.b)),
            eta_orig = f64(backward_error(P.A, xe, P.b)),
            eta_arm = f64(backward_error(As, xe, bs)),
            x1 = Float64(xe[1]), x2 = Float64(xe[2]))
end

"""
    ns_run(fmt, P, meta) -> (rows, hist_rows)

Roda todos os algoritmos de `NS_ALGS` no formato `fmt` para o problema `P`.
`meta` é uma NamedTuple com os parâmetros (família, j, solução, modo).
"""
function ns_run(fmt::Fmt, P::NSProblem, meta)
    T, fused = fmt.T, fmt.fused
    AT, bT = tofmt(T, P.A), tofmt(T, P.b)
    As, bs = exact(AT), exact(bT)
    xs = solve2_exact(As, bs)
    base = (; meta..., formato = fmt.name, delta = f64(P.delta),
            delta_arm = f64(det2(As)),
            distingue = !iszero(det2(As)),
            dados_exatos = As == P.A && bs == P.b,
            kappa = f64(cond2(P.A)), kappa_arm = f64(cond2(As)))
    rows, hist_rows = NamedTuple[], NamedTuple[]
    for alg in NS_ALGS
        tr = MagTracker()
        if alg in (:lu, :lu_pp)
            x, st, _ = gauss_solve(AT, bT; pivot = alg === :lu_pp, fused, tracker = tr)
            it, hist = 0, [x]
        else
            F = lu_generic(AT; fused, tracker = tr)
            Aref, bref = alg in (:ri_orig, :ri_orig_esc) ? (P.A, P.b) : (As, bs)
            x, st, it, hist = refine(F, bT, Aref, bref; fused,
                                     scaled = alg in (:ri_orig_esc, :ri_arm_esc))
            foreach(v -> foreach(e -> track!(tr, e), v), hist)
        end
        solved = st in (:ok, :estacionario, :ciclo, :residuo_zero, :residuo_underflow, :max_iter)
        push!(rows, (; base..., algoritmo = String(alg), status = String(st), iteracoes = it,
                     ns_metrics(x, solved, P, As, bs, xs)...,
                     mag_min = tr.lo, mag_max = tr.hi, n_fora_1_16 = tr.outside))
        if alg ∉ (:lu, :lu_pp)
            for (k, h) in enumerate(hist)
                m = ns_metrics(h, !any(isbad, h), P, As, bs, xs)
                push!(hist_rows, (; meta..., formato = fmt.name, algoritmo = String(alg),
                                  iter = k - 1, erro_orig = m.erro_orig, erro_arm = m.erro_arm,
                                  eta_orig = m.eta_orig))
            end
        end
    end
    return rows, hist_rows
end

"""
    run_near_singular(; formats, js, outdir, tag) -> (rows, hist_rows)

Roda a grade completa (famílias × soluções × modos × j × formatos) e grava
`near_singular_<tag>.csv` e `near_singular_<tag>_hist.csv` em `outdir`.
"""
function run_near_singular(; formats = FORMATS_16, js = NS_JS, outdir = joinpath(@__DIR__, "results"),
                           tag = "16", nrot = NS_NROT)
    Ts = unique(f.T for f in formats)
    rows, hist = NamedTuple[], NamedTuple[]
    configs = [(est, ang) for est in NS_STRUCTS for ang in (est === :canonica ? (0:0) : (1:nrot))]
    for (est, ang) in configs, fam in NS_FAMILIES, xk in NS_XEXACT, mode in NS_MODES, j in js
        P0 = ns_problem(fam, j, xk; estrutura = est, angulo = max(ang, 1))
        P = mode === :original ? P0 : ns_common(P0, Ts)
        meta = (estrutura = String(est), angulo = ang, familia = String(fam), j = j, solucao = String(xk), modo = String(mode),
                dados_comuns_ok = P.ok, problema_singular = P.x === nothing)
        for fmt in formats
            r, h = ns_run(fmt, P, meta)
            append!(rows, r)
            ang <= 1 && append!(hist, h)          # histórico só do 1º ângulo
        end
    end
    write_csv(joinpath(outdir, "near_singular_$(tag).csv"), rows)
    write_csv(joinpath(outdir, "near_singular_$(tag)_hist.csv"), hist)
    return rows, hist
end

# ----------------------------------------------------------------------------
# Temporização (separada das medidas de precisão)

"""
    time_near_singular(; formats, js, outdir, tag)

Tempo mediano e alocações de `gauss_solve` (sem rastreamento) em cada formato,
para alguns δ. O refinamento não é cronometrado: seu resíduo é em BigFloat.
"""
function time_near_singular(; formats = FORMATS_16, js = (4, 8),
                            outdir = joinpath(@__DIR__, "results"), tag = "16")
    rows = NamedTuple[]
    for est in NS_STRUCTS, fam in NS_FAMILIES, j in js, fmt in formats, pivot in (false, true)
        P = ns_problem(fam, j, :uns; estrutura = est)
        AT, bT = tofmt(fmt.T, P.A), tofmt(fmt.T, P.b)
        f = (A, b) -> gauss_solve(A, b; pivot, fused = fmt.fused)
        t, al = time_ns_per_call(f, AT, bT)
        push!(rows, (estrutura = String(est), familia = String(fam), j = j, formato = fmt.name,
                     algoritmo = pivot ? "lu_pp" : "lu", tempo_ns = t, alocacoes = al))
    end
    write_csv(joinpath(outdir, "near_singular_$(tag)_tempo.csv"), rows)
    return rows
end
