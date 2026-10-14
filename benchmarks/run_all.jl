# Executa os benchmarks validados e gera CSV e gráficos em results/.
#   julia --project=benchmarks benchmarks/run_all.jl [--sem-tempo]
#
# Etapas: (1) testes de validação; (2) Experimento I em 16 bits (j = 1:20) e em
# 32 bits (j = 1:32); (3) tempos, numa execução separada; (4) gráficos e resumo.
# Os experimentos II (gradiente conjugado) e III (Newton) entram aqui depois de
# validados.

include(joinpath(@__DIR__, "near_singular.jl"))
include(joinpath(@__DIR__, "plots.jl"))

Random.seed!(20261011)
const OUT = joinpath(@__DIR__, "results")
const TIMING = !("--sem-tempo" in ARGS)

println("== validação")
include(joinpath(@__DIR__, "test", "runtests.jl"))

"""
Resumo por formato. Para cada j toma-se a mediana sobre ângulos e soluções
(falha = erro infinito); reporta-se o primeiro j em que essa mediana passa de 1%,
o primeiro j com matriz armazenada exatamente singular e a mediana dos erros
direto e retroativo para j ≤ 8.
"""
function ns_summary(rows, tag)
    out = NamedTuple[]
    medj(rs, j, col) = median([isnan(r[col]) ? Inf : r[col] for r in rs if r.j == j])
    for est in NS_STRUCTS, fam in NS_FAMILIES, modo in NS_MODES, fmt in unique(r.formato for r in rows)
        sub = filter(r -> r.estrutura == String(est) && r.familia == String(fam) && r.modo == String(modo) &&
                          r.formato == fmt && r.dados_comuns_ok, rows)
        lu = filter(r -> r.algoritmo == "lu", sub)
        ri = filter(r -> r.algoritmo == "ri_orig_esc", sub)
        isempty(lu) && continue
        js = sort(unique(r.j for r in lu))
        firstj(pred) = (k = findfirst(pred, js); k === nothing ? -1 : js[k])
        push!(out, (estrutura = String(est), familia = String(fam), modo = String(modo), formato = fmt,
                    j_max = maximum(js),
                    j_singular = firstj(j -> any(r -> r.j == j && !r.distingue, lu)),
                    j_erro_lu_acima_1pc = firstj(j -> !(medj(lu, j, :erro_orig) <= 1e-2)),
                    j_erro_ri_acima_1pc = firstj(j -> !(medj(ri, j, :erro_orig) <= 1e-2)),
                    mediana_erro_lu_j_1_8 = median([r.erro_orig for r in lu if r.j <= 8]),
                    mediana_eta_lu_j_1_8 = median([r.eta_orig for r in lu if r.j <= 8])))
    end
    write_csv(joinpath(OUT, "near_singular_$(tag)_resumo.csv"), out)
    for r in out
        @printf("  %-11s %-8s %-8s %-18s j≤%2d  singular j=%3d  med. erro LU>1%% j=%3d  RI>1%% j=%3d  med.(j≤8) erro LU %.2e  η %.2e\n",
                r.estrutura, r.familia, r.modo, r.formato, r.j_max, r.j_singular, r.j_erro_lu_acima_1pc,
                r.j_erro_ri_acima_1pc, r.mediana_erro_lu_j_1_8, r.mediana_eta_lu_j_1_8)
    end
    return out
end

"""
Comparação pareada: para cada problema (estrutura, família, solução, modo, j,
ângulo) válido em que os dois formatos resolvem, compara o erro direto de A com o
de B. Reporta a fração de casos em que A erra menos, empata e erra mais, e a
mediana de log2(erro_B / erro_A) (> 0: A melhor; casos com erro nulo entram com
o piso 2^-60). Só j ≤ jmax, para ficar longe da singularidade numérica.
"""
function ns_paired(rows, tag; pairs, alg = "lu", jmax = 8)
    out = NamedTuple[]
    key(r) = (r.estrutura, r.familia, r.solucao, r.modo, r.j, r.angulo)
    for est in NS_STRUCTS, modo in NS_MODES, (fa, fb) in pairs
        sub = filter(r -> r.estrutura == String(est) && r.modo == String(modo) && r.algoritmo == alg &&
                          r.j <= jmax && r.dados_comuns_ok, rows)
        A = Dict(key(r) => r.erro_orig for r in sub if r.formato == fa)
        B = Dict(key(r) => r.erro_orig for r in sub if r.formato == fb)
        ks = [k for k in keys(A) if haskey(B, k) && !isnan(A[k]) && !isnan(B[k])]
        isempty(ks) && continue
        lr = [log2(max(B[k], 2.0^-60) / max(A[k], 2.0^-60)) for k in ks]
        push!(out, (estrutura = String(est), modo = String(modo), A = fa, B = fb, casos = length(ks),
                    A_melhor = count(>(0), lr) / length(ks), empate = count(==(0), lr) / length(ks),
                    B_melhor = count(<(0), lr) / length(ks), mediana_log2_B_sobre_A = median(lr)))
    end
    write_csv(joinpath(OUT, "near_singular_$(tag)_pareado.csv"), out)
    for r in out
        @printf("  %-11s %-8s %-18s vs %-18s n=%4d  A melhor %5.1f%%  empate %5.1f%%  B melhor %5.1f%%  med log2(eB/eA) = %+.2f\n",
                r.estrutura, r.modo, r.A, r.B, r.casos, 100r.A_melhor, 100r.empate, 100r.B_melhor,
                r.mediana_log2_B_sobre_A)
    end
    return out
end

"""
Efeito da precisão cônica: separa os casos (estrutura rotacionada, Gauss, j ≤ 8)
em que TODOS os valores intermediários do posit ficaram em [1/16, 16) — onde
Posit{16,2} tem 1 bit a mais que Float16 — dos casos em que algum saiu dessa
faixa, e compara a vantagem mediana do posit (log2(erro_float/erro_posit)).
"""
function ns_tapered(rows, tag; fnames, jmax = 8)
    out = NamedTuple[]
    key(r) = (r.familia, r.solucao, r.j, r.angulo)
    for modo in NS_MODES
        sub = filter(r -> r.estrutura == "rotacionada" && r.modo == String(modo) && r.algoritmo == "lu" &&
                          r.j <= jmax && r.dados_comuns_ok, rows)
        F = Dict(key(r) => r for r in sub if r.formato == fnames[1])
        P = Dict(key(r) => r for r in sub if r.formato == fnames[2])
        ks = [k for k in keys(P) if haskey(F, k) && !isnan(P[k].erro_orig) && !isnan(F[k].erro_orig)]
        lr(k) = log2(max(F[k].erro_orig, 2.0^-60) / max(P[k].erro_orig, 2.0^-60))
        for (faixa, pred) in (("dentro", k -> P[k].n_fora_1_16 == 0), ("fora", k -> P[k].n_fora_1_16 > 0))
            kk = filter(pred, ks)
            isempty(kk) && continue
            push!(out, (modo = String(modo), intermediarios = faixa, casos = length(kk),
                        mediana_log2_vantagem_posit = median(lr.(kk))))
        end
    end
    write_csv(joinpath(OUT, "near_singular_$(tag)_conica.csv"), out)
    for r in out
        @printf("  %-8s intermediários %-6s de [1/16,16): n=%4d  vantagem mediana do posit = %+.2f bits\n",
                r.modo, r.intermediarios, r.casos, r.mediana_log2_vantagem_posit)
    end
    return out
end

for (tag, formats, js) in (("16", FORMATS_16, 1:20), ("32", FORMATS_32, 1:32))
    println("== Experimento I ($tag bits): precisão")
    rows, _ = run_near_singular(; formats, js, outdir = OUT, tag)
    ns_summary(rows, tag)
    println("== Experimento I ($tag bits): comparação pareada (Gauss sem pivoteamento, j ≤ 8)")
    nm = [f.name for f in formats]
    ns_paired(rows, tag; pairs = [(nm[2], nm[1]), (nm[3], nm[2]), (nm[3], nm[1])])
    println("== Experimento I ($tag bits): precisão cônica (A rotacionada, Gauss, j ≤ 8)")
    ns_tapered(rows, tag; fnames = (nm[1], nm[2]))
    if TIMING
        println("== Experimento I ($tag bits): tempo")
        for r in time_near_singular(; formats, outdir = OUT, tag)
            @printf("  %-11s %-8s j=%2d %-18s %-6s %9.1f ns  %4d alocações\n", r.estrutura, r.familia, r.j, r.formato,
                    r.algoritmo, r.tempo_ns, r.alocacoes)
        end
    end
    println("== gráficos ($tag bits)")
    plots_near_singular(; dir = OUT, tag)
end
println("resultados em ", OUT)
