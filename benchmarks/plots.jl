# Gráficos dos benchmarks, gerados a partir dos CSV de results/.
# Erros exatamente nulos são desenhados no piso `ZERO_FLOOR` (rotulado "exato").

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")      # GR sem janela
using Plots
using Statistics
isdefined(@__MODULE__, :FORMATS_16) || include(joinpath(@__DIR__, "common.jl"))

const ZERO_FLOOR = 1e-8
const FMT_STYLE = Dict(
    "Float16" => (color = :gray35, ls = :dash, m = :circle),
    "Posit{16,2}" => (color = :royalblue, ls = :solid, m = :diamond),
    "Posit{16,2}+quire" => (color = :darkorange, ls = :dot, m = :utriangle),
    "Float32" => (color = :gray35, ls = :dash, m = :circle),
    "Posit{32,2}" => (color = :royalblue, ls = :solid, m = :diamond),
    "Posit{32,2}+quire" => (color = :darkorange, ls = :dot, m = :utriangle),
)

default(fontfamily = "Computer Modern", framestyle = :box, grid = true, gridalpha = 0.15,
        legendfontsize = 7, guidefontsize = 9, tickfontsize = 7, titlefontsize = 9, dpi = 150)

floor0(v) = (isnan(v) ? NaN : max(v, ZERO_FLOOR))
"Filtra linhas pelos valores dados; descarta problemas comuns descaracterizados."
sel(rows; kw...) = filter(r -> get(r, :dados_comuns_ok, "true") == "true" && all(r[k] == v for (k, v) in kw), rows)
fmtnames(rows) = unique(r[:formato] for r in rows)

"Primeiro j em que o formato arredonda 1+δ para 1."
first_nondistinct(rows, fmt) = (js = [r[:j] for r in rows if r[:formato] == fmt && r[:distingue] == "false"];
                                isempty(js) ? nothing : minimum(js))

"Quantis (q1, mediana, q3) de `v` ignorando NaN; `nothing` se vazio."
function quartis(v)
    w = sort(filter(!isnan, v))
    isempty(w) && return nothing
    return quantile(w, 0.25; sorted = true), quantile(w, 0.5; sorted = true), quantile(w, 0.75; sorted = true)
end

"""
Erro × j para um formato: mediana sobre os ângulos (estrutura rotacionada) ou o
valor único (canônica), com faixa interquartil quando houver mais de um ângulo.
Só entram os j em que a maioria dos ângulos tem solução (pivô não nulo).
"""
function series_by_j(s, col)
    js, med, lo, hi = Float64[], Float64[], Float64[], Float64[]
    for j in sort(unique(r[:j] for r in s))
        v = [floor0(r[col]) for r in s if r[:j] == j]
        count(!isnan, v) * 2 > length(v) || continue
        q = quartis(v)
        push!(js, j); push!(lo, q[1]); push!(med, q[2]); push!(hi, q[3])
    end
    return js, med, lo, hi
end

"Erro direto × log2(1/δ), um painel por (família, solução)."
function plot_error_vs_j(rows; estrutura = "canonica", modo = "original", alg = "lu", col = :erro_orig, title = "")
    panels = []
    for fam in ("diadica", "generica"), sol in ("uns", "alternado")
        sub = sel(rows; estrutura, modo, algoritmo = alg, familia = fam, solucao = sol)
        p = plot(yscale = :log10, ylims = (ZERO_FLOOR / 2, 10), xlims = (0.5, maximum(r[:j] for r in sub) + 0.5),
                 xlabel = "j = log2(1/δ)", ylabel = "erro direto relativo", title = "δ $(fam), x* $(sol)",
                 legend = :bottomright)
        for f in fmtnames(sub)
            s = sel(sub; formato = f)
            st = FMT_STYLE[f]
            js, med, lo, hi = series_by_j(s, col)
            multi = estrutura != "canonica"
            plot!(p, js, med; label = f, color = st.color, ls = st.ls, marker = st.m, ms = 3, lw = 1.5,
                  ribbon = multi ? (med .- lo, hi .- med) : nothing, fillalpha = 0.15)
            if estrutura == "canonica"
                j0 = first_nondistinct(s, f)
                j0 === nothing || vline!(p, [j0 - 0.5]; color = st.color, ls = :dashdot, lw = 0.8, label = "")
            end
        end
        hline!(p, [ZERO_FLOOR]; color = :black, lw = 0.4, label = "")
        annotate!(p, 1, ZERO_FLOOR * 2.5, text("exato", 6, :left))
        push!(panels, p)
    end
    return plot(panels...; layout = (2, 2), size = (900, 650), plot_title = title, plot_titlefontsize = 9)
end

"Erro direto × κ(A): mediana (e quartis) sobre ângulos, famílias e soluções, por j."
function plot_error_vs_kappa(rows; estrutura = "rotacionada", modo = "comum", alg = "lu",
                             us = Dict("Float16" => 2.0^-11, "Posit{16,2}" => 2.0^-12))
    sub = sel(rows; estrutura, modo, algoritmo = alg)
    p = plot(xscale = :log10, yscale = :log10, ylims = (ZERO_FLOOR / 2, 10), xlabel = "κ(A) (norma 2)",
             ylabel = "erro direto relativo", legend = :topleft, size = (650, 450),
             title = "Gauss ($alg), A $(estrutura), dados $(modo): erro × condicionamento")
    for f in fmtnames(sub)
        s = sel(sub; formato = f)
        st = FMT_STYLE[f]
        js, med, lo, hi = series_by_j(s, :erro_orig)
        ks = [median([r[:kappa] for r in s if r[:j] == j]) for j in js]
        plot!(p, ks, med; ribbon = (med .- lo, hi .- med), fillalpha = 0.15, label = f, color = st.color,
              ls = st.ls, marker = st.m, ms = 3)
    end
    ks = 10.0 .^ range(0, 7, length = 50)
    for (f, u) in sort(collect(us))
        plot!(p, ks, min.(ks .* u, 10); color = FMT_STYLE[f].color, ls = :dashdot, lw = 0.8,
              label = "κ·u ($(f), u = 2^$(Int(log2(u))))")
    end
    return p
end

"""
Refinamento iterativo: fração de casos que convergem e mediana das iterações
(entre os que convergem) × j, para cada variante. "Convergiu" = parou porque o
resíduo zerou ou o iterado ficou estacionário E o erro em relação à solução que
a variante persegue (original para ri_orig*, armazenada para ri_arm*) é ≤ `tol`.
"""
function plot_ri_summary(rows; estrutura = "rotacionada", modo = "original", tol = 2.0^-9)
    algs = ("ri_orig", "ri_orig_esc", "ri_arm", "ri_arm_esc")
    top, bot = [], []
    for alg in algs
        sub = sel(rows; estrutura, modo, algoritmo = alg)
        col = startswith(alg, "ri_orig") ? :erro_orig : :erro_arm
        pf = plot(ylims = (-0.05, 1.05), xlabel = "j = log2(1/δ)", ylabel = "fração convergida",
                  title = alg, legend = :bottomleft)
        pi = plot(xlabel = "j = log2(1/δ)", ylabel = "iterações (mediana)", title = alg, legend = :topleft,
                  ylims = (0, RI_MAXIT_PLOT + 1))
        for f in fmtnames(sub)
            s = sel(sub; formato = f)
            st = FMT_STYLE[f]
            js = sort(unique(r[:j] for r in s))
            conv(r) = r[:status] in ("residuo_zero", "estacionario") && r[col] <= tol
            frac = [count(conv, filter(r -> r[:j] == j, s)) / count(r -> r[:j] == j, s) for j in js]
            its = [(c = [r[:iteracoes] for r in s if r[:j] == j && conv(r)]; isempty(c) ? NaN : median(c)) for j in js]
            plot!(pf, js, frac; label = f, color = st.color, ls = st.ls, marker = st.m, ms = 3)
            plot!(pi, js, its; label = f, color = st.color, ls = st.ls, marker = st.m, ms = 3)
        end
        push!(top, pf); push!(bot, pi)
    end
    return plot(top..., bot...; layout = (2, 4), size = (1300, 600),
                plot_title = "Refinamento iterativo (A $(estrutura), dados $(modo)); convergiu = parou com erro ≤ 2^$(Int(log2(tol)))",
                plot_titlefontsize = 10, left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)
end
const RI_MAXIT_PLOT = 30
const NS_NROT_PLOT = 32

"Erro × iteração do refinamento para alguns j."
function plot_ri_history(hist; estrutura = "rotacionada", familia = "generica", solucao = "alternado",
                         modo = "original", js = (6, 9, 10, 11))
    panels = []
    for j in js, alg in ("ri_orig", "ri_orig_esc")
        sub = sel(hist; estrutura, familia, solucao, modo, j, algoritmo = alg)
        p = plot(yscale = :log10, ylims = (ZERO_FLOOR / 2, 10), xlabel = "iteração",
                 ylabel = "erro direto", title = "j = $j, $alg", legend = :topright)
        for f in fmtnames(sub)
            s = sort(sel(sub; formato = f), by = r -> r[:iter])
            st = FMT_STYLE[f]
            plot!(p, [r[:iter] for r in s], floor0.([r[:erro_orig] for r in s]); label = f,
                  color = st.color, ls = st.ls, marker = st.m, ms = 3)
        end
        push!(panels, p)
    end
    return plot(panels...; layout = (length(js), 2), size = (900, 260 * length(js)),
                plot_title = "Refinamento: A $(estrutura), δ $(familia), x* $(solucao), dados $(modo)", plot_titlefontsize = 10)
end

"Gera todos os gráficos do Experimento I em `outdir`."
function plots_near_singular(; dir = joinpath(@__DIR__, "results"), tag = "16")
    rows = read_csv(joinpath(dir, "near_singular_$(tag).csv"))
    hist = read_csv(joinpath(dir, "near_singular_$(tag)_hist.csv"))
    out(n) = joinpath(dir, "fig_I_$(tag)_$(n).png")
    nota = "linha traço-ponto: o formato arredonda 1+δ para 1"
    for est in ("canonica", "rotacionada")
        for modo in ("original", "comum")
            savefig(plot_error_vs_j(rows; estrutura = est, modo,
                                    title = "Gauss sem pivoteamento, A $(est), dados $(modo) — erro vs. solução de referência" *
                                            (est == "canonica" ? "\n($nota)" : "\n(mediana e quartis sobre $(NS_NROT_PLOT) ângulos)")),
                    out("$(est)_erro_vs_j_$(modo)"))
        end
        savefig(plot_error_vs_j(rows; estrutura = est, col = :erro_arm,
                                title = "Gauss sem pivoteamento, A $(est) — erro vs. problema ARMAZENADO (dados original)"),
                out("$(est)_erro_arm_vs_j"))
        savefig(plot_error_vs_j(rows; estrutura = est, alg = "ri_orig_esc",
                                title = "Refinamento (resíduo escalado), A $(est), dados originais — erro vs. problema original"),
                out("$(est)_erro_ri_vs_j"))
    end
    us = tag == "16" ? Dict("Float16" => 2.0^-11, "Posit{16,2}" => 2.0^-12) :
                       Dict("Float32" => 2.0^-24, "Posit{32,2}" => 2.0^-28)
    for modo in ("original", "comum")
        savefig(plot_error_vs_kappa(rows; modo, us), out("rotacionada_erro_vs_kappa_$(modo)"))
    end
    tol = tag == "16" ? 2.0^-9 : 2.0^-21
    for est in ("canonica", "rotacionada")
        savefig(plot_ri_summary(rows; estrutura = est, tol), out("$(est)_ri_convergencia"))
    end
    js = tag == "16" ? (6, 9, 10, 11) : (20, 23, 24, 27)
    savefig(plot_ri_history(hist; js), out("rotacionada_ri_historico"))
    return nothing
end
