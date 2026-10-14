# Benchmarks de estabilidade numérica — PositPR.jl

Compara Float16, Posit{16,2} e Posit{16,2} com quire (e as versões de 32 bits)
contra uma referência BigFloat de 256 bits.

* `FUNDAMENTACAO.md`: modelo de arredondamento, fenômenos esperados e cuidados
  para uma comparação justa.
* `RESULTADOS_I.md`: resultados e análise do Experimento I.

## Executar

```sh
julia --project=benchmarks -e 'using Pkg; Pkg.instantiate()'   # uma vez
julia --project=benchmarks benchmarks/test/runtests.jl          # validação (~10 s)
julia --project=benchmarks benchmarks/run_all.jl                # tudo (~2,5 min)
julia --project=benchmarks benchmarks/run_all.jl --sem-tempo    # sem cronometragem
```

## Arquivos

| arquivo | conteúdo |
|---|---|
| `common.jl` | formatos, conversão com um arredondamento, dados comuns, métricas em BigFloat, temporização, CSV |
| `linalg.jl` | LU genérica (Doolittle/Crout, pivoteamento opcional, `subdot` com quire) |
| `near_singular.jl` | Experimento I: problemas, refinamento misto, execução, tempos |
| `plots.jl` | gráficos a partir dos CSV |
| `run_all.jl` | validação → precisão → resumos (pareado, precisão cônica) → tempo → gráficos |
| `test/runtests.jl` | validação contra BigFloat |
| `results/` | `near_singular_{16,32}*.csv`, `fig_I_*.png`, `run_all.log` |

Os experimentos II (gradiente conjugado) e III (Newton) serão acrescentados depois
da validação do Experimento I.
