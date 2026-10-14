# Experimento I — sistemas 2×2 quase singulares: resultados

Gerado por `run_all.jl` (log completo em `results/run_all.log`). Todos os números
abaixo vêm dos CSV em `results/`; nada foi ajustado à mão.

## Desenho

| fator | níveis |
|---|---|
| estrutura | `canonica` `[1 1; 1 1+δ]`; `rotacionada` `G(θ)·A_δ`, 32 ângulos (semente fixa) |
| família de δ | `diadica` 2^-j; `generica` 2^-j·(√5−1)/2; j = 1…20 (16 bits), 1…32 (32 bits) |
| solução | `[1, 1]`, `[1, −1]` |
| dados | `original` (cada formato arredonda o problema), `comum` (mesmos dados, representáveis em todos) |
| formato | Float16, Posit{16,2}, Posit{16,2}+quire (e Float32, Posit{32,2}, Posit{32,2}+quire) |
| algoritmo | Gauss sem/com pivoteamento; refinamento misto (resíduo em BigFloat) contra dados originais ou armazenados, com/sem escala do resíduo |

Validação: `test/runtests.jl`, 2111 verificações contra BigFloat (conversão sem
arredondamento duplo, métricas, álgebra genérica, `subdot` com quire = um
arredondamento do valor exato, previsões analíticas da forma canônica).

## Observações

**O1 — Forma canônica: a eliminação é exata; só os dados importam.**
l₂₁ = 1 e u₂₂ = (1+δ) − 1 é exato, então o erro direto é o da representação de
A e b. Família diádica, x* = [1,1]: Float16 é exato até j = 9, erra 100% em
j = 10 (b₂ = 2+2^-10 não é representável) e tem pivô nulo a partir de j = 11;
Posit{16,2} é exato até j = 10, erra em j = 11, pivô nulo a partir de j = 12.
O quire produz resultados **idênticos** ao posit simples em 100% dos casos
canônicos (nada a arredondar). Pivoteamento parcial não muda nada (empate
|a₂₁| = |a₁₁|), como previsto. Em 32 bits: Float32 até j = 23, Posit{32,2} até
j = 27 (27 vs 23 bits de fração perto de 1).

**O2 — Forma rotacionada: o erro cresce com κ, abaixo de κ·u.**
`fig_I_16_rotacionada_erro_vs_kappa_*.png`: a mediana do erro acompanha κ·u com
cerca de uma ordem de grandeza de folga até o problema armazenado ficar
numericamente singular (κ ≈ 10⁴ em 16 bits).

**O3 — Comparação pareada (Gauss sem pivoteamento, A rotacionada, j ≤ 8, n = 1024 pares).**

| 16 bits | dados | A melhor | empate | B melhor | mediana log₂(e_B/e_A) |
|---|---|---|---|---|---|
| Posit vs Float16 | original | 67.5% | 7.1% | 25.4% | +1.00 |
| Posit vs Float16 | comum | 62.1% | 9.7% | 28.2% | +0.58 |
| Posit+quire vs Posit | original | 47.2% | 12.9% | 39.9% | 0.00 |
| Posit+quire vs Posit | comum | 59.8% | 12.8% | 27.4% | +0.74 |

| 32 bits | dados | A melhor | empate | B melhor | mediana log₂(e_B/e_A) |
|---|---|---|---|---|---|
| Posit vs Float32 | original | 83.6% | 2.4% | 14.0% | +4.00 |
| Posit vs Float32 | comum | 92.0% | 2.8% | 5.2% | +4.04 |
| Posit+quire vs Posit | original | 49.7% | 12.4% | 37.9% | 0.00 |
| Posit+quire vs Posit | comum | 57.2% | 12.1% | 30.7% | +0.46 |

**O4 — Precisão cônica medida.** Separando os pares (16 bits, dados comuns)
pelo rastreamento das magnitudes intermediárias do posit:
todos em [1/16, 16) → vantagem mediana do posit **+0.87 bit** (n = 501);
algum fora → **+0.19 bit** (n = 523). Com dados originais: +1.00 vs +0.94.
Em 32 bits: +4.14 vs +4.00 (Posit{32,2} ainda tem ≥ 26 bits até 2^-8, contra
23 do Float32).

**O5 — Pivô nulo e quire.** Com dois arredondamentos, u₂₂ = fl(a₂₂ − fl(l₂₁a₁₂))
pode dar exatamente 0 por cancelamento (j ≥ 11, A rotacionada: 1048 casos em
Float16, 982 em Posit). Com o quire, u₂₂ é o valor exato arredondado uma vez, e
o posit nunca arredonda um não nulo para 0: **nenhum** pivô nulo. A fatoração
segue, mas com u₂₂ dominado pelo erro dos dados.

**O6 — Refinamento iterativo em precisão mista (A rotacionada, dados originais).**
Fração de casos que param com erro ≤ 2^-9 em relação à solução perseguida:

| variante | formato | j=9 | j=10 | j=11 | j=12 | j=13 |
|---|---|---|---|---|---|---|
| ri_orig_esc | Float16 | 0.98 | 0.86 | 0.64 | 0.47 | 0.18 |
| ri_orig_esc | Posit{16,2} | 1.00 | 0.98 | 0.85 | 0.58 | 0.30 |
| ri_orig_esc | Posit{16,2}+quire | 1.00 | 0.95 | 0.84 | 0.73 | 0.44 |
| ri_arm_esc | Float16 | 0.91 | 0.81 | 0.61 | 0.41 | 0.27 |
| ri_arm_esc | Posit{16,2} | 0.92 | 0.91 | 0.80 | 0.56 | 0.28 |
| ri_arm_esc | Posit{16,2}+quire | 0.97 | 0.97 | 0.88 | 0.66 | 0.30 |

Mediana de iterações entre os que convergem, j = 10…12: Float16 4 / 7.5 / 14,
Posit 4 / 5 / 7.5, Posit+quire 4 / 5 / 7. Escalar o resíduo muda pouco aqui:
os resíduos deste problema não chegam à região de underflow do Float16. Para o
posit, a escala ajuda um pouco em j = 11…13 (por exemplo, 0.52 → 0.58 em j = 12).
Os casos com quire que não convergem terminam mais em `max_iter` e `ciclo`
(O5: não param mais em pivô nulo).

**O7 — Tempo (mediana, BenchmarkTools, solução 2×2 sem diagnóstico).**
Float16 160 ns; Posit{16,2} 385 ns (×2.4); Posit{16,2}+quire 3.9 µs (×24 o
posit; 130 alocações, acumulador BigInt). 32 bits: Float32 125 ns; Posit{32,2}
385 ns; com quire 14 µs. É emulação em software contra hardware/emulação via
Float32. Isso compara implementações, não formatos.

## Respostas às perguntas (restritas ao Experimento I)

1. **Onde Posit{16,2} erra menos que Float16?** Em todo o intervalo de κ em que
   ambos resolvem (κ ≲ 10³–10⁴), na mediana. Não é uma vitória caso a caso:
   o Float16 erra menos em ~25–28% dos pares. Além disso, o posit distingue
   1+δ de 1 por um j a mais (2^-11 contra 2^-10).
2. **Representação ou operações?** As duas coisas, em proporções diferentes.
   Com dados originais a vantagem mediana é de 1 bit e o quire não acrescenta
   nada: o erro dos dados domina. Com dados comuns sobram 0.58 bit da
   aritmética do posit e mais 0.74 bit do quire. Na forma canônica, toda a
   diferença é representação.
3. **Onde o posit vai pior?** Nos pares em que algum intermediário sai de
   [1/16, 16): a vantagem cai para 0.19 bit com dados comuns. Em ~25% dos pares,
   o Float16 erra menos. O Experimento I não tem magnitudes extremas: os dados
   são O(1) e δ só aparece como diferença. Por isso a região onde o Float16
   ganha nitidamente (|x| ≳ 2^8 ou ≲ 2^-8) quase não aparece aqui; os
   experimentos II e III é que a exercitam.
4. **Intermediários perto/longe de 1:** ver O4. O efeito é mensurável em 16 bits
   e desprezível em 32 bits nesta faixa de magnitudes.
5. **O quire reduz o erro?** Na eliminação, só quando os dados são iguais
   (+0.74 bit em 16, +0.46 em 32) e nunca na forma canônica. No refinamento,
   aumenta a fração de convergência perto da singularidade numérica (O6), em
   parte porque não gera pivô nulo (O5). Com n = 2, cada produto interno tem um
   termo: este experimento não mede o ganho de acumulação longa do quire.
6. **Troca precisão × iterações × custo:** o posit converge em menos iterações
   no refinamento (mediana 7.5 contra 14 em j = 12) e custa ×2.4 por solução.
   O quire custa ×24 a mais nesta implementação (BigInt) por ganho de fração de
   bit. Uma conclusão de custo exige uma implementação de quire com tamanho fixo.
7. **32 bits?** O padrão se mantém e fica mais nítido: +4 bits em mediana (= a
   diferença nominal de bits de fração perto de 1). O efeito da precisão cônica
   desaparece nesta faixa de magnitudes.

## Limitações

* n = 2: não há acumulação de erro em produtos internos longos, e o quire é
  testado no seu caso mais fraco.
* Não há variante "Float16 + fma": `fma(::Float16…)` do Julia passa por Float32
  e arredonda de novo. Um fma correto exigiria implementação própria.
* O modo `comum` só cobre j em que os dados comuns preservam o problema
  (j ≤ 9–14 conforme estrutura/família; colunas `dados_comuns_ok`, `j_max`).
* A PositPR lança `DivideError` em x/0 (o padrão pede NaR). Os algoritmos testam
  o pivô antes, então o resultado não é afetado.
* Tempos de uma única máquina (WSL2/Windows, Julia 1.11.5), com ruído de ±10–20%
  entre execuções.
