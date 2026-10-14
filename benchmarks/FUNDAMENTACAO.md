# Benchmarks de estabilidade numérica — fundamentação

Este documento fixa o que os benchmarks medem, por que, e quais cuidados tornam
a comparação entre `Float16`, `Posit{16,2}` e `Posit{16,2}` com quire justa.
A hipótese em teste é **neutra**: a precisão relativa maior do posit perto de 1
*pode* reduzir erros de arredondamento em alguns métodos; a precisão cônica
*pode* piorá-los quando os valores intermediários se afastam de 1. Nenhum
experimento foi desenhado para favorecer um dos formatos.

## 1. Modelo de arredondamento dos dois formatos

Para um formato com `p` bits de fração na faixa de `x`, o arredondamento ao mais
próximo satisfaz `fl(x) = x(1 + ε)`, `|ε| ≤ u(x) = 2^-(p+1)`.

* **Float16**: `p = 10` em toda a faixa normal `[2^-14, 2^16)`; abaixo disso
  (subnormais) a precisão relativa cai linearmente até 0 em `2^-24`; acima de
  `65504` há overflow para `Inf`; abaixo de `2^-25`, underflow para 0.
* **Posit{16,2}**: `p` depende da magnitude, pois o regime consome bits:

| faixa de `|x|`              | bits de fração Posit{16,2} | Float16 |
|-----------------------------|:--:|:--:|
| `[2^-4, 2^4)`               | 11 | 10 |
| `[2^-8, 2^-4)`, `[2^4, 2^8)`| 10 | 10 |
| `[2^-12, 2^-8)`, `[2^8, 2^12)` | 9 | 10 |
| `[2^-16, 2^-12)`, `[2^12, 2^16)` | 8 | 10 → 8 (subnormal abaixo de 2^-14) |
| `2^-20`                     | 7  | 4 (subnormal) |
| `2^-24`                     | 6  | 0 (menor subnormal) |

  O posit nunca vai a 0 nem a infinito (satura em `minpos = 2^-56` e
  `maxpos = 2^56`). A vantagem nominal é de **um bit, só em `[1/16, 16)`**;
  fora de `[2^-8, 2^8)` o Float16 é mais preciso, exceto na região subnormal.

Consequência direta para o experimento I: `1 + 2^-11` é representável em
Posit{16,2} e arredonda para 1 em Float16; `1 + 2^-12` arredonda para 1 nos dois.

## 2. Quire / produto interno fundido

Com o quire, `Σ aᵢbᵢ` (e `c − ab`) é calculado sem erro e arredondado **uma vez**.
O erro de um produto interno de comprimento `n` passa de `≈ n·u·Σ|aᵢbᵢ|`
(modelo padrão, pior caso) para `u·|Σ aᵢbᵢ|`. O ganho só aparece onde há
acumulação; operações isoladas `a/b` continuam com um arredondamento.
**Não** há quire para Float16 nesta comparação: o objetivo é medir o que o
padrão posit oferece; uma variante `Float16 + fma` pode ser acrescentada como
controle, mas o `fma(::Float16…)` do Julia calcula em `Float32` e arredonda de
novo, então não é um fma correto de 16 bits.

## 3. Sistemas quase singulares (Experimento I)

`A_δ = [1 1; 1 1+δ]`, `det = δ`, `κ₂(A_δ) ≈ 4/δ`. A teoria de perturbação dá

    ‖x̂ − x‖/‖x‖ ≲ κ(A) · ( ‖ΔA‖/‖A‖ + ‖Δb‖/‖b‖ ),

então há duas fontes de erro de natureza diferente:

1. **Erro de representação dos dados.** Armazenar `1+δ` já perturba `A` em
   `|ΔA| ≤ u(1)`, e o erro relativo em `δ` é `u(1)/δ`. Quando `δ < u(1)`, o
   formato não distingue `A_δ` de uma matriz singular. Esse erro não depende do
   algoritmo — por isso medimos o erro tanto em relação ao **problema original**
   quanto ao **problema armazenado** `(A_T, b_T)`.
2. **Erro de arredondamento do algoritmo.** Na eliminação, `l₂₁ = 1` e
   `u₂₂ = fl((1+δ)_T − 1)` é exato (lema de Sterbenz); o mesmo vale para
   `b₂ − b₁` quando os dados são representáveis. Assim, para a família
   `δ = 2^-j` com dados representáveis, a eliminação não comete erro algum; o
   quire não tem o que corrigir. Por isso incluímos uma família **genérica**
   `δ = 2^-j·(√5−1)/2`, cujos dados exigem arredondamento.
   Complementarmente, a estrutura **rotacionada** `G(θ)·A_δ`, `G(θ)·b`
   (32 ângulos θ ∈ (0, π/2) com semente fixa) preserva κ, det e solução, mas
   torna os dados genéricos: a eliminação passa a ter cancelamento e erro de
   arredondamento próprio, que é o que de fato distingue as aritméticas.
   Resultados são resumidos por mediana/quartis sobre os ângulos e por
   comparação **pareada** (mesmo problema, mesmo ângulo, dois formatos).
3. **Pivoteamento parcial** não muda nada aqui (`|a₂₁| = |a₁₁|`, empate ⇒ sem
   troca). Mantemos a variante para checar essa previsão e porque ela é a forma
   usada na prática.
4. **Refinamento iterativo em precisão mista**: fatoração e correção em `T`,
   resíduo `r = b − A x̂` em BigFloat. Converge se `κ(A)·u_T ≲ 1`
   (aprox. `δ ≳ 4u_T`) para o problema armazenado; contra o problema original,
   converge para a solução de `(A_T, b_T)` se o refinamento usar os dados
   armazenados e pode convergir para a solução verdadeira se usar os originais
   (o resíduo "vê" os dados exatos). O resíduo de um `x̂` bom é minúsculo:
   arredondá-lo diretamente para `T` sofre underflow no Float16 (vai a 0 ou a
   subnormal) e cai na região de poucos bits do posit. Por isso há uma variante
   **com escala** (`r/2^k` com `‖r/2^k‖∞ ∈ [1,2)`), que isola o efeito de faixa
   dinâmica do efeito de precisão.

## 4. Gradiente conjugado (Experimento II) e Newton (Experimento III)

Planejados após a validação do Experimento I. Fenômenos a observar:
* CG: perda de ortogonalidade/conjugação (atraso de convergência em aritmética
  finita, Greenbaum/Strakoš), divergência entre resíduo recursivo e verdadeiro,
  `pᵀAp ≤ 0` por arredondamento (perda de definitude), magnitudes de `rᵀr`
  descendo para fora de `[2^-8, 2^8)` — onde o posit perde bits.
* Newton para `(x−1)² − δ`: cancelamento catastrófico em `x − 1` quando a raiz
  está perto de 1; a reformulação `y = x − 1` move as grandezas para perto de
  `√δ`, longe de 1 — teste direto da precisão cônica.

## 5. Cuidados para uma comparação justa

1. **Construção em alta precisão.** Todo dado é criado em BigFloat (256 bits) e
   convertido para `T` com **um** arredondamento. Converter via `Float32` ou
   `Float64` causaria arredondamento duplo (verificado em teste).
2. **Dois experimentos de dados.** (a) *original*: cada formato recebe o
   problema arredondado por ele próprio; (b) *comum*: os dados são substituídos
   por valores exatamente representáveis em **todos** os formatos, de modo que
   a diferença entre os formatos vem só da aritmética.
   Para δ pequeno o arredondamento comum pode descaracterizar o problema
   (o det dos dados comuns se afasta de δ); casos com `|det(A_c) − δ| > δ/2`
   ficam no CSV (`dados_comuns_ok = false`) mas saem dos gráficos e resumos.
3. **Dois erros.** Erro em relação ao problema original e em relação ao
   problema armazenado; erro retroativo normwise
   `η = ‖b − Ax̂‖/(‖A‖‖x̂‖ + ‖b‖)`, calculado em BigFloat.
4. **Sem atalhos de precisão.** Dentro dos algoritmos em `T` só há operações de
   `T` (laços genéricos, sem BLAS, sem `Float64` escondido). As variantes de
   precisão mista (resíduo em BigFloat) são nomeadas explicitamente.
5. **Precisão e tempo separados.** Diagnósticos (BigFloat, rastreamento de
   magnitudes) rodam numa execução à parte; a cronometragem usa BenchmarkTools,
   com aquecimento, sobre a mesma função sem diagnóstico.
6. **Falhas são resultados.** Pivô nulo, valor não finito/NaR, estagnação e
   limite de iterações são registrados com um status, nunca descartados.
   A PositPR lança `DivideError` em `x/0` (em vez de NaR); os algoritmos testam
   o pivô antes de dividir.
7. **Reprodutibilidade.** Sementes fixas, versões no `Manifest.toml`, CSV com
   todos os parâmetros por linha.
8. **Tempo não é custo de hardware.** Float16 roda em hardware (ou emulação via
   Float32); a PositPR é emulação em software. Tempos comparam implementações,
   não formatos.
