# PositPR.jl

A pure-Julia implementation of the [Posit number system](https://posithub.org/docs/Posits4.pdf)
(the Posit Standard, 2022)
for any configuration `Posit{ps,es}` with up to 128 bits, including the
**quire** for exact sums and dot products.

All arithmetic is done with integer operations, with no floating-point
operations along the way. Every operation computes its exact result and is then
rounded **once**, to nearest with ties to even on the bit string, as the
standard requires. The suffix *PR* is a tribute to the Universidade Federal do
Paraná (UFPR), where the library was developed.

It is designed for:
- **numerical analysis** and reproducible low-precision experiments;
- **scientific computing** and exploration of custom arithmetic;
- [**integration with deep learning frameworks**](https://github.com/uzzui/PositPR-PINNs)
  such as Lux.jl and NeuralPDE.jl.

## The posit format

A posit with `ps` bits and `es` exponent bits is a bit string

```
sign | regime (variable length) | exponent (up to es bits) | fraction (remaining bits)
```

and represents `(-1)^s · useed^k · 2^e · (1 + f)` with `useed = 2^(2^es)`.
The regime is a run of identical bits that works like a "super-exponent". The
result is **tapered precision**: numbers close to 1 get the most fraction bits,
while very large or very small numbers get fewer bits and a much wider dynamic
range than an IEEE float of the same size. For example, `Posit{16,2}` has 11
fraction bits in `[1/16, 16)`, one more than `Float16`, but fewer outside
`[2^-8, 2^8)`.

There is a single zero and a single exception value, **NaR** ("not a real").
Posits never overflow to infinity or underflow to zero: results saturate at
`maxpos` and `minpos`.

## Installation

The package is not registered yet. Install it from GitHub:

```julia
] add https://github.com/LabFluid/PositPR
```

or, for development, from a local clone with `] dev path/to/PositPR`.

## Quick start

```julia
using PositPR

P = Posit{16,2}               # family; the concrete type is Posit16{2}
x = P(1.5); y = P(0.1)

x + y, x * y, x / y           # correctly rounded
sqrt(y), exp(y), sin(y)       # elementary functions (within 1 ULP)
Float64(y)                    # 0.100006103515625: the exact stored value
nextPosit(x), prevPosit(x)    # neighbours on the posit ring
isNaR(NaR(P))                 # true

sizeof(x)                     # 2: stored as a UInt16, just like a Float16
PositPR.pi(P)                 # pi rounded to Posit{16,2} (not exported)

A = P[2 1; 1 3]
A * P[1, 1]                   # generic Julia code works on posits
rand(P), randn(P)
```

Aliases select the storage word automatically: `Posit8{es}`, `Posit16{es}`,
`Posit32{es}`, `Posit64{es}`, or `PositType(ps, es)` for any size. Each posit
is stored as its bit pattern in the smallest unsigned integer that holds it
(`UInt8` … `UInt128`), so arrays of posits use the same memory as arrays of
floats of the same width.

### The quire: exact accumulation

```julia
q = Quire(P)
for (a, b) in zip(xs, ys)
    qmuladd!(q, a, b)         # q += a*b, exactly
end
Posit(q)                      # rounded once

fdp(xs, ys)                   # fused dot product
fsum(P.([1e4, 1, -1e4]))      # 1.0 (naive summation gives 0)
fma(a, b, c)                  # a*b + c with one rounding
```

The quire is a fixed-point accumulator whose least significant bit is
`minpos²`, so every product of two posits fits exactly. Here it is backed by a
`BigInt` and never overflows. Wherever the 16·ps-bit quire of the standard does
not overflow, both give the same results.

## Exported API

| area | names |
|---|---|
| types | `Posit`, `Posit8`, `Posit16`, `Posit32`, `Posit64`, `PositType`, `useed` |
| conversions | `x2Posit`, `Posit2x`, `Posit2Posit`, `PositConverter` (also `P(x)`, `Float64(p)`, `BigFloat(p)`, `convert`, `promote_rule`) |
| bit strings | `PositDecoder`, `PositEncoder`, `to_decimal`, `shift_left`, `exponent_update` |
| comparisons | `nextPosit`, `prevPosit`, `isNaR`, `NaR` (plus `==`, `<`, `abs`, `eps`, `typemax`, …) |
| quire | `Quire`, `qadd!`, `qsub!`, `qmuladd!`, `qmulsub!`, `qneg!`, `qabs!`, `qclear!`, `fdp`, `fsum` |

Arithmetic (`+ - * /`, `sqrt`, `^`, `mod`, `fma`), `exp`, `log`, `sin`, `cos`,
`tan`, `tanh`, `rand` and `randn` extend the Base functions.

## Implementation notes

- **Decoding with bit operations.** Regime, exponent and fraction are extracted
  with `leading_zeros`/`leading_ones` on the stored word.
- **One rounding point.** Every operation reduces its exact result to
  `(sign, M, P)`, meaning `±M·2^P` with an integer `M`, and passes it to a single
  function, `posit_round`, which applies round-to-nearest-even on the bit string
  and saturation. Intermediate integers use `UInt128` whenever they fit (always
  for `ps ≤ 32`) and `BigInt` otherwise.
- **Comparisons as integers.** Posits are ordered like two's-complement
  integers, so `<` and `==` are integer comparisons.
- **Elementary functions** are evaluated in a wider posit type
  (`Posit{ps+16,es}`) and rounded once to the target type.

Source files are in `src/posit/`: `type.jl`, `rounding.jl`, `operations.jl`,
`conversions.jl`, `comparisons.jl`, `encoding_decoding.jl`, `trig.jl`,
`random.jl`, `quire.jl` and `auxiliar.jl`.

## Testing

```sh
julia --project test/runtests.jl
```

The suite compares the library with an independent reference implementation
(`test/reference.jl`, `BigFloat` with 1024 bits). It is exhaustive for 8-bit
posits (all 256 patterns and all 65 025 operand pairs, for every `es` from 0 to
3) and samples 16- and 32-bit posits. All checks pass: `+ - * /` are correctly
rounded in 100% of the cases, and so are `fma`, `fdp` and `fsum`. Functions with
loops run under a watchdog process, so an infinite loop shows up as a failure
instead of hanging the suite.

`dev/stress.jl` runs 800 000 extra random operations against the oracle;
`dev/quire_demo.jl` and `dev/quire_bench.jl` illustrate the quire.

## Numerical stability benchmarks

`benchmarks/` holds a reproducible suite that compares `Float16`,
`Posit{16,2}` and `Posit{16,2}` with quire (and their 32-bit counterparts)
against a 256-bit `BigFloat` reference. It separates the error from
representing the data from the error made by the arithmetic, and records the
magnitudes of intermediate values so the effect of tapered precision can be
measured.

```sh
julia --project=benchmarks benchmarks/run_all.jl
```

The documentation is in Portuguese: see `benchmarks/README.md`,
`FUNDAMENTACAO.md` (method) and `RESULTADOS_I.md` (results for near-singular
linear systems).

## Known limitations

- `x / 0` throws `DivideError`; the standard specifies NaR.
- The quire uses `BigInt`, so it allocates and is about 10× slower than plain
  posit arithmetic. A fixed-size quire would remove both costs.
- Performance is that of a software emulation: tens of nanoseconds per
  operation, a few times slower than hardware `Float16`/`Float32`.

## Contributing

If you find a **bug**, have a **feature suggestion**, or notice something that could be improved,
please [open an issue](https://github.com/LabFluid/PositPR/issues) or submit a pull request.
Contributions, tests, and feedback are very welcome!

## References

- BOX, G. E. P.; MULLER, M. E. *A Note on the Generation of Random Normal Deviates.*  
  The Annals of Mathematical Statistics, Institute of Mathematical Statistics, v. 29, n. 2, p. 610–611, 1958.  
  DOI: [10.1214/aoms/1177706645](https://doi.org/10.1214/aoms/1177706645)

- CARMICHAEL, Z.; LANGROUDI, H. F.; KHAZANOV, C.; LILLIE, J.; GUSTAFSON, J. L.; KUDITHIPUDI, D.  
  *Deep Positron: A Deep Neural Network Using the Posit Number System.*  
  In: *2019 Design, Automation & Test in Europe Conference & Exhibition (DATE)*. 2019. p. 1421–1426.  
  DOI: [10.23919/DATE.2019.8715262](https://doi.org/10.23919/DATE.2019.8715262)

- DIXIT, V. K.; RACKAUCKAS, C. *Optimization.jl: A Unified Optimization Package.* Zenodo, Mar. 2023.  
  DOI: [10.5281/zenodo.7738525](https://doi.org/10.5281/zenodo.7738525)

- GUSTAFSON, J. L.; YONEMOTO, I. T. *Beating Floating Point at Its Own Game: Posit Arithmetic.*  
  Supercomputing Frontiers and Innovations, v. 4, n. 2, p. 71–86, 2017.  
  DOI: [10.14529/jsfi170206](https://doi.org/10.14529/jsfi170206)

- IEEE. *IEEE Standard for Floating-Point Arithmetic.* IEEE Std 754-2019 (Revision of IEEE 754-2008), 2019.  
  DOI: [10.1109/IEEESTD.2019.8766229](https://doi.org/10.1109/IEEESTD.2019.8766229)

- INNES, M.; SABA, E.; FISCHER, K.; GANDHI, D.; RUDILOSSO, M. C.; JOY, N. M.; KARMALI, T.; PAL, A.; SHAH, V.  
  *Fashionable Modelling with Flux.* CoRR, abs/1811.01457, 2018.  
  Available at: [https://arxiv.org/abs/1811.01457](https://arxiv.org/abs/1811.01457)

- INNES, M. *Flux: Elegant Machine Learning with Julia.* Journal of Open Source Software, 2018.  
  DOI: [10.21105/joss.00602](https://doi.org/10.21105/joss.00602)

- JOST, T. T. *Compilation and Optimizations for Variable Precision Floating-Point Arithmetic:  
  From Language and Libraries to Code Generation.* Université Grenoble Alpes, 2021.  
  Available at: [https://theses.hal.science/tel-03414534](https://theses.hal.science/tel-03414534)

- LU, J.; FANG, C.; XU, M.; LIN, J.; WANG, Z. *Evaluations on Deep Neural Networks Training Using Posit Number System.*  
  IEEE Transactions on Computers, v. 70, n. 2, p. 174–187, 2021.  
  DOI: [10.1109/TC.2020.2985971](https://doi.org/10.1109/TC.2020.2985971)

- LU, J.; LU, S.; WANG, Z.; FANG, C.; LIN, J.; WANG, Z.; DU, L. *Training Deep Neural Networks Using Posit Number System.*  
  arXiv: [1909.03831](https://arxiv.org/abs/1909.03831) [cs.LG], 2019.

- PAL, A. *Lux: Explicit Parameterization of Deep Neural Networks in Julia.* Zenodo, Apr. 2023.  
  DOI: [10.5281/zenodo.7808904](https://doi.org/10.5281/zenodo.7808904)

- PAL, A. *On Efficient Training & Inference of Neural Differential Equations.*  
  Massachusetts Institute of Technology, 2023.

- RACKAUCKAS, C.; MA, Y.; MARTENSEN, J.; WARNER, C.; ZUBOV, K.; SUPEKAR, R.; SKINNER, D.; RAMADHAN, A.  
  *Universal Differential Equations for Scientific Machine Learning.* arXiv preprint, 2020.  
  arXiv: [2001.04385](https://arxiv.org/abs/2001.04385)

- RAPOSO, G.; TOMÁS, P.; ROMA, N. *PositNN: Training Deep Neural Networks with Mixed Low-Precision Posit.*  
  In: *ICASSP 2021 - IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP).*  
  2021. p. 7908–7912. DOI: [10.1109/ICASSP39728.2021.9413919](https://doi.org/10.1109/ICASSP39728.2021.9413919)

- ZUBOV, K.; MCCARTHY, Z.; MA, Y.; CALISTO, F.; PAGLIARINO, V.; AZEGLIO, S.; BOTTERO, L.; LUJÁN, E.;  
  SULZER, V.; BHARAMBE, A.; VINCHHI, N.; BALAKRISHNAN, K.; UPADHYAY, D.; RACKAUCKAS, C.  
  *NeuralPDE: Automating Physics-Informed Neural Networks (PINNs) with Error Approximations.* arXiv, 2021.  
  DOI: [10.48550/ARXIV.2107.09443](https://doi.org/10.48550/ARXIV.2107.09443)

- POSIT WORKING GROUP. *Standard for Posit Arithmetic (2022).* 2022.

- USSUI, M. R. *Sistema de Ponto Flutuante com Precisão Variável Posit em Redes Neurais
  Informadas pela Física (PINNs).* Undergraduate thesis (Industrial Mathematics), UFPR, 2025.


## License

MIT © 2025 [Mateus Rebellato Ussui](http://lattes.cnpq.br/0621104339270060)  
Universidade Federal do Paraná (UFPR). See [LICENSE](LICENSE).
