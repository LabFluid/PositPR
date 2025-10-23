# PositPR.jl

A pure **Julia implementation** of the [Posit number system](https://posithub.org/docs/Posits4.pdf), providing native encoding, decoding, and arithmetic entirely performed in **Posit representation**, without relying on IEEE floating-point operations at any stage.  
Includes full type interoperability for scientific computing and machine learning applications.


---

## Overview

**PositPR.jl** implements the *Posit* numerical format, a tapered-precision alternative to IEEE 754 proposed by John L. Gustafson.  
It supports arbitrary configurations `(ps, es)` for the total bit size and exponent field, and provides a full set of arithmetic, comparison, and conversion operations.

This library is fully written in Julia and designed for:
- **Numerical analysis** and reproducible low-precision experiments  
- **Scientific computing** and custom arithmetic exploration  
- [**Integration with deep learning frameworks**](https://github.com/uzzui/PositPR-PINNs) like Lux.jl and NeuralPDE.jl  

---

## Installation

Once published:
```julia
] add PositPR
```
For development or testing:
```julia
] add https://github.com/LabFluid/PositPR
```

Then:
```julia
using PositPR
```
---

## Features

| Component | Description |
|------------|-------------|
| **Core type** | `Posit{ps, es}` with dynamic regime and exponent fields |
| **Encoding / Decoding** | Bitstring representation and reversible conversion |
| **Arithmetic** | `+`, `-`, `*`, `/`, `^`, `√` |
| **Comparisons** | Full set of relational and equality operators |
| **Conversions** | To and from `Float64`, `Float32`, and `Int` |
| **Mixed-type support** | Automatic promotion with native Julia types |
| **Mathematical functions** | Experimental: `exp`, `log`, `sqrt`, trigonometric |
| **Random generation** | Uniformly distributed Posit values |
| **Testing suite** | Comprehensive tests for all modules and edge cases |

---

## Usage

Creating Posit numbers and arithmetic operations:
```julia
T = Posit{8,2}
ONE = one(T)
ZERO = zero(T)
x = T(1.5)
y = T(0.5)

x + y     # addition
x * y     # multiplication
x / y     # division
x ^ y     # power
```

You can also create a Posit number by decoding a bit pattern:
```julia
BP = "01000000"
decoded = PositDecoder(BP, T)
```

Conversions between types:
```julia
Float64(x)
Int(x)
Posit{8,1}(x)
```

Generating random Posit numbers:
```julia
rand(T, 5)
# → Vector of random Posit values
```

## Contributing

If you find a **bug**, have a **feature suggestion**, or notice something that could be improved,  
please [open an issue](https://github.com/LabFluid/PositPR/issues) or submit a pull request.  
Contributions, tests, and feedback are very welcome!

---

## License

MIT © 2025 [Mateus Rebellato Ussui](http://lattes.cnpq.br/0621104339270060)  
Universidade Federal do Paraná (UFPR)

---

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

---


