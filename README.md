
This repository contains the scripts accompanying the article

[KM] **Decoding McEliece with a Hint -- Secret Goppa Key Parts Reveal Everything**
by _Elena Kirshanova and Alexander May_
https://eprint.iacr.org/2021/999

# Contributers

* Elena Kirshanova
* Alexander May

# Requirements

* [SageMath 9.3+](https://www.sagemath.org/)


# Description of files
Short description of the content:
* faulty_alphas.sage runs the experiments defined on p.16 to re-create Table 4 from [KM]
* run_advanced_Goppa.sage implements Algorithm 3.2 from [KM]
* run_Goppa_complete.sage implements Algorithm 3.4 from [KM]
* run_Goppa_polynomial_invertible.sage implements Algorithm 3.3 from [KM]
* run_Goppa_polynomial.sage implements Algorithm 3.3 without the assumption on the rank of H^{pub}[I]. The algorithm succeeds with probability ~0.57
* utils.sage (helper file)


# Experiments

To recreate Table 1 from [KM] run
```
sage run_Goppa_complete.sage
```
Concrete McEliece parameter set must be changed in the script.


To recreate Table 2 from [KM] run
```
sage run_advanced_Goppa.sage
```
Concrete McEliece parameter set and the parameter ell must be changed in the script.


To recreate Table 3 from [KM] run
```
sage run_advanced_polynomial_invertible.sage
```
Concrete McEliece parameter set must be changed in the script.

To recreate Table 4 from [KM] run
```
sage faulty_alphas.sage
```
Concrete McEliece parameter set and the number of faults must be changed in the script.
