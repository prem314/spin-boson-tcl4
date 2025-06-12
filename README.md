# tcl4_for_spin_boson_model
Package to generate TCL4 dynamics for generic spin boson model, for odd spectral density, based in analytical calculations.

# How to cite this work

Please cite the following works:

<pre> ```
@article{PhysRevB.111.115423,
  title = {Equivalence between the second order steady state for the spin-boson model and its quantum mean force Gibbs state},
  author = {Kumar, Prem and Athulya, K. P. and Ghosh, Sibasish},
  journal = {Phys. Rev. B},
  volume = {111},
  issue = {11},
  pages = {115423},
  numpages = {13},
  year = {2025},
  month = {Mar},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.111.115423},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.111.115423}
}
``` </pre>

To appear on arxiv soon as preprint:

<pre> ```
@article{PremAnalytic2025,
  title = {Analytic derivation and benchmarking of asymptotic TCL4 generator for general Spin-Boson Model with odd spectral density},
  author = {Kumar, Prem and Athulya, K. P. and Ghosh, Sibasish},
  journal = {arXiv preprint},
  year = {2025}
}
``` </pre>

# How to use the Supplementary material

All the main results of this work are available in `TCL4DynamicsSpinBosonResults.nb`.

## TCL General Derivation

* `TCLIntegrandCalcs.wl` starts from the expressions for the TCL2 and TCL4 generators (see Breuer & Petruccione, 2001) and derives their matrix elements in Bloch-vector form.  
* `TCL2GeneratorCalc.wl` evaluates the general TCL2 generator using those results.  
* `TCL4TripleTimeIntegral.wl` symbolically performs the triple–time integrals of the TCL4 generator matrix elements, and `TCL4OmegaIntegral.wl` then evaluates the associated frequency integrals in the long-time limit.  
* `DQDTCL4GeneratorEvaluation.wl` reproduces the DQD-dynamics figure and can also be used for any odd spectral density (see the numerical-integral caveat).  
* `O2SSCalcSpinBoson.wl` analytically evaluates the second- and fourth-order steady-state results.

## Verification of our results

* `TCL4DrudeGeneratorCalcV2.wl` and `TCL2DrudeGeneratorCalc.wl` carry out the TCL4 and TCL2 calculations for an Ohmic spectral density with a Drude cutoff.  
* `TCL4Verification.wl` cross-checks the general result against these Drude calculations and reproduces the TCL4-verification figure.  
* `TCL4PrecisionVerification.wl` repeats the check at higher numerical precision for \(F_{31}^{(4)}\) and \(F_{41}^{(4)}\).

## Benchmarking with HEOM

* `TCLVsHEOMFidelity.wl` (Mathematica) computes the TCL results and writes them to `python_objects.txt` in the same folder.  
* `tcl_vs_heom.py` (Python) reads `python_objects.txt`, runs the HEOM calculation with **QuTiP**, and compares the two approaches, reproducing the fidelity-comparison figure.

## Non-Markovianity analysis

* `SpinBosonTraceDistancePlotCode.wl` performs the non-Markovianity calculations and reproduces the trace-distance figure.

