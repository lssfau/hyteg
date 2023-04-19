# Scope

This directory contains code and job scripts for the tests and benchmarks presented in Daniel Bauer's Master thesis [*Multigrid in H(curl) on Hybrid Tetrahedral Grids*](https://www10.cs.fau.de/publications/theses/2023/Master_BauerDaniel.pdf).

# Tests and Benchmarks

Experiments from the thesis are implemented in the followings files:

| Section in thesis                                  | File                                                   |
|----------------------------------------------------|--------------------------------------------------------|
| 3.2.1 Numerical Experiment with Chebyshev Smoother | [`smoothingFactor.cpp`](smoothingFactor.cpp)           |
| 3.2.4 Numerical Experiment with Hybrid Smoother    | [`smoothingFactor.cpp`](smoothingFactor.cpp)           |
| 4.2.2 Projection                                   | [`projectionQuadrature.cpp`](projectionQuadrature.cpp) |
| 5.2 Eigenvalue Bounds for the Chebyshev Smoother   | [`chebyshevParameters.cpp`](chebyshevParameters.cpp)   |
| 5.3 Effectiveness of the Hybrid Smoother           | [`smoothingFactor.cpp`](smoothingFactor.cpp)           |
| 5.4 Grid-Independent Convergence                   | [`hIndependence.cpp`](hIndependence.cpp)               |
| 5.5 L2 Convergence Rate                            | [`L2Convergence.cpp`](L2Convergence.cpp)               |
| 5.6 Robustness to Choice of Coefficients           | [`solverConvergence.cpp`](solverConvergence.cpp)       |
| 5.7 Non-Convex Domain                              | [`solverConvergence.cpp`](solverConvergence.cpp)       |
| 5.8 Performance                                    | [`performance.cpp`](performance.cpp)                   |

Compilation and execution is documented in the jobscripts which can be found in the directory `jobs`.

All apps log progress to stdout and write results to the `output` directory.
`.dat` files contain white-space separated tabular material, while `.tex` files store key-value pairs.
