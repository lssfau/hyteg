# Composites of optimized operators for [HyTeG](https://i10git.cs.fau.de/hyteg/hyteg)

This directory contains composite operators, i.e., operators that are puzzled together using other operators.

Such operators typically correspond to block-matrices, e.g., required for saddle point problems that arise from the
discretization the Stokes-like equations.

Specifically, **all** operators in this directory solely use the generated operators from the `hyteg_operators` 
submodule in `../hyteg_operators/`.

All composites herein are in the namespace `operatorgeneration`.

Thus, the only dependencies of this library are `hyteg` and `hyteg_operators`.
