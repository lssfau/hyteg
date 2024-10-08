# Master's thesis of Michael Zikeli ''Mixed Precision Multigrid Methods on Hybrid Tetrahedral Grids''
The thesis can be found [here](https://www10.cs.fau.de/publications/theses/2024/Master_ZikeliMichael.pdf).
The defense presentation slides can be found [here](https://faubox.rrze.uni-erlangen.de/getlink/fiGij3hWL5itYYKJBKv4ZT/2024-05-29-MasterarbeitVerteidigung.pdf).

## Description
The goal of my thesis was the accurate and performant implementation of mixed precision for HyTeG. 
In scientific computing, it is common to use single or double precision for computations, however, 
lowering the precision reduces runtime and energy consumption at the cost of less reliable approximations. 
Mixed precision has been gaining attention by addressing this trade-off.

In the current version as of the 06.08.2024, mixed precision is implemented in HyTeG for the `P1` function space and tested for a Poisson problem.
With the used setup of a 3D cube made of 24 macro tets, a refinement gives an accurate solution for float16 up to level 3 (2.45e3 DoFs), 
for pure float32 up to level 9 (8.48e6 DoFs), and using mixed precision with IR, no inaccurate solution could be obtained on only one Fritz node.

However, applying lower precision did not result in a more accurate solution.
I assume this to be caused by unnecessary casting in the Chebyshev smoother, as well as bad performance for single precision vector-vector operations, as existing kernels are not optimized for precisions other than float64.

---

## Structure
* ''operators-used'' contains all operators that were generated by HOG and used for the Thesis, but are deprecated by now.
* ''scripts'' contains scripts that can be used for debugging and post-processing
* `[*].slurm.sh` are slurm scripts that have been used to submission to Fritz.
* `[*].cpp` are the applications / main files that describe a simulation. To ease the setup of different applications, this format was chosen over config files.
* `[*].h` are library like files that describe the components of a simulation that is finally described in `solvePDEOpGen.hpp`
* `[*].hpp` are standalone header libraries that might be outsourced as standalone modules to HyTeG in the future. (IR solver, PDE solver interface)

### Functionalities
#### `solvePDEOpGen.hpp` This is the main functionality for benchmarking mixed precision PDE solvers. 
This file provides an interface requiring 
* some template parameters:
    - a Problem dependent `SetupType`
    - a `quadrature order` for computing the exact values in L2-space
    - The `name of the solver` used for the benchmark. In the case of GMG or IR, the respective smoother or inner solvers are defined in another header file (_SolverRegister.h_)
    - an Operator that is used for the `smoother`
    - a `fold` of operators used to compute the `residual`
* as well as a `config file` defining the simulation.
    - This can be defined during runtime
    - All possible options and the default values can be found in `Setup.h`

The several parts of the simulation are defined by these following header files:
* `Setup.h`: Provides the definition of two different PDE problems as well as a struct describing the runtime simulation parameters. 
* `SolverRegister.h`: Describing the Hierarchy of the solvers, e.g., IR uses GMG as inner solver, which uses Chebyshev as smoother and CG as coarse grid solver...
* `GridtransformationRegister.h`: Definition and parametrization of the grid transformation functions used, i.e., coarsening, prolongation.
* `globalHeader.h`: Some arbitrary parameters that were outsourced to ease the implementation and reduce linker dependencies.

#### Newly generated solvers
* `IterativeRefinementSolver.hpp`: The idea with this solver is to precompute a residual in a higher precision (`$\bar{\text{precision}$`) and then solving the residual equation `$Ae=r$` in low precision (`$\dot{\text{precision}$`) instead of the actual problem `$Au=b$`. This can be made a standalone solver and will eventually be replaced by a mixed precision V-cycle, since there the residual is computed as well.
* `IR_VersionWithGetterFct.hpp`: This was a version where I tried to improve the solver. This will probably be discarded at some point.

#### There exist some arbitrary functionality headers that are used to gather the simulation data.
* `Tree.h`: Some `ChatGBT` generated Tree class used to gather residual and error values per simulation step and level for easier investigation of the convergence behavior.


### Applications
#### Benchmark versions used in my MT
For my MT I tested ? different benchmark options:
* `Runtime`: used to see if the performance can be increased using MP.
* `Hockey-Stick`: used to find which precision can be used where until which resolution.
  * Here two versions exist:
  - One where the simulation terminates once the residual or error converges `...-residual-convergence...`
  - One where the simulation terminates when the error falls below a given threshold `...-error-threshold...`
  The error threshold is determined by the discretization error. More information compares thesis.
* `FineTuning`: This was meant to be used to fine-tune the smoother and IR parameters, however, since the final runtime was not enough for fine-tuning, the results from this benchmark were not included in the thesis.

#### Other helpful applications
* `DoFCounter.cpp`: Standalone application that outputs the DoFs per level for the given discretization to determine the required memory size.


### Scripts
There exist execution scripts, i.e., slurm scripts as well as post-processing scripts which are located in the folder scripts.
For every benchmark executed on the `Fritz cluster`, there exists a slurm script.
There also exist a few slurm scripts that include the parameters to run the simulation with `likwid`.
Also, there is a vectorized and non-vectorized version of the runtime benchmark to see if the runtime-problem comes from a non-optimal vectorization.

More interesting are the post-processing scripts:
* `count-underflow-from-vtu.py`: I think this script is deprecated. I'm not sure though. 
* `create-timing-table.py`: This is the actual file that gathers the relevant timing values from the timing tree.
* `create-timing-table-used-for-test-timings.py`: This is a WIP file, that I created from `create-timing-table.py` to investigate some values by loosing others or having too many...
* `createResidualTablesFromJSON.py`: The residuals are gathered in a json format which is not as easy to print. This file converts the json files to a csv format.
* `memoryConsumption.py`: This script is a dictionary of all simulation parameters to return the assumed memory consumption for each configuration.
* `relate-benchmarks.py`: I wanted to know if some configurations result in a more accurate solution compared for others, so this file creates a `.dat` file with errors or residuals containing values relative to the specified file.
* `strip-none.sh`: Some files in `output_pgfplot` contain a "None" line when the simulation earlier than others do to the `createResidualTablesFromJSON.py` script, these lines are striped before printing as plot.

How all the post-processing files were executed is described in: `scriptCallHelper.txt`.

---

## TODO
The mixed precision support in HyTeG can be extended by realizing the following ToDos:
- [ ] Make IR a standalone solver
- [ ] Make a GMG solver that allows for different precision for different parts (with @Dinesh)
- [ ] Try to scale the residual by the max-norm before applying the inner solver in the IR solver to allow for lower precision (Tip from Petr)
- [ ] Model the precision for the newest version of `HOG` and `hyteg_operators` to get rid of the folder `operators-used` and allow for more recent operators.
- [ ] Extend precision support to function spaces other than `P1`.
- [ ] Fix the performance issue of single precision vector-vector operations.
- [ ] Try progressive precision FMG.
