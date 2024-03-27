# HyTeG - Overview

![](doc/logos/HYTEG_large.png)

![CI master](https://i10git.cs.fau.de/hyteg/hyteg/badges/master/pipeline.svg?&key_text=master&key_width=60)
![CI nightly](https://i10git.cs.fau.de/hyteg/hyteg/badges/nightly/pipeline.svg?&key_text=nightly&key_width=60)

## About

HyTeG (Hybrid Tetrahedral Grids) is a C++ framework for extreme-scale matrix-free finite element simulations strong 
focus on geometric multigrid.

For detailed information and references 
[have a look at the documentation](https://hyteg.pages.i10git.cs.fau.de/hyteg/index.html "HyTeG Docs").

[TOC]

## Getting started

### Quickstart

To build HyTeG, clone the [GitLab repository](https://i10git.cs.fau.de/hyteg/hyteg) via:

    $ git clone --recurse-submodules https://i10git.cs.fau.de/hyteg/hyteg.git

The option `--recurse-submodules` is **required** and will automatically initialize and clone 
[waLBerla](http://walberla.net "waLBerla homepage") and [Eigen](http://eigen.tuxfamily.org) as git submodules.

Create a build directory and invoke `cmake`:

    $ mkdir hyteg-build 
    $ cd hyteg-build
    $ cmake ../hyteg

CMake will then produce Makefiles for the included tests and applications. To build and run an application (e.g. a tutorial on isoviscous convection) invoke:

    hyteg-build $ cd tutorials/07_IsoviscousConvectionAnnulus
    hyteg-build/tutorials/07_IsoviscousConvectionAnnulus $ make
    hyteg-build/tutorials/07_IsoviscousConvectionAnnulus $ ./IsoviscousConvectionAnnulus

... or for a parallel run:

    hyteg-build/tutorials/07_IsoviscousConvectionAnnulus $ mpirun -np 4 ./IsoviscousConvectionAnnulus

### Prerequisites

Required:

* a C++17 compliant compiler (e.g. gcc, clang, Intel or MSVC)
* [CMake](https://cmake.org/ "CMake homepage") ( version >= 3.20 )

Automatically cloned via git submodules (**NO need to install/download/clone these manually**):

* [waLBerla](http://walberla.net "waLBerla homepage") for core functionalities (MPI communication, IO, logging, etc.)
* [Eigen](http://eigen.tuxfamily.org "Eigen homepage") for some linear algebra operations
* [HyTeG Operators](https://i10git.cs.fau.de/hyteg/hyteg-operators "hyteg-operators GitLab") for fast generated compute kernels
* [doxygen-awesome-css](https://jothepro.github.io/doxygen-awesome-css "doxygen-awesome-css") for awesome documentation visuals
* [mpreal](https://github.com/advanpix/mpreal "mpreal") as a wrapper for [MPFR](https://www.mpfr.org/ "MPFR homepage")

Optional:

* MPI (e.g. [OpenMPI](https://www.open-mpi.org/ "OpenMPI homepage")) for parallel runs
* [ADIOS2](https://csmd.ornl.gov/software/adios2 "ADIOS2 homepage")  for efficient parallel I/O (CMake option `-DHYTEG_BUILD_WITH_ADIOS2=yes`)
* [PETSc](https://www.mcs.anl.gov/petsc/ "PETSc homepage")  for efficient coarse grid solvers (CMake option `-DHYTEG_BUILD_WITH_PETSC=yes`)
* [Trilinos](https://trilinos.github.io/ "Trilinos homepage") for efficient coarse grid solvers (CMake option `-DHYTEG_BUILD_WITH_TRILINOS=yes`)
* [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview "ParMETIS homepage") for high-quality load balancing
* [MPFR](https://www.mpfr.org/ "MPFR homepage") for emulated floating point datatypes (CMake option `-DHYTEG_BUILD_WITH_MPFR=yes`)
* [Doxygen](https://www.doxygen.nl/ "Doxygen homepage") for building the documentation locally (version >= 1.10.0)

### Configuration options

The builds are configured with standard CMake commands and arguments (starting with `CMAKE_<...>`) and there are 
several additional configuration options from HyTeG (starting with `HYTEG_<...>`) and inherited from waLBerla
(starting with `WALBERLA_<...>`).

To pass options via the commandline, prepend `-D` to the options, like, e.g., `cmake -DHYTEG_BUILD_WITH_PETSC=yes`.

The most relevant options are listed below, with defaults in parentheses:

* `CMAKE_BUILD_TYPE` (`Release`)

  Standard CMake build types. Make sure to build in release mode when running applications and in debug for debugging
  (many asserts help debugging that are disabled in relase mode).

* `HYTEG_BUILD_WITH_PETSC` (`no`)

  Attempts to find a PETSc installation and enables fast sparse solvers (mostly relevant for multigrid coarse grid 
  solvers) and for debugging.

* `HYTEG_BUILD_WITH_TRILINOS` (`no`)

  Same as for PETSc but with less support.

* `HYTEG_BUILD_WITH_ADIOS2` (`no`)

  Finds and enables the ADIOS2 library if installed for efficient parallel I/O.

* `HYTEG_TERRANEO_MODULE` (`no`)

  Builds the shipped TerraNeo module for large-scale Geodynamics simulations. Details below.

* `HYTEG_DOWNLOAD_BOOST` (`no`)

  Downloads the C++ boost (header-only) library which is required for the TerraNeo module for instance and required if 
  not found automatically.

* `WALBERLA_OPTIMIZE_FOR_LOCALHOST` (`yes`)

  Optimizes the build for the present architecture (for instance to enable vector intrinsics if the instruction set is available).

### TerraNeo module

TerraNeo is a module of HyTeG that is providing functionality for running mantle convection models from Geodynamics. 
As this is a specialised application, the module is not built by default. 
In order to compile the corresponding sources, tests and apps (re)run CMake with the following option

    -DHYTEG_TERRANEO_MODULE=yes

Please note that the module depends on [Boost](https://www.boost.org/) library, which is a header-only library.

CMake will search for installed Boost libraries. Should these not be found, you can tell it to download them also by setting

    -DHYTEG_DOWNLOAD_BOOST=yes

## Documentation

Our [documentation page](https://hyteg.pages.i10git.cs.fau.de/hyteg/index.html "HyTeG Docs") 
provides additional documentation beyond this README, such as tutorials, API reference, etc.

  _**Consult the tutorial programs under `tutorials/` and the [generated documentation](Tutorials.html) to get started with
the software.**_

See also the publications below for more in-depth discussion of applications, scalability, and results computed with HyTeG.

For an overview of the TerraNeo project, refer to [the corresponding site](http://terraneo.fau.de).

## Publications

If you are interested in more background information or are looking for a way to cite us, a list of articles is found 
below. _If in doubt, cite the first article._

Overview:

* Kohl, N., Thönnes, D., Drzisga, D., Bartuschat, D., & Rüde, U. (2019).
  _The HyTeG finite-element software framework for scalable multigrid solvers_.
  International Journal of Parallel, Emergent and Distributed Systems.
  [10.1080/17445760.2018.1506453](https://doi.org/10.1080/17445760.2018.1506453)

Finite element data structures:

* Kohl, N., Bauer, D., Böhm, F., & Rüde, U. (2024). 
  _Fundamental data structures for matrix-free finite elements on hybrid tetrahedral grids_. 
  International Journal of Parallel, Emergent and Distributed Systems.
  [10.1080/17445760.2023.2266875](https://doi.org/10.1080/17445760.2023.2266875)

Multigrid for Stokes

* Kohl, N., & Rüde, U. (2022). 
  _Textbook efficiency: massively parallel matrix-free multigrid for the Stokes system_. 
  SIAM Journal on Scientific Computing.
  [10.1137/20M1376005](https://doi.org/10.1137/20M1376005)

Eulerian-Lagrangian methods

* Kohl, N., Mohr, M., Eibl, S., & Rüde, U. (2022). 
  A Massively Parallel Eulerian-Lagrangian Method for Advection-Dominated Transport in Viscous Fluids. 
  SIAM Journal on Scientific Computing.
  [10.1137/21M1402510](https://doi.org/10.1137/21M1402510)

Performance engineering

* Thönnes, D., & Rüde, U. (2023). 
  Model-Based Performance Analysis of the HyTeG Finite Element Framework. 
  In Proceedings of the Platform for Advanced Scientific Computing Conference.
  [10.1145/3592979.3593422](https://doi.org/10.1145/3592979.3593422)

TerraNeo

* Bauer, S., Bunge, H. P., Drzisga, D., Ghelichkhan, S., Huber, M., Kohl, N., Mohr, M., Rüde, U., Thönnes, D., & Wohlmuth, B. I. (2020).
  TerraNeo — Mantle Convection Beyond a Trillion Degrees of Freedom. 
  In Software for Exascale Computing-SPPEXA 2016-2019. 
  Springer International Publishing.
  [10.1007/978-3-030-47956-5_19](https://doi.org/10.1007/978-3-030-47956-5_19)

Surrogates

* Drzisga, D., Wagner, A., & Wohlmuth, B. (2023). 
  A matrix-free ILU realization based on surrogates. 
  SIAM Journal on Scientific Computing.
  [10.1137/22M1529415](https://doi.org/10.1137/22M1529415)

## Contributing

To contribute, you need an account for this GitLab instance. Please contact Dominik Thönnes `dominik.thoennes@fau.de`
for details. 

### Code Style

To keep our code consistently formatted, we use [ClangFormat](https://clang.llvm.org/docs/ClangFormat.html).
There is a `.clang-format` file located in the root directory where all the formatting rules are documented.
The rules can be automatically applied by using:

    clang-format -i $FileName

### Merge requests

A merge request (MR) can be in three different states:

1. Draft: Work in progress, not ready for review or merging.
   When opening a new MR, mark it as draft (by prefixing the title with "Draft: " or using the UI) unless it is already ready for review.
   Please formulate a descriptive title and provide a summary of your changes in the description.
2. Ready for review: Finished from the author's point of view.
   Remove the "Draft" flag to mark the MR as ready and ask a maintainer for review.
   *Every MR must have at least one approving review before it can be merged.*
3. Ready for merging: All work on this MR is done, it has been approved, all discussions are resolved and the CI passed.
   Thanks for your contribution to HyTeG!
   After making sure that "Delete source branch" is checked, you may press that Merge button.
   Please be aware that the CI checks the HEAD of your branch, not the result of the merge.
   In case the master and your branch diverged, it might be a good idea to rebase on the current master or merge the master into your branch so that the CI can check the combination of both changes.

### CCache

Due to the large amount of generated files it can be advisable to activate ccache.
To do so use the CMake setting
    
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache

See also [this StackOverflow answer](https://stackoverflow.com/a/37828605).

## Contact

Nils Kohl `nils.kohl@lmu.de` or Dominik Thönnes `dominik.thoennes@fau.de`.

## License

This project is licensed under the [GNU General Public License v3.0 or later](https://www.gnu.org/licenses/gpl-3.0.html).

## Acknowledgements

Artwork by Manuel Weimann.
