# TinyHHG C++

TinyHHG (HHG stands for hierarchical hybrid grids) is a C++ framework for large scale high performance finite element simulations based on (but not limited to) geometric multigrid.


### Documentation

Access the [main documentation](https://i10git.cs.fau.de/terraneo/hhg_doku "Main Documentation") (work in progress) or the [doxygen documentation](http://terraneo.pages.walberla.net/tinyhhg/index.html "TinyHHG Doxygen").


### Dependencies

The framework is built on top of the core of the [waLBerla framework](http://walberla.net "waLBerla homepage") and therefore requires the its source code to be cloned.

Required:

* the waLBerla source code ([waLBerla GitLab repository](https://i10git.cs.fau.de/walberla/walberla "waLBerla repository"))
* [Boost](https://www.boost.org "boost homepage")
* [CMake](https://cmake.org/ "CMake homepage")
* a C++14 compliant compiler (e.g. gcc, clang, Intel or MSVC)

Optional:

* MPI (e.g. [OpenMPI](https://www.open-mpi.org/ "OpenMPI homepage")) for parallel runs
* [PETSc](https://www.mcs.anl.gov/petsc/ "PETSc homepage") for efficient coarse grid solvers
* [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview "ParMETIS homepage") for high-quality load balancing
* VTK visualization software (e.g. [ParaView](https://www.paraview.org/ "ParaView homepage"))


### Build instructions

To build TinyHHG, clone the TinyHHG and the waLBerla source code:

    $ git clone --recurse-submodules https://i10git.cs.fau.de/terraneo/tinyhhg

`--recurse-submodules` will automatically initialize and clone walberla as a submodule.


    $ mkdir tinyhhg-build 
    $ cd tinyhhg-build
    $ cmake ../tinyhhg

CMake will then produce Makefiles for the included tests and applications. To build and run an application (e.g. a stabilized stokes solver) invoke:

    tinyhhg-build $ cd apps
    tinyhhg-build/apps $ make stokes_stab
    tinyhhg-build/apps $ ./stokes_stab

... or for a parallel run:

    tinyhhg-build/apps $ mpirun -np 4 ./stokes_stab


