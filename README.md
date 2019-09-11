# HyTeG

HyTeG (Hybrid Tetrahedral Grids) is a C++ framework for large scale high performance finite element simulations based on (but not limited to) geometric multigrid.


### Build instructions

To build HyTeG, clone the HyTeG and the waLBerla source code:

    $ git clone --recurse-submodules https://i10git.cs.fau.de/hyteg/hyteg.git

`--recurse-submodules` will automatically initialize and clone walberla as a submodule.


    $ mkdir hyteg-build 
    $ cd hyteg-build
    $ cmake ../hyteg

CMake will then produce Makefiles for the included tests and applications. To build and run an application (e.g. a multigrid benchmark setting) invoke:

    hyteg-build $ cd apps/MultigridStudies
    hyteg-build/apps/MultigridStudies $ make
    hyteg-build/apps/MultigridStudies $ ./MultigridStudies

... or for a parallel run:

    hyteg-build/apps/MultigridStudies $ mpirun -np 4 ./MultigridStudies


### Documentation

The [Doxygen documentation](http://hyteg.pages.walberla.net/hyteg/index.html "HyTeG Doxygen") provides some basic tutorials for example applications.

If you are interested in more background information you may either have
a look at 

* our article [The HyTeG finite-element software framework for scalable multigrid solvers](https://www.tandfonline.com/doi/abs/10.1080/17445760.2018.1506453) - please cite this if you use the software

  ```
  @article{doi:10.1080/17445760.2018.1506453,
  author = {Nils Kohl and Dominik Thönnes and Daniel Drzisga and Dominik 
  Bartuschat and Ulrich Rüde},
  title = {The {HyTeG} finite-element software framework for scalable multigrid solvers},
  journal = {International Journal of Parallel, Emergent and Distributed Systems},
  volume = {34},
  number = {5},
  pages = {477-496},
  year  = {2019},
  publisher = {Taylor & Francis},
  doi = {10.1080/17445760.2018.1506453}
  }
  ```
  
* the [TerraNeo web page](http://terraneo.fau.de) providing information and publications regarding the related research project


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


### Notes

#### CCache

Due to the large amount of generated files it is advisable to activate ccache.
To do so use the CMake setting
    
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache

See also [this StackOverflow answer](https://stackoverflow.com/a/37828605).