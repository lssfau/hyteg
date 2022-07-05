# ![HyTeG](doc/logos/HYTEG_large.png)

HyTeG (Hybrid Tetrahedral Grids) is a C++ framework for large scale high performance finite element simulations based on
(but not limited to) geometric multigrid.


## Build instructions

To build HyTeG, clone via:

    $ git clone --recurse-submodules https://i10git.cs.fau.de/hyteg/hyteg.git

The option `--recurse-submodules` is **required** and will automatically initialize and clone 
[waLBerla](http://walberla.net "waLBerla homepage") as a submodule.

Create a build directory and invoke cmake:

    $ mkdir hyteg-build 
    $ cd hyteg-build
    $ cmake ../hyteg

CMake will then produce Makefiles for the included tests and applications. To build and run an application (e.g. a 
multigrid benchmark setting) invoke:

    hyteg-build $ cd apps/MultigridStudies
    hyteg-build/apps/MultigridStudies $ make
    hyteg-build/apps/MultigridStudies $ ./MultigridStudies

... or for a parallel run:

    hyteg-build/apps/MultigridStudies $ mpirun -np 4 ./MultigridStudies

## Modules

### TerraNeo

TerraNeo is a module of HyTeG that is providing functionality for running mantle convection models from Geodynamics. As this is a specialised application, the module is not build by default. In order to compile the corresponding sources, tests and apps (re)run CMake with the following option

    -DHYTEG_TERRANEO_MODULE=yes

Please note that the module depends on [CGAL (The Computational Geometry Algorithms Library)](https://www.cgal.org/), which in turn has some dependencies on the [Boost](https://www.boost.org/) libraries. Both are header-only libraries.

CMake will download CGAL once as part of the configuration process. If you already have CGAL installed on your system, you can alternatively provide the path to that installation using

    -DCGAL_INCLUDE_PATH=<path to your local installation>

CMake will search for installed Boost libraries. Should these not be found, you can tell it to download them also by setting

    -DHYTEG_DOWNLOAD_BOOST=yes

## Documentation

The [Doxygen documentation](https://hyteg.pages.i10git.cs.fau.de/hyteg/index.html "HyTeG Doxygen") provides some basic 
tutorials for example applications.

If you are interested in more background information you may either have
a look at 

* our article [The HyTeG finite-element software framework for scalable multigrid solvers](https://www.tandfonline.com/doi/abs/10.1080/17445760.2018.1506453) - please cite this if you use the software


  ```
  @article{doi:10.1080/17445760.2018.1506453,
    author = {Nils Kohl and Dominik Thönnes and Daniel Drzisga and Dominik Bartuschat and Ulrich Rüde},
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
  
* the [TerraNeo web page](http://terraneo.fau.de) providing information and publications regarding the related research 
  project
  
* our article [TerraNeo—Mantle Convection Beyond a Trillion Degrees of Freedom](https://doi.org/10.1007/978-3-030-47956-5_19)
  summarizing recent achievements during the TerraNeo project


## Dependencies

The framework is built on top of the core of the [waLBerla framework](http://walberla.net "waLBerla homepage").
Its repository is included via git submodule. So just clone with 

    $ git clone --recurse-submodules https://i10git.cs.fau.de/hyteg/hyteg.git

as written above, to set up waLBerla automatically.

Required:

* [CMake](https://cmake.org/ "CMake homepage")
* a C++17 compliant compiler (e.g. gcc, clang, Intel or MSVC)

Optional:

* [Eigen](http://eigen.tuxfamily.org "Eigen homepage") for some linear algebra operations
  
  Eigen is, (like waLBerla) automatically cloned as a git submodule.
  
  CMake will automatically find the Eigen submodule, there is no need to specify a path
  or to download Eigen at all.

* MPI (e.g. [OpenMPI](https://www.open-mpi.org/ "OpenMPI homepage")) for parallel runs
* [PETSc](https://www.mcs.anl.gov/petsc/ "PETSc homepage") and/or [Trilinos](https://trilinos.github.io/ "Trilinos homepage") for efficient coarse grid solvers
* [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview "ParMETIS homepage") for high-quality load balancing


## Notes

### Code Style

To keep our code consistently formatted, we use [ClangFormat](https://clang.llvm.org/docs/ClangFormat.html).
There is a `.clang-format` file located in the root directory where all the formatting rules are documented.
The rules can be automatically applied by using:

    clang-format -i $FileName

### CCache

Due to the large amount of generated files it is advisable to activate ccache.
To do so use the CMake setting
    
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache

See also [this StackOverflow answer](https://stackoverflow.com/a/37828605).

## Acknowledgements

Artwork by Manuel Weimann.
