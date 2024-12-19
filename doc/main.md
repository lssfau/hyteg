[TOC]
# Hybrid Tetrahedral Grids (HyTeG)

This page serves as an entry point to our doxygen documentation. Before
working with the software this page should be read entirely. We try to
give an short but complete overview of the features from a developer 
point of view. For more theoretical background see our publications
especially [HyTeG](https://doi.org/10.1080/17445760.2018.1506453)

## Install

### Requirements
The only external dependency for our framework is 
[waLBerla](https://www.walberla.net), which is a multi-physics framework
also developed at our group.
We use a lot of core functionality and included it therefore as a submodule.

Additionally we require a compiler that supports C++17.

MPI is required if the application should be executed
in parallel. For example using [Open MPI](https://www.open-mpi.org/).

### Obtain the software
To build HyTeG, clone the HyTeG and the waLBerla source code:

    git clone --recurse-submodules https://i10git.cs.fau.de/hyteg/hyteg.git

`--recurse-submodules` will automatically initialize and clone waLBerla 
as a submodule. Alternatively this can be achieved by:

    cd hyteg
    git submodule init
    git submodule update

Please be aware that submodule are not automatically updated when 
switching branches or checkout out commits. One can however achieve this
behaviour in the git config

    git config --global submodule.recurse true
    
### Installing the software

We use CMake as our build system. These commands inside the hyteg repository 
build our whole framework:

    mkdir hyteg-build 
    cd hyteg-build
    cmake ../hyteg
    make
    
### Customizing the build

CMake build as customize by providing the `cmake` command with addtional
variables using `-D...`. For example to change the compiler to clang one
uses:

    cmake ../hyteg -DCMAKE_CXX_COMPILER=clang++
    
By default CMake uses the `Release` mode which adds `-O3 -DNDEBUG` as 
compiler flags. Other options are e.g. `Debug` or `RelWithDebInfo` which
can be set using e.g.

    cmake ../hyteg -DCMAKE_BUILD_TYPE=Debug
    
A kind of "GUI" to edit the CMake options is available by using

    ccmake ../hyteg
    
Additional options are available via waLBerla.
MPI is enabled by default but can be disabled using:

    -DHYTEG_BUILD_WITH_MPI=OFF
    
    
### Additional components
### TerraNeo

TerraNeo is a module of HyTeG that is providing functionality for running
mantle convection models from Geodynamics. As this is a specialised application,
the module is not build by default. In order to compile the corresponding
sources, tests and apps (re)run CMake with the following option

    -DHYTEG_TERRANEO_MODULE=yes

Please note that the module depends on
[CGAL (The Computational Geometry Algorithms Library)](https://www.cgal.org/),
which in turn has some dependencies on the [Boost](https://www.boost.org/)
libraries. Both are header-only libraries.

CMake will download CGAL once as part of the configuration process. If you
already have CGAL installed on your system, you can alternatively provide the
path to that installation using

    -DCGAL_INCLUDE_PATH=<path to your local installation>

CMake will search for installed Boost libraries. Should these not be found, you
can tell it to download them also by setting

    -DHYTEG_DOWNLOAD_BOOST=yes

#### PETSc

We implement an interface to [PETSc](https://www.mcs.anl.gov/petsc/),
which can be used as a preconditioner or solver for suitable applications.

One has to include the wrapper `PETScWrapper.hpp` and enable PETSc in CMake:

    -DHYTEG_BUILD_WITH_PETSC=ON

#### Eigen

The [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) template
library is also available through `EigenWrapper.hpp`. Eigen is automatically
cloned as a submodule. So you do not need to download it. Instead, simply set

    -DHYTEG_BUILD_WITTH_EIGEN=ON

#### LIKWID

[LIKWID](https://github.com/RRZE-HPC/likwid) is a tool suite to measure
performance counters and can also be used to instrument the code and
create specific regions which can be analyzed instead of analyzing the
whole application.  
These markers are available through LikwidWrapper.hpp and are enabled
using:

    -DHYTEG_BUILD_WITH_LIKWID=ON

## Framework overview

Here an overview of the structure of the HyTeG source code is given. 
For details about the functionality see the corresponding doxygen pages.

### `tutorials`

The best entry point to get to know our framework are the tutorials.
These are also linked in the Doxygen documentation: [Tutorial](pages.html)

### `src`

This folder contains the actual source code of the HyTeG framework.
Each subfolder in this directory is treated as a module by cmake. Meaning
that these are compiled into library files and linked accordingly.

### `apps`

Contains some example applications

### `tests`

All the tests for the HyTeG framework are located here. The structure
resembles the same as `src`. All of these tests are executed by our
automated continuous integration pipeline.
