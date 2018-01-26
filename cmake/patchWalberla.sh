#!/bin/bash
## this script patches walberla such that most of the modules are removed from the all target
## it is not very nice
sed -i 's/set ( excludedModules  )/set ( excludedModules "pe" "blockforest" "boundary" "communication" "cuda" "pde" "domain_decomposition" "fft" "field" "gather" "geometry" "gui" "lbm" "mesh" "pe_coupling" "postprocessing" "python_coupling" "simd" "stencil" "timeloop" "vtk")/' src/CMakeLists.txt
sed -i '/add_subdirectory( tests EXCLUDE_FROM_ALL )/d' CMakeLists.txt
sed -i '/add_subdirectory ( tests )/d' CMakeLists.txt
sed -i '/add_subdirectory ( apps )/d' CMakeLists.txt
