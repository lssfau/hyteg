add_executable       ( FA.01_GeometricMultigrid FA.01_GeometricMultigrid.cpp )
target_link_libraries( FA.01_GeometricMultigrid hyteg walberla::core constant_stencil_operator )

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( ../../../data/meshes/annulus.msh)

if(CMAKE_BUILD_TYPE MATCHES "Release")
   waLBerla_execute_test( NO_MODULE_LABEL NAME FA.01_GeometricMultigrid )
endif()