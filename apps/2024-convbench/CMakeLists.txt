if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

# include(CMakePrintHelpers)
# cmake_print_variables(Python3_VERSION)
# cmake_print_variables(Python3_INCLUDE_DIRS)

# include_directories(${Python3_INCLUDE_DIRS})

waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *sh )

waLBerla_add_executable( NAME 2DStokesFSConvTest
        FILES  2DStokesFSConvTest.cpp
        DEPENDS hyteg mixed_operator)

waLBerla_add_executable( NAME 3DStokesFSConvTest
        FILES  3DStokesFSConvTest.cpp
        DEPENDS hyteg mixed_operator)

waLBerla_add_executable( NAME TestDivKGrad
        FILES  TestDivKGrad.cpp
        DEPENDS hyteg mixed_operator)

waLBerla_add_executable( NAME CompressibleAnnulusStatic
        FILES  CompressibleAnnulusStatic.cpp
        DEPENDS hyteg mixed_operator opgen-composites-stokes opgen-terraneo)

waLBerla_add_executable( NAME 2DStokesFSConvRediscTest
        FILES  2DStokesFSConvRediscTest.cpp
        DEPENDS hyteg mixed_operator opgen-composites-stokes opgen-k_mass opgen-terraneo)

waLBerla_add_executable( NAME SphericalShellBenchRotation
        FILES  SphericalShellBenchRotation.cpp
        DEPENDS hyteg terraneo mixed_operator convection_particles opgen-composites-stokes opgen-div_k_grad opgen-advection opgen-supg_diffusion opgen-supg_advection opgen-full_stokes opgen-terraneo opgen-mass opgen-k_mass opgen-epsilon opgen-divergence opgen-gradient)
