if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

add_executable       ( smoothingFactor smoothingFactor.cpp common.hpp common.cpp )
target_link_libraries( smoothingFactor hyteg mixed_operator )

add_executable       ( massForm massForm.cpp
            forms/N1E1Form_mass_B_fp32.hpp    forms/N1E1Form_mass_B_fp32.cpp
            forms/N1E1Form_mass_B_fp64.hpp    forms/N1E1Form_mass_B_fp64.cpp
            forms/N1E1Form_mass_BInv_fp32.hpp forms/N1E1Form_mass_BInv_fp32.cpp
            forms/N1E1Form_mass_BInv_fp64.hpp forms/N1E1Form_mass_BInv_fp64.cpp
            forms/N1E1Form_mass_LU_fp32.hpp   forms/N1E1Form_mass_LU_fp32.cpp
            forms/N1E1Form_mass_LU_fp64.hpp   forms/N1E1Form_mass_LU_fp64.cpp )
target_link_libraries( massForm hyteg)

add_executable       ( projectionQuadrature projectionQuadrature.cpp )
target_link_libraries( projectionQuadrature hyteg )

add_executable       ( chebyshevParameters chebyshevParameters.cpp common.hpp common.cpp )
target_link_libraries( chebyshevParameters hyteg mixed_operator )

add_executable       ( hIndependence hIndependence.cpp common.hpp common.cpp )
target_link_libraries( hIndependence hyteg mixed_operator )

add_executable       ( L2Convergence L2Convergence.cpp common.hpp common.cpp )
target_link_libraries( L2Convergence hyteg mixed_operator )

add_executable       ( solverConvergence solverConvergence.cpp common.hpp common.cpp )
target_link_libraries( solverConvergence hyteg mixed_operator )

add_executable       ( performance performance.cpp common.hpp common.cpp )
target_link_libraries( performance hyteg mixed_operator )
