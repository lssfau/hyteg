if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_add_executable(NAME smoothingFactor
      FILES smoothingFactor.cpp common.hpp common.cpp
      DEPENDS hyteg mixed_operator)

waLBerla_add_executable(NAME massForm
      FILES massForm.cpp
            forms/N1E1Form_mass_B_fp32.hpp    forms/N1E1Form_mass_B_fp32.cpp
            forms/N1E1Form_mass_B_fp64.hpp    forms/N1E1Form_mass_B_fp64.cpp
            forms/N1E1Form_mass_BInv_fp32.hpp forms/N1E1Form_mass_BInv_fp32.cpp
            forms/N1E1Form_mass_BInv_fp64.hpp forms/N1E1Form_mass_BInv_fp64.cpp
            forms/N1E1Form_mass_LU_fp32.hpp   forms/N1E1Form_mass_LU_fp32.cpp
            forms/N1E1Form_mass_LU_fp64.hpp   forms/N1E1Form_mass_LU_fp64.cpp
      DEPENDS hyteg)

waLBerla_add_executable(NAME projectionQuadrature
      FILES projectionQuadrature.cpp
      DEPENDS hyteg)

waLBerla_add_executable(NAME chebyshevParameters
      FILES chebyshevParameters.cpp common.hpp common.cpp
      DEPENDS hyteg mixed_operator)

waLBerla_add_executable(NAME hIndependence
      FILES hIndependence.cpp common.hpp common.cpp
      DEPENDS hyteg mixed_operator)

waLBerla_add_executable(NAME L2Convergence
      FILES L2Convergence.cpp common.hpp common.cpp
      DEPENDS hyteg mixed_operator)

waLBerla_add_executable(NAME solverConvergence
      FILES solverConvergence.cpp common.hpp common.cpp
      DEPENDS hyteg mixed_operator)

waLBerla_add_executable(NAME performance
      FILES performance.cpp common.hpp common.cpp
      DEPENDS hyteg mixed_operator)
