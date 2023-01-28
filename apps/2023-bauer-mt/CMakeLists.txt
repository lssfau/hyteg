if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_add_executable(NAME smoothingFactor
      FILES smoothingFactor.cpp common.hpp common.cpp KeyValueStore.hpp Table.hpp
      DEPENDS hyteg)

waLBerla_add_executable(NAME hIndependence
      FILES hIndependence.cpp common.hpp common.cpp KeyValueStore.hpp Table.hpp
      DEPENDS hyteg)

waLBerla_add_executable(NAME L2Convergence
      FILES L2Convergence.cpp common.hpp common.cpp KeyValueStore.hpp Table.hpp
      DEPENDS hyteg)

waLBerla_add_executable(NAME solverConvergence
      FILES solverConvergence.cpp common.hpp common.cpp KeyValueStore.hpp Table.hpp
      DEPENDS hyteg)
