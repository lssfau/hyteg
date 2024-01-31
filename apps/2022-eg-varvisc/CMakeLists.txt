

if (HYTEG_BUILD_WITH_PETSC)
waLBerla_add_executable( NAME MultiSinkerBenchmark
        FILES MultiSinkerBenchmark.cpp
      DEPENDS hyteg core mixedOperator)

waLBerla_add_executable( NAME SolViBenchmark
        FILES SolViBenchmark.cpp
      DEPENDS hyteg core mixedOperator)

waLBerla_add_executable( NAME EGSmoothViscTests
        FILES EGSmoothViscTests.cpp
      DEPENDS hyteg core mixedOperator)

waLBerla_add_executable( NAME EGJumpTests.cpp
        FILES EGJumpTests.cpp
      DEPENDS hyteg core mixedOperator)
endif ()
