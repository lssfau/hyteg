
waLBerla_add_executable( NAME MultiSinkerBenchmark
        FILES MultiSinkerBenchmark.cpp
        DEPENDS hyteg )

waLBerla_add_executable( NAME SolViBenchmark
        FILES SolViBenchmark.cpp
        DEPENDS hyteg )

waLBerla_add_executable( NAME SmoothViscTests
        FILES EGSmoothViscTests.cpp
        DEPENDS hyteg )

