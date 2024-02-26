add_library(opgen-composites-stokes-divergence INTERFACE
        P2ToP1DivergenceOperator.hpp
        )

target_link_libraries(opgen-composites-stokes-divergence INTERFACE hyteg opgen-divergence)