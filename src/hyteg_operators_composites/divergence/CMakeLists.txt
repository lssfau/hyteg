add_library(opgen-composites-divergence INTERFACE
        P2ToP1DivergenceOperator.hpp
        )

target_link_libraries(opgen-composites-divergence INTERFACE hyteg opgen-divergence)