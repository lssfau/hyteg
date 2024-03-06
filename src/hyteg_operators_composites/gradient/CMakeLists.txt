add_library(opgen-composites-gradient INTERFACE
        P1ToP2GradientOperator.hpp
        )

target_link_libraries(opgen-composites-gradient INTERFACE hyteg opgen-gradient)