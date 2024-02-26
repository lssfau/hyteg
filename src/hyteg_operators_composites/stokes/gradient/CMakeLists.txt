add_library(opgen-composites-stokes-gradient INTERFACE
        P1ToP2GradientOperator.hpp
        )

target_link_libraries(opgen-composites-stokes-gradient INTERFACE hyteg opgen-gradient)