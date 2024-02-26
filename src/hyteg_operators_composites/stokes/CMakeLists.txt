add_subdirectory(divergence)
add_subdirectory(gradient)
add_subdirectory(viscousblock)

add_library(opgen-composites-stokes INTERFACE
        P2P1StokesConstantOperator.hpp
        P2P1StokesEpsilonOperator.hpp
        )

target_link_libraries(opgen-composites-stokes INTERFACE hyteg opgen-composites-stokes-divergence opgen-composites-stokes-gradient opgen-composites-stokes-viscousblock)
