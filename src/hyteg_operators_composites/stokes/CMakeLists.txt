add_library(opgen-composites-stokes INTERFACE
        P2P1StokesConstantOperator.hpp
        P2P1StokesEpsilonOperator.hpp
        P2P1StokesFullOperator.hpp
        )

target_link_libraries(opgen-composites-stokes INTERFACE hyteg opgen-composites-divergence opgen-composites-gradient opgen-composites-viscousblock)
