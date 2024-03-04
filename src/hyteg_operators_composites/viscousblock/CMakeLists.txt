add_library(opgen-composites-viscousblock INTERFACE
        P2ViscousBlockEpsilonOperator.hpp
        P2ViscousBlockFullOperator.hpp
        P2ViscousBlockLaplaceOperator.hpp
        )

target_link_libraries(opgen-composites-viscousblock INTERFACE hyteg opgen-diffusion opgen-epsilon opgen-full_stokes)