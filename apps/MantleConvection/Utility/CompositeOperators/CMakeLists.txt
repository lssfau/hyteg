add_library(MC-composites INTERFACE
        P1ToP2GradientOperator.hpp
        P2RhoGMassOperator.hpp
        P2ToP1DivergenceOperator.hpp
        P2toP1GradRhoRhoOperator.hpp
        P2ToP1GradRhoRhoDivergenceOperator.hpp
        )

target_link_libraries(MC-composites INTERFACE hyteg MC-divergence MC-gradient MC-grad_rho_rho MC-rho_g_mass MC-grad_rho_rho_divergence)