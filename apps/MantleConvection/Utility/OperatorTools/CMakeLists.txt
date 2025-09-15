target_sources( hyteg
    PRIVATE
    ABlockProjectionWrapper.hpp
    AdditiveOperator.hpp
    ApproximateSchurComplementOperator.hpp
    BBlockProjectionWrapper.hpp
    BlockPreconditioners.hpp
    BTBlockProjectionWrapper.hpp
    CombinedABlockOperator.hpp
    StabilisationProjectionWrapper.hpp
    SchurOperator.hpp
    OptimisationTypedefs.hpp
    OperatorUpdater.hpp
    BFBT.hpp
)

add_library(MC-operators INTERFACE)

target_link_libraries(MC-operators INTERFACE hyteg opgen-adiabatic_heating opgen-supg_adiabatic_heating opgen-advection opgen-supg_advection opgen-diffusion_inv_rho opgen-supg_diffusion opgen-div_k_grad opgen-divergence opgen-full_stokes opgen-grad_rho_rho)