target_sources( hyteg
    PRIVATE
    P2FormDivKGrad.hpp
    QuadratureRules.hpp
    ShapeFunctionMacros.hpp
    P2FormDivKGrad.cpp
    P2FormLaplace.hpp     
)
add_subdirectory( kernels )
