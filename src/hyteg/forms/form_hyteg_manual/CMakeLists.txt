target_sources( hyteg
    PRIVATE
    P2FormDivKGrad.hpp
    QuadratureRules.hpp
    ShapeFunctionMacros.hpp
    P2FormDivKGrad.cpp
    P2FormLaplace.hpp
    N1E1FormCurlCurl.hpp
    N1E1FormMass.hpp
)
add_subdirectory( kernels )
