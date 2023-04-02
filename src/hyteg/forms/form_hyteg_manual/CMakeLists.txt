target_sources( hyteg
    PRIVATE
    P2FormDivKGrad.hpp
    P2FormDivKGrad.cpp
    P2FormLaplace.hpp
    N1E1FormCurlCurl.hpp
    N1E1FormMass.hpp
    QuadratureRules.hpp
    ShapeFunctionMacros.hpp
    SphericalElementFormMass.hpp
)
add_subdirectory( kernels )
