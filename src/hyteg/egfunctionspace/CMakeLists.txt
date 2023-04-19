target_sources( hyteg
    PRIVATE
    EGFunction.hpp
    EGFunction.cpp
    EGBasis.hpp
    EGOperators.hpp
    EGConstEpsilonForm.hpp
    EGEpsilonForm.hpp
    EGVectorLaplaceForm.hpp
    EGNIPGVectorLaplaceForm.hpp
    EGIIPGVectorLaplaceForm.hpp
    EGDivForm.hpp
    EGDivtForm.hpp
        EGOperatorsNitscheBC.hpp
    EGVectorLaplaceFormNitscheBC.hpp
    EGConvTestUtils.hpp
)
