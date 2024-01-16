target_sources( hyteg
    PRIVATE
    EdgeDoFMacroCell.hpp
    EdgeDoFOrientation.hpp
    EdgeDoFAdditivePackInfo.hpp
    EdgeDoFProjectNormalOperator.cpp
    EdgeDoFPackInfo.hpp
    EdgeDoFFunction.hpp
    EdgeDoFAdditivePackInfo.cpp
    EdgeDoFProjectNormalOperator.hpp
    EdgeDoFPackInfo.cpp
    EdgeDoFFunction.cpp
    EdgeDoFIndexing.hpp
    EdgeDoFPetsc.hpp
    EdgeDoFMacroEdge.hpp
    EdgeDoFOperatorTypeDefs.hpp
    EdgeDoFMacroFace.hpp     
)
add_subdirectory( freeslip )
add_subdirectory( generatedKernels )
