target_sources( hyteg
    PRIVATE
    EdgeDoFMacroCell.hpp
    EdgeDoFOrientation.hpp
    EdgeDoFOperator.cpp
    EdgeDoFAdditivePackInfo.hpp
    EdgeDoFProjectNormalOperator.cpp
    EdgeDoFPackInfo.hpp
    EdgeDoFFunction.hpp
    EdgeDoFOperator.hpp
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
