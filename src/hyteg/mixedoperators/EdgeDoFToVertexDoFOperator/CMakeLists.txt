target_sources( hyteg
    PRIVATE
    EdgeDoFToVertexDoFPetsc.hpp
    EdgeDoFToVertexDoFOperator.hpp
    EdgeDoFToVertexDoFOperator.cpp
    EdgeDoFToVertexDoFApply.hpp     
)
add_subdirectory( generatedKernels )
