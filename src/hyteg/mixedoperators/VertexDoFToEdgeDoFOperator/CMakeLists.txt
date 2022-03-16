target_sources( hyteg
    PRIVATE
    VertexDoFToEdgeDoFOperator.hpp
    VertexDoFToEdgeDoFApply.hpp
    VertexDoFToEdgeDoFPetsc.hpp
    VertexDoFToEdgeDoFOperator.cpp     
)
add_subdirectory( generatedKernels )
