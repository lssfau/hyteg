target_sources( hyteg
    PRIVATE
    VertexDoFMacroCell.hpp
    P1VariableOperator.hpp
    P1ProjectNormalOperator.hpp
    P1RotationOperator.hpp
    P1Operator.hpp
    VertexDoFFunction.cpp
    VertexDoFAdditivePackInfo.hpp
    VertexDoFIndexing.hpp
    P1Elements.hpp
    VertexDoFMemory.hpp
    P1HelperFunctions.hpp
    VertexDoFFunction.hpp
    VertexDoFMacroEdge.hpp
    VertexDoFPackInfo.hpp
    VertexDoFIndexing.cpp
    P1Petsc.hpp
    P1SurrogateOperator.hpp
    P1SurrogateOperator_optimized.hpp
    P1SurrogateOperator_optimized.cpp
    VertexDoFMacroFace.hpp
    P1VectorFunction.hpp
    P1ProjectNormalOperator.cpp
    P1RotationOperator.cpp
    P1Function.hpp
    VertexDoFMacroVertex.hpp
)
add_subdirectory( variablestencil )
add_subdirectory( freeslip )
add_subdirectory( generatedKernels )
add_subdirectory( polynomial )
