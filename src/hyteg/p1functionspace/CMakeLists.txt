target_sources( hyteg
    PRIVATE
    VertexDoFMacroCell.hpp
    P1VariableOperator.hpp
    P1ProjectNormalOperator.hpp
    P1Operator.hpp
    P1DataHandling.hpp
    VertexDoFFunction.cpp
    VertexDoFAdditivePackInfo.hpp
    VertexDoFIndexing.hpp
    P1Elements.hpp
    VertexDoFMemory.hpp
    P1HelperFunctions.hpp
    P1ConstantOperator.hpp
    VertexDoFFunction.hpp
    VertexDoFMacroEdge.hpp
    VertexDoFPackInfo.hpp
    VertexDoFIndexing.cpp
    P1Petsc.hpp
    P1ConstantOperator.cpp
    P1SurrogateOperator.hpp
    VertexDoFMacroFace.hpp
    P1VectorFunction.hpp
    P1ProjectNormalOperator.cpp
    P1DataHandling.cpp
    P1Function.hpp
    VertexDoFMacroVertex.hpp     
)
add_subdirectory( variablestencil )
add_subdirectory( freeslip )
add_subdirectory( generatedKernels )
add_subdirectory( polynomial )
