target_sources( hyteg
    PRIVATE
    P2Elements.hpp
    P2Elements3D.hpp
    P2Function.cpp
    P2Function.hpp
    P2FunctionApplyOperator.hpp
    P2MacroCell.cpp
    P2MacroCell.hpp
    P2MacroFace.cpp
    P2MacroFace.hpp
    P2MacroVertex.cpp
    P2MacroVertex.hpp     
    P2Multigrid.hpp
    P2PlusBubbleFunction.cpp
    P2PlusBubbleFunction.hpp
    P2PlusBubbleVectorFunction.hpp
    P2ProjectNormalOperator.cpp
    P2ProjectNormalOperator.hpp
    P2RotationOperator.cpp
    P2RotationOperator.hpp
    P2SurrogateOperator.hpp
    P2TransferOperators.hpp
    P2VariableOperator.hpp
    P2VectorFunction.hpp
)
add_subdirectory( variablestencil )
add_subdirectory( polynomial )
