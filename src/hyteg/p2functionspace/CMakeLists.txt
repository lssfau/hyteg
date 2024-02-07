target_sources( hyteg
    PRIVATE
    P2MacroVertex.cpp
    P2Elements.hpp
    P2VariableOperator.hpp
    P2Function.cpp
    P2SurrogateOperator.hpp
    P2TransferOperators.hpp
    P2ProjectNormalOperator.cpp
    P2MacroFace.cpp
    P2Elements3D.hpp
    P2MacroCell.hpp
    P2Function.hpp
    P2ProjectNormalOperator.hpp
    P2MacroFace.hpp
    P2VectorFunction.hpp
    P2MacroCell.cpp
    P2Multigrid.hpp
    P2MacroVertex.hpp     
    P2FunctionApplyOperator.hpp
)
add_subdirectory( variablestencil )
add_subdirectory( polynomial )
