target_sources( hyteg
    PRIVATE
    MixedDummyOperators.hpp
    P0ScalarToP1VectorOperator.hpp
    P0ToP1Operator.hpp
    P1ToP0Operator.hpp
    P1ToP2SurrogateOperator.hpp
    P1ToP2VariableOperator.hpp
    P1VectorToP0ScalarOperator.hpp
    P0ScalarToP1VectorOperator.hpp
    P1ToDGOperator.hpp
    DGToP1Operator.hpp
    P2ToP1SurrogateOperator.hpp
    P2ToP1VariableOperator.hpp
)
add_subdirectory( variablestencil )
add_subdirectory( VertexDoFToEdgeDoFOperator )
add_subdirectory( EdgeDoFToVertexDoFOperator )
add_subdirectory( polynomial )
