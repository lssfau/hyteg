target_sources( hyteg
    PRIVATE
    P2ToP1VariableOperator.hpp
    P2ToP1Operator.hpp
    P1ToP2SurrogateOperator.hpp
    P2ToP1SurrogateOperator.hpp
    P1ToP2VariableOperator.hpp
    MixedDummyOperators.hpp
    P1ToP2Operator.hpp
    P1ToP0Operator.hpp
    P0ToP1Operator.hpp
    P1VectorToP0ScalarOperator.hpp
    P0ScalarToP1VectorOperator.hpp
)
add_subdirectory( variablestencil )
add_subdirectory( VertexDoFToEdgeDoFOperator )
add_subdirectory( EdgeDoFToVertexDoFOperator )
add_subdirectory( polynomial )
