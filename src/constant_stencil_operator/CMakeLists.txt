add_library( constant_stencil_operator)
target_link_libraries( constant_stencil_operator PUBLIC hyteg PRIVATE core)

target_sources( constant_stencil_operator
      PRIVATE
      EdgeDoFOperator.hpp
      EdgeDoFOperator.cpp
      P1ConstantOperator.cpp
      P1ConstantOperator.hpp
      P1EpsilonStokesOperator.hpp
      P1Transport.hpp
      P2ConstantOperator.cpp
      P2ConstantOperator.hpp
      P2ConstantEpsilonOperator.hpp
      P2ConstantFullViscousOperator.hpp
      P2MacroEdge.cpp
      P2MacroEdge.hpp
)
add_subdirectory( EdgeDoFGeneratedKernels )
add_subdirectory( P1generatedKernels )
add_subdirectory( P2generatedKernels )