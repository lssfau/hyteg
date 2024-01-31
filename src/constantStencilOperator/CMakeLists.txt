add_library(constantStencilOperator)
target_link_libraries(constantStencilOperator PUBLIC hyteg PRIVATE core)

target_sources(constantStencilOperator
      PRIVATE
      P1ConstantOperator.cpp
      P1ConstantOperator.hpp
      P1EpsilonStokesOperator.hpp
      P1Transport.hpp
      P2ConstantOperator.cpp
      P2ConstantOperator.hpp
      P2ConstantEpsilonOperator.hpp
      P2ConstantFullViscousOperator.hpp
)