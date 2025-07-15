target_sources( hyteg
  PRIVATE
  CoefficientQuadratureAveraging.hpp
  DGForm.hpp
  Form.hpp
  N1E1Form.hpp
  N1E1LinearCombinationForm.hpp
  P0Form.cpp
  P0Form.hpp
  P1Form.cpp
  P1Form.hpp
  P1LinearCombinationForm.hpp
  P1RowSumForm.hpp
  P1WrapperForm.hpp
  P1ZeroForm.hpp
  P2Form.cpp
  P2Form.hpp
  P2LinearCombinationForm.hpp
  P2RowSumForm.hpp
)
add_subdirectory( form_hyteg_base )
add_subdirectory( form_hyteg_dg )
add_subdirectory( form_hyteg_manual )
add_subdirectory( form_fenics_base )
add_subdirectory( form_hyteg_generated )
add_subdirectory( form_fenics_generated )

