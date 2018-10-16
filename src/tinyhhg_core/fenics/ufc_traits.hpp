#pragma once

#include "tinyhhg_core/types/matrix.hpp"

class p1_tet_diffusion_cell_integral_0_otherwise;
class p1_tet_div_tet_cell_integral_0_otherwise;
class p1_tet_div_tet_cell_integral_1_otherwise;
class p1_tet_div_tet_cell_integral_2_otherwise;
class p1_tet_divt_tet_cell_integral_0_otherwise;
class p1_tet_divt_tet_cell_integral_1_otherwise;
class p1_tet_divt_tet_cell_integral_2_otherwise;
class p1_tet_mass_cell_integral_0_otherwise;
class p1_tet_pspg_tet_cell_integral_0_otherwise;

class p2_tet_diffusion_cell_integral_0_otherwise;

namespace hhg {
namespace fenics {

class UndefinedAssembly;
class NoAssemble;
class Dummy10x10Assembly;

template< typename UFCOperator >
struct UFCTrait;

#define HHG_CREATE_UFCOPERATOR_SPECIALIZATION( UFC_TYPE, LOCAL_STIFFNESS_MATRIX_TYPE ) \
   template<> \
   struct UFCTrait< UFC_TYPE > \
   { \
      typedef LOCAL_STIFFNESS_MATRIX_TYPE LocalStiffnessMatrix_T; \
   }


// P1 / Tet traits

HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p1_tet_diffusion_cell_integral_0_otherwise, Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p1_tet_div_tet_cell_integral_0_otherwise,   Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p1_tet_div_tet_cell_integral_1_otherwise,   Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p1_tet_div_tet_cell_integral_2_otherwise,   Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p1_tet_divt_tet_cell_integral_0_otherwise,  Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p1_tet_divt_tet_cell_integral_1_otherwise,  Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p1_tet_divt_tet_cell_integral_2_otherwise,  Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p1_tet_mass_cell_integral_0_otherwise,      Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p1_tet_pspg_tet_cell_integral_0_otherwise,  Matrix4r );


// P2 / Tet traits

HHG_CREATE_UFCOPERATOR_SPECIALIZATION( p2_tet_diffusion_cell_integral_0_otherwise, Matrix10r );

// other

HHG_CREATE_UFCOPERATOR_SPECIALIZATION( UndefinedAssembly,  Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( NoAssemble,         Matrix4r );
HHG_CREATE_UFCOPERATOR_SPECIALIZATION( Dummy10x10Assembly, Matrix10r );


#undef HHG_CREATE_UFCOPERATOR_SPECIALIZATION


}
}