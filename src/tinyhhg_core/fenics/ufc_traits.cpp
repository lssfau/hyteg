#include "ufc_traits.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_pspg.h"
#include "tinyhhg_core/p1functionspace/generated/p1_stokes_epsilon.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_div_tet.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_divt_tet.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_pspg_tet.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_diffusion.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace hhg {
namespace fenics {

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

}
}