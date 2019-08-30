#include "ufc_traits.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/forms/form_fenics_generated/p1_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p1_div.h"
#include "hyteg/forms/form_fenics_generated/p1_divt.h"
#include "hyteg/forms/form_fenics_generated/p1_mass.h"
#include "hyteg/forms/form_fenics_generated/p1_pspg.h"
#include "hyteg/forms/form_fenics_generated/p1_stokes_epsilon.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_div_tet.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_divt_tet.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_mass.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_pspg_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_diffusion.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace hyteg {
namespace fenics {

//template struct UFCTrait< p1_tet_mass_cell_integral_0_otherwise >;
//template struct UFCTrait< p1_tet_diffusion_cell_integral_0_otherwise >;
//template struct UFCTrait< p1_tet_div_tet_cell_integral_0_otherwise >;
//template struct UFCTrait< p1_tet_div_tet_cell_integral_1_otherwise >;
//template struct UFCTrait< p1_tet_div_tet_cell_integral_2_otherwise >;
//template struct UFCTrait< p1_tet_divt_tet_cell_integral_0_otherwise >;
//template struct UFCTrait< p1_tet_divt_tet_cell_integral_1_otherwise >;
//template struct UFCTrait< p1_tet_divt_tet_cell_integral_2_otherwise >;
//template struct UFCTrait< p1_tet_pspg_tet_cell_integral_0_otherwise >;
//
//template struct UFCTrait< p2_tet_diffusion_cell_integral_0_otherwise >;
//
//template struct UFCTrait< UndefinedAssembly >;
//template struct UFCTrait< NoAssemble >;
//template struct UFCTrait< Dummy10x10Assembly >;

} // namespace fenics
} // namespace hyteg
