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
class p2_tet_mass_cell_integral_0_otherwise;
class p2_tet_pspg_tet_cell_integral_0_otherwise;
class p2_tet_div_tet_cell_integral_0_otherwise;
class p2_tet_div_tet_cell_integral_1_otherwise;
class p2_tet_div_tet_cell_integral_2_otherwise;
class p2_tet_divt_tet_cell_integral_0_otherwise;
class p2_tet_divt_tet_cell_integral_1_otherwise;
class p2_tet_divt_tet_cell_integral_2_otherwise;

class p2_to_p1_tet_div_tet_cell_integral_0_otherwise;
class p2_to_p1_tet_div_tet_cell_integral_1_otherwise;
class p2_to_p1_tet_div_tet_cell_integral_2_otherwise;

class p1_to_p2_tet_divt_tet_cell_integral_0_otherwise;
class p1_to_p2_tet_divt_tet_cell_integral_1_otherwise;
class p1_to_p2_tet_divt_tet_cell_integral_2_otherwise;

namespace hyteg {
namespace fenics {

class UndefinedAssembly;
class NoAssemble;
class Dummy10x10Assembly;

template< typename UFCOperator >
struct UFCTrait;

template<> struct UFCTrait< p1_tet_diffusion_cell_integral_0_otherwise > { typedef Matrix4r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_tet_div_tet_cell_integral_0_otherwise > { typedef Matrix4r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_tet_div_tet_cell_integral_1_otherwise > { typedef Matrix4r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_tet_div_tet_cell_integral_2_otherwise > { typedef Matrix4r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_tet_divt_tet_cell_integral_0_otherwise > { typedef Matrix4r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_tet_divt_tet_cell_integral_1_otherwise > { typedef Matrix4r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_tet_divt_tet_cell_integral_2_otherwise > { typedef Matrix4r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_tet_mass_cell_integral_0_otherwise > { typedef Matrix4r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_tet_pspg_tet_cell_integral_0_otherwise > { typedef Matrix4r LocalStiffnessMatrix_T; };

template<> struct UFCTrait< p2_tet_diffusion_cell_integral_0_otherwise > { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_tet_mass_cell_integral_0_otherwise >      { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_tet_pspg_tet_cell_integral_0_otherwise >  { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_tet_div_tet_cell_integral_0_otherwise >   { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_tet_div_tet_cell_integral_1_otherwise >   { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_tet_div_tet_cell_integral_2_otherwise >   { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_tet_divt_tet_cell_integral_0_otherwise >  { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_tet_divt_tet_cell_integral_1_otherwise >  { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_tet_divt_tet_cell_integral_2_otherwise >  { typedef Matrix10r LocalStiffnessMatrix_T; };

template<> struct UFCTrait< p2_to_p1_tet_div_tet_cell_integral_0_otherwise > { typedef Matrixr<4, 10> LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_to_p1_tet_div_tet_cell_integral_1_otherwise > { typedef Matrixr<4, 10> LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p2_to_p1_tet_div_tet_cell_integral_2_otherwise > { typedef Matrixr<4, 10> LocalStiffnessMatrix_T; };

template<> struct UFCTrait< p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > { typedef Matrixr<10, 4> LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > { typedef Matrixr<10, 4> LocalStiffnessMatrix_T; };
template<> struct UFCTrait< p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > { typedef Matrixr<10, 4> LocalStiffnessMatrix_T; };

template<> struct UFCTrait< UndefinedAssembly > { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< NoAssemble > { typedef Matrix10r LocalStiffnessMatrix_T; };
template<> struct UFCTrait< Dummy10x10Assembly > { typedef Matrix10r LocalStiffnessMatrix_T; };

}
}