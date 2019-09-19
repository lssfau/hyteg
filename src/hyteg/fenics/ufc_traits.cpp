/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
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
