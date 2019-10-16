/*
 * Copyright (c) 2019 Nils Kohl, Dominik Thoennes.
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

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "all.hpp"
#include "core/Macros.h"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace vertexdof {
namespace macrocell {
namespace generated {

void sor_3D_macrocell_P1_colored(double * RESTRICT _data_p1CellDst_group_0, double const * RESTRICT const _data_p1CellDst_group_0_const, double * RESTRICT _data_p1CellDst_group_1, double const * RESTRICT const _data_p1CellDst_group_1_const, double * RESTRICT _data_p1CellDst_group_2, double const * RESTRICT const _data_p1CellDst_group_2_const, double * RESTRICT _data_p1CellDst_group_3, double const * RESTRICT const _data_p1CellDst_group_3_const, double * RESTRICT _data_p1CellDst_group_4, double const * RESTRICT const _data_p1CellDst_group_4_const, double * RESTRICT _data_p1CellDst_group_5, double const * RESTRICT const _data_p1CellDst_group_5_const, double * RESTRICT _data_p1CellDst_group_6, double const * RESTRICT const _data_p1CellDst_group_6_const, double * RESTRICT _data_p1CellDst_group_7, double const * RESTRICT const _data_p1CellDst_group_7_const, double const * RESTRICT const _data_p1CellRhs_group_0_const, double const * RESTRICT const _data_p1CellRhs_group_1_const, double const * RESTRICT const _data_p1CellRhs_group_2_const, double const * RESTRICT const _data_p1CellRhs_group_3_const, double const * RESTRICT const _data_p1CellRhs_group_4_const, double const * RESTRICT const _data_p1CellRhs_group_5_const, double const * RESTRICT const _data_p1CellRhs_group_6_const, double const * RESTRICT const _data_p1CellRhs_group_7_const, int64_t color_group, int32_t level, std::map< hyteg::indexing::IndexIncrement, double > p1CellStencil, double relax);

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg