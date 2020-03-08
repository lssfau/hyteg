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
#include "core/Macros.h"
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

void sor_2D_macroface_P2_update_edgedofs(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceRhs_X, double * RESTRICT _data_edgeFaceRhs_XY, double * RESTRICT _data_edgeFaceRhs_Y, double const * RESTRICT const _data_edge_stencil_at_edge_x, double const * RESTRICT const _data_edge_stencil_at_edge_xy, double const * RESTRICT const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * RESTRICT const _data_vertex_stencil_at_edge_x, double const * RESTRICT const _data_vertex_stencil_at_edge_xy, double const * RESTRICT const _data_vertex_stencil_at_edge_y, int level, double relax);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg