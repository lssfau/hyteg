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
#include "core/DataTypes.h"
#include "core/Macros.h"
#define RESTRICT WALBERLA_RESTRICT

using walberla::real_t;

namespace hyteg {
namespace P2 {
namespace macrocell {
namespace generated {

void restrict_3D_macrocell_P2_update_edgedofs_level_1_to_0(real_t * RESTRICT _data_edgeCoarseDst_X, real_t * RESTRICT _data_edgeCoarseDst_XY, real_t * RESTRICT _data_edgeCoarseDst_XZ, real_t * RESTRICT _data_edgeCoarseDst_Y, real_t * RESTRICT _data_edgeCoarseDst_YZ, real_t * RESTRICT _data_edgeCoarseDst_Z, real_t const * RESTRICT const _data_edgeFineSrc_X, real_t const * RESTRICT const _data_edgeFineSrc_XY, real_t const * RESTRICT const _data_edgeFineSrc_XYZ, real_t const * RESTRICT const _data_edgeFineSrc_XZ, real_t const * RESTRICT const _data_edgeFineSrc_Y, real_t const * RESTRICT const _data_edgeFineSrc_YZ, real_t const * RESTRICT const _data_edgeFineSrc_Z, real_t const * RESTRICT const _data_vertexFineSrc, int coarse_level, double num_neighbor_cells_edge0, double num_neighbor_cells_edge1, double num_neighbor_cells_edge2, double num_neighbor_cells_edge3, double num_neighbor_cells_edge4, double num_neighbor_cells_edge5, double num_neighbor_cells_face0, double num_neighbor_cells_face1, double num_neighbor_cells_face2, double num_neighbor_cells_face3);

} // namespace generated
} // namespace macrocell
} // namespace P2
} // namespace hyteg