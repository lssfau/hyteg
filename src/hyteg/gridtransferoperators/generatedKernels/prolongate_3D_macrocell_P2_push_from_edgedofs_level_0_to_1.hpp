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
namespace macrocell {
namespace generated {

void prolongate_3D_macrocell_P2_push_from_edgedofs_level_0_to_1(double const * RESTRICT const _data_edgeCoarseSrc_X, double const * RESTRICT const _data_edgeCoarseSrc_XY, double const * RESTRICT const _data_edgeCoarseSrc_XZ, double const * RESTRICT const _data_edgeCoarseSrc_Y, double const * RESTRICT const _data_edgeCoarseSrc_YZ, double const * RESTRICT const _data_edgeCoarseSrc_Z, double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_XYZ, double * RESTRICT _data_edgeFineDst_XZ, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_edgeFineDst_YZ, double * RESTRICT _data_edgeFineDst_Z, double * RESTRICT _data_vertexFineDst, int coarse_level, double num_neighbor_cells_edge0, double num_neighbor_cells_edge1, double num_neighbor_cells_edge2, double num_neighbor_cells_edge3, double num_neighbor_cells_edge4, double num_neighbor_cells_edge5, double num_neighbor_cells_face0, double num_neighbor_cells_face1, double num_neighbor_cells_face2, double num_neighbor_cells_face3);

} // namespace generated
} // namespace macrocell
} // namespace P2
} // namespace hyteg