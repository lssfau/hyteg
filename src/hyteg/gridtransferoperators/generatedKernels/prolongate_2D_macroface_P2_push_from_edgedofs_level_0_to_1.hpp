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
namespace P2 {
namespace macroface {
namespace generated {

void prolongate_2D_macroface_P2_push_from_edgedofs_level_0_to_1(double * RESTRICT _data_edgeCoarseSrc_X, double * RESTRICT _data_edgeCoarseSrc_XY, double * RESTRICT _data_edgeCoarseSrc_Y, double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexFineDst, int32_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg