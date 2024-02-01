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
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#include <map>
#include "core/Macros.h"
#define RESTRICT WALBERLA_RESTRICT
#include "sor_3D_macroedge_P2_update_vertexdofs_impl.hpp"

namespace hyteg {
namespace P2 {
namespace macroedge {
namespace generated {

void sor_3D_macroedge_P2_update_vertexdofs(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil);

} // namespace generated
} // namespace macroedge
} // namespace P2
} // namespace hyteg