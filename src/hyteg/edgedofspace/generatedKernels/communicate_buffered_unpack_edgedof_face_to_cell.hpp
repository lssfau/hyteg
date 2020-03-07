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
#include "communicate_buffered_unpack_edgedof_face_to_cell_impl.hpp"

namespace hyteg {
namespace edgedof {
namespace comm {
namespace generated {

void communicate_buffered_unpack_edgedof_face_to_cell(double * RESTRICT _data_edge_cell_dst_X, double * RESTRICT _data_edge_cell_dst_XY, double * RESTRICT _data_edge_cell_dst_XZ, double * RESTRICT _data_edge_cell_dst_Y, double * RESTRICT _data_edge_cell_dst_YZ, double * RESTRICT _data_edge_cell_dst_Z, double const * RESTRICT const _data_recv_buffer, int level, int neighbor_cell_local_vertex_id_0, int neighbor_cell_local_vertex_id_1, int neighbor_cell_local_vertex_id_2, int recv_buffer_first_element_idx);

} // namespace generated
} // namespace comm
} // namespace edgedof
} // namespace hyteg