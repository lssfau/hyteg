/*
 * Copyright (c) 2019-2023 Nils Kohl, Dominik Thoennes, Michael Zikeli.
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
#include "communicate_directly_vertexdof_cell_to_face_impl.hpp"

namespace hyteg {
namespace vertexdof {
namespace comm {
namespace generated {

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level, int neighbor_cell_local_vertex_id_0, int neighbor_cell_local_vertex_id_1, int neighbor_cell_local_vertex_id_2);

} // namespace generated
} // namespace comm
} // namespace vertexdof
} // namespace hyteg