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

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

template < typename ValueType >
void prolongate_2D_macroface_P1_push_additive( ValueType const* RESTRICT const _data_vertexCoarseSrc,
                                               ValueType* RESTRICT             _data_vertexFineDst,
                                               int                             coarse_level,
                                               ValueType                       num_neighbor_faces_edge0,
                                               ValueType                       num_neighbor_faces_edge1,
                                               ValueType                       num_neighbor_faces_edge2,
                                               ValueType                       num_neighbor_faces_vertex0,
                                               ValueType                       num_neighbor_faces_vertex1,
                                               ValueType                       num_neighbor_faces_vertex2 );

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg