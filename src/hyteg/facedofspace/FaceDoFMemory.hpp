/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#pragma once

#include "hyteg/Levelinfo.hpp"
#include "core/debug/CheckFunctions.h"

namespace hyteg {

using walberla::uint_t;

/// Vertex Memory layout
/// the vertex memory has two entries for each adjacent face,
/// where the first entry is the Gray Face DoF and the Second one is the Blue Face DoF
/// the Gray Face DoF is owned by the Vertex and the Blue Face DoF is owned by the Face
inline uint_t faceDoFMacroVertexFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  return primitive.getNumNeighborFaces() * 2;
}

inline uint_t faceDoFMacroEdgeFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  const size_t num_cell_dofs = primitive.getNumNeighborFaces() * ( 2 * levelinfo::num_microedges_per_edge( level ) - 1 );
  return num_cell_dofs;
}

inline uint_t faceDoFMacroFaceFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( primitive );
  return levelinfo::num_microfaces_per_face( level );
}


}
