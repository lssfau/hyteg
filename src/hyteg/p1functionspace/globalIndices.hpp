/*
 * Copyright (c) 2025 Benjamin Mann
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

#include <hyteg/Levelinfo.hpp>
#include <hyteg/volumedofspace/CellDoFIndexing.hpp>
#include <hyteg/volumedofspace/FaceDoFIndexing.hpp>

namespace hyteg {
namespace p1 {
// alternative to vertexdof::getVertexDoFDataIndicesFromMicroCell() to improve performance
static inline void getGlobalIndices3D( const celldof::CellType  cType,
                                const uint_t             lvl,
                                const indexing::Index&   microElement,
                                std::array< uint_t, 4 >& globalDofIndices )
{
   // get local indices of micro vertices
   std::array< indexing::Index, 4 > v = celldof::macrocell::getMicroVerticesFromMicroCell( microElement, cType );
   // convert to global indexing
   const auto n = levelinfo::num_microvertices_per_edge( lvl );
   for ( uint_t k = 0; k < 4; ++k )
   {
      globalDofIndices[k] = indexing::macroCellIndex( n, v[k].x(), v[k].y(), v[k].z() );
   }
}

// alternative to vertexdof::getVertexDoFDataIndicesFromMicroFace() to improve performance
static inline void getGlobalIndices2D( const facedof::FaceType  fType,
                                const uint_t             lvl,
                                const indexing::Index&   microElement,
                                std::array< uint_t, 3 >& globalDofIndices )
{
   // get indices of micro vertices
   std::array< indexing::Index, 3 > v = facedof::macroface::getMicroVerticesFromMicroFace( microElement, fType );
   // convert to global indexing
   const auto n = levelinfo::num_microvertices_per_edge( lvl );
   for ( uint_t k = 0; k < 3; ++k )
   {
      globalDofIndices[k] = indexing::macroFaceIndex( n, v[k].x(), v[k].y() );
   }
}

} // namespace p1

} // namespace hyteg
