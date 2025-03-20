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
#include <hyteg/polynomial/elementwise/data.hpp>
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

/* convert local tet indices to local cube indices
   indices are encodes as i = x + 2y + 4z = x + y<<1 + z<<2
   where x,y,z ∈ {0,1}
*/
static constexpr surrogate::ElementTypeWiseData< std::array< uint_t, 4 >, 3 > cubeIndices{
    std::array< uint_t, 4 >{ 0b000, 0b001, 0b010, 0b100 }, // WHITE UP
    std::array< uint_t, 4 >{ 0b001, 0b011, 0b010, 0b101 }, // BLUE UP
    std::array< uint_t, 4 >{ 0b001, 0b010, 0b101, 0b100 }, // GREEN UP
    std::array< uint_t, 4 >{ 0b011, 0b111, 0b110, 0b101 }, // WHITE DOWN
    std::array< uint_t, 4 >{ 0b101, 0b110, 0b100, 0b010 }, // BLUE DOWN
    std::array< uint_t, 4 >{ 0b010, 0b011, 0b101, 0b110 }  // GREEN DOWN
};

/**
 * @brief Computes the global vertices of the micro-cube containing the six tets corresponding to the given index
 * @param lvl The refinement level of the grid. Determines the number of micro-vertices per edge.
 * @param microElement The micro-element specified by its x, y, and z indices.
 * @param globalDofIndices Array to store the global indices of the cube's vertices.
 */
static inline void
    getGlobalCubeIndices3D( const uint_t lvl, const indexing::Index& microElement, std::array< uint_t, 8 >& globalDofIndices )
{
   const auto x = microElement.x();
   const auto y = microElement.y();
   const auto z = microElement.z();
   const auto n = levelinfo::num_microvertices_per_edge( lvl );

   // local indices of micro cube
   // indices are encodes as i = dx + 2dy + 4dz = dx + dy<<1 + dz<<2 where dx,dy ∈ {0,1}
   globalDofIndices[0b000] = indexing::macroCellIndex( n, x + 0, y + 0, z + 0 );
   globalDofIndices[0b001] = indexing::macroCellIndex( n, x + 1, y + 0, z + 0 );
   globalDofIndices[0b010] = indexing::macroCellIndex( n, x + 0, y + 1, z + 0 );
   globalDofIndices[0b011] = indexing::macroCellIndex( n, x + 1, y + 1, z + 0 );
   globalDofIndices[0b100] = indexing::macroCellIndex( n, x + 0, y + 0, z + 1 );
   globalDofIndices[0b101] = indexing::macroCellIndex( n, x + 1, y + 0, z + 1 );
   globalDofIndices[0b110] = indexing::macroCellIndex( n, x + 0, y + 1, z + 1 );
   globalDofIndices[0b111] = indexing::macroCellIndex( n, x + 1, y + 1, z + 1 );
}

} // namespace p1
} // namespace hyteg
