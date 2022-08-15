/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include <array>

#include "core/Abort.h"

#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"

namespace hyteg {
namespace n1e1 {

using celldof::CellType;
using indexing::Index;

namespace macrocell {

/// Returns an array of the four logical micro-vertex-indices that span the micro-cell of the given indices and cell type.
/// This function is an improved version of \sa hyteg::celldofspace::macrocell::getMicroVerticesFromMicroCell that ensures
/// consistent orientation of edges.
/// An edge is defined to be oriented from the vertex with the lower index to the vertex with the higher index.
/// The returned array is sorted in such a way, that micro-edges have the same orientation as the corresponding macro-edges.
/// Furthermore, the orientation of an edge is independent of the adjacent micro-cell that is used to obtain the edge.
/// The orienation of edges is crucial for certain FE spaces, e.g. Nédélec spaces with edge DoFs.
inline std::array< Index, 4 > getMicroVerticesFromMicroCell( const Index& microCellIndex, const CellType& microCellType )
{
   const idx_t cellX = microCellIndex.x();
   const idx_t cellY = microCellIndex.y();
   const idx_t cellZ = microCellIndex.z();

   switch ( microCellType )
   {
   case CellType::WHITE_UP:
      return std::array< Index, 4 >( { { Index( cellX, cellY, cellZ ),
                                         Index( cellX + 1, cellY, cellZ ),
                                         Index( cellX, cellY + 1, cellZ ),
                                         Index( cellX, cellY, cellZ + 1 ) } } );
   case CellType::BLUE_UP:
      return std::array< Index, 4 >( { { Index( cellX + 1, cellY, cellZ ),
                                         Index( cellX, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY, cellZ + 1 ) } } );
   case CellType::GREEN_UP:
      return std::array< Index, 4 >( { { Index( cellX + 1, cellY, cellZ ),
                                         Index( cellX, cellY + 1, cellZ ),
                                         Index( cellX, cellY, cellZ + 1 ),
                                         Index( cellX + 1, cellY, cellZ + 1 ) } } );
   case CellType::WHITE_DOWN:
      return std::array< Index, 4 >( { { Index( cellX + 1, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY, cellZ + 1 ),
                                         Index( cellX, cellY + 1, cellZ + 1 ),
                                         Index( cellX + 1, cellY + 1, cellZ + 1 ) } } );
   case CellType::BLUE_DOWN:
      return std::array< Index, 4 >( { { Index( cellX, cellY + 1, cellZ ),
                                         Index( cellX, cellY, cellZ + 1 ),
                                         Index( cellX + 1, cellY, cellZ + 1 ),
                                         Index( cellX, cellY + 1, cellZ + 1 ) } } );
   case CellType::GREEN_DOWN:
      return std::array< Index, 4 >( { { Index( cellX, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY, cellZ + 1 ),
                                         Index( cellX, cellY + 1, cellZ + 1 ) } } );
   default:
      WALBERLA_ABORT( "Not implemented for this cell type." );
      break;
   }
   return std::array< Index, 4 >();
}

} // namespace macrocell

/// Return data indices for the edge dofs of a given micro-cell in a macro-cell
inline void getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( const indexing::Index&   microCellIndex,
                                                              const celldof::CellType& cellType,
                                                              const uint_t             level,
                                                              std::array< uint_t, 6 >& edgeDoFIndices )
{
   // get indices of micro-vertices forming the micro-cell
   std::array< indexing::Index, 4 > verts = n1e1::macrocell::getMicroVerticesFromMicroCell( microCellIndex, cellType );

   // loop over micro-vertex pairs
   std::array< std::pair< uint_t, uint_t >, 6 > pairs = { std::make_pair( 2, 3 ),
                                                          std::make_pair( 1, 3 ),
                                                          std::make_pair( 1, 2 ),
                                                          std::make_pair( 0, 3 ),
                                                          std::make_pair( 0, 2 ),
                                                          std::make_pair( 0, 1 ) };

   for ( uint_t k = 0; k < 6; ++k )
   {
      // generate IndexIncrements for vertices in the current pair
      uint_t n1 = pairs[k].first;
      uint_t n2 = pairs[k].second;

      indexing::IndexIncrement vertexIdx1( (int) verts[n1].x(), (int) verts[n1].y(), (int) verts[n1].z() );
      indexing::IndexIncrement vertexIdx2( (int) verts[n2].x(), (int) verts[n2].y(), (int) verts[n2].z() );

      // get edge orientation
      edgedof::EdgeDoFOrientation orientation = edgedof::calcEdgeDoFOrientation( vertexIdx1, vertexIdx2 );

      // get index for this edgeDoF
      indexing::IndexIncrement eIdx = edgedof::calcEdgeDoFIndex( vertexIdx1, vertexIdx2 );

      // finally get the data index of the edgeDoF
      edgeDoFIndices[k] = edgedof::macrocell::index( level, eIdx.x(), eIdx.y(), eIdx.z(), orientation );
   }
}

} // namespace n1e1
} // namespace hyteg
