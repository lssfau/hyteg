/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/facedofspace_old/FaceDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"


namespace hyteg {
namespace dgfunction {
namespace macroedge {

template < typename ValueType >
inline void projectP1( const uint_t&                                               Level,
                       Edge&                                                       edge,
                       const std::shared_ptr< PrimitiveStorage >&                  storage,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId,
                       UpdateType                                                  updateType )
{
   auto src = edge.getData( srcId )->getPointer( Level );
   auto dst = edge.getData( dstId )->getPointer( Level );

   size_t    rowsize = levelinfo::num_microvertices_per_edge( Level );
   ValueType tmp;

   // first face (south)
   {
      for ( uint_t i = 1; i < rowsize - 2; ++i )
      {
         tmp = 1.0 / 3.0 *
               ( src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                 src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] +
                 src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_SE )] );

         if ( updateType == Replace )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )] += tmp;
         }
      }
   }

   // second face (north)
   if ( edge.getNumNeighborFaces() == 2 )
   {
      for ( uint_t i = 1; i < rowsize - 2; ++i )
      {
         tmp = 1.0 / 3.0 *
               ( src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                 src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] +
                 src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_N )] );

         if ( updateType == Replace )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )] += tmp;
         }
      }
   }
}

} // namespace macroedge
} // namespace dgfunction
} //namespace hyteg
