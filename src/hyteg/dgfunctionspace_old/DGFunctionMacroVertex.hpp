/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

namespace hyteg {
namespace dgfunction {
namespace macrovertex {

template < typename ValueType >
inline void projectP1( const uint_t&                                                 Level,
                       Vertex&                                                       vertex,
                       const std::shared_ptr< PrimitiveStorage >&                    storage,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId,
                       UpdateType                                                    updateType )
{
   auto src = vertex.getData( srcId )->getPointer( Level );
   auto dst = vertex.getData( dstId )->getPointer( Level );

   ValueType tmp;

   for ( auto faceIt : vertex.neighborFaces() )
   {
      Face* face = storage->getFace( faceIt );

      uint_t localFaceId = vertex.face_index( face->getID() );

      uint_t faceMemoryIndex = 2 * localFaceId;

      std::vector< PrimitiveID > adjEdgeIds = face->adjacent_edges( vertex.getID() );
      std::vector< Edge* >       adjEdges;
      adjEdges.push_back( storage->getEdge( adjEdgeIds[0] ) );
      adjEdges.push_back( storage->getEdge( adjEdgeIds[1] ) );

      uint_t p1EdgeId0 = vertex.edge_index( adjEdgeIds[0] ) + 1;
      uint_t p1EdgeId1 = vertex.edge_index( adjEdgeIds[1] ) + 1;

      tmp = 1.0 / 3.0 * ( src[0] + src[p1EdgeId0] + src[p1EdgeId1] );

      if ( updateType == Replace )
      {
         dst[faceMemoryIndex] = tmp;
      }
      else if ( updateType == Add )
      {
         dst[faceMemoryIndex] += tmp;
      }
   }
}

} // namespace macrovertex
} // namespace dgfunction
} //namespace hyteg
