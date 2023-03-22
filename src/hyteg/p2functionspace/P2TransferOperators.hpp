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

#include "core/DataTypes.h"

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitives/Vertex.hpp"

namespace hyteg {
namespace P2 {

using walberla::uint_t;

namespace macroface {

template < typename ValueType >
inline void prolongateP1ToP2( const uint_t&                                               Level,
                              const Face&                                                 face,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Face >& p2VertexDoFMemoryID,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Face >& p2EdgeDoFMemoryID,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Face >& p1VertexDoFMemoryID )
{
   auto p2Vertices = face.getData( p2VertexDoFMemoryID )->getPointer( Level );
   auto p2Edges    = face.getData( p2EdgeDoFMemoryID )->getPointer( Level );

   const auto p1Vertices = face.getData( p1VertexDoFMemoryID )->getPointer( Level );

   for ( const auto& it : vertexdof::macroface::Iterator( Level, 0 ) )
   {
      const uint_t idx = vertexdof::macroface::indexFromVertex( Level, it.x(), it.y(), stencilDirection::VERTEX_C );
      p2Vertices[idx]  = p1Vertices[idx];
   }

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      p2Edges[edgedof::macroface::indexFromHorizontalEdge( Level, it.x(), it.y(), stencilDirection::EDGE_HO_C )] =
          ValueType( 0.5 ) *
          ( p1Vertices[vertexdof::macroface::indexFromHorizontalEdge( Level, it.x(), it.y(), stencilDirection::VERTEX_W )] +
            p1Vertices[vertexdof::macroface::indexFromHorizontalEdge( Level, it.x(), it.y(), stencilDirection::VERTEX_E )] );

      p2Edges[edgedof::macroface::indexFromDiagonalEdge( Level, it.x(), it.y(), stencilDirection::EDGE_DI_C )] =
          ValueType( 0.5 ) *
          ( p1Vertices[vertexdof::macroface::indexFromDiagonalEdge( Level, it.x(), it.y(), stencilDirection::VERTEX_NW )] +
            p1Vertices[vertexdof::macroface::indexFromDiagonalEdge( Level, it.x(), it.y(), stencilDirection::VERTEX_SE )] );

      p2Edges[edgedof::macroface::indexFromVerticalEdge( Level, it.x(), it.y(), stencilDirection::EDGE_VE_C )] =
          ValueType( 0.5 ) *
          ( p1Vertices[vertexdof::macroface::indexFromVerticalEdge( Level, it.x(), it.y(), stencilDirection::VERTEX_N )] +
            p1Vertices[vertexdof::macroface::indexFromVerticalEdge( Level, it.x(), it.y(), stencilDirection::VERTEX_S )] );
   }
}

template < typename ValueType >
inline void restrictP2ToP1( const uint_t&                                               Level,
                            const Face&                                                 face,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Face >& p2VertexDoFMemoryID,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Face >& p2EdgeDoFMemoryID,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Face >& p1VertexDoFMemoryID )
{
   const auto p2Vertices = face.getData( p2VertexDoFMemoryID )->getPointer( Level );
   const auto p2Edges    = face.getData( p2EdgeDoFMemoryID )->getPointer( Level );

   auto p1Vertices = face.getData( p1VertexDoFMemoryID )->getPointer( Level );

   for ( const auto& it : vertexdof::macroface::Iterator( Level, 1 ) )
   {
      ValueType    tmp;
      const uint_t idx = vertexdof::macroface::indexFromVertex( Level, it.x(), it.y(), stencilDirection::VERTEX_C );
      tmp              = p2Vertices[idx];

      const std::array< stencilDirection, 6 > directVertexDoFNeighbors = { {
          stencilDirection::EDGE_HO_W,
          stencilDirection::EDGE_HO_E,
          stencilDirection::EDGE_DI_NW,
          stencilDirection::EDGE_DI_SE,
          stencilDirection::EDGE_VE_N,
          stencilDirection::EDGE_VE_S,
      } };

      for ( const auto& neighbor : directVertexDoFNeighbors )
      {
         const uint_t neighborEdgeIdx = edgedof::macroface::indexFromVertex( Level, it.x(), it.y(), neighbor );
         tmp += ValueType( 0.5 ) * p2Edges[neighborEdgeIdx];
      }

      p1Vertices[idx] = tmp;
   }
}

} // namespace macroface

namespace macroedge {

template < typename ValueType >
inline void prolongateP1ToP2( const uint_t&                                               Level,
                              const Edge&                                                 edge,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& p2VertexDoFMemoryID,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& p2EdgeDoFMemoryID,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& p1VertexDoFMemoryID )
{
   auto p2Vertices = edge.getData( p2VertexDoFMemoryID )->getPointer( Level );
   auto p2Edges    = edge.getData( p2EdgeDoFMemoryID )->getPointer( Level );

   const auto p1Vertices = edge.getData( p1VertexDoFMemoryID )->getPointer( Level );

   for ( const auto& it : edgedof::macroedge::Iterator( Level, 0 ) )
   {
      p2Edges[edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C )] = 0;
   }

   for ( const auto& it : vertexdof::macroedge::Iterator( Level, 1 ) )
   {
      const uint_t    idx               = vertexdof::macroedge::indexFromVertex( Level, it.x(), stencilDirection::VERTEX_C );
      const ValueType p1VertexValue     = p1Vertices[idx];
      const ValueType p1VertexValueHalf = ValueType( 0.5 ) * p1VertexValue;

      p2Vertices[idx] = p1VertexValue;

      for ( const auto& neighbor : edgedof::macroedge::neighborsOnEdgeFromVertex )
      {
         const uint_t neighborEdgeIdx = edgedof::macroedge::indexFromVertex( Level, it.x(), neighbor );
         p2Edges[neighborEdgeIdx] += p1VertexValueHalf;
      }
   }
}

template < typename ValueType >
inline void restrictP2ToP1( const uint_t&                                               Level,
                            const Edge&                                                 edge,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& p2VertexDoFMemoryID,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& p2EdgeDoFMemoryID,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& p1VertexDoFMemoryID )
{
   const auto p2Vertices = edge.getData( p2VertexDoFMemoryID )->getPointer( Level );
   const auto p2Edges    = edge.getData( p2EdgeDoFMemoryID )->getPointer( Level );

   auto p1Vertices = edge.getData( p1VertexDoFMemoryID )->getPointer( Level );

   for ( const auto& it : vertexdof::macroedge::Iterator( Level, 1 ) )
   {
      ValueType    tmp;
      const uint_t idx = vertexdof::macroedge::indexFromVertex( Level, it.x(), stencilDirection::VERTEX_C );
      tmp              = p2Vertices[idx];

      tmp += ValueType( 0.5 ) * p2Edges[edgedof::macroedge::indexFromVertex( Level, it.x(), stencilDirection::EDGE_HO_W )];
      tmp += ValueType( 0.5 ) * p2Edges[edgedof::macroedge::indexFromVertex( Level, it.x(), stencilDirection::EDGE_HO_E )];

      tmp += ValueType( 0.5 ) * p2Edges[edgedof::macroedge::indexFromVertex( Level, it.x(), stencilDirection::EDGE_VE_S )];
      tmp += ValueType( 0.5 ) * p2Edges[edgedof::macroedge::indexFromVertex( Level, it.x(), stencilDirection::EDGE_DI_SE )];

      if ( edge.getNumNeighborFaces() == 2 )
      {
         tmp += ValueType( 0.5 ) * p2Edges[edgedof::macroedge::indexFromVertex( Level, it.x(), stencilDirection::EDGE_VE_N )];
         tmp += ValueType( 0.5 ) * p2Edges[edgedof::macroedge::indexFromVertex( Level, it.x(), stencilDirection::EDGE_DI_NW )];
      }

      p1Vertices[idx] = tmp;
   }
}

} // namespace macroedge

namespace macrovertex {

template < typename ValueType >
inline void prolongateP1ToP2( const uint_t&                                                 Level,
                              const Vertex&                                                 vertex,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& p2VertexDoFMemoryID,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& p2EdgeDoFMemoryID,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& p1VertexDoFMemoryID )
{
   WALBERLA_UNUSED( p2EdgeDoFMemoryID );

   auto       p2Vertices = vertex.getData( p2VertexDoFMemoryID )->getPointer( Level );
   const auto p1Vertices = vertex.getData( p1VertexDoFMemoryID )->getPointer( Level );

   p2Vertices[0] = p1Vertices[0];
}

template < typename ValueType >
inline void restrictP2ToP1( const uint_t&                                                 Level,
                            const Vertex&                                                 vertex,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& p2VertexDoFMemoryID,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& p2EdgeDoFMemoryID,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& p1VertexDoFMemoryID )
{
   const auto p2Vertices = vertex.getData( p2VertexDoFMemoryID )->getPointer( Level );
   const auto p2Edges    = vertex.getData( p2EdgeDoFMemoryID )->getPointer( Level );

   auto p1Vertices = vertex.getData( p1VertexDoFMemoryID )->getPointer( Level );

   p1Vertices[0] = p2Vertices[0];

   const uint_t numNeighborEdges = vertex.getNumNeighborEdges();

   for ( uint_t i = 0; i < numNeighborEdges; i++ )
   {
      p1Vertices[0] += ValueType( 0.5 ) * p2Edges[i];
   }
}

} // namespace macrovertex

} // namespace P2
} // namespace hyteg
