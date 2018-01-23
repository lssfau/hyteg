
#pragma once

#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "core/DataTypes.h"

namespace hhg {
namespace P2 {

using walberla::uint_t;

namespace macroface {

template< typename ValueType, uint_t Level >
inline void prolongateP1ToP2Tmpl( const Face & face,
                                  const PrimitiveDataID< FunctionMemory< ValueType >, Face > & p2VertexDoFMemoryID,
                                  const PrimitiveDataID< FunctionMemory< ValueType >, Face > & p2EdgeDoFMemoryID,
                                  const PrimitiveDataID< FunctionMemory< ValueType >, Face > & p1VertexDoFMemoryID )
{
        auto p2Vertices = face.getData( p2VertexDoFMemoryID )->getPointer( Level );
        auto p2Edges    = face.getData( p2EdgeDoFMemoryID   )->getPointer( Level );

  const auto p1Vertices = face.getData( p1VertexDoFMemoryID )->getPointer( Level );

  for ( const auto & it : vertexdof::macroface::Iterator( Level, 0 ) )
  {
    const uint_t idx = vertexdof::macroface::indexFromVertex< Level >( it.col(), it.row(), stencilDirection::VERTEX_C );
    p2Vertices[ idx ] = p1Vertices[ idx ];
  }

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    p2Edges[ edgedof::macroface::indexFromHorizontalEdge< Level >( it.col(), it.row(), stencilDirection::EDGE_HO_C ) ] =
        0.5 * (   p1Vertices[ vertexdof::macroface::indexFromHorizontalEdge< Level >( it.col(), it.row(), stencilDirection::VERTEX_W ) ]
                + p1Vertices[ vertexdof::macroface::indexFromHorizontalEdge< Level >( it.col(), it.row(), stencilDirection::VERTEX_E ) ] );

    p2Edges[ edgedof::macroface::indexFromDiagonalEdge< Level >( it.col(), it.row(), stencilDirection::EDGE_DI_C ) ] =
        0.5 * (   p1Vertices[ vertexdof::macroface::indexFromDiagonalEdge< Level >( it.col(), it.row(), stencilDirection::VERTEX_NW ) ]
                + p1Vertices[ vertexdof::macroface::indexFromDiagonalEdge< Level >( it.col(), it.row(), stencilDirection::VERTEX_SE ) ] );

    p2Edges[ edgedof::macroface::indexFromVerticalEdge< Level >( it.col(), it.row(), stencilDirection::EDGE_VE_C ) ] =
        0.5 * (   p1Vertices[ vertexdof::macroface::indexFromVerticalEdge< Level >( it.col(), it.row(), stencilDirection::VERTEX_N ) ]
                + p1Vertices[ vertexdof::macroface::indexFromVerticalEdge< Level >( it.col(), it.row(), stencilDirection::VERTEX_S ) ] );
  }

}

SPECIALIZE_WITH_VALUETYPE(void, prolongateP1ToP2Tmpl, prolongateP1ToP2)

template< typename ValueType, uint_t Level >
inline void restrictP2ToP1Tmpl( const Face & face,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Face > & p2VertexDoFMemoryID,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Face > & p2EdgeDoFMemoryID,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Face > & p1VertexDoFMemoryID )
{
  const auto p2Vertices = face.getData( p2VertexDoFMemoryID )->getPointer( Level );
  const auto p2Edges    = face.getData( p2EdgeDoFMemoryID   )->getPointer( Level );

        auto p1Vertices = face.getData( p1VertexDoFMemoryID )->getPointer( Level );

  for ( const auto & it : vertexdof::macroface::Iterator( Level, 1 ) )
  {
    ValueType tmp;
    const uint_t idx = vertexdof::macroface::indexFromVertex< Level >( it.col(), it.row(), stencilDirection::VERTEX_C );
    tmp = p2Vertices[ idx ];

    const std::array< stencilDirection, 6 > directVertexDoFNeighbors = {{
                                                                          stencilDirection::EDGE_HO_W,
                                                                          stencilDirection::EDGE_HO_E,
                                                                          stencilDirection::EDGE_DI_NW,
                                                                          stencilDirection::EDGE_DI_SE,
                                                                          stencilDirection::EDGE_VE_N,
                                                                          stencilDirection::EDGE_VE_S,
                                                                       }};

    for ( const auto & neighbor : directVertexDoFNeighbors )
    {
      const uint_t neighborEdgeIdx = edgedof::macroface::indexFromVertex< Level >( it.col(), it.row(), neighbor );
      tmp += 0.5 * p2Edges[ neighborEdgeIdx ];
    }

    p1Vertices[ idx ] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, restrictP2ToP1Tmpl, restrictP2ToP1)

} // namespace macroface

namespace macroedge {

template< typename ValueType, uint_t Level >
inline void prolongateP1ToP2Tmpl( const Edge & edge,
                                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & p2VertexDoFMemoryID,
                                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & p2EdgeDoFMemoryID,
                                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & p1VertexDoFMemoryID )
{
        auto p2Vertices = edge.getData( p2VertexDoFMemoryID )->getPointer( Level );
        auto p2Edges    = edge.getData( p2EdgeDoFMemoryID   )->getPointer( Level );

  const auto p1Vertices = edge.getData( p1VertexDoFMemoryID )->getPointer( Level );

  for ( const auto & it : edgedof::macroedge::Iterator( Level, 0 ) )
  {
    p2Edges[ edgedof::macroedge::indexFromHorizontalEdge< Level >( it.col(), stencilDirection::EDGE_HO_C ) ] = 0;
  }

  for ( const auto & it : vertexdof::macroedge::Iterator( Level, 1 ) )
  {
    const uint_t idx = vertexdof::macroedge::indexFromVertex< Level >( it.col(), stencilDirection::VERTEX_C );
    const ValueType p1VertexValue     = p1Vertices[ idx ];
    const ValueType p1VertexValueHalf = 0.5 * p1VertexValue;

    p2Vertices[ idx ] = p1VertexValue;

    for ( const auto & neighbor : edgedof::macroedge::neighborsOnEdgeFromVertex )
    {
      const uint_t neighborEdgeIdx = edgedof::macroedge::indexFromVertex< Level >( it.col(), neighbor );
      p2Edges[ neighborEdgeIdx ] += p1VertexValueHalf;
    }
  }
}

SPECIALIZE_WITH_VALUETYPE(void, prolongateP1ToP2Tmpl, prolongateP1ToP2)

template< typename ValueType, uint_t Level >
inline void restrictP2ToP1Tmpl( const Edge & edge,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & p2VertexDoFMemoryID,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & p2EdgeDoFMemoryID,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & p1VertexDoFMemoryID )
{
  const auto p2Vertices = edge.getData( p2VertexDoFMemoryID )->getPointer( Level );
  const auto p2Edges    = edge.getData( p2EdgeDoFMemoryID   )->getPointer( Level );

        auto p1Vertices = edge.getData( p1VertexDoFMemoryID )->getPointer( Level );

  for ( const auto & it : vertexdof::macroedge::Iterator( Level, 1 ) )
  {
    ValueType tmp;
    const uint_t idx = vertexdof::macroedge::indexFromVertex< Level >( it.col(), stencilDirection::VERTEX_C );
    tmp = p2Vertices[ idx ];

    tmp += 0.5 * p2Edges[ edgedof::macroedge::indexFromVertex< Level >( it.col(), stencilDirection::EDGE_HO_W  ) ];
    tmp += 0.5 * p2Edges[ edgedof::macroedge::indexFromVertex< Level >( it.col(), stencilDirection::EDGE_HO_E )  ];

    tmp += 0.5 * p2Edges[ edgedof::macroedge::indexFromVertex< Level >( it.col(), stencilDirection::EDGE_VE_S  ) ];
    tmp += 0.5 * p2Edges[ edgedof::macroedge::indexFromVertex< Level >( it.col(), stencilDirection::EDGE_DI_SE ) ];

    if ( edge.getNumNeighborFaces() == 2 )
    {
      tmp += 0.5 * p2Edges[ edgedof::macroedge::indexFromVertex< Level >( it.col(), stencilDirection::EDGE_VE_N  ) ];
      tmp += 0.5 * p2Edges[ edgedof::macroedge::indexFromVertex< Level >( it.col(), stencilDirection::EDGE_DI_NW ) ];
    }

    p1Vertices[ idx ] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, restrictP2ToP1Tmpl, restrictP2ToP1)

} // namespace macroedge

namespace macrovertex {

template< typename ValueType, uint_t Level >
inline void prolongateP1ToP2Tmpl( const Vertex & vertex,
                                  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex > & p2VertexDoFMemoryID,
                                  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex > & p2EdgeDoFMemoryID,
                                  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex > & p1VertexDoFMemoryID )
{
  WALBERLA_UNUSED( p2EdgeDoFMemoryID );

        auto p2Vertices = vertex.getData( p2VertexDoFMemoryID )->getPointer( Level );
  const auto p1Vertices = vertex.getData( p1VertexDoFMemoryID )->getPointer( Level );

  p2Vertices[ 0 ] = p1Vertices[ 0 ];
}

SPECIALIZE_WITH_VALUETYPE(void, prolongateP1ToP2Tmpl, prolongateP1ToP2)

template< typename ValueType, uint_t Level >
inline void restrictP2ToP1Tmpl( const Vertex & vertex,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Vertex > & p2VertexDoFMemoryID,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Vertex > & p2EdgeDoFMemoryID,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Vertex > & p1VertexDoFMemoryID )
{
  const auto p2Vertices = vertex.getData( p2VertexDoFMemoryID )->getPointer( Level );
  const auto p2Edges    = vertex.getData( p2EdgeDoFMemoryID   )->getPointer( Level );

        auto p1Vertices = vertex.getData( p1VertexDoFMemoryID )->getPointer( Level );

  p1Vertices[ 0 ] = p2Vertices[ 0 ];

  const uint_t numNeighborEdges = vertex.getNumNeighborEdges();

  for ( uint_t i = 0; i < numNeighborEdges; i++ )
  {
    p1Vertices[ 0 ] += 0.5 * p2Edges[ i ];
  }

  WALBERLA_LOG_DEVEL( p1Vertices[0] );
}

SPECIALIZE_WITH_VALUETYPE(void, restrictP2ToP1Tmpl, restrictP2ToP1)

} // namespace macrovertex


} // namespace P2
}
