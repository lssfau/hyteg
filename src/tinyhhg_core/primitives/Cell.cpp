
#include "tinyhhg_core/primitives/Cell.hpp"
#include "core/mpi/BufferDataTypeExtensions.h"

namespace hyteg {

Cell::Cell( const PrimitiveID & primitiveID,
            const std::vector< PrimitiveID > & vertexIDs,
            const std::vector< PrimitiveID > & edgeIDs,
            const std::vector< PrimitiveID > & faceIDs,
            const std::array< Point3D, 4 >   & coordinates,
            const std::array< std::map< uint_t, uint_t >, 6 > & edgeLocalVertexToCellLocalVertexMaps,
            const std::array< std::map< uint_t, uint_t >, 4 > & faceLocalVertexToCellLocalVertexMaps ) :
  Primitive( primitiveID ), coordinates_( coordinates ),
  edgeLocalVertexToCellLocalVertexMaps_( edgeLocalVertexToCellLocalVertexMaps ),
  faceLocalVertexToCellLocalVertexMaps_( faceLocalVertexToCellLocalVertexMaps )
{
  WALBERLA_ASSERT_EQUAL( vertexIDs.size(), 4, "Only tetrahedron cells are supported (number of vertices mismatches)." );
  WALBERLA_ASSERT_EQUAL( edgeIDs.size(),   6, "Only tetrahedron cells are supported (number of edges mismatches)." );
  WALBERLA_ASSERT_EQUAL( faceIDs.size(),   4, "Only tetrahedron cells are supported (number of faces mismatches)." );

  WALBERLA_ASSERT_EQUAL( edgeLocalVertexToCellLocalVertexMaps[0].size(), 2 );
  WALBERLA_ASSERT_EQUAL( edgeLocalVertexToCellLocalVertexMaps[1].size(), 2 );
  WALBERLA_ASSERT_EQUAL( edgeLocalVertexToCellLocalVertexMaps[2].size(), 2 );
  WALBERLA_ASSERT_EQUAL( edgeLocalVertexToCellLocalVertexMaps[3].size(), 2 );
  WALBERLA_ASSERT_EQUAL( edgeLocalVertexToCellLocalVertexMaps[4].size(), 2 );
  WALBERLA_ASSERT_EQUAL( edgeLocalVertexToCellLocalVertexMaps[5].size(), 2 );

  WALBERLA_ASSERT_EQUAL( faceLocalVertexToCellLocalVertexMaps[0].size(), 3 );
  WALBERLA_ASSERT_EQUAL( faceLocalVertexToCellLocalVertexMaps[1].size(), 3 );
  WALBERLA_ASSERT_EQUAL( faceLocalVertexToCellLocalVertexMaps[2].size(), 3 );
  WALBERLA_ASSERT_EQUAL( faceLocalVertexToCellLocalVertexMaps[3].size(), 3 );

  neighborVertices_.assign( vertexIDs.begin(), vertexIDs.end() );
  neighborEdges_.assign   ( edgeIDs.begin(), edgeIDs.end() );
  neighborFaces_.assign   ( faceIDs.begin(), faceIDs.end() );
}

uint_t Cell::getLocalFaceID( const PrimitiveID & faceID ) const
{
  WALBERLA_ASSERT( neighborPrimitiveExists( faceID ) );
  for ( uint_t localFaceID = 0; localFaceID < 4; localFaceID++ )
  {
    if ( neighborFaces_[ localFaceID ] == faceID )
    {
      return localFaceID;
    }
  }
  return std::numeric_limits< uint_t >::max();
}

uint_t Cell::getLocalEdgeID( const PrimitiveID & edgeID ) const
{
  WALBERLA_ASSERT( neighborPrimitiveExists( edgeID ) );
  for ( uint_t localEdgeID = 0; localEdgeID < 6; localEdgeID++ )
  {
    if ( neighborEdges_[ localEdgeID ] == edgeID )
    {
      return localEdgeID;
    }
  }
  return std::numeric_limits< uint_t >::max();
}

uint_t Cell::getLocalVertexID( const PrimitiveID & vertexID ) const
{
  WALBERLA_ASSERT( neighborPrimitiveExists( vertexID ) );
  for ( uint_t localVertexID = 0; localVertexID < getNumNeighborVertices(); localVertexID++ )
  {
    if ( neighborVertices_[ localVertexID ] == vertexID )
    {
      return localVertexID;
    }
  }
  return std::numeric_limits< uint_t >::max();
}

void Cell::serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << coordinates_;
  sendBuffer << edgeLocalVertexToCellLocalVertexMaps_;
  sendBuffer << faceLocalVertexToCellLocalVertexMaps_;
}

void Cell::deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> coordinates_;
  recvBuffer >> edgeLocalVertexToCellLocalVertexMaps_;
  recvBuffer >> faceLocalVertexToCellLocalVertexMaps_;
}

}
