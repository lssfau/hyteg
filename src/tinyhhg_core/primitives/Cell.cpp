
#include "tinyhhg_core/primitives/Cell.hpp"

namespace hhg {

Cell::Cell( const PrimitiveID & primitiveID,
            const std::vector< PrimitiveID > & vertexIDs,
            const std::vector< PrimitiveID > & edgeIDs,
            const std::vector< PrimitiveID > & faceIDs,
            const std::array< Point3D, 4 >   & coordinates,
            const std::array< std::map< uint_t, uint_t >, 4 > & cellLocalVertexToFaceLocalVertexMaps ) :
  Primitive( primitiveID ), coordinates_( coordinates ), faceLocalVertexToCellLocalVertexMaps_( cellLocalVertexToFaceLocalVertexMaps )
{
  WALBERLA_ASSERT_EQUAL( vertexIDs.size(), 4, "Only tetrahedron cells are supported (number of vertices mismatches)." );
  WALBERLA_ASSERT_EQUAL( edgeIDs.size(),   6, "Only tetrahedron cells are supported (number of edges mismatches)." );
  WALBERLA_ASSERT_EQUAL( faceIDs.size(),   4, "Only tetrahedron cells are supported (number of faces mismatches)." );

  WALBERLA_ASSERT_EQUAL( cellLocalVertexToFaceLocalVertexMaps[0].size(), 3 );
  WALBERLA_ASSERT_EQUAL( cellLocalVertexToFaceLocalVertexMaps[1].size(), 3 );
  WALBERLA_ASSERT_EQUAL( cellLocalVertexToFaceLocalVertexMaps[2].size(), 3 );
  WALBERLA_ASSERT_EQUAL( cellLocalVertexToFaceLocalVertexMaps[3].size(), 3 );

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

}
