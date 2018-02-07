
#include "tinyhhg_core/primitives/Cell.hpp"

namespace hhg {

Cell::Cell( const PrimitiveID & primitiveID,
            const std::vector< PrimitiveID > & vertexIDs,
            const std::vector< PrimitiveID > & edgeIDs,
            const std::vector< PrimitiveID > & faceIDs,
            const std::array< Point3D, 4 >   & coordinates ) :
  Primitive( primitiveID ), coordinates_( coordinates )
{
  WALBERLA_ASSERT_EQUAL( vertexIDs.size(), 4, "Only tetrahedron cells are supported (number of vertices mismatches)." );
  WALBERLA_ASSERT_EQUAL( edgeIDs.size(),   6, "Only tetrahedron cells are supported (number of edges mismatches)." );
  WALBERLA_ASSERT_EQUAL( faceIDs.size(), 4, "Only tetrahedron cells are supported (number of faces mismatches)." );

  neighborVertices_.assign( vertexIDs.begin(), vertexIDs.end() );
  neighborEdges_.assign   ( edgeIDs.begin(), edgeIDs.end() );
  neighborFaces_.assign   ( faceIDs.begin(), faceIDs.end() );
}

}
