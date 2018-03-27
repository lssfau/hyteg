
#include "tinyhhg_core/mesh/MeshInfo.hpp"

#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"

#include <array>
#include <vector>

namespace hhg {


void MeshInfo::addEdge( const Edge & edge )
{
  WALBERLA_CHECK_UNEQUAL( edge.getVertices()[0], edge.getVertices()[1], "[Mesh] Mesh contains edge with zero length." );

  std::array< IDType, 2 > sortedVertexIDs;

  if ( edge.getVertices()[0] < edge.getVertices()[1] )
  {
    sortedVertexIDs = edge.getVertices();
  }
  else if ( edge.getVertices()[1] < edge.getVertices()[0] )
  {
    sortedVertexIDs[0] = edge.getVertices()[1];
    sortedVertexIDs[1] = edge.getVertices()[0];
  }

  if (edges_.count( sortedVertexIDs ) == 0)
  {
    edges_[ sortedVertexIDs ] = Edge( sortedVertexIDs, edge.getDoFType() );
  }
}


void MeshInfo::addFace( const Face & face )
{
  std::vector< IDType > sortedVertexIDs = face.getVertices();
  std::sort( sortedVertexIDs.begin(), sortedVertexIDs.end() );

  WALBERLA_CHECK_EQUAL( std::set< IDType >( sortedVertexIDs.begin(), sortedVertexIDs.end() ).size(), sortedVertexIDs.size(),
                        "[Mesh] Mesh contains face with duplicate vertices." );

  if ( faces_.count( sortedVertexIDs ) == 0 )
  {
    faces_[ sortedVertexIDs ] = Face( sortedVertexIDs, face.getDoFType() );
  }

}

} // namespace hhg
