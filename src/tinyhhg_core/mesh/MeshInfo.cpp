
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


void MeshInfo:: deriveEdges()
{

  MeshInfo::FaceContainer faces = this->getFaces();
  MeshInfo::VertexContainer verts = this->getVertices();

  DoFType edgeType = Inner;

  for ( const auto & it : faces )
    {
      // extract the three nodes of the face
      std::vector<IDType> fNode = it.second.getVertices();

      // determine their position w.r.t. the boundary
      std::vector<DoFType> dT( 3 );
      dT[0] = verts.find( fNode[0] )->second.getDoFType();
      dT[1] = verts.find( fNode[1] )->second.getDoFType();
      dT[2] = verts.find( fNode[2] )->second.getDoFType();

      // set the three edges of triangle, edge is on boundary, if both
      // its vertices are
      edgeType = ( dT[0] == DirichletBoundary &&
                   dT[1] == DirichletBoundary ) ? DirichletBoundary : Inner;
      this->addEdge( Edge( { fNode[0], fNode[1] }, edgeType ) );

      edgeType = ( dT[0] == DirichletBoundary &&
                   dT[2] == DirichletBoundary ) ? DirichletBoundary : Inner;
      this->addEdge( Edge( { fNode[0], fNode[2] }, edgeType ) );

      edgeType = ( dT[1] == DirichletBoundary &&
                   dT[2] == DirichletBoundary ) ? DirichletBoundary : Inner;
      this->addEdge( Edge( { fNode[1], fNode[2] }, edgeType ) );
    }
}

} // namespace hhg
