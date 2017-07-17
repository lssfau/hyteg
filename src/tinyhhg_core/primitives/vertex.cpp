#include "vertex.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include <core/mpi/MPIManager.h>

namespace hhg
{

using walberla::uint_c;

Vertex::Vertex(size_t _id, const Point3D& _coords)
  : Primitive( PrimitiveStorage(0, SetupPrimitiveStorage( MeshInfo::emptyMeshInfo(), uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ))),
				PrimitiveID( _id ) ),
				id(_id), rank(id % uint_c(walberla::mpi::MPIManager::instance()->numProcesses())),
				type(Inner), coords(_coords)
{
}

Vertex::Vertex( PrimitiveStorage & storage, const SetupPrimitiveStorage & setupStorage, const PrimitiveID & primitiveID ) :
  Primitive( storage, primitiveID ), id( primitiveID.getID() ),
  rank( setupStorage.getTargetRank( primitiveID ) ),
  type( Inner ), coords( setupStorage.getVertex( primitiveID )->getCoordinates() )
{
  const SetupVertex * setupVertex = setupStorage.getVertex( primitiveID );

  WALBERLA_ASSERT_EQUAL( setupVertex->getNumLowerDimNeighbors(), 0 );

  for ( auto higherDimNeighbor  = setupVertex->beginHigherDimNeighbors();
		         higherDimNeighbor != setupVertex->endHigherDimNeighbors();
		         higherDimNeighbor++ )
  {
	  WALBERLA_ASSERT( setupStorage.edgeExists( *higherDimNeighbor ) );
	  higherDimNeighbors_[ higherDimNeighbor->getID() ] = setupStorage.getTargetRank( *higherDimNeighbor );
  }
}

void Vertex::addEdge(Edge* edge)
{
  edges.push_back(edge);
}

void Vertex::addFace(Face* face)
{
  faces.push_back(face);
}

size_t Vertex::edge_index(const Edge& edge) const
{
  auto it = std::find(edges.begin(),edges.end(),&edge);
  return uint_c((it - edges.begin()));
}

size_t Vertex::face_index(const Face& face) const
{
  auto it = std::find(faces.begin(), faces.end(), &face);
  return uint_c((it - faces.begin()));
}

std::ostream& operator<<(std::ostream &os, const hhg::Vertex &vertex)
{
  return os << "Vertex { id = " << vertex.id << "; "
            << "coords = [" << vertex.coords[0] << ", " << vertex.coords[1] << ", " << vertex.coords[2] << "]; "
            << "}";
}

}
