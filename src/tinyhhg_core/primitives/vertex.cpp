#include "vertex.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include <core/mpi/MPIManager.h>

namespace hhg
{

using walberla::uint_c;

Vertex::Vertex(size_t _id, const Point3D& _coords)
  : Primitive( PrimitiveStorage(0, SetupPrimitiveStorage( MeshInfo::emptyMeshInfo(), uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ))),
				SetupVertex(SetupPrimitiveStorage( MeshInfo::emptyMeshInfo(), uint_c( walberla::mpi::MPIManager::instance()->numProcesses() )), _id, Point3D()) ),
				id(_id), rank(id % uint_c(walberla::mpi::MPIManager::instance()->numProcesses())),
				type(Inner), coords(_coords)
{
}

Vertex::Vertex( PrimitiveStorage & storage, const SetupVertex & setupVertex ) :
  Primitive( storage, setupVertex ), id( setupVertex.getPrimitiveID().getID() ),
  rank( setupVertex.getPrimitiveID().getID() % uint_c(walberla::mpi::MPIManager::instance()->numProcesses()) ),
  type( Inner ), coords( setupVertex.getCoordinates() )
{
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

std::ostream& operator<<(std::ostream &os, const hhg::Vertex &vertex)
{
  return os << "Vertex { id = " << vertex.id << "; "
            << "coords = [" << vertex.coords[0] << ", " << vertex.coords[1] << ", " << vertex.coords[2] << "]; "
            << "}";
}

}
