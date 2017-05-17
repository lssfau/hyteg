#include "vertex.hpp"

#include <core/mpi/MPIManager.h>

namespace hhg
{

using walberla::uint_c;

Vertex::Vertex(size_t _id, const Point3D& _coords)
  : id(_id), rank(id % uint_c(walberla::mpi::MPIManager::instance()->numProcesses())), coords(_coords)
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
