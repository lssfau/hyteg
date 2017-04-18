#include "vertex.hpp"

#include <core/mpi/MPIManager.h>

namespace hhg
{

Vertex::Vertex(size_t _id, const Point3D& _coords)
  : id(_id), rank(id % walberla::mpi::MPIManager::instance()->numProcesses()), coords(_coords)
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
  size_t i = 1;

  for (Edge* e : edges)
  {
    if (&edge == e)
    {
      return i;
    }
    ++i;
  }
}

std::ostream& operator<<(std::ostream &os, const hhg::Vertex &vertex)
{
  return os << "Vertex { id = " << vertex.id << "; "
            << "coords = [" << vertex.coords[0] << ", " << vertex.coords[1] << ", " << vertex.coords[2] << "]; "
            << "}";
}

}
