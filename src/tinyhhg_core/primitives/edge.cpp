#include "edge.hpp"
#include "vertex.hpp"

#include <core/mpi/MPIManager.h>

#include <fmt/format.h>

namespace hhg
{

Edge::Edge(size_t _id, size_t _type, Vertex* _v0, Vertex* _v1)
  : id(_id), rank(id % walberla::mpi::MPIManager::instance()->numProcesses()), type(_type), v0(_v0), v1(_v1)
{
  direction = v1->coords - v0->coords;
  length = direction.norm();
  tangent = direction / length;
  normal_2d = Point3D({tangent[1], -tangent[0], 0.0});

  // fmt::print("direction = {}\n", direction);
  // fmt::print("length = {}\n", length);
  // fmt::print("tangent = {}\n", tangent);
  // fmt::print("normal_2d = {}\n", normal_2d);
}

void Edge::addFace(Face* face)
{
  faces.push_back(face);
}

size_t Edge::vertex_index(const Vertex& vertex) const
{
  if (&vertex == v0)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

size_t Edge::face_index(const Face& face) const
{
  if (&face == faces[0])
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

std::ostream& operator<<(std::ostream &os, const hhg::Edge &edge)
{
  return os << "Edge { id = " << edge.id << "; "
            << "type = " << edge.type << "; "
            << "v0 = " << edge.v0->id << "; "
            << "v1 = " << edge.v1->id << "; }";
}


void EdgeStencilMemory::free()
{
	for (auto el : data)
	{
		delete[] el.second;
	}
	data.clear();
}

void EdgeP1Memory::free()
{
	for (auto el : data)
	{
		delete[] el.second;
	}
	data.clear();
}

}
