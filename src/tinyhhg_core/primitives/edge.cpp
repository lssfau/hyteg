#include "edge.hpp"
#include "vertex.hpp"

#include <fmt/format.h>

#include <core/mpi/MPIManager.h>

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

namespace hhg
{

using walberla::uint_c;

Edge::Edge(size_t _id, DoFType _type, Vertex* _v0, Vertex* _v1)
  : Primitive( PrimitiveID( _id ) ),
	  id(_id), rank(id % uint_c(walberla::mpi::MPIManager::instance()->numProcesses())), type(_type), v0(_v0), v1(_v1)
{

  coords = {{v0->coords, v1->coords}};
  direction = v1->coords - v0->coords;
  length = direction.norm();
  tangent = direction / length;
  const std::array<walberla::real_t,3> init{{tangent[1], -tangent[0], 0.0}};
  normal_2d = Point3D(init);

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
  else if (&vertex == v1)
  {
    return 1;
  }
  else
  {
    return 2;
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

Vertex* Edge::get_opposite_vertex(const Vertex& vertex) const
{
  if (&vertex == v0)
  {
    return v1;
  }
  else if (&vertex == v1)
  {
    return v0;
  }
  else
  {
    return NULL;
  }
}

std::ostream& operator<<(std::ostream &os, const hhg::Edge &edge)
{
  return os << "Edge { id = " << edge.id << "; "
            << "type = " << edge.type << "; "
            << "v0 = " << edge.v0->id << "; "
            << "v1 = " << edge.v1->id << "; }";
}

}
