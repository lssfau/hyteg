#include <fmt/ostream.h>

#include "face.hpp"
#include "edge.hpp"
#include "vertex.hpp"
#include "tinyhhg_core/types/flags.hpp"
#include "tinyhhg_core/math.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include <core/mpi/MPIManager.h>

#include <cstddef>

namespace hhg
{
using walberla::uint_c;

Face::Face(size_t _id, Edge* _edges[3])
  : Primitive( PrimitiveStorage( 0, SetupPrimitiveStorage( MeshInfo::emptyMeshInfo(), uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ))),
	       SetupFace(_id, 0, 0, 0, std::array< int, 3 >(), std::array< Point3D, 3 >()) ), id(_id), rank(id % uint_c(walberla::mpi::MPIManager::instance()->numProcesses())), type(Inner)
{
  for (size_t i=0; i < 3; ++i)
  {
    edges[i] = _edges[i];
  }

  Vertex* v0_0 = edges[0]->v0;
  Vertex* v0_1 = edges[0]->v1;
  Vertex* v1_0 = edges[1]->v0;
  Vertex* v1_1 = edges[1]->v1;
  Vertex* v2_0 = edges[2]->v0;
  Vertex* v2_1 = edges[2]->v1;

  if (v0_1 == v1_0 && v1_1 == v2_0 && v2_1 == v0_0)
  {
    edge_orientation = {{1, 1, 1}};
  }
  else if (v0_1 == v1_0 && v1_1 == v2_1 && v2_0 == v0_0)
  {
    edge_orientation = {{1, 1, -1}};
  }
  else if (v0_1 == v1_1 && v1_0 == v2_0 && v2_1 == v0_0)
  {
    edge_orientation = {{1, -1, 1}};
  }
  else if (v0_1 == v1_1 && v1_0 == v2_1 && v2_0 == v0_0)
  {
    edge_orientation = {{1, -1, -1}};
  }
  else if (v0_0 == v1_0 && v1_1 == v2_0 && v2_1 == v0_1)
  {
    edge_orientation = {{-1, 1, 1}};
  }
  else if (v0_0 == v1_0 && v1_1 == v2_1 && v2_0 == v0_1)
  {
    edge_orientation = {{-1, 1, -1}};
  }
  else if (v0_0 == v1_1 && v1_0 == v2_0 && v2_1 == v0_1)
  {
    edge_orientation = {{-1, -1, 1}};
  }
  else if (v0_0 == v1_1 && v1_0 == v2_1 && v2_0 == v0_1)
  {
    edge_orientation = {{-1, -1, -1}};
  }

  if (edge_orientation[0] == 1)
  {
    vertices.push_back(v0_0);
    vertices.push_back(v0_1);
  }
  else
  {
    vertices.push_back(v0_1);
    vertices.push_back(v0_0);
  }

  if (edge_orientation[1] == 1)
  {
    vertices.push_back(v1_1);
  }
  else
  {
    vertices.push_back(v1_0);
  }

  coords = {{vertices[0]->coords, vertices[1]->coords, vertices[2]->coords}};

  std::array<Point3D, 2> B({{coords[1]-coords[0], coords[2] - coords[0]}});
  area = std::abs(0.5 * math::det2(B));
}

Face::Face( PrimitiveStorage & storage, const SetupFace & setupFace )
  : Primitive( storage, setupFace ), id( setupFace.getPrimitiveID().getID() ),
    rank(setupFace.getPrimitiveID().getID() % uint_c(walberla::mpi::MPIManager::instance()->numProcesses())),
    type(Inner)
{
  edges[0] = storage.getEdge( setupFace.getEdgeID0() );
  edges[1] = storage.getEdge( setupFace.getEdgeID1() );
  edges[2] = storage.getEdge( setupFace.getEdgeID2() );

  edge_orientation = setupFace.getEdgeOrientation();
  coords = setupFace.getCoordinates();

  std::array<Point3D, 2> B({{coords[1]-coords[0], coords[2] - coords[0]}});
  area = std::abs(0.5 * math::det2(B));
}

size_t Face::vertex_index(const Vertex& vertex) const
{
  for (size_t i = 0; i < 3; ++i)
  {
    if (&vertex == vertices[i])
    {
      return i;
    }
  }

  return std::numeric_limits<std::size_t>::max();
}

size_t Face::edge_index(const Edge& edge) const
{
  for (size_t i = 0; i < 3; ++i)
  {
    if (&edge == edges[i])
    {
      return i;
    }
  }

  return std::numeric_limits<std::size_t>::max();
}

std::vector<Edge*> Face::adjacent_edges(const Vertex& vertex) const
{
  std::vector<Edge*> e;

  for (size_t i = 0; i < 3; ++i)
  {
    if (edges[i]->vertex_index(vertex) != 2)
    {
      e.push_back(edges[i]);
    }
  }

  return e;
}

Vertex* Face::get_vertex_opposite_to_edge(const Edge& edge) const
{
  std::vector<Vertex*> v(vertices);

  std::remove(v.begin(), v.end(), edge.v0);
  std::remove(v.begin(), v.end(), edge.v1);

  return v[0];
}

std::ostream& operator<<(std::ostream &os, const hhg::Face &face)
{
  return os << "Face { id = " << face.id << "; "
            << "edges[0] = " << face.edges[0]->id << "; "
            << "edges[1] = " << face.edges[1]->id << "; "
            << "edges[2] = " << face.edges[2]->id << "; "
            << "}";
}

}
