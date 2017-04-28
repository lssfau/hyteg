#include "edge.hpp"
#include "face.hpp"
#include "vertex.hpp"
#include "tinyhhg_core/types/flags.hpp"
#include "tinyhhg_core/math.hpp"

#include <core/mpi/MPIManager.h>

#include <fmt/ostream.h>
#include <cstddef>

namespace hhg
{

Face::Face(size_t _id, Edge* _edges[3])
  : id(_id), rank(id % walberla::mpi::MPIManager::instance()->numProcesses()), type(Inner)
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
    edge_orientation = {1, 1, 1};
  }
  else if (v0_1 == v1_0 && v1_1 == v2_1 && v2_0 == v0_0)
  {
    edge_orientation = {1, 1, -1};
  }
  else if (v0_1 == v1_1 && v1_0 == v2_0 && v2_1 == v0_0)
  {
    edge_orientation = {1, -1, 1};
  }
  else if (v0_1 == v1_1 && v1_0 == v2_1 && v2_0 == v0_0)
  {
    edge_orientation = {1, -1, -1};
  }
  else if (v0_0 == v1_0 && v1_1 == v2_0 && v2_1 == v0_1)
  {
    edge_orientation = {-1, 1, 1};
  }
  else if (v0_0 == v1_0 && v1_1 == v2_1 && v2_0 == v0_1)
  {
    edge_orientation = {-1, 1, -1};
  }
  else if (v0_0 == v1_1 && v1_0 == v2_0 && v2_1 == v0_1)
  {
    edge_orientation = {-1, -1, 1};
  }
  else if (v0_0 == v1_1 && v1_0 == v2_1 && v2_0 == v0_1)
  {
    edge_orientation = {-1, -1, -1};
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

  coords = { vertices[0]->coords, vertices[1]->coords, vertices[2]->coords };

  std::array<Point3D, 2> B({coords[1]-coords[0], coords[2] - coords[0]});
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

  return -1;
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

  return -1;
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

std::ostream& operator<<(std::ostream &os, const hhg::Face &face)
{
  return os << "Face { id = " << face.id << "; "
            << "edges[0] = " << face.edges[0]->id << "; "
            << "edges[1] = " << face.edges[1]->id << "; "
            << "edges[2] = " << face.edges[2]->id << "; "
            << "}";
}

}