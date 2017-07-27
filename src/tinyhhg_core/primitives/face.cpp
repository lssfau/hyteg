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

uint_t Face::vertex_index(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 3);

  for (size_t i = 0; i < 3; ++i)
  {
    if (vertex == neighborVertices_[i].getID())
    {
      return i;
    }
  }

  WALBERLA_ASSERT(false, "Face::vertex_index: Vertex does not belong to face")
  return std::numeric_limits<std::size_t>::max();
}

uint_t Face::edge_index(const PrimitiveID& edge) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborEdges(), 3);

  for (size_t i = 0; i < 3; ++i)
  {
    if (edge.getID() == neighborEdges_[i].getID())
    {
      return i;
    }
  }

  WALBERLA_ASSERT(false, "Face::edge_index: Edge does not belong to face")
  return std::numeric_limits<std::size_t>::max();
}

std::vector<PrimitiveID> Face::adjacent_edges(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborEdges(), 3);

  std::vector<PrimitiveID> e;

  if (vertex_index(vertex) == 0) {
    e.push_back(neighborEdges_[0]);
    e.push_back(neighborEdges_[2]);
  } else if (vertex_index(vertex) == 1) {
    e.push_back(neighborEdges_[0]);
    e.push_back(neighborEdges_[1]);
  } else if (vertex_index(vertex) == 2) {
    e.push_back(neighborEdges_[1]);
    e.push_back(neighborEdges_[2]);
  }

  WALBERLA_ASSERT_EQUAL(e.size(), 2);
  return e;
}

PrimitiveID Face::get_vertex_opposite_to_edge(const PrimitiveID& edge) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 3);
  WALBERLA_ASSERT_EQUAL(getNumNeighborEdges(), 3);

  if (edge_index(edge) == 0) {
    return neighborVertices_[2];
  } else if (edge_index(edge) == 1) {
    return neighborVertices_[0];
  } else if (edge_index(edge) == 2) {
    return neighborVertices_[1];
  } else {
    WALBERLA_ASSERT(false, "Face::get_vertex_opposite_to_edge: Edge does not belong to face");
  }
}

std::ostream& operator<<(std::ostream &os, const hhg::Face &face)
{
  return os << "Face { id = " << face.getID().getID() << "; "
            << "neighborEdges_[0] = " << face.neighborEdges_[0].getID() << "; "
            << "neighborEdges_[1] = " << face.neighborEdges_[1].getID() << "; "
            << "neighborEdges_[2] = " << face.neighborEdges_[2].getID() << "; "
            << "}";
}

}
