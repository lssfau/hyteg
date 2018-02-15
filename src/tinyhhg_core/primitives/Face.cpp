#include "Face.hpp"
#include "Edge.hpp"
#include "Vertex.hpp"
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

  WALBERLA_ASSERT(false, "Face::vertex_index: Vertex does not belong to face");
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

  WALBERLA_ASSERT(false, "Face::edge_index: Edge does not belong to face");
  return std::numeric_limits<std::size_t>::max();
}

uint_t Face::cell_index( const PrimitiveID & cell ) const
{
  WALBERLA_ASSERT_LESS_EQUAL( getNumNeighborCells(), 2 );

  for (size_t i = 0; i < 2; ++i)
  {
    if ( cell.getID() == neighborCells_[i].getID() )
    {
      return i;
    }
  }

  WALBERLA_ASSERT(false, "Face::cell_index: Cell does not belong to face");
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
  }

  WALBERLA_ABORT("Face::get_vertex_opposite_to_edge: Edge does not belong to face");
}

PrimitiveID Face::get_edge_between_vertices(const PrimitiveID& v0, const PrimitiveID& v1) const
{
  std::vector<PrimitiveID> edges_v0 = adjacent_edges(v0);
  std::vector<PrimitiveID> edges_v1 = adjacent_edges(v1);

  if (edges_v0[0] == edges_v1[0]) {
    return edges_v0[0];
  }

  if (edges_v0[0] == edges_v1[1]) {
    return edges_v0[0];
  }

  if (edges_v0[1] == edges_v1[0]) {
    return edges_v0[1];
  }

  if (edges_v0[1] == edges_v1[1]) {
    return edges_v0[1];
  }

  WALBERLA_ABORT("Face::get_edge_between_vertices: Vertex v1 does not belong to face");
}

std::ostream& operator<<(std::ostream &os, const hhg::Face &face)
{
  return os << "Face { id = " << face.getID().getID() << "; "
            << "neighborEdges_[0] = " << face.neighborEdges_[0].getID() << "; "
            << "neighborEdges_[1] = " << face.neighborEdges_[1].getID() << "; "
            << "neighborEdges_[2] = " << face.neighborEdges_[2].getID() << "; "
            << "}";
}

void Face::serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << type;
  sendBuffer << area;
  sendBuffer << edge_orientation[0];
  sendBuffer << edge_orientation[1];
  sendBuffer << edge_orientation[2];
  sendBuffer << coords[0];
  sendBuffer << coords[1];
  sendBuffer << coords[2];
}

void Face::deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> type;
  recvBuffer >> area;
  recvBuffer >> edge_orientation[0];
  recvBuffer >> edge_orientation[1];
  recvBuffer >> edge_orientation[2];
  recvBuffer >> coords[0];
  recvBuffer >> coords[1];
  recvBuffer >> coords[2];
}

}



