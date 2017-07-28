#include "edge.hpp"
#include "vertex.hpp"

#include <fmt/format.h>

#include <core/mpi/MPIManager.h>

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

namespace hhg
{

using walberla::uint_c;

uint_t Edge::vertex_index(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 2)

  if (vertex.getID() == neighborVertices_[0].getID())
  {
    return 0;
  }
  else if (vertex.getID() == neighborVertices_[1].getID())
  {
    return 1;
  }
  else
  {
    WALBERLA_ASSERT(false, "Edge::vertex_index: Vertex does not belong to edge");
    return std::numeric_limits<std::size_t>::max();
  }
}

uint_t Edge::face_index(const PrimitiveID& face) const
{
  WALBERLA_ASSERT_GREATER_EQUAL(getNumNeighborFaces(), 1);
  WALBERLA_ASSERT_LESS_EQUAL(getNumNeighborFaces(), 2);

  if (face.getID() == neighborFaces_[0].getID()) {
    return 0;
  } else if(getNumNeighborFaces() == 2 && face.getID() == neighborFaces_[1].getID()) {
    return 1;
  } else {
    WALBERLA_ASSERT(false, "Edge::face_index: Face does not belong to edge");
    return std::numeric_limits<uint_t>::max();
  }
}

PrimitiveID Edge::get_opposite_vertex(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 2)

  if (vertex.getID() == neighborVertices_[0].getID())
  {
    return neighborVertices_[0];
  }
  else if (vertex.getID() == neighborVertices_[1].getID())
  {
    return neighborVertices_[1];
  }

  WALBERLA_ABORT("Edge::get_opposite_vertex: Vertex does not belong to edge");
}

std::ostream& operator<<(std::ostream &os, const hhg::Edge &edge)
{
  return os << "Edge { id = " << edge.getID().getID() << "; "
            << "type = " << edge.type << "; "
            << "neighborVertices_[0] = " << edge.neighborVertices_[0].getID() << "; "
            << "neighborVertices_[1] = " << edge.neighborVertices_[1].getID() << "; }";
}

}
