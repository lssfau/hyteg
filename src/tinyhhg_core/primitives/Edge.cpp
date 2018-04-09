#include "Edge.hpp"
#include "Vertex.hpp"

#include <core/mpi/MPIManager.h>

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

namespace hhg
{

using walberla::uint_c;

uint_t Edge::vertex_index(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 2);

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

  for ( uint_t localFaceID = 0; localFaceID < getNumNeighborFaces(); localFaceID++ )
  {
    if ( face.getID() == neighborFaces_[ localFaceID ].getID() )
    {
      return localFaceID;
    }
  }

  WALBERLA_ASSERT(false, "Edge::face_index: Face does not belong to edge");
  return std::numeric_limits<std::size_t>::max();
}

PrimitiveID Edge::get_opposite_vertex(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 2);

  if (vertex.getID() == neighborVertices_[0].getID())
  {
    return neighborVertices_[1];
  }
  else if (vertex.getID() == neighborVertices_[1].getID())
  {
    return neighborVertices_[0];
  }

  WALBERLA_ABORT("Edge::get_opposite_vertex: Vertex does not belong to edge");
}

bool Edge::opposite_face_exists(const PrimitiveID& face) const
{
  if (face.getID() == neighborFaces_[0].getID()) {
    if (getNumNeighborFaces() == 2) {
      return true;
    } else {
      return false;
    }
  }

  if (getNumNeighborFaces() == 2 && face.getID() == neighborFaces_[1].getID())
  {
    return true;
  }

  WALBERLA_ABORT("Edge::opposite_face_exists: Face does not belong to edge");
}

PrimitiveID Edge::get_opposite_face(const PrimitiveID& face) const
{
  if (face.getID() == neighborFaces_[0].getID()) {
    if (getNumNeighborFaces() == 2) {
      return neighborFaces_[1];
    } else {
      WALBERLA_ABORT("Edge::get_opposite_face: Requesting face that does not exist");
    }
  }

  if (getNumNeighborFaces() == 2 && face.getID() == neighborFaces_[1].getID())
  {
    return neighborFaces_[0];
  }

  WALBERLA_ABORT("Edge::get_opposite_face: Face does not belong to edge");
}

bool Edge::onBoundary() const {
  return testFlag(dofType_, hhg::Boundary);
}

std::ostream& operator<<(std::ostream &os, const hhg::Edge &edge)
{
  return os << "Edge { id = " << edge.getID().getID() << "; "
            << "type = " << edge.dofType_ << "; "
            << "neighborVertices_[0] = " << edge.neighborVertices_[0].getID() << "; "
            << "neighborVertices_[1] = " << edge.neighborVertices_[1].getID() << "; }";
}

void Edge::serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << dofType_;
  sendBuffer << coordinates_[0];
  sendBuffer << coordinates_[1];
  sendBuffer << direction_;
  sendBuffer << length_;
  sendBuffer << tangent_;
  sendBuffer << normal2D_;
}

void Edge::deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> dofType_;
  recvBuffer >> coordinates_[0];
  recvBuffer >> coordinates_[1];
  recvBuffer >> direction_;
  recvBuffer >> length_;
  recvBuffer >> tangent_;
  recvBuffer >> normal2D_;
}

}
