#include "Vertex.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include <core/mpi/MPIManager.h>

namespace hyteg {

using walberla::uint_c;

Vertex::Vertex( const PrimitiveID & primitiveID, const Point3D & coordinates )
  : Primitive( primitiveID ), dofType_(Inner), coordinates_( coordinates )
{}

uint_t Vertex::edge_index(const PrimitiveID& edge) const
{
  auto it = std::find(neighborEdges().begin(), neighborEdges().end(), edge);
  return uint_c((it - neighborEdges().begin()));
}

uint_t Vertex::face_index(const PrimitiveID& face) const
{
  auto it = std::find(neighborFaces().begin(), neighborFaces().end(), face);
  return uint_c((it - neighborFaces().begin()));
}

uint_t Vertex::cell_index(const PrimitiveID& cell) const
{
  auto it = std::find(neighborCells().begin(), neighborCells().end(), cell);
  return uint_c((it - neighborCells().begin()));
}

std::ostream& operator<<(std::ostream &os, const hyteg::Vertex &vertex)
{
  return os << "Vertex { id = " << vertex.getID().getID() << "; "
            << "coords = [" << vertex.coordinates_[0] << ", " << vertex.coordinates_[1] << ", " << vertex.coordinates_[2] << "]; "
            << "}";
}

void Vertex::serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << dofType_;
  sendBuffer << coordinates_;
}

void Vertex::deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> dofType_;
  recvBuffer >> coordinates_;
}

}
