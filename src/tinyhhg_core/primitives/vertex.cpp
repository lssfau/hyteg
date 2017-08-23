#include "vertex.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include <core/mpi/MPIManager.h>

namespace hhg
{

using walberla::uint_c;

Vertex::Vertex( const PrimitiveID & primitiveID, const Point3D & coordinates )
  : Primitive( primitiveID ), dofType_(Inner), coords( coordinates )
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

std::ostream& operator<<(std::ostream &os, const hhg::Vertex &vertex)
{
  return os << "Vertex { id = " << vertex.getID().getID() << "; "
            << "coords = [" << vertex.coords[0] << ", " << vertex.coords[1] << ", " << vertex.coords[2] << "]; "
            << "}";
}

void Vertex::serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << dofType_;
  sendBuffer << coords;
}

void Vertex::deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> dofType_;
  recvBuffer >> coords;
}

}
