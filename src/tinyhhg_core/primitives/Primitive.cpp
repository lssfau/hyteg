
#include <tinyhhg_core/primitives/Primitive.hpp>
#include <tinyhhg_core/primitivestorage/PrimitiveStorage.hpp>

#include <core/mpi/BufferDataTypeExtensions.h>

namespace hhg {

void Primitive::getNeighborPrimitives( std::vector< PrimitiveID > & neighborPrimitives ) const
{
  getNeighborVertices( neighborPrimitives );

  std::vector< PrimitiveID > someNeighbors;

  getNeighborEdges( someNeighbors );
  neighborPrimitives.insert( neighborPrimitives.end(), someNeighbors.begin(), someNeighbors.end() );

  getNeighborFaces( someNeighbors );
  neighborPrimitives.insert( neighborPrimitives.end(), someNeighbors.begin(), someNeighbors.end() );
}

void Primitive::serialize( walberla::mpi::SendBuffer & sendBuffer ) const
{
  serializePrimitive( sendBuffer );
  serializeSubclass( sendBuffer );
}

void Primitive::deserialize( walberla::mpi::RecvBuffer & recvBuffer )
{
  deserializePrimitive( recvBuffer );
  deserializeSubclass( recvBuffer );
}

void Primitive::serializePrimitive( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << primitiveID_;
  sendBuffer << neighborVertices_;
  sendBuffer << neighborEdges_;
  sendBuffer << neighborFaces_;
}

void Primitive::deserializePrimitive( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> primitiveID_;
  recvBuffer >> neighborVertices_;
  recvBuffer >> neighborEdges_;
  recvBuffer >> neighborFaces_;
}

} // namespace hhg

