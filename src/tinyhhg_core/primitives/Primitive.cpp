
#include <tinyhhg_core/primitives/Primitive.hpp>
#include <tinyhhg_core/primitivestorage/PrimitiveStorage.hpp>

#include <core/mpi/BufferDataTypeExtensions.h>

namespace hhg {

void Primitive::getNeighborPrimitives( std::vector< PrimitiveID > & neighborPrimitives ) const
{
  getNeighborVertices( neighborPrimitives );

  std::vector< PrimitiveID > neighborEdges;
  std::vector< PrimitiveID > neighborFaces;

  getNeighborEdges( neighborEdges );
  getNeighborFaces( neighborFaces );

  neighborPrimitives.insert( neighborPrimitives.end(), neighborEdges.begin(), neighborEdges.end() );
  neighborPrimitives.insert( neighborPrimitives.end(), neighborFaces.begin(), neighborFaces.end() );

}

void Primitive::serialize ( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << primitiveID_;
  sendBuffer << neighborVertices_;
  sendBuffer << neighborEdges_;
  sendBuffer << neighborFaces_;
}

void Primitive::deserialize ( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> primitiveID_;
  recvBuffer >> neighborVertices_;
  recvBuffer >> neighborEdges_;
  recvBuffer >> neighborFaces_;
}


} // namespace hhg

