
#include <tinyhhg_core/primitives/Primitive.hpp>
#include <tinyhhg_core/primitivestorage/PrimitiveStorage.hpp>

#include <core/mpi/BufferDataTypeExtensions.h>

#include <algorithm>

namespace hhg {

bool Primitive::neighborPrimitiveExists( const PrimitiveID & primitiveID ) const
{
  std::vector< PrimitiveID > neighborIDs;
  getNeighborPrimitives( neighborIDs );
  return std::find( neighborIDs.begin(), neighborIDs.end(), primitiveID ) != neighborIDs.end();
}

void Primitive::getNeighborPrimitives( std::vector< PrimitiveID > & neighborPrimitives ) const
{
  getNeighborVertices( neighborPrimitives );

  std::vector< PrimitiveID > someNeighbors;

  getNeighborEdges( someNeighbors );
  neighborPrimitives.insert( neighborPrimitives.end(), someNeighbors.begin(), someNeighbors.end() );

  getNeighborFaces( someNeighbors );
  neighborPrimitives.insert( neighborPrimitives.end(), someNeighbors.begin(), someNeighbors.end() );

  getNeighborCells( someNeighbors );
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
  sendBuffer << neighborCells_;
  geometryMap_->serialize(sendBuffer);
}

void Primitive::deserializePrimitive( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> primitiveID_;
  recvBuffer >> neighborVertices_;
  recvBuffer >> neighborEdges_;
  recvBuffer >> neighborFaces_;
  recvBuffer >> neighborCells_;
  geometryMap_ = GeometryMap::deserialize(recvBuffer);
}

void Primitive::setGeometryMap(const std::shared_ptr<GeometryMap>& newMap) {
  geometryMap_ = newMap;
}

const std::shared_ptr<GeometryMap>& Primitive::getGeometryMap() const {
  return geometryMap_;
}

} // namespace hhg

