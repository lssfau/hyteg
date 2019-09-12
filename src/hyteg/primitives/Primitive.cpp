/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <hyteg/primitives/Primitive.hpp>
#include <hyteg/primitivestorage/PrimitiveStorage.hpp>

#include <core/mpi/BufferDataTypeExtensions.h>

#include <algorithm>

namespace hyteg {

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
  sendBuffer << meshBoundaryFlag_;
  sendBuffer << neighborVertices_;
  sendBuffer << neighborEdges_;
  sendBuffer << neighborFaces_;
  sendBuffer << neighborCells_;
  GeometryMap::serialize(geometryMap_, sendBuffer);
}

void Primitive::deserializePrimitive( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> primitiveID_;
  recvBuffer >> meshBoundaryFlag_;
  recvBuffer >> neighborVertices_;
  recvBuffer >> neighborEdges_;
  recvBuffer >> neighborFaces_;
  recvBuffer >> neighborCells_;
  geometryMap_ = GeometryMap::deserialize(recvBuffer);
}

const std::shared_ptr<GeometryMap>& Primitive::getGeometryMap() const {
  return geometryMap_;
}

} // namespace hyteg

