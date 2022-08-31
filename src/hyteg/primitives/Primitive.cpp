/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include <algorithm>
#include <core/mpi/BufferDataTypeExtensions.h>
#include <hyteg/primitives/Primitive.hpp>
#include <hyteg/primitivestorage/PrimitiveStorage.hpp>

namespace hyteg {

bool Primitive::neighborPrimitiveExists( const PrimitiveID& primitiveID ) const
{
   std::vector< PrimitiveID > neighborIDs;
   getNeighborPrimitives( neighborIDs );
   return std::find( neighborIDs.begin(), neighborIDs.end(), primitiveID ) != neighborIDs.end();
}

void Primitive::getNeighborPrimitives( std::vector< PrimitiveID >& neighborPrimitives ) const
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

void Primitive::serialize( walberla::mpi::SendBuffer& sendBuffer ) const
{
   serializePrimitive( sendBuffer );
   serializeSubclass( sendBuffer );
}

void Primitive::deserialize( walberla::mpi::RecvBuffer& recvBuffer )
{
   deserializePrimitive( recvBuffer );
   deserializeSubclass( recvBuffer );
}

void Primitive::serializePrimitive( walberla::mpi::SendBuffer& sendBuffer ) const
{
   sendBuffer << primitiveID_;
   sendBuffer << meshBoundaryFlag_;
   sendBuffer << neighborVertices_;
   sendBuffer << neighborEdges_;
   sendBuffer << neighborFaces_;
   sendBuffer << neighborCells_;
   sendBuffer << childVertices_;
   sendBuffer << childEdges_;
   sendBuffer << childFaces_;
   sendBuffer << childCells_;
   sendBuffer << parent_;
   GeometryMap::serialize( geometryMap_, sendBuffer );
}

void Primitive::deserializePrimitive( walberla::mpi::RecvBuffer& recvBuffer )
{
   recvBuffer >> primitiveID_;
   recvBuffer >> meshBoundaryFlag_;
   recvBuffer >> neighborVertices_;
   recvBuffer >> neighborEdges_;
   recvBuffer >> neighborFaces_;
   recvBuffer >> neighborCells_;
   recvBuffer >> childVertices_;
   recvBuffer >> childEdges_;
   recvBuffer >> childFaces_;
   recvBuffer >> childCells_;
   recvBuffer >> parent_;
   geometryMap_ = GeometryMap::deserialize( recvBuffer );
}

const std::shared_ptr< GeometryMap >& Primitive::getGeometryMap() const
{
   return geometryMap_;
}

bool Primitive::hasChildren() const
{
   return !( childVertices_.empty() && childEdges_.empty() && childFaces_.empty() && childCells_.empty() );
}

std::vector< PrimitiveID > Primitive::childVertices() const
{
   return childVertices_;
}

std::vector< PrimitiveID > Primitive::childEdges() const
{
   return childEdges_;
}

std::vector< PrimitiveID > Primitive::childFaces() const
{
   return childFaces_;
}

std::vector< PrimitiveID > Primitive::childCells() const
{
   return childCells_;
}

void Primitive::addChildVertices( const std::vector< PrimitiveID >& pids )
{
   childVertices_.insert( childVertices_.end(), pids.begin(), pids.end() );
}

void Primitive::addChildEdges( const std::vector< PrimitiveID >& pids )
{
   childEdges_.insert( childEdges_.end(), pids.begin(), pids.end() );
}

void Primitive::addChildFaces( const std::vector< PrimitiveID >& pids )
{
   childFaces_.insert( childFaces_.end(), pids.begin(), pids.end() );
}

void Primitive::addChildCells( const std::vector< PrimitiveID >& pids )
{
   childCells_.insert( childCells_.end(), pids.begin(), pids.end() );
}

void Primitive::removeChildren( const std::vector< PrimitiveID >& pids )
{
   for ( const auto& pid : pids )
   {
      std::remove( childVertices_.begin(), childVertices_.end(), pid );
      std::remove( childEdges_.begin(), childEdges_.end(), pid );
      std::remove( childFaces_.begin(), childFaces_.end(), pid );
      std::remove( childCells_.begin(), childCells_.end(), pid );
   }
}

bool Primitive::hasParent() const
{
   return parent_.has_value();
}

PrimitiveID Primitive::parent() const
{
   WALBERLA_CHECK( hasParent(), "Requested parent that does not exist." );
   return parent_.value();
}

void Primitive::setParent( const PrimitiveID& pid )
{
   parent_ = pid;
}

void Primitive::clearParent()
{
   parent_.reset();
}

} // namespace hyteg
