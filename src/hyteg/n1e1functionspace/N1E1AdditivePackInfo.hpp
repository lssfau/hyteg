/*
 * Copyright (c) 2022 Daniel Bauer.
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
#pragma once

#include "core/DataTypes.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/communication/DoFSpacePackInfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFAdditivePackInfo.hpp"

namespace hyteg {

using walberla::uint_t;

template < typename ValueType >
class FunctionMemory;
template < typename DataType, typename PrimitiveType >
class PrimitiveDataID;

class Vertex;
class Edge;
class Face;
class Cell;
class PrimitiveStorage;
class PrimitiveID;

namespace n1e1 {

template < typename ValueType >
class N1E1AdditivePackInfo : public communication::DoFSpacePackInfo< ValueType >
{
 public:
   N1E1AdditivePackInfo( uint_t                                                 level,
                         PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                         PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                         PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                         PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                         std::weak_ptr< PrimitiveStorage >                      storage );

   void packVertexForEdge( const Vertex* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;
   void unpackEdgeFromVertex( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;
   void communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const override;

   void packEdgeForVertex( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;
   void unpackVertexFromEdge( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;
   void communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const override;

   void packEdgeForFace( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;
   void unpackFaceFromEdge( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;
   void communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const override;

   void packFaceForEdge( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;
   void unpackEdgeFromFace( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;
   void communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const override;

   void packFaceForVertex( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;
   void unpackVertexFromFace( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;
   void communicateLocalFaceToVertex( const Face* sender, Vertex* receiver ) const override;

   void packCellForFace( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;
   inline void unpackFaceFromCell( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override
   {
      dofPackInfo_.unpackFaceFromCell( receiver, sender, buffer );
   }
   void communicateLocalCellToFace( const Cell* sender, Face* receiver ) const override;

   void packCellForEdge( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;
   inline void unpackEdgeFromCell( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override
   {
      dofPackInfo_.unpackEdgeFromCell( receiver, sender, buffer );
   }
   void communicateLocalCellToEdge( const Cell* sender, Edge* receiver ) const override;

   void packCellForVertex( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;
   void unpackVertexFromCell( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;
   void communicateLocalCellToVertex( const Cell* sender, Vertex* receiver ) const override;

 private:
   using communication::DoFSpacePackInfo< ValueType >::level_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDVertex_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDEdge_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDFace_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDCell_;
   using communication::DoFSpacePackInfo< ValueType >::storage_;

   EdgeDoFAdditivePackInfo< ValueType > dofPackInfo_;
};

} // namespace n1e1
} //namespace hyteg
