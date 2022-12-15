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
#include "hyteg/edgedofspace/EdgeDoFPackInfo.hpp"

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
class N1E1PackInfo : public communication::DoFSpacePackInfo< ValueType >
{
 public:
   N1E1PackInfo( uint_t                                                 level,
                 PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                 PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                 PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                 PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                 std::weak_ptr< PrimitiveStorage >                      storage );

   inline void
       packVertexForEdge( const Vertex* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override
   {
      dofPackInfo_.packVertexForEdge( sender, receiver, buffer );
   }
   inline void unpackEdgeFromVertex( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override
   {
      dofPackInfo_.unpackEdgeFromVertex( receiver, sender, buffer );
   }
   inline void communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const override
   {
      dofPackInfo_.communicateLocalVertexToEdge( sender, receiver );
   }

   inline void
       packEdgeForVertex( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override
   {
      dofPackInfo_.packEdgeForVertex( sender, receiver, buffer );
   }
   inline void
       unpackVertexFromEdge( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override
   {
      dofPackInfo_.unpackVertexFromEdge( receiver, sender, buffer );
   }
   inline void communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const override
   {
      dofPackInfo_.communicateLocalEdgeToVertex( sender, receiver );
   }

   inline void
       packEdgeForFace( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override
   {
      dofPackInfo_.packEdgeForFace( sender, receiver, buffer );
   }
   inline void unpackFaceFromEdge( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override
   {
      dofPackInfo_.unpackFaceFromEdge( receiver, sender, buffer );
   }
   inline void communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const override
   {
      dofPackInfo_.communicateLocalEdgeToFace( sender, receiver );
   }

   inline void packFaceForEdge( const Face*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const override
   {
      WALBERLA_ABORT( "Macro-face to macro-edge communication not implemented!" );
   }
   inline void unpackEdgeFromFace( Edge*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const override
   {
      WALBERLA_ABORT( "Macro-face to macro-edge communication not implemented!" );
   }
   inline void communicateLocalFaceToEdge( const Face*, Edge* ) const override
   {
      WALBERLA_ABORT( "Macro-face to macro-edge communication not implemented!" );
   }

   inline void
       packFaceForCell( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override
   {
      dofPackInfo_.packFaceForCell( sender, receiver, buffer );
   }
   void unpackCellFromFace( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;
   void communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const override;

   inline void
       packVertexForCell( const Vertex* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override
   {
      dofPackInfo_.packVertexForCell( sender, receiver, buffer );
   }
   inline void unpackCellFromVertex( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override
   {
      dofPackInfo_.unpackCellFromVertex( receiver, sender, buffer );
   }
   inline void communicateLocalVertexToCell( const Vertex* sender, Cell* receiver ) const override
   {
      dofPackInfo_.communicateLocalVertexToCell( sender, receiver );
   }

   inline void
       packEdgeForCell( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override
   {
      dofPackInfo_.packEdgeForCell( sender, receiver, buffer );
   }
   void unpackCellFromEdge( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;
   void communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const override;

 private:
   using communication::DoFSpacePackInfo< ValueType >::level_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDVertex_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDEdge_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDFace_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDCell_;
   using communication::DoFSpacePackInfo< ValueType >::storage_;

   EdgeDoFPackInfo< ValueType > dofPackInfo_;
};

} // namespace n1e1
} //namespace hyteg
