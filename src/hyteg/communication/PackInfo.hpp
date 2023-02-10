/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "core/Abort.h"

#include "hyteg/PrimitiveID.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {
// namespace containing function for communication between primitives
namespace communication {

enum class PackType
{
   PACK = 0,
   UNPACK = 1,
   DIRECT = 2
};

/// /brief Abstract class for pack and unpack functions of primitives
///
/// This abstract class provides the necessary functions to handle communication between the different primitives
/// (Vertex, Edge, Face, Cell). For each pair of primitives which need to communicate,
/// there must be a pack, unpack and a communicateLocal function
class PackInfo
{
 public:
   virtual ~PackInfo() {}

   /// @name Vertex to Edge
   ///@{
   /// pack data from Vertex into SendBuffer for Edge
   virtual void
       packVertexForEdge( const Vertex* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const = 0;
   /// unpack data into Edge from RecvBuffer from Vertex
   virtual void unpackEdgeFromVertex( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const = 0;
   /// transfer data locally (meaning both are on the same MPI Process) from Vertex to Edge
   virtual void communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const = 0;
   ///@}

   /// @name Edge to Vertex
   ///@{
   virtual void packEdgeForVertex( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const = 0;

   virtual void unpackVertexFromEdge( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const = 0;

   virtual void communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const = 0;
   ///@}

   /// @name Edge to Face
   ///@{
   virtual void packEdgeForFace( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const = 0;

   virtual void unpackFaceFromEdge( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const = 0;

   virtual void communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const = 0;
   ///@}

   /// @name Face to Edge
   ///@{
   virtual void packFaceForEdge( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const = 0;

   virtual void unpackEdgeFromFace( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const = 0;

   virtual void communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const = 0;
   ///@}

   /// @name Face to Cell
   ///@{
   virtual void packFaceForCell( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-face to macro-cell communication not implemented!" );
   }

   virtual void unpackCellFromFace( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-face to macro-cell communication not implemented!" );
   }

   virtual void communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
   {
      WALBERLA_ABORT( "Macro-face to macro-cell communication not implemented!" );
   }
   ///@}

   /// @name Face to Vertex
   ///@{
   virtual void packFaceForVertex( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-face to macro-vertex communication not implemented!" );
   }

   virtual void unpackVertexFromFace( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-face to macro-vertex communication not implemented!" );
   }

   virtual void communicateLocalFaceToVertex( const Face* sender, Vertex* receiver ) const
   {
      WALBERLA_ABORT( "Macro-face to macro-vertex communication not implemented!" );
   }
   ///@}

   /// @name Cell to Face
   ///@{
   virtual void packCellForFace( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-face communication not implemented!" );
   }

   virtual void unpackFaceFromCell( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-face communication not implemented!" );
   }

   virtual void communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-face communication not implemented!" );
   }
   ///@}

   /// @name Cell to Edge
   ///@{
   virtual void packCellForEdge( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-edge communication not implemented!" );
   }

   virtual void unpackEdgeFromCell( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-edge communication not implemented!" );
   }

   virtual void communicateLocalCellToEdge( const Cell* sender, Edge* receiver ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-edge communication not implemented!" );
   }
   ///@}

   /// @name Cell to Vertex
   ///@{
   virtual void packCellForVertex( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-vertex communication not implemented!" );
   }

   virtual void unpackVertexFromCell( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-vertex communication not implemented!" );
   }

   virtual void communicateLocalCellToVertex( const Cell* sender, Vertex* receiver ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-vertex communication not implemented!" );
   }
   ///@}

   /// @name Vertex to Cell
   ///@{
   virtual void packVertexForCell( const Vertex* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-vertex to macro-cell communication not implemented!" );
   }

   virtual void unpackCellFromVertex( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-vertex to macro-cell communication not implemented!" );
   }

   virtual void communicateLocalVertexToCell( const Vertex* sender, Cell* receiver ) const
   {
      WALBERLA_ABORT( "Macro-vertex to macro-cell communication not implemented!" );
   }
   ///@}

   /// @name Edge to Cell
   ///@{
   virtual void packEdgeForCell( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-edge to macro-cell communication not implemented!" );
   }

   virtual void unpackCellFromEdge( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-edge to macro-cell communication not implemented!" );
   }

   virtual void communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const
   {
      WALBERLA_ABORT( "Macro-edge to macro-cell communication not implemented!" );
   }
   ///@}

   /// @name Face to Face
   ///@{
   virtual void packFaceForFace( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-face to macro-face communication not implemented!" );
   }

   virtual void unpackFaceFromFace( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-face to macro-face communication not implemented!" );
   }

   virtual void communicateLocalFaceToFace( const Face* sender, Face* receiver ) const
   {
      WALBERLA_ABORT( "Macro-face to macro-face communication not implemented!" );
   }
   ///@}

   /// @name Cell to Cell
   ///@{
   virtual void packCellForCell( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-cell communication not implemented!" );
   }

   virtual void unpackCellFromCell( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-cell communication not implemented!" );
   }

   virtual void communicateLocalCellToCell( const Cell* sender, Cell* receiver ) const
   {
      WALBERLA_ABORT( "Macro-cell to macro-cell communication not implemented!" );
   }
   ///@}

   /// @name Generic communication methods
   ///@{
   template < typename SenderType, typename ReceiverType >
   void pack( const SenderType* sender, const PrimitiveID& receiverID, walberla::mpi::SendBuffer& sendBuffer ) const
   {
      static_assert( sizeof( SenderType ) == 0 /* always false */, "Invalid primitive type" );
   }

   template < typename SenderType, typename ReceiverType >
   void unpack( ReceiverType* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
   {
      static_assert( sizeof( SenderType ) == 0 /* always false */, "Invalid primitive type" );
   }

   template < typename SenderType, typename ReceiverType >
   void communicateLocal( const SenderType* sender, ReceiverType* receiver ) const
   {
      static_assert( sizeof( SenderType ) == 0 /* always false */, "Invalid primitive type" );
   }
   ///@}
};

template <>
inline void PackInfo::pack< Vertex, Edge >( const Vertex*              sender,
                                            const PrimitiveID&         receiverID,
                                            walberla::mpi::SendBuffer& sendBuffer ) const
{
   packVertexForEdge( sender, receiverID, sendBuffer );
}
template <>
inline void PackInfo::pack< Edge, Vertex >( const Edge*                sender,
                                            const PrimitiveID&         receiverID,
                                            walberla::mpi::SendBuffer& sendBuffer ) const
{
   packEdgeForVertex( sender, receiverID, sendBuffer );
}
template <>
inline void
    PackInfo::pack< Edge, Face >( const Edge* sender, const PrimitiveID& receiverID, walberla::mpi::SendBuffer& sendBuffer ) const
{
   packEdgeForFace( sender, receiverID, sendBuffer );
}
template <>
inline void
    PackInfo::pack< Face, Edge >( const Face* sender, const PrimitiveID& receiverID, walberla::mpi::SendBuffer& sendBuffer ) const
{
   packFaceForEdge( sender, receiverID, sendBuffer );
}
template <>
inline void
    PackInfo::pack< Face, Cell >( const Face* sender, const PrimitiveID& receiverID, walberla::mpi::SendBuffer& sendBuffer ) const
{
   packFaceForCell( sender, receiverID, sendBuffer );
}
template <>
inline void PackInfo::pack< Face, Vertex >( const Face*                sender,
                                            const PrimitiveID&         receiverID,
                                            walberla::mpi::SendBuffer& sendBuffer ) const
{
   packFaceForVertex( sender, receiverID, sendBuffer );
}
template <>
inline void
    PackInfo::pack< Cell, Face >( const Cell* sender, const PrimitiveID& receiverID, walberla::mpi::SendBuffer& sendBuffer ) const
{
   packCellForFace( sender, receiverID, sendBuffer );
}
template <>
inline void
    PackInfo::pack< Cell, Edge >( const Cell* sender, const PrimitiveID& receiverID, walberla::mpi::SendBuffer& sendBuffer ) const
{
   packCellForEdge( sender, receiverID, sendBuffer );
}
template <>
inline void PackInfo::pack< Cell, Vertex >( const Cell*                sender,
                                            const PrimitiveID&         receiverID,
                                            walberla::mpi::SendBuffer& sendBuffer ) const
{
   packCellForVertex( sender, receiverID, sendBuffer );
}
template <>
inline void PackInfo::pack< Vertex, Cell >( const Vertex*              sender,
                                            const PrimitiveID&         receiverID,
                                            walberla::mpi::SendBuffer& sendBuffer ) const
{
   packVertexForCell( sender, receiverID, sendBuffer );
}
template <>
inline void
    PackInfo::pack< Edge, Cell >( const Edge* sender, const PrimitiveID& receiverID, walberla::mpi::SendBuffer& sendBuffer ) const
{
   packEdgeForCell( sender, receiverID, sendBuffer );
}

template <>
inline void
    PackInfo::pack< Face, Face >( const Face* sender, const PrimitiveID& receiverID, walberla::mpi::SendBuffer& sendBuffer ) const
{
   packFaceForFace( sender, receiverID, sendBuffer );
}

template <>
inline void
    PackInfo::pack< Cell, Cell >( const Cell* sender, const PrimitiveID& receiverID, walberla::mpi::SendBuffer& sendBuffer ) const
{
   packCellForCell( sender, receiverID, sendBuffer );
}

template <>
inline void
    PackInfo::unpack< Vertex, Edge >( Edge* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackEdgeFromVertex( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Edge, Vertex >( Vertex* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackVertexFromEdge( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Edge, Face >( Face* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackFaceFromEdge( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Face, Edge >( Edge* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackEdgeFromFace( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Face, Cell >( Cell* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackCellFromFace( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Face, Vertex >( Vertex* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackVertexFromFace( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Cell, Face >( Face* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackFaceFromCell( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Cell, Edge >( Edge* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackEdgeFromCell( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Cell, Vertex >( Vertex* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackVertexFromCell( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Vertex, Cell >( Cell* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackCellFromVertex( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Edge, Cell >( Cell* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackCellFromEdge( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Face, Face >( Face* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackFaceFromFace( receiver, senderID, recvBuffer );
}
template <>
inline void
    PackInfo::unpack< Cell, Cell >( Cell* receiver, const PrimitiveID& senderID, walberla::mpi::RecvBuffer& recvBuffer ) const
{
   unpackCellFromCell( receiver, senderID, recvBuffer );
}

template <>
inline void PackInfo::communicateLocal< Vertex, Edge >( const Vertex* sender, Edge* receiver ) const
{
   communicateLocalVertexToEdge( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Edge, Vertex >( const Edge* sender, Vertex* receiver ) const
{
   communicateLocalEdgeToVertex( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Edge, Face >( const Edge* sender, Face* receiver ) const
{
   communicateLocalEdgeToFace( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Face, Edge >( const Face* sender, Edge* receiver ) const
{
   communicateLocalFaceToEdge( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Face, Cell >( const Face* sender, Cell* receiver ) const
{
   communicateLocalFaceToCell( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Face, Vertex >( const Face* sender, Vertex* receiver ) const
{
   communicateLocalFaceToVertex( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Cell, Face >( const Cell* sender, Face* receiver ) const
{
   communicateLocalCellToFace( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Cell, Edge >( const Cell* sender, Edge* receiver ) const
{
   communicateLocalCellToEdge( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Cell, Vertex >( const Cell* sender, Vertex* receiver ) const
{
   communicateLocalCellToVertex( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Vertex, Cell >( const Vertex* sender, Cell* receiver ) const
{
   communicateLocalVertexToCell( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Edge, Cell >( const Edge* sender, Cell* receiver ) const
{
   communicateLocalEdgeToCell( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Face, Face >( const Face* sender, Face* receiver ) const
{
   communicateLocalFaceToFace( sender, receiver );
}
template <>
inline void PackInfo::communicateLocal< Cell, Cell >( const Cell* sender, Cell* receiver ) const
{
   communicateLocalCellToCell( sender, receiver );
}

} // namespace communication
} // namespace hyteg
