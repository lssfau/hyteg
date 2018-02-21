#pragma once

#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/primitives/all.hpp"

#include "core/Abort.h"

namespace hhg {
/// namespace containing function for communication between primitives
namespace communication {

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
  virtual void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const = 0;
  /// unpack data into Edge from RecvBuffer from Vertex
  virtual void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const = 0;
  /// transfer data locally (meaning both are on the same MPI Process) from Vertex to Edge
  virtual void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) const = 0;
  ///@}

  /// @name Edge to Vertex
  ///@{
  virtual void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const = 0;

  virtual void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const = 0;

  virtual void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) const = 0;
  ///@}

  /// @name Edge to Face
  ///@{
  virtual void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const = 0;

  virtual void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const = 0;

  virtual void communicateLocalEdgeToFace(const Edge *sender, Face *receiver) const = 0;
  ///@}

  /// @name Face to Edge
  ///@{
  virtual void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const = 0;

  virtual void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const = 0;

  virtual void communicateLocalFaceToEdge(const Face *sender, Edge *receiver) const = 0;
  ///@}

  /// @name Face to Cell
  ///@{
  virtual void packFaceForCell(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
  {
    WALBERLA_ABORT( "Macro-face to macro-cell communication not implemented!" );
  }

  virtual void unpackCellFromFace(Cell *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
  {
    WALBERLA_ABORT( "Macro-face to macro-cell communication not implemented!" );
  }

  virtual void communicateLocalFaceToCell(const Face *sender, Cell *receiver) const
  {
    WALBERLA_ABORT( "Macro-face to macro-cell communication not implemented!" );
  }
  ///@}

  /// @name Cell to Face
  ///@{
  virtual void packCellForFace(const Cell *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
  {
    WALBERLA_ABORT( "Macro-cell to macro-face communication not implemented!" );
  }

  virtual void unpackFaceFromCell(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
  {
    WALBERLA_ABORT( "Macro-cell to macro-face communication not implemented!" );
  }

  virtual void communicateLocalCellToFace(const Cell *sender, Face *receiver) const
  {
    WALBERLA_ABORT( "Macro-cell to macro-face communication not implemented!" );
  }
  ///@}

  /// @name Generic communication methods
  ///@{
  template< typename SenderType, typename ReceiverType >
  void pack  ( const SenderType * sender,     const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer ) const { static_assert( sizeof( SenderType ) == 0 /* always false */, "Invalid primitive type" ); }

  template< typename SenderType, typename ReceiverType >
  void unpack(       ReceiverType * receiver, const PrimitiveID & senderID,   walberla::mpi::RecvBuffer & recvBuffer ) const { static_assert( sizeof( SenderType ) == 0 /* always false */, "Invalid primitive type" ); }

  template< typename SenderType, typename ReceiverType >
  void communicateLocal( const SenderType * sender, ReceiverType * receiver ) const { static_assert( sizeof( SenderType ) == 0 /* always false */, "Invalid primitive type" ); }
  ///@}

};

template<>
inline void PackInfo::pack< Vertex, Edge > ( const Vertex * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer ) const
{
  packVertexForEdge( sender, receiverID, sendBuffer );
}
template<>
inline void PackInfo::pack< Edge, Vertex > ( const Edge * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer ) const
{
  packEdgeForVertex( sender, receiverID, sendBuffer );
}
template<>
inline void PackInfo::pack< Edge, Face > ( const Edge * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer ) const
{
  packEdgeForFace( sender, receiverID, sendBuffer );
}
template<>
inline void PackInfo::pack< Face, Edge > ( const Face * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer ) const
{
  packFaceForEdge( sender, receiverID, sendBuffer );
}
template<>
inline void PackInfo::pack< Face, Cell > ( const Face * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer ) const
{
  packFaceForCell( sender, receiverID, sendBuffer );
}
template<>
inline void PackInfo::pack< Cell, Face > ( const Cell * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer ) const
{
  packCellForFace( sender, receiverID, sendBuffer );
}


template<>
inline void PackInfo::unpack< Vertex, Edge > ( Edge * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer ) const
{
  unpackEdgeFromVertex( receiver, senderID, recvBuffer );
}
template<>
inline void PackInfo::unpack< Edge, Vertex > ( Vertex * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer ) const
{
  unpackVertexFromEdge( receiver, senderID, recvBuffer );
}
template<>
inline void PackInfo::unpack< Edge, Face > ( Face * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer ) const
{
  unpackFaceFromEdge( receiver, senderID, recvBuffer );
}
template<>
inline void PackInfo::unpack< Face, Edge > ( Edge * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer ) const
{
  unpackEdgeFromFace( receiver, senderID, recvBuffer );
}
template<>
inline void PackInfo::unpack< Face, Cell > ( Cell * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer ) const
{
  unpackCellFromFace( receiver, senderID, recvBuffer );
}
template<>
inline void PackInfo::unpack< Cell, Face > ( Face * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer ) const
{
  unpackFaceFromCell( receiver, senderID, recvBuffer );
}


template<>
inline void PackInfo::communicateLocal< Vertex, Edge >( const Vertex * sender, Edge * receiver ) const
{
  communicateLocalVertexToEdge( sender, receiver );
}
template<>
inline void PackInfo::communicateLocal< Edge, Vertex >( const Edge * sender, Vertex * receiver ) const
{
  communicateLocalEdgeToVertex( sender, receiver );
}
template<>
inline void PackInfo::communicateLocal< Edge, Face >( const Edge * sender, Face * receiver ) const
{
  communicateLocalEdgeToFace( sender, receiver );
}
template<>
inline void PackInfo::communicateLocal< Face, Edge >( const Face * sender, Edge * receiver ) const
{
  communicateLocalFaceToEdge( sender, receiver );
}
template<>
inline void PackInfo::communicateLocal< Face, Cell >( const Face * sender, Cell * receiver ) const
{
  communicateLocalFaceToCell( sender, receiver );
}
template<>
inline void PackInfo::communicateLocal< Cell, Face >( const Cell * sender, Face * receiver ) const
{
  communicateLocalCellToFace( sender, receiver );
}


} // namespace communication
} // namespace hhg
