#pragma once

#include <tinyhhg_core/primitives/Edge.hpp>
#include "tinyhhg_core/primitiveid.hpp"

namespace hhg {
/// namespace containing function for communication between primitives
namespace communication {
// WIP
///
/// /brief Abstract class for pack and unpack functions of primitives
///
/// This abstract class provides the necessary functions to handle communication between the different primitives
/// (Vertex, Edge, Face). For each pair of primitives which need to communicate,
/// there must be a pack, unpack and a communicateLocal function
///
class PackInfo {
public:

    virtual ~PackInfo() {}

  /// @name Vertex to Edge
  ///@{
  /// pack data from Vertex into SendBuffer for Edge
  virtual void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) = 0;
  /// unpack data into Edge from RecvBuffer from Vertex
  virtual void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) = 0;
  /// transfer data locally (meaning both are on the same MPI Process) from Vertex to Edge
  virtual void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) = 0;
  ///@}

  /// @name Edge to Vertex
  ///@{
  virtual void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) = 0;

  virtual void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) = 0;

  virtual void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) = 0;
  ///@}

  /// @name Edge to Face
  ///@{
  virtual void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) = 0;

  virtual void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) = 0;

  virtual void communicateLocalEdgeToFace(const Edge *sender, Face *receiver) = 0;
  ///@}

  /// @name Face to Edge
  ///@{
  virtual void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) = 0;

  virtual void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) = 0;

  virtual void communicateLocalFaceToEdge(const Face *sender, Edge *receiver) = 0;
  ///@}

  /// @name Generic communication methods
  ///@{
  template< typename SenderType, typename ReceiverType >
  void pack  ( const SenderType * sender,     const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer ) { static_assert( sizeof( SenderType ) == 0 /* always false */, "Invalid primitive type" ); }

  template< typename SenderType, typename ReceiverType >
  void unpack(       ReceiverType * receiver, const PrimitiveID & senderID,   walberla::mpi::RecvBuffer & recvBuffer ) { static_assert( sizeof( SenderType ) == 0 /* always false */, "Invalid primitive type" ); }

  template< typename SenderType, typename ReceiverType >
  void communicateLocal( const SenderType * sender, ReceiverType * receiver ) { static_assert( sizeof( SenderType ) == 0 /* always false */, "Invalid primitive type" ); }
  ///@}

};

template<>
inline void PackInfo::pack< Vertex, Edge > ( const Vertex * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer )
{
  packVertexForEdge( sender, receiverID, sendBuffer );
}
template<>
inline void PackInfo::pack< Edge, Vertex > ( const Edge * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer )
{
  packEdgeForVertex( sender, receiverID, sendBuffer );
}
template<>
inline void PackInfo::pack< Edge, Face > ( const Edge * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer )
{
  packEdgeForFace( sender, receiverID, sendBuffer );
}
template<>
inline void PackInfo::pack< Face, Edge > ( const Face * sender, const PrimitiveID & receiverID, walberla::mpi::SendBuffer & sendBuffer )
{
  packFaceForEdge( sender, receiverID, sendBuffer );
}


template<>
inline void PackInfo::unpack< Vertex, Edge > ( Edge * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer )
{
  unpackEdgeFromVertex( receiver, senderID, recvBuffer );
}
template<>
inline void PackInfo::unpack< Edge, Vertex > ( Vertex * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer )
{
  unpackVertexFromEdge( receiver, senderID, recvBuffer );
}
template<>
inline void PackInfo::unpack< Edge, Face > ( Face * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer )
{
  unpackFaceFromEdge( receiver, senderID, recvBuffer );
}
template<>
inline void PackInfo::unpack< Face, Edge > ( Edge * receiver, const PrimitiveID & senderID, walberla::mpi::RecvBuffer & recvBuffer )
{
  unpackEdgeFromFace( receiver, senderID, recvBuffer );
}


template<>
inline void PackInfo::communicateLocal< Vertex, Edge >( const Vertex * sender, Edge * receiver )
{
  communicateLocalVertexToEdge( sender, receiver );
}
template<>
inline void PackInfo::communicateLocal< Edge, Vertex >( const Edge * sender, Vertex * receiver )
{
  communicateLocalEdgeToVertex( sender, receiver );
}
template<>
inline void PackInfo::communicateLocal< Edge, Face >( const Edge * sender, Face * receiver )
{
  communicateLocalEdgeToFace( sender, receiver );
}
template<>
inline void PackInfo::communicateLocal< Face, Edge >( const Face * sender, Edge * receiver )
{
  communicateLocalFaceToEdge( sender, receiver );
}



} // namespace communication
} // namespace hhg
