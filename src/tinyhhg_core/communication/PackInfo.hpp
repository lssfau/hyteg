#pragma once

#include <tinyhhg_core/primitives/edge.hpp>
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
/// unpack data into Vertex from RecvBuffer from Edge
  virtual void unpackVertexFromEdge(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) = 0;
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

};

} // namespace communication
} // namespace hhg
