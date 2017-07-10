#pragma once

#include "tinyhhg_core/communication/PackInfo.hpp"
#include "p1bubblememory.hpp"

namespace hhg {
namespace communication {

class P1BubblePackInfo : public PackInfo {

public:
  P1BubblePackInfo(uint_t level,
                   PrimitiveDataID<VertexP1BubbleFunctionMemory, Vertex> dataIDVertex,
                   PrimitiveDataID<EdgeP1BubbleFunctionMemory, Edge> dataIDEdge,
                   PrimitiveDataID<FaceP1BubbleFunctionMemory, Face> dataIDFace,
                   std::weak_ptr<PrimitiveStorage> storage)
          : level_(level),
            dataIDVertex_(dataIDVertex),
            dataIDEdge_(dataIDEdge),
            dataIDFace_(dataIDFace),
            storage_(storage){

  }
  virtual void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer);

  virtual void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer);

  virtual void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver);

  virtual void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer);

  virtual void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer);

  virtual void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver);

  virtual void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer);

  virtual void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer);

  virtual void communicateLocalEdgeToFace(const Edge *sender, Face *receiver);

  virtual void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer);

  virtual void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer);

  virtual void communicateLocalFaceToEdge(const Face *sender, Edge *receiver);


private:
  uint_t level_;
  PrimitiveDataID<VertexP1BubbleFunctionMemory, Vertex> dataIDVertex_;
  PrimitiveDataID<EdgeP1BubbleFunctionMemory, Edge> dataIDEdge_;
  PrimitiveDataID<FaceP1BubbleFunctionMemory, Face> dataIDFace_;
  std::weak_ptr<hhg::PrimitiveStorage> storage_;
};

} //namespace communication
} //namespace hhg