#pragma once

#include "tinyhhg_core/communication/PackInfo.hpp"
#include "BubbleMemory.hpp"

namespace hhg {

class BubblePackInfo : public communication::PackInfo {

public:
  BubblePackInfo(uint_t level,
                   PrimitiveDataID<VertexBubbleFunctionMemory, Vertex> dataIDVertex,
                   PrimitiveDataID<EdgeBubbleFunctionMemory, Edge> dataIDEdge,
                   PrimitiveDataID<FaceBubbleFunctionMemory, Face> dataIDFace,
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
  PrimitiveDataID<VertexBubbleFunctionMemory, Vertex> dataIDVertex_;
  PrimitiveDataID<EdgeBubbleFunctionMemory, Edge> dataIDEdge_;
  PrimitiveDataID<FaceBubbleFunctionMemory, Face> dataIDFace_;
  std::weak_ptr<hhg::PrimitiveStorage> storage_;
};

} //namespace hhg