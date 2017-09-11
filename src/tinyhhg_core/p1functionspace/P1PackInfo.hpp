#pragma once

#include "tinyhhg_core/communication/PackInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"

namespace hhg {

class P1PackInfo : public communication::PackInfo {

public:
  P1PackInfo(uint_t level,
                   PrimitiveDataID<VertexP1FunctionMemory< real_t >, Vertex> dataIDVertex,
                   PrimitiveDataID<EdgeP1FunctionMemory< real_t >, Edge> dataIDEdge,
                   PrimitiveDataID<FaceP1FunctionMemory< real_t >, Face> dataIDFace,
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
  PrimitiveDataID<VertexP1FunctionMemory< real_t >, Vertex> dataIDVertex_;
  PrimitiveDataID<EdgeP1FunctionMemory< real_t >, Edge> dataIDEdge_;
  PrimitiveDataID<FaceP1FunctionMemory< real_t >, Face> dataIDFace_;
  std::weak_ptr<hhg::PrimitiveStorage> storage_;
};

} //namespace hhg
