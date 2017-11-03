#pragma once
#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"

namespace hhg {

using walberla::uint_t;

template<typename ValueType>
class EdgeDoFPackInfo : public communication::DoFSpacePackInfo<ValueType> {

public:
  EdgeDoFPackInfo(uint_t level,
                  PrimitiveDataID<FunctionMemory<ValueType>, Vertex> dataIDVertex,
                  PrimitiveDataID<FunctionMemory<ValueType>, Edge> dataIDEdge,
                  PrimitiveDataID<FunctionMemory<ValueType>, Face> dataIDFace,
                  std::weak_ptr<PrimitiveStorage> storage)
    : communication::DoFSpacePackInfo<ValueType>(level, dataIDVertex, dataIDEdge, dataIDFace, storage) {

  }

  void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) override;

  void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) override;

  void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  void communicateLocalEdgeToFace(const Edge *sender, Face *receiver) override;

  void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  void communicateLocalFaceToEdge(const Face *sender, Edge *receiver) override;

private:
  using communication::DoFSpacePackInfo<ValueType>::level_;
  using communication::DoFSpacePackInfo<ValueType>::dataIDVertex_;
  using communication::DoFSpacePackInfo<ValueType>::dataIDEdge_;
  using communication::DoFSpacePackInfo<ValueType>::dataIDFace_;
  using communication::DoFSpacePackInfo<ValueType>::storage_;

};

} //namespace hhg