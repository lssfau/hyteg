#pragma once

#include "tinyhhg_core/communication/PackInfo.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"

namespace hhg{


template< typename ValueType >
class DGPackInfo : public communication::PackInfo {

public:
  DGPackInfo(uint_t level,
                 PrimitiveDataID<FunctionMemory< ValueType >, Vertex> dataIDVertex,
                 PrimitiveDataID<FunctionMemory< ValueType >, Edge> dataIDEdge,
                 PrimitiveDataID<FunctionMemory< ValueType >, Face> dataIDFace,
                 std::weak_ptr<PrimitiveStorage> storage)
      : level_(level),
        dataIDVertex_(dataIDVertex),
        dataIDEdge_(dataIDEdge),
        dataIDFace_(dataIDFace),
        storage_(storage){

  }
  virtual void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  virtual void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  virtual void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) override;

  virtual void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  virtual void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  virtual void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) override;

  virtual void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  virtual void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  virtual void communicateLocalEdgeToFace(const Edge *sender, Face *receiver) override;

  virtual void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  virtual void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  virtual void communicateLocalFaceToEdge(const Face *sender, Edge *receiver) override;


private:
  uint_t level_;
  PrimitiveDataID<FunctionMemory< ValueType >, Vertex> dataIDVertex_;
  PrimitiveDataID<FunctionMemory< ValueType >, Edge> dataIDEdge_;
  PrimitiveDataID<FunctionMemory< ValueType >, Face> dataIDFace_;
  std::weak_ptr<hhg::PrimitiveStorage> storage_;
};

template< typename ValueType >
void DGPackInfo< ValueType >::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {

}

template< typename ValueType >
void DGPackInfo< ValueType >::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) {

}


}
