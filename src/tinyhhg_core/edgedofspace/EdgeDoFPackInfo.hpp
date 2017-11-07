#pragma once

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/levelinfo.hpp"

namespace hhg {

namespace {
SPECIALIZE(uint_t,indexing::edgedof::macroedge::indexFromHorizontalEdge,edgeIndexFromHorizontalEdge)
SPECIALIZE(uint_t,indexing::edgedof::macroface::indexFromHorizontalEdge,faceIndexFromHorizontalEdge)
}

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

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {

}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {

}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) {

}


template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {

}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {

}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) {

}

template< typename ValueType >
void EdgeDoFPackInfo< ValueType >::packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );
  for(uint_t i = 0; i < levelinfo::num_microedges_per_edge( level_ ); ++i){
    buffer << edgeData[edgeIndexFromHorizontalEdge(level_,i,stencilDirection::EDGE_HO_C)];
  }
}

template< typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  using hhg::indexing::edgedof::macroface::BorderIterator;
  using indexing::edgedof::macroface::indexFromHorizontalEdge;
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender);
  indexing::FaceBorderDirection faceDir = indexing::getFaceBorderDirection(edgeIndexOnFace,receiver->edge_orientation[edgeIndexOnFace]);
  for(const auto& it : BorderIterator(level_,faceDir,0)){
    if(edgeIndexOnFace == 0) {
      buffer >> faceData[faceIndexFromHorizontalEdge(level_, it.col(), it.row(), stencilDirection::EDGE_HO_C)];
    } else if(edgeIndexOnFace == 1){
      buffer >> faceData[faceIndexFromHorizontalEdge(level_, it.col(), it.row(), stencilDirection::EDGE_DI_N)];
    } else if(edgeIndexOnFace == 2){
    buffer >> faceData[faceIndexFromHorizontalEdge(level_, it.col(), it.row(), stencilDirection::EDGE_VE_NW)];
    } else {
      WALBERLA_ABORT("Wrong edgeIndexOnFace")
    }
  }
}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) {

}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {

}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {

}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) {

}



} //namespace hhg