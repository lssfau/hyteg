#pragma once

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/levelinfo.hpp"

namespace hhg {

namespace {
SPECIALIZE(uint_t,indexing::edgedof::macroedge::indexFromHorizontalEdge,edgeIndexFromHorizontalEdge)
SPECIALIZE(uint_t,indexing::edgedof::macroedge::indexFromVertex,edgeIndexFromVertex)
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
  using hhg::indexing::edgedof::macroface::BorderIterator;
  using indexing::edgedof::macroface::indexFromHorizontalEdge;
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender->getID());
  indexing::FaceBorderDirection faceDir = indexing::getFaceBorderDirection(edgeIndexOnFace,receiver->edge_orientation[edgeIndexOnFace]);
  uint_t indexOnEdge = 0;
  for(const auto& it : BorderIterator(level_,faceDir,0)){
    if(edgeIndexOnFace == 0) {
      faceData[faceIndexFromHorizontalEdge(level_, it.col(), it.row(), stencilDirection::EDGE_HO_C)] =
        edgeData[edgeIndexFromHorizontalEdge(level_,indexOnEdge,stencilDirection::EDGE_HO_C)];
    } else if(edgeIndexOnFace == 1){
      faceData[faceIndexFromHorizontalEdge(level_, it.col(), it.row(), stencilDirection::EDGE_DI_N)] =
        edgeData[edgeIndexFromHorizontalEdge(level_,indexOnEdge,stencilDirection::EDGE_HO_C)];
    } else if(edgeIndexOnFace == 2){
      faceData[faceIndexFromHorizontalEdge(level_, it.col(), it.row(), stencilDirection::EDGE_VE_NW)] =
        edgeData[edgeIndexFromHorizontalEdge(level_,indexOnEdge,stencilDirection::EDGE_HO_C)];
    } else {
      WALBERLA_ABORT("Wrong edgeIndexOnFace")
    }
    ++indexOnEdge;
  }
}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using hhg::indexing::edgedof::macroface::BorderIterator;
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = sender->edge_index(receiver);
  indexing::FaceBorderDirection faceBorderDir = indexing::getFaceBorderDirection(edgeIndexOnFace,sender->edge_orientation[edgeIndexOnFace]);
  stencilDirection faceDirOne;
  stencilDirection faceDirTwo;
  stencilDirection faceDirThree;
  if(edgeIndexOnFace == 0) {
    faceDirOne = stencilDirection::EDGE_HO_C;
    faceDirTwo = stencilDirection::EDGE_DI_N;
    faceDirThree = stencilDirection::EDGE_VE_NW;
  } else if(edgeIndexOnFace == 1){
    faceDirOne = stencilDirection::EDGE_DI_N;
    faceDirTwo = stencilDirection::EDGE_VE_NW;
    faceDirThree = stencilDirection::EDGE_HO_C;
  } else if(edgeIndexOnFace == 2){
    faceDirOne = stencilDirection::EDGE_VE_NW;
    faceDirTwo = stencilDirection::EDGE_HO_C;
    faceDirThree = stencilDirection::EDGE_DI_N;
  } else {
    WALBERLA_ABORT("Wrong edgeIndexOnFace")
  }
  for(const auto& it : BorderIterator(level_,faceBorderDir,1)){
      buffer << faceData[faceIndexFromHorizontalEdge(level_, it.col(), it.row(), faceDirOne)];
  }
  for(const auto& it : BorderIterator(level_,faceBorderDir,0)){
    buffer << faceData[faceIndexFromHorizontalEdge(level_, it.col(), it.row(), faceDirTwo)];
  }
}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );
  uint_t faceIdOnEdge = receiver->face_index(sender);
  stencilDirection dirHorizontal = faceIdOnEdge == 0 ? stencilDirection::EDGE_HO_SE : stencilDirection::EDGE_HO_NW;
  /// first edge is south edge by convention
  for (uint_t i = 1; i < levelinfo::num_microvertices_per_edge(level_) -1; ++i) {
    buffer >> edgeData[edgeIndexFromVertex(level_, i, dirHorizontal)];
  }
  stencilDirection dirDiagonal = faceIdOnEdge == 0 ? stencilDirection::EDGE_DI_SE : stencilDirection::EDGE_DI_NW;
  for (uint_t i = 1; i < levelinfo::num_microvertices_per_edge(level_); ++i) {
    buffer >> edgeData[edgeIndexFromVertex(level_, i, dirDiagonal)];
  }
}

template< typename ValueType>
void EdgeDoFPackInfo< ValueType >::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) {
  using hhg::indexing::edgedof::macroface::BorderIterator;
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = sender->edge_index(receiver->getID());
  indexing::FaceBorderDirection faceBorderDir = indexing::getFaceBorderDirection(edgeIndexOnFace,sender->edge_orientation[edgeIndexOnFace]);
  stencilDirection faceDir;
  if(edgeIndexOnFace == 0) {
    faceDir = stencilDirection::EDGE_HO_C;
  } else if(edgeIndexOnFace == 1){
    faceDir = stencilDirection::EDGE_DI_N;
  } else if(edgeIndexOnFace == 2){
    faceDir = stencilDirection::EDGE_VE_NW;
  } else {
    WALBERLA_ABORT("Wrong edgeIndexOnFace")
  }

  ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );
  uint_t faceIdOnEdge = receiver->face_index(sender->getID());
  stencilDirection edgeDir = faceIdOnEdge == 0 ? stencilDirection::EDGE_HO_SE : stencilDirection::EDGE_HO_NW;

  uint_t indexOnEdge = 1;
  for(const auto& it : BorderIterator(level_,faceBorderDir,1)){
    edgeData[edgeIndexFromVertex(level_, indexOnEdge, edgeDir)] =
        faceData[faceIndexFromHorizontalEdge(level_, it.col(), it.row(), faceDir)];
    ++indexOnEdge;
  }
}



} //namespace hhg