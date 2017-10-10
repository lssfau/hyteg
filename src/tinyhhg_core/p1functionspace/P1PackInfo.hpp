#pragma once

#include "tinyhhg_core/communication/PackInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"

#include "tinyhhg_core/p1functionspace/P1PackInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Edge.hpp"
#include "tinyhhg_core/p1functionspace/P1EdgeIndex.hpp"
#include "tinyhhg_core/p1functionspace/P1Face.hpp"
#include "tinyhhg_core/p1functionspace/P1FaceIndex.hpp"

namespace hhg {

template< typename ValueType >
class P1PackInfo : public communication::PackInfo {

public:
  P1PackInfo(uint_t level,
                   PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> dataIDVertex,
                   PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> dataIDEdge,
                   PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> dataIDFace,
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
  PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> dataIDVertex_;
  PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> dataIDEdge_;
  PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> dataIDFace_;
  std::weak_ptr<hhg::PrimitiveStorage> storage_;
};


/// @name Vertex to Edge
///@{

template< typename ValueType >
void P1PackInfo< ValueType >::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer){
  WALBERLA_UNUSED(receiver);
  ValueType *vertexData = sender->getData(dataIDVertex_)->getPointer( level_ );
  buffer << vertexData[0];
}

template< typename ValueType >
void P1PackInfo< ValueType >::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer){
  using namespace P1Edge::EdgeCoordsVertex;
  ValueType* edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  //position in edge memory
  uint_t pos;
  if(receiver->vertex_index(sender) == 0){
    pos = 0;
  } else if(receiver->vertex_index(sender) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 1;
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << sender.getID() << " is not in Edge: " << receiver)
  }
  buffer >> edgeData[edge_index(level_,pos,VERTEX_C)];
}

template< typename ValueType >
void P1PackInfo< ValueType >::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver){
  using namespace P1Edge::EdgeCoordsVertex;
  ValueType *vertexData = sender->getData(dataIDVertex_)->getPointer( level_ );
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  uint_t pos;
  if(receiver->vertex_index(sender->getID()) == 0){
    pos = 0;
  } else if(receiver->vertex_index(sender->getID()) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 1;
  } else {
    WALBERLA_LOG_WARNING("Vertex: " << sender << " is not in Edge: " << receiver)
  }
  edgeData[edge_index(level_,pos,VERTEX_C)] = vertexData[0];
}

///@}
/// @name Edge to Vertex
///@{

template< typename ValueType >
void P1PackInfo< ValueType >::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using namespace P1Edge::EdgeCoordsVertex;
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  uint_t vertexIdOnEdge = sender->vertex_index(receiver);
  //the last element would be the vertex itself so we have to send the next one
  if(vertexIdOnEdge == 0){
    buffer << edgeData[edge_index(level_,1u,VERTEX_C)];
  } else if(vertexIdOnEdge == 1){
    buffer << edgeData[edge_index(level_,levelinfo::num_microvertices_per_edge(level_)-2,VERTEX_C)];
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver.getID() << " is not in Edge: " << sender);
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  ValueType *vertexData = receiver->getData(dataIDVertex_)->getPointer( level_ );
  uint_t edgeIdOnVertex = receiver->edge_index(sender);
  buffer >> vertexData[edgeIdOnVertex + 1];
}

template< typename ValueType >
void P1PackInfo< ValueType >::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) {
  using namespace P1Edge::EdgeCoordsVertex;
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *vertexData = receiver->getData(dataIDVertex_)->getPointer( level_ );
  uint_t vertexIdOnEdge = sender->vertex_index(receiver->getID());
  uint_t edgeIdOnVertex = receiver->edge_index(sender->getID());
  //the last element would be the vertex itself so we have to send the next one
  if(vertexIdOnEdge == 0){
    uint_t idx = edge_index(level_,1u,VERTEX_C);
    vertexData[edgeIdOnVertex+1] = edgeData[idx];
  } else if(vertexIdOnEdge == 1){
    uint_t idx = edge_index(level_,levelinfo::num_microvertices_per_edge(level_)-2,VERTEX_C);
    vertexData[edgeIdOnVertex+1] = edgeData[idx];
  } else {
    WALBERLA_LOG_WARNING("Vertex: " << receiver << " is not contained in Edge: " << sender);
  }
}

///@}
/// @name Edge to Face
///@{

template< typename ValueType >
void P1PackInfo< ValueType >::packEdgeForFace(const Edge *sender, const PrimitiveID &/*receiver*/, walberla::mpi::SendBuffer &buffer) {
  using namespace P1Edge::EdgeCoordsVertex;
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  uint_t v_perEdge = levelinfo::num_microvertices_per_edge(level_);

  for (uint_t i = 0; i < v_perEdge; ++i) {
    buffer << edgeData[edge_index(level_,i,VERTEX_C)];
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  using namespace hhg::P1Face;
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender);
  for(auto it = indexIterator(edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace], VERTEX, level_);
      it != indexIterator(); ++it){
    buffer >> faceData[*it];
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) {
  using namespace hhg::P1Face;
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender->getID());
  uint_t idx = 0;
  for(auto it = indexIterator(edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace], VERTEX, level_);
      it != indexIterator(); ++it)
  {
    faceData[*it] = edgeData[idx];
    idx++;
  }
}

///@}
/// @name Face to Edge
///@{

template< typename ValueType >
void P1PackInfo< ValueType >::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using namespace hhg::P1Face;
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = sender->edge_index(receiver);
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], VERTEX_INNER, level_);
      it != indexIterator(); ++it){
    buffer << faceData[*it];
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  uint_t v_perEdge = levelinfo::num_microvertices_per_edge(level_);
  uint_t edgeIdOnFace = receiver->face_index(sender);
  P1Edge::EdgeCoordsVertex::DirVertex dir;
  //the first face is the south face and the second the north face
  if(edgeIdOnFace == 0)
  {
    dir = P1Edge::EdgeCoordsVertex::VERTEX_SE;
  }
  else if(edgeIdOnFace == 1)
  {
    dir = P1Edge::EdgeCoordsVertex::VERTEX_N;
  }
  for (uint_t i = 0; i < v_perEdge -1; ++i)
  {
    buffer >> edgeData[P1Edge::EdgeCoordsVertex::edge_index(level_,i,dir)];
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) {
  using namespace hhg::P1Face;
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t faceIdOnEdge = receiver->face_index(sender->getID());
  P1Edge::EdgeCoordsVertex::DirVertex dir;
  uint_t edgeIdOnFace = sender->edge_index(receiver->getID());
  //the first face is the south face and the second the north face
  if(faceIdOnEdge == 0)
  {
    dir = P1Edge::EdgeCoordsVertex::VERTEX_SE;
  }
  else if(faceIdOnEdge == 1)
  {
    dir = P1Edge::EdgeCoordsVertex::VERTEX_N;
  }
  uint_t idx = 0;
  for(auto it = indexIterator(edgeIdOnFace, sender->edge_orientation[edgeIdOnFace], VERTEX_INNER, level_);
      it != indexIterator(); ++it) {
    edgeData[hhg::P1Edge::EdgeCoordsVertex::edge_index(level_,idx,dir)] = faceData[*it];
    idx++;
  }
}


} //namespace hhg
