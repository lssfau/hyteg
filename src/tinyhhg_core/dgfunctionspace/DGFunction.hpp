#pragma once

#include "tinyhhg_core/Function.hpp"
#include "DGPackInfo.hpp"
#include "DGMemory.hpp"
#include "DGDataHandling.hpp"
#include "DGVertex.hpp"
#include "DGEdge.hpp"
#include "DGFace.hpp"

namespace hhg {

template<typename ValueType>
class DGFunction : public Function<DGFunction<ValueType> >
{
public:

  DGFunction(const std::string &name, const std::shared_ptr<PrimitiveStorage> &storage, uint_t minLevel,
             uint_t maxLevel) :
    Function<DGFunction<ValueType> >(name, storage, minLevel, maxLevel) {
    auto vertexDGFunctionMemoryDataHandling =
      std::make_shared<VertexDGFunctionMemoryDataHandling<ValueType> >(minLevel, maxLevel);
    auto edgeDGFunctionMemoryDataHandling =
      std::make_shared<EdgeDGFunctionMemoryDataHandling<ValueType> >(minLevel, maxLevel);
    auto faceDGFunctionMemoryDataHandling =
      std::make_shared<FaceDGFunctionMemoryDataHandling<ValueType> >(minLevel, maxLevel);


    storage->addFaceData(faceDataID_, faceDGFunctionMemoryDataHandling, name);
    storage->addEdgeData(edgeDataID_, edgeDGFunctionMemoryDataHandling, name);
    storage->addVertexData(vertexDataID_, vertexDGFunctionMemoryDataHandling, name);
    for (uint_t level = minLevel; level <= maxLevel; ++level) {
      //communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
      communicators_[level]->addPackInfo(
        std::make_shared<DGPackInfo<ValueType> >(level, vertexDataID_, edgeDataID_, faceDataID_, storage_));
    }
  }

  const PrimitiveDataID<FunctionMemory<ValueType>, Vertex> &getVertexDataID() const { return vertexDataID_; }

  const PrimitiveDataID<FunctionMemory<ValueType>, Edge> &getEdgeDataID() const { return edgeDataID_; }

  const PrimitiveDataID<FunctionMemory<ValueType>, Face> &getFaceDataID() const { return faceDataID_; }

  void interpolate(std::function<ValueType(const Point3D&, const std::vector<ValueType>&)>& expr, const std::vector<DGFunction<ValueType>*> srcFunctions, uint_t level, DoFType flag = All);

private:

  using Function<DGFunction<ValueType> >::storage_;
  using Function<DGFunction<ValueType> >::communicators_;

  PrimitiveDataID<FunctionMemory<ValueType>, Vertex> vertexDataID_;
  PrimitiveDataID<FunctionMemory<ValueType>, Edge> edgeDataID_;
  PrimitiveDataID<FunctionMemory<ValueType>, Face> faceDataID_;

  inline void interpolate_impl(std::function<ValueType(const Point3D&)>& expr, uint_t level, DoFType flag = All) override;

  inline void assign_impl(const std::vector<ValueType> scalars,
                          const std::vector<DGFunction< ValueType >*> functions,
                          uint_t level,
                          DoFType flag = All) override;

  inline void add_impl(const std::vector<ValueType> scalars,
                       const std::vector<DGFunction< ValueType >*> functions,
                       uint_t level,
                       DoFType flag = All) override;

  inline real_t dot_impl(DGFunction< ValueType >& rhs, uint_t level, DoFType flag = All) override;

  inline void prolongate_impl(uint_t level, DoFType flag = All) override;

  inline void prolongateQuadratic_impl(uint_t level, DoFType flag = All) override;

  inline void restrict_impl(uint_t level, DoFType flag = All) override;

  inline void enumerate_impl(uint_t level, uint_t& num) override;

};

template< typename ValueType >
void DGFunction< ValueType >::add_impl(const std::vector<ValueType> scalars, const std::vector<DGFunction<ValueType> *> functions, uint_t level,
                          DoFType flag) {

  // Collect all source IDs in a vector
//  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> srcVertexIDs;
//  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>>   srcEdgeIDs;
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>>   srcFaceIDs;

  for (auto& function : functions)
  {
//    srcVertexIDs.push_back(function->vertexDataID_);
//    srcEdgeIDs.push_back(function->edgeDataID_);
    srcFaceIDs.push_back(function->faceDataID_);
  }

  for (auto &it : storage_->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      DGFace::add< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
    }
  }

}

template< typename ValueType >
void DGFunction< ValueType >::interpolate(std::function<ValueType(const Point3D &, const std::vector<ValueType>&)> &expr,
                                          const std::vector<DGFunction<ValueType>*> srcFunctions,
                                          uint_t level,
                                          DoFType flag) {

  this->startTiming( "Interpolate");

  // Collect all source IDs in a vector
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> srcVertexIDs;
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>>   srcEdgeIDs;
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>>   srcFaceIDs;

  for (auto& function : srcFunctions)
  {
    srcVertexIDs.push_back(function->vertexDataID_);
    srcEdgeIDs.push_back(function->edgeDataID_);
    srcFaceIDs.push_back(function->faceDataID_);
  }

  for (auto &it : storage_->getVertices()) {
    Vertex &vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag)) {
      DGVertex::interpolate< ValueType >(level, vertex, vertexDataID_, srcVertexIDs, expr, storage_);
    }
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();

  for (auto &it : storage_->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.getDoFType(), flag)) {
      DGEdge::interpolate< ValueType >(level, edge, edgeDataID_, srcEdgeIDs, expr, storage_);
    }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for (auto &it : storage_->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      DGFace::interpolate< ValueType >(level, face, faceDataID_, srcFaceIDs, expr);
    }
  }

  communicators_[level]->template endCommunication<Edge, Face>();

  this->stopTiming( "Interpolate");

}

template< typename ValueType >
inline void DGFunction< ValueType >::interpolate_impl(std::function<ValueType(const Point3D&)>& expr, uint_t level, DoFType flag) {
  std::function<ValueType(const Point3D&, const std::vector<ValueType>&)> tmpExpr = [&] (const Point3D& x, const std::vector<ValueType>&) {
    return expr(x);
  };
  WALBERLA_ABORT("this does not work currently use other interpolate")
  //interpolate(tmpExpr, {}, level, flag);
}

template< typename ValueType >
void DGFunction< ValueType >::assign_impl(const std::vector<ValueType> scalars,
                                          const std::vector<DGFunction<ValueType> *> functions,
                                          uint_t level,
                                          DoFType flag) {
  // Collect all source IDs in a vector
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> srcVertexIDs;
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>>   srcEdgeIDs;
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>>   srcFaceIDs;

  for (auto& function : functions)
  {
    srcVertexIDs.push_back(function->vertexDataID_);
    srcEdgeIDs.push_back(function->edgeDataID_);
    srcFaceIDs.push_back(function->faceDataID_);
  }

  for (auto &it : storage_->getVertices()) {
    Vertex &vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag)) {
      DGVertex::assign< ValueType >(level, vertex, scalars, srcVertexIDs, vertexDataID_);
    }
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();

  for (auto &it : storage_->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.getDoFType(), flag)) {
      DGEdge::assign< ValueType >(level, edge, scalars, srcEdgeIDs, edgeDataID_);
    }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for (auto &it : storage_->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      DGFace::assign< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
    }
  }

  communicators_[level]->template endCommunication<Edge, Face>();
}

template< typename ValueType >
real_t DGFunction< ValueType >::dot_impl(DGFunction<ValueType> &rhs,
                                         uint_t level,
                                         DoFType flag) {
  return 0;
}

template< typename ValueType >
void DGFunction< ValueType >::prolongate_impl(uint_t level, DoFType flag) {

}

template< typename ValueType >
void DGFunction< ValueType >::prolongateQuadratic_impl(uint_t level, DoFType flag) {

}

template< typename ValueType >
void DGFunction< ValueType >::restrict_impl(uint_t level, DoFType flag) {

}

template< typename ValueType >
void DGFunction< ValueType >::enumerate_impl(uint_t level, uint_t &num) {
  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;
    DGVertex::enumerate(vertex,vertexDataID_,level,num);
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();
  communicators_[level]->template endCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;
    DGEdge::enumerate< ValueType >(level, edge, edgeDataID_, num);
  }

  communicators_[level]->template startCommunication<Edge, Face>();
  communicators_[level]->template endCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;
    DGFace::enumerate< ValueType >(level, face, faceDataID_, num);
  }

  communicators_[level]->template startCommunication<Face, Edge>();
  communicators_[level]->template endCommunication<Face, Edge>();

  communicators_[level]->template startCommunication<Edge, Vertex>();
  communicators_[level]->template endCommunication<Edge, Vertex>();
}

}
