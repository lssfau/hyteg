#pragma once

#include "tinyhhg_core/Function.hpp"
#include "DGPackInfo.hpp"
#include "DGMemory.hpp"
#include "DGDataHandling.hpp"
#include "DGVertex.hpp"
#include "DGEdge.hpp"
#include "DGFace.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"

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
        std::make_shared<DGPackInfo<ValueType> >(level, vertexDataID_, edgeDataID_, faceDataID_, this->getStorage()));
    }
  }

  inline void interpolate(std::function<ValueType(const Point3D &)>& expr,
                          uint_t level,
                          DoFType flag = All);

  inline void interpolateExtended(std::function<ValueType(const Point3D&, const std::vector<ValueType>&)>& expr,
                                  const std::vector<DGFunction<ValueType>*> srcFunctions,
                                  uint_t level,
                                  DoFType flag = All);

  inline void assign(const std::vector<ValueType> scalars,
                     const std::vector<DGFunction< ValueType >*> functions,
                     uint_t level,
                     DoFType flag = All);

  inline void add(const std::vector<ValueType> scalars,
                  const std::vector<DGFunction< ValueType >*> functions,
                  uint_t level,
                  DoFType flag = All);

  const PrimitiveDataID<FunctionMemory<ValueType>, Vertex> &getVertexDataID() const { return vertexDataID_; }

  const PrimitiveDataID<FunctionMemory<ValueType>, Edge> &getEdgeDataID() const { return edgeDataID_; }

  const PrimitiveDataID<FunctionMemory<ValueType>, Face> &getFaceDataID() const { return faceDataID_; }

  void projectP1(P1Function< real_t >& src, uint_t level, DoFType flag, UpdateType updateType = Replace);

private:

  using Function<DGFunction<ValueType> >::communicators_;

  PrimitiveDataID<FunctionMemory<ValueType>, Vertex> vertexDataID_;
  PrimitiveDataID<FunctionMemory<ValueType>, Edge> edgeDataID_;
  PrimitiveDataID<FunctionMemory<ValueType>, Face> faceDataID_;



  inline void enumerate_impl(uint_t level, uint_t& num) override;

};

template< typename ValueType >
void DGFunction< ValueType >::add(const std::vector<ValueType> scalars, const std::vector<DGFunction<ValueType> *> functions, uint_t level,
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

  for (auto &it : this->getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      DGFace::add< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
    }
  }

}

template< typename ValueType >
inline void DGFunction< ValueType >::interpolate(std::function< ValueType( const Point3D& ) >& expr,
                                                 uint_t level, DoFType flag)
{
  std::function< ValueType(const Point3D&,const std::vector<ValueType>&)> exprExtended = [&expr](const hhg::Point3D& x, const std::vector<ValueType>&) {
      return expr(x);
  };
  interpolateExtended( exprExtended, {}, level, flag );
}

template< typename ValueType >
void DGFunction< ValueType >::interpolateExtended(std::function<ValueType(const Point3D &, const std::vector<ValueType>&)> &expr,
                                                  const std::vector<DGFunction<ValueType>*> srcFunctions,
                                                  uint_t level,
                                                  DoFType flag)
{
  this->startTiming( "Interpolate" );
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

  for (auto &it : this->getStorage()->getVertices()) {
    Vertex &vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag)) {
      DGVertex::interpolate< ValueType >(level, vertex, vertexDataID_, srcVertexIDs, expr, this->getStorage());
    }
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();

  for (auto &it : this->getStorage()->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.getDoFType(), flag)) {
      DGEdge::interpolate< ValueType >(level, edge, edgeDataID_, srcEdgeIDs, expr, this->getStorage());
    }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for (auto &it : this->getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      DGFace::interpolate< ValueType >(level, face, faceDataID_, srcFaceIDs, expr);
    }
  }

  communicators_[level]->template endCommunication<Edge, Face>();
  this->stopTiming( "Interpolate" );
}

template< typename ValueType >
void DGFunction< ValueType >::assign(const std::vector<ValueType> scalars,
                                     const std::vector<DGFunction<ValueType> *> functions,
                                     uint_t level,
                                     DoFType flag)
{
  this->startTiming( "Assign" );
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

  for (auto &it : this->getStorage()->getVertices()) {
    Vertex &vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag)) {
      DGVertex::assign< ValueType >(level, vertex, scalars, srcVertexIDs, vertexDataID_);
    }
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();

  for (auto &it : this->getStorage()->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.getDoFType(), flag)) {
      DGEdge::assign< ValueType >(level, edge, scalars, srcEdgeIDs, edgeDataID_);
    }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for (auto &it : this->getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      DGFace::assign< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
    }
  }

  communicators_[level]->template endCommunication<Edge, Face>();
  this->stopTiming( "Assign" );
}


template< typename ValueType >
void DGFunction< ValueType >::enumerate_impl(uint_t level, uint_t &num)
{
  this->startTiming( "Enumerate" );
  for (auto& it : this->getStorage()->getVertices()) {
    Vertex& vertex = *it.second;
    DGVertex::enumerate(vertex,vertexDataID_,level,num);
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();
  communicators_[level]->template endCommunication<Vertex, Edge>();

  for (auto& it : this->getStorage()->getEdges()) {
    Edge& edge = *it.second;
    DGEdge::enumerate< ValueType >(level, edge, edgeDataID_, num);
  }

  communicators_[level]->template startCommunication<Edge, Face>();
  communicators_[level]->template endCommunication<Edge, Face>();

  for (auto& it : this->getStorage()->getFaces()) {
    Face& face = *it.second;
    DGFace::enumerate< ValueType >(level, face, faceDataID_, num);
  }

  communicators_[level]->template startCommunication<Face, Edge>();
  communicators_[level]->template endCommunication<Face, Edge>();

  communicators_[level]->template startCommunication<Edge, Vertex>();
  communicators_[level]->template endCommunication<Edge, Vertex>();
  this->stopTiming( "Enumerate" );
}

template< typename ValueType >
void DGFunction< ValueType >::projectP1(P1Function< real_t >& src, uint_t level, DoFType flag, UpdateType updateType)
{
  this->startTiming( "projectP1" );

  src.getCommunicator(level)->template startCommunication<Edge, Vertex>();
  src.getCommunicator(level)->template startCommunication<Face, Edge>();
  src.getCommunicator(level)->template endCommunication<Edge, Vertex>();

  for (auto& it : this->getStorage()->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      DGVertex::projectP1< real_t >(level, vertex, this->getStorage(), src.getVertexDataID(), this->getVertexDataID(), updateType);
    }
  }

  this->getCommunicator(level)->template startCommunication<Vertex, Edge>();

  src.getCommunicator(level)->template endCommunication<Face, Edge>();

  for (auto& it : this->getStorage()->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      DGEdge::projectP1< real_t >(level, edge, this->getStorage(), src.getEdgeDataID(), this->getEdgeDataID(), updateType);
    }
  }

  this->getCommunicator(level)->template endCommunication<Vertex, Edge>();

  this->getCommunicator(level)->template startCommunication<Edge, Face>();

  for (auto& it : this->getStorage()->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      DGFace::projectP1<real_t>(level, face, this->getStorage(), src.getFaceDataID(), this->getFaceDataID(), updateType);
    }
  }

  this->getCommunicator(level)->template endCommunication<Edge, Face>();

  this->stopTiming( "projectP1" );
}

}
