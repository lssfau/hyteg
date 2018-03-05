#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/bubblefunctionspace/BubbleMemory.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFace.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleDataHandling.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubblePackInfo.hpp"

namespace hhg {

template< typename ValueType >
class BubbleFunction : public Function< BubbleFunction< ValueType > >
{
public:

  BubbleFunction( const std::string &name, const std::shared_ptr< PrimitiveStorage > &storage, uint_t minLevel, uint_t maxLevel ) :
      Function< BubbleFunction< ValueType > >( name, storage, minLevel, maxLevel )
  {
    auto faceBubbleFunctionMemoryDataHandling = std::make_shared< FaceBubbleFunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
    auto edgeBubbleFunctionMemoryDataHandling = std::make_shared< EdgeBubbleFunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
    auto vertexBubbleFunctionMemoryDataHandling = std::make_shared< VertexBubbleFunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
    storage->addFaceData( faceDataID_, faceBubbleFunctionMemoryDataHandling, name );
    storage->addEdgeData( edgeDataID_, edgeBubbleFunctionMemoryDataHandling, name );
    storage->addVertexData( vertexDataID_, vertexBubbleFunctionMemoryDataHandling, name );
    for ( uint_t level = minLevel; level <= maxLevel; ++level )
    {
      communicators_[level]->addPackInfo( std::make_shared< BubblePackInfo< ValueType > >( level, vertexDataID_, edgeDataID_, faceDataID_, this->getStorage() ) );
    }
  }

  void assign(const std::vector<ValueType> scalars,
              const std::vector<BubbleFunction<ValueType> *> functions,
              size_t level,
              DoFType flag = All);

  void add(const std::vector<ValueType> scalars,
           const std::vector<BubbleFunction<ValueType> *> functions,
           size_t level,
           DoFType flag = All);

  real_t dot(BubbleFunction<ValueType> &rhs, size_t level, DoFType flag = All);

  const PrimitiveDataID<VertexBubbleFunctionMemory< ValueType >, Vertex> &getVertexDataID() const { return vertexDataID_; }

  const PrimitiveDataID<EdgeBubbleFunctionMemory< ValueType >, Edge> &getEdgeDataID() const { return edgeDataID_; }

  const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &getFaceDataID() const { return faceDataID_; }

 private:

  using Function< BubbleFunction< ValueType > >::communicators_;

  void enumerate_impl(size_t level, uint_t& num);

  PrimitiveDataID<VertexBubbleFunctionMemory< ValueType >, Vertex> vertexDataID_;
  PrimitiveDataID<EdgeBubbleFunctionMemory< ValueType >, Edge> edgeDataID_;
  PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> faceDataID_;
};


template< typename ValueType >
void BubbleFunction< ValueType >::assign(const std::vector< ValueType > scalars,
                            const std::vector<BubbleFunction<ValueType> *> functions,
                            size_t level,
                            DoFType flag)
{
    this->startTiming( "Assign" );
    // Collect all source IDs in a vector
    std::vector<PrimitiveDataID<VertexBubbleFunctionMemory< ValueType >, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeBubbleFunctionMemory< ValueType >, Edge>> srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face>> srcFaceIDs;

    for (auto &function : functions) {
      srcVertexIDs.push_back(function->vertexDataID_);
      srcEdgeIDs.push_back(function->edgeDataID_);
      srcFaceIDs.push_back(function->faceDataID_);
    }

    for (auto &it : this->getStorage()->getFaces()) {
      Face &face = *it.second;

      if (testFlag(face.type, flag)) {
        BubbleFace::assign< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
      }
    }
    this->stopTiming( "Assign" );
}

template< typename ValueType >
void BubbleFunction< ValueType >::add(const std::vector< ValueType > scalars,
                         const std::vector<BubbleFunction< ValueType > *> functions,
                         size_t level,
                         DoFType flag)
{
    this->startTiming( "Add" );
    // Collect all source IDs in a vector
    std::vector<PrimitiveDataID<VertexBubbleFunctionMemory< ValueType >, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeBubbleFunctionMemory< ValueType >, Edge>> srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face>> srcFaceIDs;

    for (auto &function : functions) {
      srcVertexIDs.push_back(function->vertexDataID_);
      srcEdgeIDs.push_back(function->edgeDataID_);
      srcFaceIDs.push_back(function->faceDataID_);
    }

    for (auto &it : this->getStorage()->getFaces()) {
      Face &face = *it.second;

      if (testFlag(face.type, flag)) {
        BubbleFace::add< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
      }
    }
  this->stopTiming( "Add" );
}

template< typename ValueType >
real_t BubbleFunction< ValueType >::dot(BubbleFunction< ValueType > &rhs, uint_t level, DoFType flag)
{
  this->startTiming( "Dot" );
  real_t scalarProduct = 0.0;

  for (auto &it : this->getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      scalarProduct += BubbleFace::dot< ValueType >(level, face, faceDataID_, rhs.faceDataID_);
    }
  }

  walberla::mpi::allReduceInplace(scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm());
  this->stopTiming( "Dot" );
  return scalarProduct;
}

template< typename ValueType >
void BubbleFunction< ValueType >::enumerate_impl(size_t level, uint_t& num)
{
  this->startTiming( "Enumerate" );
  for (auto &it : this->getStorage()->getFaces()) {
    Face &face = *it.second;

    BubbleFace::enumerate< ValueType >(level, face, faceDataID_, num);
  }
  communicators_[level]->template startCommunication<Face, Edge>();
  communicators_[level]->template endCommunication<Face, Edge>();

  communicators_[level]->template startCommunication<Edge, Vertex>();
  communicators_[level]->template endCommunication<Edge, Vertex>();
  this->stopTiming( "Enumerate" );
}

}
