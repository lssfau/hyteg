#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/bubblefunctionspace/BubbleMemory.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFace.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleDataHandling.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubblePackInfo.hpp"

namespace hhg {

template< typename ValueType >
class BubbleFunction : public Function< BubbleFunction< ValueType > > {
 public:

  using Function< BubbleFunction< ValueType > >::storage_;
  using Function< BubbleFunction< ValueType > >::communicators_;

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
      communicators_[level]->setLocalCommunicationMode( communication::BufferedCommunicator::BUFFERED_MPI );
      communicators_[level]->addPackInfo( std::make_shared< BubblePackInfo< ValueType > >( level, vertexDataID_, edgeDataID_, faceDataID_, storage_ ) );
    }
  }

  const PrimitiveDataID<VertexBubbleFunctionMemory< ValueType >, Vertex> &getVertexDataID() const { return vertexDataID_; }

  const PrimitiveDataID<EdgeBubbleFunctionMemory< ValueType >, Edge> &getEdgeDataID() const { return edgeDataID_; }

  const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &getFaceDataID() const { return faceDataID_; }

 private:

  /// Interpolates a given expression to a P1Function
  void interpolate_impl(std::function<ValueType(const Point3D &)> &expr, uint_t level, DoFType flag = All);

  void assign_impl(const std::vector<ValueType> scalars,
              const std::vector<BubbleFunction<ValueType> *> functions,
              size_t level,
              DoFType flag = All);

  void add_impl(const std::vector<ValueType> scalars,
           const std::vector<BubbleFunction<ValueType> *> functions,
           size_t level,
           DoFType flag = All);

  real_t dot_impl(BubbleFunction<ValueType> &rhs, size_t level, DoFType flag = All);

  void prolongate_impl(size_t level, DoFType flag = All) {
    WALBERLA_UNUSED( level );
    WALBERLA_UNUSED( flag );
    WALBERLA_ABORT("Bubble prolongate not implemented");
  }

  void prolongateQuadratic_impl(size_t level, DoFType flag = All) {
    WALBERLA_UNUSED( level );
    WALBERLA_UNUSED( flag );
    WALBERLA_ABORT("Bubble prolongate quadratic not implemented");
  }

  void restrict_impl(size_t level, DoFType flag = All) {
    WALBERLA_UNUSED( level );
    WALBERLA_UNUSED( flag );
    WALBERLA_ABORT("Bubble restrict not implemented");
  }

  void enumerate_impl(size_t level, uint_t& num);

  PrimitiveDataID<VertexBubbleFunctionMemory< ValueType >, Vertex> vertexDataID_;
  PrimitiveDataID<EdgeBubbleFunctionMemory< ValueType >, Edge> edgeDataID_;
  PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> faceDataID_;
};

template< typename ValueType >
void BubbleFunction< ValueType >::interpolate_impl(std::function<ValueType(const hhg::Point3D &)> &expr, uint_t level, DoFType flag) {
  // TODO: implement Bubble interpolation. It is not required for Dirichlet only interpolation in most of the apps
  WALBERLA_UNUSED( expr );
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( flag );
  WALBERLA_ASSERT(false, "BubbleFunction::interpolate is not implemented!");
}

template< typename ValueType >
void BubbleFunction< ValueType >::assign_impl(const std::vector< ValueType > scalars,
                            const std::vector<BubbleFunction<ValueType> *> functions,
                            size_t level,
                            DoFType flag) {
    // Collect all source IDs in a vector
    std::vector<PrimitiveDataID<VertexBubbleFunctionMemory< ValueType >, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeBubbleFunctionMemory< ValueType >, Edge>> srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face>> srcFaceIDs;

    for (auto &function : functions) {
      srcVertexIDs.push_back(function->vertexDataID_);
      srcEdgeIDs.push_back(function->edgeDataID_);
      srcFaceIDs.push_back(function->faceDataID_);
    }

    for (auto &it : storage_->getFaces()) {
      Face &face = *it.second;

      if (testFlag(face.type, flag)) {
        BubbleFace::assign< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
      }
    }
}

template< typename ValueType >
void BubbleFunction< ValueType >::add_impl(const std::vector< ValueType > scalars,
                         const std::vector<BubbleFunction< ValueType > *> functions,
                         size_t level,
                         DoFType flag) {
    // Collect all source IDs in a vector
    std::vector<PrimitiveDataID<VertexBubbleFunctionMemory< ValueType >, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeBubbleFunctionMemory< ValueType >, Edge>> srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face>> srcFaceIDs;

    for (auto &function : functions) {
      srcVertexIDs.push_back(function->vertexDataID_);
      srcEdgeIDs.push_back(function->edgeDataID_);
      srcFaceIDs.push_back(function->faceDataID_);
    }

    for (auto &it : storage_->getFaces()) {
      Face &face = *it.second;

      if (testFlag(face.type, flag)) {
        BubbleFace::add< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
      }
    }
}

template< typename ValueType >
real_t BubbleFunction< ValueType >::dot_impl(BubbleFunction< ValueType > &rhs, uint_t level, DoFType flag) {
  real_t scalarProduct = 0.0;

  for (auto &it : storage_->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      scalarProduct += BubbleFace::dot< ValueType >(level, face, faceDataID_, rhs.faceDataID_);
    }
  }

  walberla::mpi::allReduceInplace(scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm());

  return scalarProduct;
}

template< typename ValueType >
void BubbleFunction< ValueType >::enumerate_impl(size_t level, uint_t& num)
{
  for (auto &it : storage_->getFaces()) {
    Face &face = *it.second;

    BubbleFace::enumerate< ValueType >(level, face, faceDataID_, num);
  }
  communicators_[level]->template startCommunication<Face, Edge>();
  communicators_[level]->template endCommunication<Face, Edge>();

  communicators_[level]->template startCommunication<Edge, Vertex>();
  communicators_[level]->template endCommunication<Edge, Vertex>();
}

}