#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/bubblefunctionspace/BubbleMemory.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFace.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleDataHandling.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubblePackInfo.hpp"

namespace hhg {


class BubbleFunction : public Function< BubbleFunction > {
 public:
  BubbleFunction( const std::string &name, const std::shared_ptr< PrimitiveStorage > &storage, uint_t minLevel, uint_t maxLevel ) :
      Function( name, storage, minLevel, maxLevel )
  {
    auto faceBubbleFunctionMemoryDataHandling = std::make_shared< FaceBubbleFunctionMemoryDataHandling >( minLevel, maxLevel );
    auto edgeBubbleFunctionMemoryDataHandling = std::make_shared< EdgeBubbleFunctionMemoryDataHandling >( minLevel, maxLevel );
    auto vertexBubbleFunctionMemoryDataHandling = std::make_shared< VertexBubbleFunctionMemoryDataHandling >( minLevel, maxLevel );
    storage->addFaceData( faceDataID_, faceBubbleFunctionMemoryDataHandling, name );
    storage->addEdgeData( edgeDataID_, edgeBubbleFunctionMemoryDataHandling, name );
    storage->addVertexData( vertexDataID_, vertexBubbleFunctionMemoryDataHandling, name );
    for ( uint_t level = minLevel; level <= maxLevel; ++level )
    {
      communicators_[level]->setLocalCommunicationMode( communication::BufferedCommunicator::BUFFERED_MPI );
      communicators_[level]->addPackInfo( std::make_shared< BubblePackInfo >( level, vertexDataID_, edgeDataID_, faceDataID_, storage_ ) );
    }
  }

  const PrimitiveDataID<VertexBubbleFunctionMemory< real_t >, Vertex> &getVertexDataID() const { return vertexDataID_; }

  const PrimitiveDataID<EdgeBubbleFunctionMemory< real_t >, Edge> &getEdgeDataID() const { return edgeDataID_; }

  const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &getFaceDataID() const { return faceDataID_; }

 private:

  /// Interpolates a given expression to a P1Function
  void interpolate_impl(std::function<real_t(const Point3D &)> &expr, uint_t level, DoFType flag = All);

  void assign_impl(const std::vector<walberla::real_t> scalars,
              const std::vector<BubbleFunction *> functions,
              size_t level,
              DoFType flag = All);

  void add_impl(const std::vector<walberla::real_t> scalars,
           const std::vector<BubbleFunction *> functions,
           size_t level,
           DoFType flag = All);

  real_t dot_impl(BubbleFunction &rhs, size_t level, DoFType flag = All);

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

  PrimitiveDataID<VertexBubbleFunctionMemory< real_t >, Vertex> vertexDataID_;
  PrimitiveDataID<EdgeBubbleFunctionMemory< real_t >, Edge> edgeDataID_;
  PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> faceDataID_;
};

void BubbleFunction::interpolate_impl(std::function<real_t(const hhg::Point3D &)> &expr, uint_t level, DoFType flag) {
  // TODO: implement Bubble interpolation. It is not required for Dirichlet only interpolation in most of the apps
  WALBERLA_UNUSED( expr );
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( flag );
  WALBERLA_ASSERT(false, "BubbleFunction::interpolate is not implemented!");
}

void BubbleFunction::assign_impl(const std::vector<walberla::real_t> scalars,
                            const std::vector<BubbleFunction *> functions,
                            size_t level,
                            DoFType flag) {
    // Collect all source IDs in a vector
    std::vector<PrimitiveDataID<VertexBubbleFunctionMemory< real_t >, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeBubbleFunctionMemory< real_t >, Edge>> srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face>> srcFaceIDs;

    for (auto &function : functions) {
      srcVertexIDs.push_back(function->vertexDataID_);
      srcEdgeIDs.push_back(function->edgeDataID_);
      srcFaceIDs.push_back(function->faceDataID_);
    }

    for (auto &it : storage_->getFaces()) {
      Face &face = *it.second;

      if (testFlag(face.type, flag)) {
        BubbleFace::assign(level, face, scalars, srcFaceIDs, faceDataID_);
      }
    }
}

void BubbleFunction::add_impl(const std::vector<walberla::real_t> scalars,
                         const std::vector<BubbleFunction *> functions,
                         size_t level,
                         DoFType flag) {
    // Collect all source IDs in a vector
    std::vector<PrimitiveDataID<VertexBubbleFunctionMemory< real_t >, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeBubbleFunctionMemory< real_t >, Edge>> srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face>> srcFaceIDs;

    for (auto &function : functions) {
      srcVertexIDs.push_back(function->vertexDataID_);
      srcEdgeIDs.push_back(function->edgeDataID_);
      srcFaceIDs.push_back(function->faceDataID_);
    }

    for (auto &it : storage_->getFaces()) {
      Face &face = *it.second;

      if (testFlag(face.type, flag)) {
        BubbleFace::add(level, face, scalars, srcFaceIDs, faceDataID_);
      }
    }
}

real_t BubbleFunction::dot_impl(BubbleFunction &rhs, uint_t level, DoFType flag) {
  real_t scalarProduct = 0.0;

  for (auto &it : storage_->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      scalarProduct += BubbleFace::dot(level, face, faceDataID_, rhs.faceDataID_);
    }
  }

  walberla::mpi::allReduceInplace(scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm());

  return scalarProduct;
}

void BubbleFunction::enumerate_impl(size_t level, uint_t& num)
{
  for (auto &it : storage_->getFaces()) {
    Face &face = *it.second;

    BubbleFace::enumerate(level, face, faceDataID_, num);
  }
  communicators_[level]->startCommunication<Face, Edge>();
  communicators_[level]->endCommunication<Face, Edge>();

  communicators_[level]->startCommunication<Edge, Vertex>();
  communicators_[level]->endCommunication<Edge, Vertex>();
}

}
