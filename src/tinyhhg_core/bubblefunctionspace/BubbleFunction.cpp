#include "BubbleFunction.hpp"
#include "BubbleDataHandling.hpp"
#include "BubbleFace.hpp"
#include "BubblePackInfo.hpp"

namespace hhg {

BubbleFunction::BubbleFunction(const std::string &name,
                               const std::shared_ptr<PrimitiveStorage> &storage,
                               uint_t minLevel,
                               uint_t maxLevel)
    : Function(name, storage, minLevel, maxLevel) {
  auto faceBubbleFunctionMemoryDataHandling = std::make_shared< FaceBubbleFunctionMemoryDataHandling >(minLevel, maxLevel);
  auto edgeBubbleFunctionMemoryDataHandling = std::make_shared< EdgeBubbleFunctionMemoryDataHandling >(minLevel, maxLevel);
  auto vertexBubbleFunctionMemoryDataHandling = std::make_shared< VertexBubbleFunctionMemoryDataHandling >(minLevel, maxLevel);
  storage->addFaceData(faceDataID_, faceBubbleFunctionMemoryDataHandling, name);
  storage->addEdgeData(edgeDataID_, edgeBubbleFunctionMemoryDataHandling, name);
  storage->addVertexData(vertexDataID_, vertexBubbleFunctionMemoryDataHandling, name);
  for (uint_t level = minLevel; level <= maxLevel; ++level) {
    communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
    communicators_[level]->addPackInfo(std::make_shared<BubblePackInfo>(level,
                                                                    vertexDataID_,
                                                                    edgeDataID_,
                                                                    faceDataID_,
                                                                    storage_));
  }
}

BubbleFunction::~BubbleFunction() {
  //TODO implement!
}

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
    std::vector<PrimitiveDataID<VertexBubbleFunctionMemory, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeBubbleFunctionMemory, Edge>> srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceBubbleFunctionMemory, Face>> srcFaceIDs;

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
    std::vector<PrimitiveDataID<VertexBubbleFunctionMemory, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeBubbleFunctionMemory, Edge>> srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceBubbleFunctionMemory, Face>> srcFaceIDs;

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
