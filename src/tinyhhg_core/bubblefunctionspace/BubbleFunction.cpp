#include "BubbleFunction.hpp"
#include "BubbleDataHandling.hpp"
#include "BubbleEdge.hpp"
#include "BubbleFace.hpp"
#include "BubblePackInfo.hpp"
#include "BubbleVertex.hpp"

namespace hhg {

BubbleFunction::BubbleFunction(const std::string &name,
                               const std::shared_ptr<PrimitiveStorage> &storage,
                               uint_t minLevel,
                               uint_t maxLevel)
    : Function(name, storage, minLevel, maxLevel) {
  FaceBubbleFunctionMemoryDataHandling faceBubbleFunctionMemoryDataHandling(minLevel, maxLevel);
  EdgeBubbleFunctionMemoryDataHandling edgeBubbleFunctionMemoryDataHandling(minLevel, maxLevel);
  VertexBubbleFunctionMemoryDataHandling vertexBubbleFunctionMemoryDataHandling(minLevel, maxLevel);
  faceDataID_ = storage->addFaceData(faceBubbleFunctionMemoryDataHandling, name);
  edgeDataID_ = storage->addEdgeData(edgeBubbleFunctionMemoryDataHandling, name);
  vertexDataID_ = storage->addVertexData(vertexBubbleFunctionMemoryDataHandling, name);
  for (uint_t level = minLevel; level <= maxLevel; ++level) {
    //    communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
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

void BubbleFunction::interpolate(std::function<real_t(const hhg::Point3D &)> &expr, uint_t level, DoFType flag) {
  // TODO: implement Bubble interpolation. It is not required for Dirichlet only interpolation in most of the apps
  WALBERLA_ASSERT(false, "BubbleFunction::interpolate is not implemented!")
}

void BubbleFunction::assign(const std::vector<walberla::real_t> scalars,
                            const std::vector<BubbleFunction *> functions,
                            size_t level,
                            DoFType flag) {
  //  // Collect all source IDs in a vector
  //  std::vector<PrimitiveDataID<VertexP1FunctionMemory, Vertex>> srcVertexIDs;
  //  std::vector<PrimitiveDataID<EdgeP1FunctionMemory, Edge>> srcEdgeIDs;
  //  std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>> srcFaceIDs;
  //
  //  for (auto &function : functions) {
  //    srcVertexIDs.push_back(function->vertexDataID_);
  //    srcEdgeIDs.push_back(function->edgeDataID_);
  //    srcFaceIDs.push_back(function->faceDataID_);
  //  }
  //
  //  for (auto &it : storage_->getVertices()) {
  //    Vertex &vertex = *it.second;
  //
  //    if (testFlag(vertex.type, flag)) {
  //      P1Vertex::assign(vertex, scalars, srcVertexIDs, vertexDataID_, level);
  //    }
  //  }
  //
  //  communicators_[level]->startCommunication<Vertex, Edge>();
  //
  //  for (auto &it : storage_->getEdges()) {
  //    Edge &edge = *it.second;
  //
  //    if (testFlag(edge.type, flag)) {
  //      P1Edge::assign(edge, scalars, srcEdgeIDs, edgeDataID_, level);
  //    }
  //  }
  //
  //  communicators_[level]->endCommunication<Vertex, Edge>();
  //  communicators_[level]->startCommunication<Edge, Face>();
  //
  //  for (auto &it : storage_->getFaces()) {
  //    Face &face = *it.second;
  //
  //    if (testFlag(face.type, flag)) {
  //      P1Face::assign(face, scalars, srcFaceIDs, faceDataID_, level);
  //    }
  //  }
  //
  //  communicators_[level]->endCommunication<Edge, Face>();
}

void BubbleFunction::add(const std::vector<walberla::real_t> scalars,
                         const std::vector<BubbleFunction *> functions,
                         size_t level,
                         DoFType flag) {
  //  // Collect all source IDs in a vector
  //  std::vector<PrimitiveDataID<VertexP1FunctionMemory, Vertex>> srcVertexIDs;
  //  std::vector<PrimitiveDataID<EdgeP1FunctionMemory, Edge>> srcEdgeIDs;
  //  std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>> srcFaceIDs;
  //
  //  for (auto &function : functions) {
  //    srcVertexIDs.push_back(function->vertexDataID_);
  //    srcEdgeIDs.push_back(function->edgeDataID_);
  //    srcFaceIDs.push_back(function->faceDataID_);
  //  }
  //
  //  for (auto &it : storage_->getVertices()) {
  //    Vertex &vertex = *it.second;
  //
  //    if (testFlag(vertex.type, flag)) {
  //      P1Vertex::add(vertex, scalars, srcVertexIDs, vertexDataID_, level);
  //    }
  //  }
  //
  //  communicators_[level]->startCommunication<Vertex, Edge>();
  //
  //  for (auto &it : storage_->getEdges()) {
  //    Edge &edge = *it.second;
  //
  //    if (testFlag(edge.type, flag)) {
  //      P1Edge::add(edge, scalars, srcEdgeIDs, edgeDataID_, level);
  //    }
  //  }
  //
  //  communicators_[level]->endCommunication<Vertex, Edge>();
  //  communicators_[level]->startCommunication<Edge, Face>();
  //
  //  for (auto &it : storage_->getFaces()) {
  //    Face &face = *it.second;
  //
  //    if (testFlag(face.type, flag)) {
  //      P1Face::add(face, scalars, srcFaceIDs, faceDataID_, level);
  //    }
  //  }
  //
  //  communicators_[level]->endCommunication<Edge, Face>();
}

real_t BubbleFunction::dot(BubbleFunction &rhs, uint_t level, DoFType flag) {
  real_t scalarProduct = 0.0;

  for (auto &it : storage_->getVertices()) {
    Vertex &vertex = *it.second;

    if (testFlag(vertex.type, flag)) {
      scalarProduct += BubbleVertex::dot(vertex, vertexDataID_, rhs.vertexDataID_, level);
    }
  }

  for (auto &it : storage_->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.type, flag)) {
      scalarProduct += BubbleEdge::dot(edge, edgeDataID_, rhs.edgeDataID_, level);
    }
  }

  for (
    auto &it :
      storage_->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      scalarProduct += BubbleFace::dot(face, faceDataID_, rhs.faceDataID_, level);
    }
  }

  walberla::mpi::allReduceInplace(scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm());

  return scalarProduct;
}
}
