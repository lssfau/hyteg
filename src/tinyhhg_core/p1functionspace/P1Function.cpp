#include "P1Function.hpp"
#include "P1DataHandling.hpp"
#include "p1vertex.hpp"
#include "p1edge.hpp"
#include "p1face.hpp"
#include "P1PackInfo.hpp"

namespace hhg {

P1Function::P1Function(const std::string& name, const std::shared_ptr< PrimitiveStorage > & storage, uint_t minLevel, uint_t maxLevel)
    : Function(name, storage, minLevel, maxLevel)
{
    FaceP1FunctionMemoryDataHandling faceP1FunctionMemoryDataHandling(minLevel, maxLevel);
    EdgeP1FunctionMemoryDataHandling edgeP1FunctionMemoryDataHandling(minLevel, maxLevel);
    VertexP1FunctionMemoryDataHandling vertexP1FunctionMemoryDataHandling(minLevel, maxLevel);
    faceDataID_ = storage->addFaceData(faceP1FunctionMemoryDataHandling, name);
    edgeDataID_ = storage->addEdgeData(edgeP1FunctionMemoryDataHandling, name);
    vertexDataID_ = storage->addVertexData(vertexP1FunctionMemoryDataHandling, name);
  for(uint_t i = minLevel; i <= maxLevel; ++i){
    communicators_[i]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
    communicators_[i]->addPackInfo(std::make_shared<P1PackInfo>(i,vertexDataID_,edgeDataID_,faceDataID_,storage_));
  }
}

P1Function::~P1Function()
{
    //TODO implement!
}

void P1Function::interpolate(std::function<real_t(const hhg::Point3D&)>& expr, uint_t level, DoFType flag)
{
    for (auto& it : storage_->getVertices()) {
        Vertex& vertex = *it.second;

        if (testFlag(vertex.type, flag)) {
            P1Vertex::interpolate(vertex, vertexDataID_, expr, level);
        }
    }

    communicators_[level]->startCommunication<Vertex, Edge>();

    for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        if (testFlag(edge.type, flag)) {
            P1Edge::interpolate(edge, edgeDataID_, expr, level);
        }
    }

    communicators_[level]->endCommunication<Vertex, Edge>();
    communicators_[level]->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        if (testFlag(face.type, flag)) {
            P1Face::interpolate(level, face, faceDataID_, expr);
        }
    }

    communicators_[level]->endCommunication<Edge, Face>();
}

void P1Function::assign(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag)
{
    // Collect all source IDs in a vector
    std::vector<PrimitiveDataID<VertexP1FunctionMemory, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeP1FunctionMemory, Edge>>     srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>>     srcFaceIDs;

    for (auto& function : functions)
    {
        srcVertexIDs.push_back(function->vertexDataID_);
        srcEdgeIDs.push_back(function->edgeDataID_);
        srcFaceIDs.push_back(function->faceDataID_);
    }

    for (auto& it : storage_->getVertices()) {
        Vertex& vertex = *it.second;

        if (testFlag(vertex.type, flag)) {
            P1Vertex::assign(vertex, scalars, srcVertexIDs, vertexDataID_, level);
        }
    }

    communicators_[level]->startCommunication<Vertex, Edge>();

    for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        if (testFlag(edge.type, flag)) {
            P1Edge::assign(edge, scalars, srcEdgeIDs, edgeDataID_, level);
        }
    }

    communicators_[level]->endCommunication<Vertex, Edge>();
    communicators_[level]->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        if (testFlag(face.type, flag)) {
            P1Face::assign(face, scalars, srcFaceIDs, faceDataID_, level);
        }
    }

    communicators_[level]->endCommunication<Edge, Face>();
}

void P1Function::add(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag)
{
  // Collect all source IDs in a vector
  std::vector<PrimitiveDataID<VertexP1FunctionMemory, Vertex>> srcVertexIDs;
  std::vector<PrimitiveDataID<EdgeP1FunctionMemory, Edge>>     srcEdgeIDs;
  std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>>     srcFaceIDs;

  for (auto& function : functions)
  {
      srcVertexIDs.push_back(function->vertexDataID_);
      srcEdgeIDs.push_back(function->edgeDataID_);
      srcFaceIDs.push_back(function->faceDataID_);
  }

  for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.type, flag)) {
          P1Vertex::add(vertex, scalars, srcVertexIDs, vertexDataID_, level);
      }
  }

  communicators_[level]->startCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.type, flag)) {
          P1Edge::add(edge, scalars, srcEdgeIDs, edgeDataID_, level);
      }
  }

  communicators_[level]->endCommunication<Vertex, Edge>();
  communicators_[level]->startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag)) {
          P1Face::add(face, scalars, srcFaceIDs, faceDataID_, level);
      }
  }

  communicators_[level]->endCommunication<Edge, Face>();
}

real_t P1Function::dot(P1Function& rhs, size_t level, DoFType flag)
{
  real_t scalarProduct = 0.0;

  for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.type, flag)) {
        scalarProduct += P1Vertex::dot(vertex, vertexDataID_, rhs.vertexDataID_, level);
      }
  }

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.type, flag)) {
        scalarProduct += P1Edge::dot(edge, edgeDataID_, rhs.edgeDataID_, level);
      }
  }

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag)) {
        scalarProduct += P1Face::dot(face, faceDataID_, rhs.faceDataID_, level);
      }
  }

  walberla::mpi::allReduceInplace( scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );

  return scalarProduct;
}

void P1Function::prolongate(size_t sourceLevel, DoFType flag)
{
  const size_t destinationLevel = sourceLevel + 1;

  for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.type, flag))
      {
        P1Vertex::prolongate(vertex, vertexDataID_, sourceLevel);
      }
  }

  communicators_[destinationLevel]->startCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.type, flag))
      {
        P1Edge::prolongate(edge, edgeDataID_, sourceLevel);
      }
  }

  communicators_[destinationLevel]->endCommunication<Vertex, Edge>();
  communicators_[destinationLevel]->startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::prolongate(sourceLevel, face, faceDataID_);
      }
  }

  communicators_[destinationLevel]->endCommunication<Edge, Face>();
}

void P1Function::prolongateQuadratic(size_t level, DoFType flag){

}

void P1Function::restrict(size_t level, DoFType flag)
{

}
}
