#include "P1Function.hpp"
#include "P1DataHandling.hpp"
#include "P1Vertex.hpp"
#include "P1Edge.hpp"
#include "P1Face.hpp"
#include "P1PackInfo.hpp"

namespace hhg {

P1Function::P1Function(const std::string& name, const std::shared_ptr< PrimitiveStorage > & storage, uint_t minLevel, uint_t maxLevel)
    : Function(name, storage, minLevel, maxLevel)
{
    auto faceP1FunctionMemoryDataHandling = std::make_shared< FaceP1FunctionMemoryDataHandling >(minLevel, maxLevel);
    auto edgeP1FunctionMemoryDataHandling = std::make_shared< EdgeP1FunctionMemoryDataHandling >(minLevel, maxLevel);
    auto vertexP1FunctionMemoryDataHandling = std::make_shared< VertexP1FunctionMemoryDataHandling >(minLevel, maxLevel);
    storage->addFaceData(faceDataID_, faceP1FunctionMemoryDataHandling, name);
    storage->addEdgeData(edgeDataID_, edgeP1FunctionMemoryDataHandling, name);
    storage->addVertexData(vertexDataID_, vertexP1FunctionMemoryDataHandling, name);
  for(uint_t level = minLevel; level <= maxLevel; ++level){
//    communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
    communicators_[level]->addPackInfo(std::make_shared<P1PackInfo>(level,vertexDataID_,edgeDataID_,faceDataID_,storage_));
  }
}

P1Function::~P1Function()
{
    //TODO implement!
}

void P1Function::interpolate_impl(std::function<real_t(const hhg::Point3D&)>& expr, uint_t level, DoFType flag)
{
    for (auto& it : storage_->getVertices()) {
        Vertex& vertex = *it.second;

        if (testFlag(vertex.getDoFType(), flag)) {
            P1Vertex::interpolate(vertex, vertexDataID_, expr, level);
        }
    }

    communicators_[level]->startCommunication<Vertex, Edge>();

    for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        if (testFlag(edge.getDoFType(), flag)) {
            P1Edge::interpolate(level, edge, edgeDataID_, expr);
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

void P1Function::assign_impl(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag)
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

        if (testFlag(vertex.getDoFType(), flag)) {
            P1Vertex::assign(vertex, scalars, srcVertexIDs, vertexDataID_, level);
        }
    }

    communicators_[level]->startCommunication<Vertex, Edge>();

    for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        if (testFlag(edge.getDoFType(), flag)) {
            P1Edge::assign(level, edge, scalars, srcEdgeIDs, edgeDataID_);
        }
    }

    communicators_[level]->endCommunication<Vertex, Edge>();
    communicators_[level]->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        if (testFlag(face.type, flag)) {
            P1Face::assign(level, face, scalars, srcFaceIDs, faceDataID_);
        }
    }

    communicators_[level]->endCommunication<Edge, Face>();
}

void P1Function::add_impl(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag)
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

      if (testFlag(vertex.getDoFType(), flag)) {
          P1Vertex::add(vertex, scalars, srcVertexIDs, vertexDataID_, level);
      }
  }

  communicators_[level]->startCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag)) {
          P1Edge::add(level, edge, scalars, srcEdgeIDs, edgeDataID_);
      }
  }

  communicators_[level]->endCommunication<Vertex, Edge>();
  communicators_[level]->startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag)) {
          P1Face::add(level, face, scalars, srcFaceIDs, faceDataID_);
      }
  }

  communicators_[level]->endCommunication<Edge, Face>();
}

real_t P1Function::dot_impl(P1Function& rhs, size_t level, DoFType flag)
{
  real_t scalarProduct = 0.0;

  for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag)) {
        scalarProduct += P1Vertex::dot(vertex, vertexDataID_, rhs.vertexDataID_, level);
      }
  }

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag)) {
        scalarProduct += P1Edge::dot(level, edge, edgeDataID_, rhs.edgeDataID_);
      }
  }

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag)) {
        scalarProduct += P1Face::dot(level, face, faceDataID_, rhs.faceDataID_);
      }
  }

  walberla::mpi::allReduceInplace( scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );

  return scalarProduct;
}

void P1Function::prolongate_impl(size_t sourceLevel, DoFType flag)
{
  const size_t destinationLevel = sourceLevel + 1;

  for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        P1Vertex::prolongate(vertex, vertexDataID_, sourceLevel);
      }
  }

  communicators_[destinationLevel]->startCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::prolongate(sourceLevel, edge, edgeDataID_);
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

void P1Function::prolongateQuadratic_impl(size_t sourceLevel, DoFType flag)
{
  const size_t destinationLevel = sourceLevel + 1;

  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      P1Vertex::prolongateQuadratic(vertex, vertexDataID_, sourceLevel);
    }
  }

  communicators_[destinationLevel]->startCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      P1Edge::prolongateQuadratic(sourceLevel, edge, edgeDataID_);
    }
  }

  communicators_[destinationLevel]->endCommunication<Vertex, Edge>();
  communicators_[destinationLevel]->startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      P1Face::prolongateQuadratic(sourceLevel, face, faceDataID_);
    }
  }

  communicators_[destinationLevel]->endCommunication<Edge, Face>();
}

void P1Function::restrict_impl(size_t sourceLevel, DoFType flag)
{
  const size_t destinationLevel = sourceLevel - 1;

  // start pulling vertex halos
  communicators_[sourceLevel]->startCommunication<Edge, Vertex>();

  // start pulling edge halos
  communicators_[sourceLevel]->startCommunication<Face, Edge>();

  // end pulling vertex halos
  communicators_[sourceLevel]->endCommunication<Edge, Vertex>();

  for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        P1Vertex::restrict(vertex, vertexDataID_, sourceLevel);
      }
  }

  communicators_[destinationLevel]->startCommunication<Vertex, Edge>();

  // end pulling edge halos
  communicators_[sourceLevel]->endCommunication<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::restrict(sourceLevel, edge, edgeDataID_);
      }
  }

  communicators_[destinationLevel]->endCommunication<Vertex, Edge>();

  communicators_[destinationLevel]->startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::restrict(sourceLevel, face, faceDataID_);
      }
  }

  communicators_[destinationLevel]->endCommunication<Edge, Face>();

}

void P1Function::enumerate_impl(uint_t level, uint_t& num)
{
  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;
    P1Vertex::enumerate(vertex, vertexDataID_, level, num);
  }

  communicators_[level]->startCommunication<Vertex, Edge>();
  communicators_[level]->endCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;
    P1Edge::enumerate(level, edge, edgeDataID_, num);
  }

  communicators_[level]->startCommunication<Edge, Face>();
  communicators_[level]->endCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;
    P1Face::enumerate(face, faceDataID_, level, num);
  }

  //communicators_[level]->startCommunication<Face, Edge>();
  //communicators_[level]->endCommunication<Face, Edge>();

  //communicators_[level]->startCommunication<Edge, Vertex>();
  //communicators_[level]->endCommunication<Edge, Vertex>();
}


void P1Function::createVectorFromFunction_impl(P1Function &numerator,Vec &vec, uint_t level,DoFType flag)
{
  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      P1Vertex::createVectorFromFunction(vertex, vertexDataID_, numerator.getVertexDataID(), vec, level);
    }
  }


  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      P1Edge::createVectorFromFunction(level, edge, edgeDataID_, numerator.getEdgeDataID(), vec);
    }
  }


  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      P1Face::createVectorFromFunction(level, face, faceDataID_, numerator.getFaceDataID(), vec);
    }
  }
}


void P1Function::createFunctionFromVector_impl(P1Function &numerator,Vec &vec, uint_t level,DoFType flag)
{
  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      P1Vertex::createFunctionFromVector(vertex, vertexDataID_, numerator.getVertexDataID(), vec, level);
    }
  }

  communicators_[level]->startCommunication<Vertex, Edge>();
  communicators_[level]->endCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      P1Edge::createFunctionFromVector(level, edge, edgeDataID_, numerator.getEdgeDataID(), vec);
    }
  }

  communicators_[level]->startCommunication<Edge, Face>();
  communicators_[level]->endCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      P1Face::createFunctionFromVector(level, face, faceDataID_, numerator.getFaceDataID(), vec);
    }
  }
}

}
