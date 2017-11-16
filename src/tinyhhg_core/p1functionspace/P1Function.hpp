#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/communication/BufferedCommunication.hpp"

#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "tinyhhg_core/p1functionspace/P1Vertex.hpp"
#include "tinyhhg_core/p1functionspace/P1Edge.hpp"
#include "tinyhhg_core/p1functionspace/P1Face.hpp"
#include "tinyhhg_core/p1functionspace/P1DataHandling.hpp"
#include "tinyhhg_core/p1functionspace/P1PackInfo.hpp"

namespace hhg {



template< typename ValueType >
class P1Function : public Function< P1Function< ValueType > >
{
public:

  P1Function( const std::string& name, const std::shared_ptr< PrimitiveStorage > & storage, uint_t minLevel, uint_t maxLevel ) :
      Function< P1Function< ValueType > >( name, storage, minLevel, maxLevel )
  {
    auto faceP1FunctionMemoryDataHandling = std::make_shared< FaceP1FunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
    auto edgeP1FunctionMemoryDataHandling = std::make_shared< EdgeP1FunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
    auto vertexP1FunctionMemoryDataHandling = std::make_shared< VertexP1FunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
    storage->addFaceData( faceDataID_, faceP1FunctionMemoryDataHandling, name );
    storage->addEdgeData( edgeDataID_, edgeP1FunctionMemoryDataHandling, name );
    storage->addVertexData( vertexDataID_, vertexP1FunctionMemoryDataHandling, name );
    for ( uint_t level = minLevel; level <= maxLevel; ++level )
    {
      communicators_[level]->addPackInfo( std::make_shared< P1PackInfo< ValueType > >( level, vertexDataID_, edgeDataID_, faceDataID_, storage_ ) );
    }
  }

  const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &getVertexDataID() const { return vertexDataID_; }

  const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &getEdgeDataID() const { return edgeDataID_; }

  const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &getFaceDataID() const { return faceDataID_; }

  // TODO: split this function into impl
  inline void integrateDG(DGFunction< ValueType >& rhs, P1Function< ValueType >& rhsP1, uint_t level, DoFType flag);

  // TODO: write more general version
  inline real_t getMaxValue(uint_t level);

private:

  using Function< P1Function< ValueType > >::storage_;
  using Function< P1Function< ValueType > >::communicators_;

  /// Interpolates a given expression to a P1Function
  inline void interpolate_impl(std::function<ValueType(const Point3D&)>& expr, uint_t level, DoFType flag = All);

  inline void assign_impl(const std::vector<ValueType> scalars, const std::vector<P1Function< ValueType >*> functions, uint_t level, DoFType flag = All);

  inline void add_impl(const std::vector<ValueType> scalars, const std::vector<P1Function< ValueType >*> functions, uint_t level, DoFType flag = All);

  inline real_t dot_impl(P1Function< ValueType >& rhs, uint_t level, DoFType flag = All);

  inline void prolongate_impl(uint_t sourceLevel, DoFType flag = All);

  inline void prolongateQuadratic_impl(uint_t sourceLevel, DoFType flag = All);

  inline void restrict_impl(uint_t sourceLevel, DoFType flag = All);

  inline void enumerate_impl(uint_t level, uint_t& num);

  PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> vertexDataID_;
  PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> edgeDataID_;
  PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> faceDataID_;
};

template< typename ValueType >
inline void P1Function< ValueType >::interpolate_impl(std::function< ValueType(const hhg::Point3D&) > & expr, uint_t level, DoFType flag)
{
    for (auto& it : storage_->getVertices()) {
        Vertex& vertex = *it.second;

        if (testFlag(vertex.getDoFType(), flag)) {
            P1Vertex::interpolate(vertex, vertexDataID_, expr, level);
        }
    }

    communicators_[level]->template startCommunication<Vertex, Edge>();

    for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        if (testFlag(edge.getDoFType(), flag)) {
            P1Edge::interpolate< ValueType >(level, edge, edgeDataID_, expr);
        }
    }

    communicators_[level]->template endCommunication<Vertex, Edge>();
    communicators_[level]->template startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        if (testFlag(face.type, flag)) {
            P1Face::interpolate< ValueType >(level, face, faceDataID_, expr);
        }
    }

    communicators_[level]->template endCommunication<Edge, Face>();
}

template< typename ValueType >
inline void P1Function< ValueType >::assign_impl(const std::vector<ValueType> scalars, const std::vector<P1Function< ValueType >*> functions, size_t level, DoFType flag)
{
    // Collect all source IDs in a vector
    std::vector<PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex>> srcVertexIDs;
    std::vector<PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge>>     srcEdgeIDs;
    std::vector<PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>>     srcFaceIDs;

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

    communicators_[level]->template startCommunication<Vertex, Edge>();

    for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        if (testFlag(edge.getDoFType(), flag)) {
            P1Edge::assign< ValueType >(level, edge, scalars, srcEdgeIDs, edgeDataID_);
        }
    }

    communicators_[level]->template endCommunication<Vertex, Edge>();
    communicators_[level]->template startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        if (testFlag(face.type, flag)) {
            P1Face::assign< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
        }
    }

    communicators_[level]->template endCommunication<Edge, Face>();
}

template< typename ValueType >
inline void P1Function< ValueType >::add_impl(const std::vector<ValueType> scalars, const std::vector<P1Function< ValueType >*> functions, size_t level, DoFType flag)
{
  // Collect all source IDs in a vector
  std::vector<PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex>> srcVertexIDs;
  std::vector<PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge>>     srcEdgeIDs;
  std::vector<PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>>     srcFaceIDs;

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

  communicators_[level]->template startCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag)) {
          P1Edge::add< ValueType >(level, edge, scalars, srcEdgeIDs, edgeDataID_);
      }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag)) {
          P1Face::add< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
      }
  }

  communicators_[level]->template endCommunication<Edge, Face>();
}

template< typename ValueType >
inline real_t P1Function< ValueType >::dot_impl(P1Function< ValueType >& rhs, size_t level, DoFType flag)
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
        scalarProduct += P1Edge::dot< ValueType >(level, edge, edgeDataID_, rhs.edgeDataID_);
      }
  }

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag)) {
        scalarProduct += P1Face::dot< ValueType >(level, face, faceDataID_, rhs.faceDataID_);
      }
  }

  walberla::mpi::allReduceInplace( scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );

  return scalarProduct;
}

template< typename ValueType >
inline void P1Function< ValueType >::prolongate_impl(size_t sourceLevel, DoFType flag)
{
  const size_t destinationLevel = sourceLevel + 1;

  for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        P1Vertex::prolongate(vertex, vertexDataID_, sourceLevel);
      }
  }

  communicators_[destinationLevel]->template startCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::prolongate< ValueType >(sourceLevel, edge, edgeDataID_);
      }
  }

  communicators_[destinationLevel]->template endCommunication<Vertex, Edge>();
  communicators_[destinationLevel]->template startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::prolongate< ValueType >(sourceLevel, face, faceDataID_);
      }
  }

  communicators_[destinationLevel]->template endCommunication<Edge, Face>();
}

template< typename ValueType >
inline void P1Function< ValueType >::prolongateQuadratic_impl(size_t sourceLevel, DoFType flag)
{
  const size_t destinationLevel = sourceLevel + 1;

  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      P1Vertex::prolongateQuadratic(vertex, vertexDataID_, sourceLevel);
    }
  }

  communicators_[destinationLevel]->template startCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      P1Edge::prolongateQuadratic< ValueType >(sourceLevel, edge, edgeDataID_);
    }
  }

  communicators_[destinationLevel]->template endCommunication<Vertex, Edge>();
  communicators_[destinationLevel]->template startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      P1Face::prolongateQuadratic< ValueType >(sourceLevel, face, faceDataID_);
    }
  }

  communicators_[destinationLevel]->template endCommunication<Edge, Face>();
}

template< typename ValueType >
inline void P1Function< ValueType >::restrict_impl(size_t sourceLevel, DoFType flag)
{
  const size_t destinationLevel = sourceLevel - 1;

  // start pulling vertex halos
  communicators_[sourceLevel]->template startCommunication<Edge, Vertex>();

  // start pulling edge halos
  communicators_[sourceLevel]->template startCommunication<Face, Edge>();

  // end pulling vertex halos
  communicators_[sourceLevel]->template endCommunication<Edge, Vertex>();

  for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        P1Vertex::restrict(vertex, vertexDataID_, sourceLevel);
      }
  }

  communicators_[destinationLevel]->template startCommunication<Vertex, Edge>();

  // end pulling edge halos
  communicators_[sourceLevel]->template endCommunication<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::restrict< ValueType >(sourceLevel, edge, edgeDataID_);
      }
  }

  communicators_[destinationLevel]->template endCommunication<Vertex, Edge>();

  communicators_[destinationLevel]->template startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::restrict< ValueType >(sourceLevel, face, faceDataID_);
      }
  }

  communicators_[destinationLevel]->template endCommunication<Edge, Face>();

}

template< typename ValueType >
inline void P1Function< ValueType >::enumerate_impl(uint_t level, uint_t& num)
{
  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;
    P1Vertex::enumerate(level, vertex, vertexDataID_, num);
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();
  communicators_[level]->template endCommunication<Vertex, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;
    P1Edge::enumerate< ValueType >(level, edge, edgeDataID_, num);
  }

  communicators_[level]->template startCommunication<Edge, Face>();
  communicators_[level]->template endCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;
    P1Face::enumerate< ValueType >(level, face, faceDataID_, num);
  }

  communicators_[level]->template startCommunication<Face, Edge>();
  communicators_[level]->template endCommunication<Face, Edge>();

  communicators_[level]->template startCommunication<Edge, Vertex>();
  communicators_[level]->template endCommunication<Edge, Vertex>();
}

template< typename ValueType >
inline void P1Function< ValueType >::integrateDG(DGFunction< ValueType >& rhs, P1Function< ValueType >& rhsP1, uint_t level, DoFType flag)
{
  this->startTiming( "integrateDG" );

  rhsP1.getCommunicator(level)->template startCommunication<Edge, Vertex>();
  rhsP1.getCommunicator(level)->template startCommunication<Face, Edge>();

  rhs.getCommunicator(level)->template startCommunication<Face, Edge>();
  rhs.getCommunicator(level)->template endCommunication<Face, Edge>();

  rhs.getCommunicator(level)->template startCommunication<Edge, Vertex>();
  rhs.getCommunicator(level)->template endCommunication<Edge, Vertex>();

  rhsP1.getCommunicator(level)->template endCommunication<Edge, Vertex>();

  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag)) {
      P1Vertex::integrateDG< ValueType >(vertex, storage_, rhs.getVertexDataID(), rhsP1.getVertexDataID(), vertexDataID_, level);
    }
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();
  rhsP1.getCommunicator(level)->template endCommunication<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag)) {
      P1Edge::integrateDG< ValueType >(level, edge, storage_, rhs.getEdgeDataID(), rhsP1.getEdgeDataID(), edgeDataID_);
    }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag)) {
      P1Face::integrateDG< ValueType >(level, face, rhs.getFaceDataID(), rhsP1.getFaceDataID(), faceDataID_);
    }
  }

  communicators_[level]->template endCommunication<Edge, Face>();

  this->stopTiming( "integrateDG" );
}

inline void projectMean(P1Function<real_t>& pressure, hhg::P1Function<real_t>& tmp, uint_t level) {

  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) {
    return 1.0;
  };

  tmp.interpolate(ones, level);

  real_t numGlobalVertices = tmp.dot(tmp, level, hhg::All);
  real_t mean = pressure.dot(tmp, level, hhg::All);

  pressure.assign({1.0, -mean/numGlobalVertices}, {&pressure, &tmp}, level, hhg::All);
}

template< typename ValueType >
inline real_t P1Function< ValueType >::getMaxValue(uint_t level)
{
  communicators_[level]->template startCommunication<Vertex, Edge>();
  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();
  communicators_[level]->template endCommunication<Edge, Face>();

  real_t localMax = std::numeric_limits<real_t>::min();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;
    localMax = std::max(localMax, P1Face::getMaxValue< ValueType >(level, face, faceDataID_));
  }

  real_t globalMax = walberla::mpi::allReduce(localMax, walberla::mpi::MAX);

  return globalMax;
}

}
