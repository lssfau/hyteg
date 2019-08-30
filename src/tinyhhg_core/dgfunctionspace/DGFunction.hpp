#pragma once

#include "tinyhhg_core/Function.hpp"
#include "DGPackInfo.hpp"
#include "DGMemory.hpp"
#include "DGDataHandling.hpp"
#include "DGVertex.hpp"
#include "DGEdge.hpp"
#include "DGFace.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/boundary/BoundaryConditions.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"

namespace hyteg {

template<typename ValueType>
class DGFunction : public Function<DGFunction<ValueType> >
{
public:

  DGFunction(const std::string &name, const std::shared_ptr<PrimitiveStorage> &storage, uint_t minLevel, uint_t maxLevel ) :
    DGFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create012BC() )
  {}

  DGFunction(const std::string &name, const std::shared_ptr<PrimitiveStorage> &storage, uint_t minLevel,
             uint_t maxLevel, BoundaryCondition boundaryCondition ) :
    Function<DGFunction<ValueType> >(name, storage, minLevel, maxLevel), boundaryCondition_( boundaryCondition )
  {
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
  inline void enumerate( uint_t level, ValueType offset );

  inline void enumerate( uint_t level );

  inline real_t getMaxValue( const uint_t level, DoFType flag = All );
  inline real_t getMinValue( const uint_t level, DoFType flag = All );
  inline real_t getMaxMagnitude( const uint_t level, DoFType flag = All );

  const PrimitiveDataID<FunctionMemory<ValueType>, Vertex> &getVertexDataID() const { return vertexDataID_; }

  const PrimitiveDataID<FunctionMemory<ValueType>, Edge> &getEdgeDataID() const { return edgeDataID_; }

  const PrimitiveDataID<FunctionMemory<ValueType>, Face> &getFaceDataID() const { return faceDataID_; }

  void projectP1(P1Function< real_t >& src, uint_t level, DoFType flag, UpdateType updateType = Replace);

  inline BoundaryCondition getBoundaryCondition() const { return boundaryCondition_; }

  template< typename SenderType, typename ReceiverType >
  inline void startCommunication( const uint_t & level ) const { communicators_.at( level )->template startCommunication< SenderType, ReceiverType >(); }

  template< typename SenderType, typename ReceiverType >
  inline void endCommunication( const uint_t & level ) const { communicators_.at( level )->template endCommunication< SenderType, ReceiverType >(); }

  template< typename SenderType, typename ReceiverType >
  inline void communicate( const uint_t & level ) const { communicators_.at( level )->template communicate< SenderType, ReceiverType >(); }

  inline void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode & localCommunicationMode )
  {
    for ( auto & communicator : communicators_ )
    {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
    }
  }

private:

  using Function<DGFunction<ValueType> >::communicators_;

  PrimitiveDataID<FunctionMemory<ValueType>, Vertex> vertexDataID_;
  PrimitiveDataID<FunctionMemory<ValueType>, Edge> edgeDataID_;
  PrimitiveDataID<FunctionMemory<ValueType>, Face> faceDataID_;

  BoundaryCondition boundaryCondition_;

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

    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag))
    {
      DGFace::add< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
    }
  }

}

template< typename ValueType >
inline void DGFunction< ValueType >::interpolate(std::function< ValueType( const Point3D& ) >& expr,
                                                 uint_t level, DoFType flag)
{
  std::function< ValueType(const Point3D&,const std::vector<ValueType>&)> exprExtended = [&expr](const hyteg::Point3D& x, const std::vector<ValueType>&) {
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

    const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if (testFlag(vertexBC, flag)) {
      DGVertex::interpolate< ValueType >(level, vertex, vertexDataID_, srcVertexIDs, expr, this->getStorage());
    }
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();

  for (auto &it : this->getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag)) {
      DGEdge::interpolate< ValueType >(level, edge, edgeDataID_, srcEdgeIDs, expr, this->getStorage());
    }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for (auto &it : this->getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag)) {
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

    const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if (testFlag(vertexBC, flag)) {
      DGVertex::assign< ValueType >(level, vertex, scalars, srcVertexIDs, vertexDataID_);
    }
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();

  for (auto &it : this->getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag)) {
      DGEdge::assign< ValueType >(level, edge, scalars, srcEdgeIDs, edgeDataID_);
    }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for (auto &it : this->getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag)) {
      DGFace::assign< ValueType >(level, face, scalars, srcFaceIDs, faceDataID_);
    }
  }

  communicators_[level]->template endCommunication<Edge, Face>();
  this->stopTiming( "Assign" );
}

template< typename ValueType >
void DGFunction< ValueType >::enumerate(uint_t level) {
   enumerate( level, static_cast< ValueType >(0) );
}

template< typename ValueType >
void DGFunction< ValueType >::enumerate(uint_t level, ValueType offset)
{
  this->startTiming( "Enumerate" );

  uint_t counter = hyteg::numberOfLocalDoFs< VertexDoFFunctionTag >( *( this->getStorage() ), level );

  std::vector< uint_t > dofs_per_rank = walberla::mpi::allGather( counter );

  ValueType startOnRank = offset;

  for( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
  {
    startOnRank += dofs_per_rank[i];
  }
  for (auto& it : this->getStorage()->getVertices()) {
    Vertex& vertex = *it.second;
    DGVertex::enumerate(vertex,vertexDataID_,level,startOnRank);
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();
  communicators_[level]->template endCommunication<Vertex, Edge>();

  for (auto& it : this->getStorage()->getEdges()) {
    Edge& edge = *it.second;
    DGEdge::enumerate< ValueType >(level, edge, edgeDataID_, startOnRank);
  }

  communicators_[level]->template startCommunication<Edge, Face>();
  communicators_[level]->template endCommunication<Edge, Face>();

  for (auto& it : this->getStorage()->getFaces()) {
    Face& face = *it.second;
    DGFace::enumerate< ValueType >(level, face, faceDataID_, startOnRank);
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

  src.startCommunication<Edge, Vertex>( level );
  src.startCommunication<Face, Edge>( level );
  src.endCommunication<Edge, Vertex>( level );

  for (auto& it : this->getStorage()->getVertices())
  {
    Vertex& vertex = *it.second;

    const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if (testFlag(vertexBC, flag))
    {
      DGVertex::projectP1< real_t >(level, vertex, this->getStorage(), src.getVertexDataID(), this->getVertexDataID(), updateType);
    }
  }

  startCommunication<Vertex, Edge>( level );

  src.endCommunication<Face, Edge>( level );

  for (auto& it : this->getStorage()->getEdges())
  {
    Edge& edge = *it.second;

    const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

    if (testFlag(edgeBC, flag))
    {
      DGEdge::projectP1< real_t >(level, edge, this->getStorage(), src.getEdgeDataID(), this->getEdgeDataID(), updateType);
    }
  }

  endCommunication<Vertex, Edge>( level );

  startCommunication<Edge, Face>( level );

  for (auto& it : this->getStorage()->getFaces()) {
    Face& face = *it.second;

    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag))
    {
      DGFace::projectP1<real_t>(level, face, this->getStorage(), src.getFaceDataID(), this->getFaceDataID(), updateType);
    }
  }

  endCommunication<Edge, Face>( level );

  this->stopTiming( "projectP1" );
}

template< typename ValueType >
real_t DGFunction< ValueType >::getMaxValue( const uint_t level, DoFType flag ) {

  real_t localMax = -std::numeric_limits< ValueType >::max();

  for( auto& it : this->getStorage()->getFaces() ) {
    Face& face = *it.second;
    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if( testFlag( faceBC, flag ) )
    {
      localMax = std::max( localMax, DGFace::getMaxValue< ValueType >( level, face, faceDataID_ ));
    }
  }

  walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );
  return localMax;
}

template< typename ValueType >
real_t DGFunction< ValueType >::getMinValue( const uint_t level, DoFType flag ) {

  real_t localMin = std::numeric_limits< ValueType >::min();

  for( auto& it : this->getStorage()->getFaces() ) {
    Face& face = *it.second;
    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if( testFlag( faceBC, flag ) )
    {
      localMin = std::min( localMin, DGFace::getMinValue< ValueType >( level, face, faceDataID_ ));
    }
  }

  walberla::mpi::allReduceInplace( localMin, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() );
  return localMin;
}

template< typename ValueType >
real_t DGFunction< ValueType >::getMaxMagnitude( const uint_t level, DoFType flag ) {

  real_t localMax = real_t(0.0);

  for( auto& it : this->getStorage()->getFaces() ) {
    Face& face = *it.second;
    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if( testFlag( faceBC, flag ) )
    {
      localMax = std::max( localMax, DGFace::getMaxMagnitude< ValueType >( level, face, faceDataID_ ));
    }
  }

  walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );
  return localMax;
}

}
