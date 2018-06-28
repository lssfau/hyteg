#pragma once

#include "core/mpi/all.h"

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/communication/BufferedCommunication.hpp"

#include "EdgeDoFMacroFace.hpp"
#include "EdgeDoFMacroEdge.hpp"
#include "EdgeDoFPackInfo.hpp"

#include "tinyhhg_core/boundary/BoundaryConditions.hpp"

namespace hhg {

namespace edgedof {
///@name Size Functions
///@{

inline uint_t edgeDoFMacroVertexFunctionMemorySize( const uint_t &level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  return 2 * primitive.getNumNeighborEdges();
}

inline uint_t edgeDoFMacroEdgeFunctionMemorySize( const uint_t &level, const Primitive & primitive )
{
  return levelinfo::num_microedges_per_edge( level ) + primitive.getNumNeighborFaces() * ( 3 * ( levelinfo::num_microedges_per_edge( level ) ) - 1 );
}

inline uint_t edgeDoFMacroFaceFunctionMemorySize( const uint_t &level, const Primitive & primitive )
{
  WALBERLA_UNUSED( primitive );
  return 3 * ( ( ( levelinfo::num_microedges_per_edge( level ) + 1 ) * levelinfo::num_microedges_per_edge( level ) ) / 2 );
}

///@}

}// namespace edgedof

template< typename ValueType >
class EdgeDoFFunction : public Function< EdgeDoFFunction< ValueType > >
{
public:

  EdgeDoFFunction( const std::string & name, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & minLevel, const uint_t & maxLevel ) :
    EdgeDoFFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create012BC() )
  {}

  EdgeDoFFunction( const std::string & name, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & minLevel, const uint_t & maxLevel, const BoundaryCondition & boundaryCondition ) :
      Function< EdgeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel ), boundaryCondition_( boundaryCondition )
  {
    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Vertex >> vertexDataHandling =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Vertex >>(minLevel, maxLevel, edgedof::edgeDoFMacroVertexFunctionMemorySize);

    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Edge >> edgeDataHandling   =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Edge   >>(minLevel, maxLevel, edgedof::edgeDoFMacroEdgeFunctionMemorySize);

    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Face >> faceDataHandling   =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Face   >>(minLevel, maxLevel, edgedof::edgeDoFMacroFaceFunctionMemorySize);


    storage->addVertexData( vertexDataID_, vertexDataHandling, name );
    storage->addEdgeData(   edgeDataID_,   edgeDataHandling,   name );
    storage->addFaceData(   faceDataID_,   faceDataHandling,   name );

    for (uint_t level = minLevel; level <= maxLevel; ++level) {
      //communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
      communicators_[level]->addPackInfo(
        std::make_shared<EdgeDoFPackInfo<ValueType> >(level, vertexDataID_, edgeDataID_, faceDataID_, this->getStorage()));
    }
  }

  inline void
  assign( const std::vector< ValueType > scalars, const std::vector< EdgeDoFFunction< ValueType >* > functions,
          uint_t level, DoFType flag = All );

  inline void
  add( const std::vector< ValueType > scalars, const std::vector< EdgeDoFFunction< ValueType >* > functions,
       uint_t level, DoFType flag = All );

  /// Interpolates a given expression to a EdgeDoFFunction

  inline void
  interpolate( const std::function< ValueType( const Point3D & ) >& expr,
                          uint_t level, DoFType flag = All);

  inline void
  interpolateExtended( const std::function<ValueType(const Point3D &, const std::vector<ValueType>&)> &expr,
                       const std::vector<EdgeDoFFunction<ValueType>*> srcFunctions,
                      uint_t level,
                      DoFType flag = All);

  inline real_t
  dot( EdgeDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All );

  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex>   & getVertexDataID() const { return vertexDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Edge>   & getEdgeDataID()   const { return edgeDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Face>   & getFaceDataID()   const { return faceDataID_; }


  inline real_t getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true );

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

    using Function< EdgeDoFFunction< ValueType > >::communicators_;

    inline void
    enumerate_impl( uint_t level, uint_t& num );

    PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
    PrimitiveDataID< FunctionMemory< ValueType >, Edge > edgeDataID_;
    PrimitiveDataID< FunctionMemory< ValueType >, Face > faceDataID_;

    BoundaryCondition boundaryCondition_;
};

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::interpolate(const std::function< ValueType( const Point3D& ) >& expr,
                                                      uint_t level, DoFType flag)
{
  std::function< ValueType(const Point3D&,const std::vector<ValueType>&)> exprExtended = [&expr](const hhg::Point3D& x, const std::vector<ValueType>&) {
      return expr(x);
  };
  interpolateExtended( exprExtended, {}, level, flag );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::interpolateExtended(const std::function<ValueType(const Point3D &, const std::vector<ValueType>&)> &expr,
                                                              const std::vector<EdgeDoFFunction<ValueType>*> srcFunctions,
                                                              uint_t level,
                                                              DoFType flag)
{
  this->startTiming( "Interpolate" );
  // Collect all source IDs in a vector
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>>   srcEdgeIDs;
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>>   srcFaceIDs;

  for (auto& function : srcFunctions)
  {
    srcEdgeIDs.push_back(function->edgeDataID_);
    srcFaceIDs.push_back(function->faceDataID_);
  }

  for ( auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      edgedof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
    }
  }

  communicators_[ level ]->template startCommunication< Edge, Face >();

  for ( auto & it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      edgedof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
    }
  }

  communicators_[ level ]->template endCommunication< Edge, Face >();
  this->stopTiming( "Interpolate" );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::assign(const std::vector<ValueType> scalars, const std::vector<EdgeDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
  this->startTiming( "Assign" );
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge >>     srcEdgeIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face >>     srcFaceIDs;

  for ( auto& function : functions )
  {
      srcEdgeIDs.push_back( function->edgeDataID_);
      srcFaceIDs.push_back( function->faceDataID_ );
  }

  for ( auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      edgedof::macroedge::assign< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
    }
  }

  communicators_[ level ]->template startCommunication< Edge, Face >();

  for ( auto & it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      edgedof::macroface::assign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
    }
  }

  communicators_[ level ]->template endCommunication< Edge, Face >();
  this->stopTiming( "Assign" );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::add(const std::vector<ValueType> scalars, const std::vector<EdgeDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
  this->startTiming( "Add" );
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge >>     srcEdgeIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face >>     srcFaceIDs;

  for ( auto& function : functions )
  {
      srcEdgeIDs.push_back( function->edgeDataID_);
      srcFaceIDs.push_back( function->faceDataID_ );
  }

  for ( auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      edgedof::macroedge::add< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
    }
  }

  communicators_[ level ]->template startCommunication< Edge, Face >();

  for ( auto & it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      edgedof::macroface::add< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
    }
  }

  communicators_[ level ]->template endCommunication< Edge, Face >();
  this->stopTiming( "Add" );
}

template< typename ValueType >
inline real_t EdgeDoFFunction< ValueType >::dot(EdgeDoFFunction< ValueType >& rhs, size_t level, DoFType flag)
{
  this->startTiming( "Dot" );
  real_t scalarProduct =  0.0 ;

  for ( auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      scalarProduct += edgedof::macroedge::dot< ValueType >( level, edge, edgeDataID_, rhs.edgeDataID_ );
    }
  }

  for ( auto & it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      scalarProduct += edgedof::macroface::dot< ValueType >( level, face, faceDataID_, rhs.faceDataID_ );
    }
  }

  walberla::mpi::allReduceInplace( scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
  this->stopTiming( "Dot" );
  return scalarProduct;
}


template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::enumerate_impl(uint_t level, uint_t& num)
{
  this->startTiming( "Enumerate" );
  for (auto& it : this->getStorage()->getEdges()) {
    Edge& edge = *it.second;
    edgedof::macroedge::enumerate< ValueType >(level, edge, edgeDataID_, num);
  }

  communicators_[level]->template startCommunication<Edge, Face>();


  for (auto& it : this->getStorage()->getFaces()) {
    Face& face = *it.second;
    edgedof::macroface::enumerate< ValueType >(level, face, faceDataID_, num);
  }

  communicators_[level]->template endCommunication<Edge, Face>();

  communicators_[level]->template startCommunication<Face, Edge>();
  communicators_[level]->template endCommunication<Face, Edge>();

  communicators_[level]->template startCommunication<Edge, Vertex>();
  communicators_[level]->template endCommunication<Edge, Vertex>();
  this->stopTiming( "Enumerate" );
}


template< typename ValueType >
inline real_t EdgeDoFFunction< ValueType >::getMaxMagnitude( uint_t level, DoFType flag, bool mpiReduce )
{
  real_t localMax = real_t(0.0);

  for( auto& it : this->getStorage()->getEdges() )
  {
    Edge &edge = *it.second;
    const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if ( testFlag( edgeBC, flag ) )
    {
      localMax = std::max( localMax, edgedof::macroedge::getMaxMagnitude< ValueType >( level, edge, edgeDataID_ ));
    }
  }

  for( auto& it : this->getStorage()->getFaces() )
  {
    Face &face = *it.second;
    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if ( testFlag( faceBC, flag ) )
    {
      localMax = std::max( localMax, edgedof::macroface::getMaxMagnitude< ValueType >( level, face, faceDataID_ ));
    }
  }

  if( mpiReduce )
  {
    walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );
  }

  return localMax;
}


}// namespace hhg
