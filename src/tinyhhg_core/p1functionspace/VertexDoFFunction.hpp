#pragma once

#include <tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp>

#include <tinyhhg_core/p1functionspace/VertexDoFPackInfo.hpp>
#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/communication/BufferedCommunication.hpp"

#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/p1functionspace/P1DataHandling.hpp"

#include "tinyhhg_core/boundary/BoundaryConditions.hpp"

namespace hhg {
namespace vertexdof {

template< typename ValueType >
class VertexDoFFunction : public Function< VertexDoFFunction< ValueType > >
{
public:

  VertexDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage > & storage, uint_t minLevel, uint_t maxLevel ) :
    VertexDoFFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create012BC() )
  {}

  VertexDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage > & storage, uint_t minLevel, uint_t maxLevel, BoundaryCondition boundaryCondition ) :
      Function< VertexDoFFunction< ValueType > >( name, storage, minLevel, maxLevel ), boundaryCondition_( boundaryCondition )
  {
    auto cellVertexDoFFunctionMemoryDataHandling   = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Cell   > >( minLevel, maxLevel, vertexDoFMacroCellFunctionMemorySize );
    auto faceVertexDoFFunctionMemoryDataHandling   = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face   > >( minLevel, maxLevel, vertexDoFMacroFaceFunctionMemorySize );
    auto edgeVertexDoFFunctionMemoryDataHandling   = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Edge   > >( minLevel, maxLevel, vertexDoFMacroEdgeFunctionMemorySize   );
    auto vertexVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Vertex > >( minLevel, maxLevel, vertexDoFMacroVertexFunctionMemorySize );

    storage->addCellData( cellDataID_, cellVertexDoFFunctionMemoryDataHandling, name );
    storage->addFaceData( faceDataID_, faceVertexDoFFunctionMemoryDataHandling, name );
    storage->addEdgeData( edgeDataID_, edgeVertexDoFFunctionMemoryDataHandling, name );
    storage->addVertexData( vertexDataID_, vertexVertexDoFFunctionMemoryDataHandling, name );

    for ( uint_t level = minLevel; level <= maxLevel; ++level )
    {
      communicators_[level]->addPackInfo( std::make_shared< VertexDoFPackInfo< ValueType > >( level, vertexDataID_, edgeDataID_, faceDataID_, cellDataID_, this->getStorage() ) );
    }
  }

  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex> & getVertexDataID() const { return vertexDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >, Edge>   & getEdgeDataID()   const { return edgeDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >, Face>   & getFaceDataID()   const { return faceDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >, Cell>   & getCellDataID()   const { return cellDataID_; }

  inline void assign(const std::vector<ValueType> scalars, const std::vector<VertexDoFFunction< ValueType >*> functions, uint_t level, DoFType flag = All);

  inline void add(const std::vector<ValueType> scalars, const std::vector<VertexDoFFunction< ValueType >*> functions, uint_t level, DoFType flag = All);

  inline real_t dot(VertexDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All);

  inline void prolongate(uint_t sourceLevel, DoFType flag = All);

  inline void prolongateQuadratic(uint_t sourceLevel, DoFType flag = All);

  inline void restrict(uint_t sourceLevel, DoFType flag = All);

  inline void integrateDG(DGFunction< ValueType >& rhs, VertexDoFFunction< ValueType >& rhsP1, uint_t level, DoFType flag);

  /// Interpolates a given expression to a VertexDoFFunction
  inline void interpolate(std::function< ValueType( const Point3D & ) >& expr,
                          uint_t level, DoFType flag = All);

  inline void interpolateExtended(std::function< ValueType( const Point3D&, const std::vector<ValueType>& ) >& expr,
                                  const std::vector<VertexDoFFunction*> srcFunctions,
                                  uint_t level, DoFType flag = All);

  // TODO: write more general version(s)
  inline real_t getMaxValue( uint_t level, DoFType flag = All );
  inline real_t getMinValue( uint_t level, DoFType flag = All );
  inline real_t getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true );

  inline uint_t getNumLocalDoFs ( const uint_t & level ) const;
  inline uint_t getNumGlobalDoFs( const uint_t & level ) const;

  inline BoundaryCondition getBoundaryCondition() const { return boundaryCondition_; }

  template< typename SenderType, typename ReceiverType >
  inline void startCommunication( const uint_t & level ) const { communicators_.at( level )->startCommunication< SenderType, ReceiverType >(); }

  template< typename SenderType, typename ReceiverType >
  inline void endCommunication( const uint_t & level ) const { communicators_.at( level )->endCommunication< SenderType, ReceiverType >(); }

  template< typename SenderType, typename ReceiverType >
  inline void communicate( const uint_t & level ) const { communicators_.at( level )->communicate< SenderType, ReceiverType >(); }

  inline void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode & localCommunicationMode )
  {
    for ( auto & communicator : communicators_ )
    {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
    }
  }

private:

  using Function< VertexDoFFunction< ValueType > >::communicators_;

  inline void enumerate_impl(uint_t level, uint_t& num);

  PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
  PrimitiveDataID< FunctionMemory< ValueType >, Edge >   edgeDataID_;
  PrimitiveDataID< FunctionMemory< ValueType >, Face >   faceDataID_;
  PrimitiveDataID< FunctionMemory< ValueType >, Cell >   cellDataID_;

  BoundaryCondition boundaryCondition_;
};

template< typename ValueType >
inline void VertexDoFFunction< ValueType >::interpolate(std::function< ValueType( const Point3D& ) >& expr,
                                                        uint_t level, DoFType flag)
{
  std::function< ValueType(const Point3D&,const std::vector<ValueType>&)> exprExtended = [&expr](const hhg::Point3D& x, const std::vector<ValueType>&) {
      return expr(x);
  };
  interpolateExtended( exprExtended, {}, level, flag );
}

template< typename ValueType >
inline void VertexDoFFunction< ValueType >::interpolateExtended(std::function< ValueType( const Point3D&, const std::vector<ValueType>& ) >& expr,
                                                                const std::vector<VertexDoFFunction*> srcFunctions,
                                                                uint_t level, DoFType flag)
{
  this->startTiming( "Interpolate" );
  // Collect all source IDs in a vector
  std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
  std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
  std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
  std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

  for ( const auto & function : srcFunctions )
  {
    srcVertexIDs.push_back(function->vertexDataID_);
    srcEdgeIDs.push_back(function->edgeDataID_);
    srcFaceIDs.push_back(function->faceDataID_);
    srcCellIDs.push_back(function->cellDataID_);
  }

  for ( const auto & it : this->getStorage()->getVertices() )
  {
    Vertex & vertex = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macrovertex::interpolate( vertex, vertexDataID_, srcVertexIDs, expr, level );
    }
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();

  for ( const auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
    }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for ( auto& it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
    }
  }

  communicators_[level]->template endCommunication<Edge, Face>();

  for ( const auto & it : this->getStorage()->getCells() )
  {
    Cell & cell = *it.second;

    if ( testFlag(  boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag  ) )
    {
      vertexdof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
    }

  }
  this->stopTiming( "Interpolate" );
}

template< typename ValueType >
inline void VertexDoFFunction< ValueType >::assign(const std::vector<ValueType> scalars, const std::vector<VertexDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
  this->startTiming( "Assign" );
    // Collect all source IDs in a vector
    std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
    std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge > >     srcEdgeIDs;
    std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face > >     srcFaceIDs;
    std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Cell > >     srcCellIDs;

    for ( const auto & function : functions )
    {
      srcVertexIDs.push_back( function->vertexDataID_ );
      srcEdgeIDs.push_back( function->edgeDataID_ );
      srcFaceIDs.push_back( function->faceDataID_ );
      srcCellIDs.push_back( function->cellDataID_ );
    }

    for ( const auto & it : this->getStorage()->getVertices() )
    {
      Vertex & vertex = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
        vertexdof::macrovertex::assign< ValueType >( vertex, scalars, srcVertexIDs, vertexDataID_, level );
      }
    }

    communicators_[level]->template startCommunication< Vertex, Edge >();

    for ( const auto & it : this->getStorage()->getEdges() )
    {
      Edge & edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
        vertexdof::macroedge::assign< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
    }

    communicators_[level]->template endCommunication< Vertex, Edge >();
    communicators_[level]->template startCommunication< Edge, Face >();

    for ( const auto & it : this->getStorage()->getFaces() )
    {
      Face & face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
        vertexdof::macroface::assign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
      }
    }

    communicators_[level]->template endCommunication<Edge, Face>();

    for ( const auto & it : this->getStorage()->getCells() )
    {
      Cell & cell = *it.second;
      if ( testFlag(  boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag  ) )
      {
        vertexdof::macrocell::assign< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
      }
    }
  this->stopTiming( "Assign" );
}

template< typename ValueType >
inline void VertexDoFFunction< ValueType >::add(const std::vector<ValueType> scalars, const std::vector<VertexDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
  this->startTiming( "Add" );
  // Collect all source IDs in a vector
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge > >     srcEdgeIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face > >     srcFaceIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Cell > >     srcCellIDs;

  for ( auto& function : functions )
  {
    srcVertexIDs.push_back( function->vertexDataID_ );
    srcEdgeIDs.push_back( function->edgeDataID_ );
    srcFaceIDs.push_back( function->faceDataID_ );
    srcCellIDs.push_back( function->cellDataID_ );
  }

  for ( const auto & it : this->getStorage()->getVertices() )
  {
    Vertex & vertex = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macrovertex::add( vertex, scalars, srcVertexIDs, vertexDataID_, level );
    }
  }

  communicators_[level]->template startCommunication< Vertex, Edge >();

  for ( const auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macroedge::add< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
    }
  }

  communicators_[level]->template endCommunication< Vertex, Edge >();
  communicators_[level]->template startCommunication< Edge, Face >();

  for ( const auto & it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macroface::add< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
    }
  }

  communicators_[level]->template endCommunication< Edge, Face >();

  for ( const auto & it : this->getStorage()->getCells() )
  {
    Cell & cell = *it.second;
    if ( testFlag(  boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag  ) )
    {
      vertexdof::macrocell::add< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
    }
  }
  this->stopTiming( "Add" );
}

template< typename ValueType >
inline real_t VertexDoFFunction< ValueType >::dot(VertexDoFFunction< ValueType >& rhs, size_t level, DoFType flag)
{
  this->startTiming( "Dot" );
  real_t scalarProduct = 0.0;

  for ( const auto & it : this->getStorage()->getVertices() )
  {
    Vertex& vertex = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
    {
      scalarProduct += vertexdof::macrovertex::dot( vertex, vertexDataID_, rhs.vertexDataID_, level );
    }
  }

  for ( const auto & it : this->getStorage()->getEdges() )
  {
    Edge& edge = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      scalarProduct += vertexdof::macroedge::dot< ValueType >( level, edge, edgeDataID_, rhs.edgeDataID_ );
    }
  }

  for ( const auto & it : this->getStorage()->getFaces() )
  {
    Face& face = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      scalarProduct += vertexdof::macroface::dot< ValueType >( level, face, faceDataID_, rhs.faceDataID_ );
    }
  }

  for ( const auto & it : this->getStorage()->getCells() )
  {
    Cell& cell = *it.second;
    if ( testFlag(  boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag  ) )
    {
      scalarProduct += vertexdof::macrocell::dot< ValueType >( level, cell, cellDataID_, rhs.cellDataID_ );
    }
  }

  walberla::mpi::allReduceInplace( scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );

  this->stopTiming( "Dot" );
  return scalarProduct;
}

template< typename ValueType >
inline void VertexDoFFunction< ValueType >::prolongate(size_t sourceLevel, DoFType flag)
{
  this->startTiming( "Prolongate" );
  const size_t destinationLevel = sourceLevel + 1;

  for (auto& it : this->getStorage()->getVertices()) {
      Vertex& vertex = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
        vertexdof::macrovertex::prolongate(vertex, vertexDataID_, sourceLevel);
      }
  }

  communicators_[destinationLevel]->template startCommunication<Vertex, Edge>();

  for (auto& it : this->getStorage()->getEdges()) {
      Edge& edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
        vertexdof::macroedge::prolongate< ValueType >(sourceLevel, edge, edgeDataID_);
      }
  }

  communicators_[destinationLevel]->template endCommunication<Vertex, Edge>();
  communicators_[destinationLevel]->template startCommunication<Edge, Face>();

  for (auto& it : this->getStorage()->getFaces()) {
      Face& face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
        vertexdof::macroface::prolongate< ValueType >(sourceLevel, face, faceDataID_);
      }
  }

  communicators_[destinationLevel]->template endCommunication<Edge, Face>();
  this->stopTiming( "Prolongate" );
}

template< typename ValueType >
inline void VertexDoFFunction< ValueType >::prolongateQuadratic(size_t sourceLevel, DoFType flag)
{
  this->startTiming( "Prolongate Quadratic" );
  const size_t destinationLevel = sourceLevel + 1;

  for (auto& it : this->getStorage()->getVertices()) {
    Vertex& vertex = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macrovertex::prolongateQuadratic(vertex, vertexDataID_, sourceLevel);
    }
  }

  communicators_[destinationLevel]->template startCommunication<Vertex, Edge>();

  for (auto& it : this->getStorage()->getEdges()) {
    Edge& edge = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macroedge::prolongateQuadratic< ValueType >(sourceLevel, edge, edgeDataID_);
    }
  }

  communicators_[destinationLevel]->template endCommunication<Vertex, Edge>();
  communicators_[destinationLevel]->template startCommunication<Edge, Face>();

  for (auto& it : this->getStorage()->getFaces()) {
    Face& face = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macroface::prolongateQuadratic< ValueType >(sourceLevel, face, faceDataID_);
    }
  }

  communicators_[destinationLevel]->template endCommunication<Edge, Face>();
  this->stopTiming( "Prolongate Quadratic" );
}

template< typename ValueType >
inline void VertexDoFFunction< ValueType >::restrict(size_t sourceLevel, DoFType flag)
{
  this->startTiming( "Restrict" );
  const size_t destinationLevel = sourceLevel - 1;

  // start pulling vertex halos
  communicators_[sourceLevel]->template startCommunication<Edge, Vertex>();

  // start pulling edge halos
  communicators_[sourceLevel]->template startCommunication<Face, Edge>();

  // end pulling vertex halos
  communicators_[sourceLevel]->template endCommunication<Edge, Vertex>();

  for (auto& it : this->getStorage()->getVertices()) {
      Vertex& vertex = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
        vertexdof::macrovertex::restrict(vertex, vertexDataID_, sourceLevel);
      }
  }

  communicators_[destinationLevel]->template startCommunication<Vertex, Edge>();

  // end pulling edge halos
  communicators_[sourceLevel]->template endCommunication<Face, Edge>();

  for (auto& it : this->getStorage()->getEdges()) {
      Edge& edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
        vertexdof::macroedge::restrict< ValueType >(sourceLevel, edge, edgeDataID_);
      }
  }

  communicators_[destinationLevel]->template endCommunication<Vertex, Edge>();

  communicators_[destinationLevel]->template startCommunication<Edge, Face>();

  for (auto& it : this->getStorage()->getFaces()) {
      Face& face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
        vertexdof::macroface::restrict< ValueType >(sourceLevel, face, faceDataID_);
      }
  }

  communicators_[destinationLevel]->template endCommunication<Edge, Face>();
  this->stopTiming( "Restrict" );
}

template< typename ValueType >
inline void VertexDoFFunction< ValueType >::enumerate_impl(uint_t level, uint_t& num)
{
  this->startTiming( "Enumerate" );
  for (auto& it : this->getStorage()->getVertices()) {
    Vertex& vertex = *it.second;
    vertexdof::macrovertex::enumerate(level, vertex, vertexDataID_, num);
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();
  communicators_[level]->template endCommunication<Vertex, Edge>();

  for (auto& it : this->getStorage()->getEdges()) {
    Edge& edge = *it.second;
    vertexdof::macroedge::enumerate< ValueType >(level, edge, edgeDataID_, num);
  }

  communicators_[level]->template startCommunication<Edge, Face>();
  communicators_[level]->template endCommunication<Edge, Face>();

  for (auto& it : this->getStorage()->getFaces()) {
    Face& face = *it.second;
    vertexdof::macroface::enumerate< ValueType >(level, face, faceDataID_, num);
  }

  communicators_[level]->template startCommunication<Face, Edge>();
  communicators_[level]->template endCommunication<Face, Edge>();

  communicators_[level]->template startCommunication<Edge, Vertex>();
  communicators_[level]->template endCommunication<Edge, Vertex>();
  this->stopTiming( "Enumerate" );
}

template< typename ValueType >
inline void VertexDoFFunction< ValueType >::integrateDG(DGFunction< ValueType >& rhs, VertexDoFFunction< ValueType >& rhsP1, uint_t level, DoFType flag)
{
  this->startTiming( "integrateDG" );

  rhsP1.startCommunication<Edge, Vertex>( level );
  rhsP1.startCommunication<Face, Edge>( level );

  rhs.startCommunication<Face, Edge>( level );
  rhs.endCommunication<Face, Edge>( level );

  rhs.startCommunication<Edge, Vertex>( level );
  rhs.endCommunication<Edge, Vertex>( level );

  rhsP1.endCommunication<Edge, Vertex>( level );

  for (auto& it : this->getStorage()->getVertices()) {
    Vertex& vertex = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macrovertex::integrateDG< ValueType >(vertex, this->getStorage(), rhs.getVertexDataID(), rhsP1.getVertexDataID(), vertexDataID_, level);
    }
  }

  communicators_[level]->template startCommunication<Vertex, Edge>();
  rhsP1.endCommunication<Face, Edge>( level );

  for (auto& it : this->getStorage()->getEdges()) {
    Edge& edge = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macroedge::integrateDG< ValueType >(level, edge, this->getStorage(), rhs.getEdgeDataID(), rhsP1.getEdgeDataID(), edgeDataID_);
    }
  }

  communicators_[level]->template endCommunication<Vertex, Edge>();
  communicators_[level]->template startCommunication<Edge, Face>();

  for (auto& it : this->getStorage()->getFaces()) {
    Face& face = *it.second;

    if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      vertexdof::macroface::integrateDG< ValueType >(level, face, rhs.getFaceDataID(), rhsP1.getFaceDataID(), faceDataID_);
    }
  }

  communicators_[level]->template endCommunication<Edge, Face>();

  this->stopTiming( "integrateDG" );
}

inline void projectMean(VertexDoFFunction<real_t>& pressure, VertexDoFFunction<real_t>& tmp, uint_t level) {

  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) {
    return 1.0;
  };

  tmp.interpolate(ones, level);

  real_t numGlobalVertices = tmp.dot(tmp, level, hhg::All);
  real_t mean = pressure.dot(tmp, level, hhg::All);

  pressure.assign({1.0, -mean/numGlobalVertices}, {&pressure, &tmp}, level, hhg::All);
}


template< typename ValueType >
inline real_t VertexDoFFunction< ValueType >::getMaxValue( uint_t level, DoFType flag )
{
  real_t localMax = - std::numeric_limits<real_t>::max();

  for( auto& it : this->getStorage()->getFaces() ) {
    Face& face = *it.second;
    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if ( testFlag( faceBC, flag ) ) {
      localMax = std::max(localMax, vertexdof::macroface::getMaxValue< ValueType >( level, face, faceDataID_ ));
    }
  }

  for( auto& it : this->getStorage()->getEdges() ) {
    Edge& edge = *it.second;
    const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if ( testFlag( edgeBC, flag ) ) {
      localMax = std::max(localMax, vertexdof::macroedge::getMaxValue< ValueType >( level, edge, edgeDataID_ ));
    }
  }

  for( auto& it : this->getStorage()->getVertices() ) {
    Vertex& vertex = *it.second;
    const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if ( testFlag( vertexBC, flag ) ) {
      localMax = std::max(localMax, vertexdof::macrovertex::getMaxValue< ValueType >( level, vertex, vertexDataID_ ));
    }
  }

  real_t globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );

  return globalMax;
}


template< typename ValueType >
inline real_t VertexDoFFunction< ValueType >::getMaxMagnitude( uint_t level, DoFType flag, bool mpiReduce )
{
  real_t localMax = real_t(0.0);

  for( auto& it : this->getStorage()->getFaces() ) {
    Face& face = *it.second;
    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if ( testFlag( faceBC, flag ) ) {
      localMax = std::max(localMax, vertexdof::macroface::getMaxMagnitude< ValueType >( level, face, faceDataID_ ));
    }
  }

  for( auto& it : this->getStorage()->getEdges() ) {
    Edge& edge = *it.second;
    const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if ( testFlag( edgeBC, flag ) ) {
      localMax = std::max(localMax, vertexdof::macroedge::getMaxMagnitude< ValueType >( level, edge, edgeDataID_ ));
    }
  }

  for( auto& it : this->getStorage()->getVertices() ) {
    Vertex& vertex = *it.second;
    const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if ( testFlag( vertexBC, flag ) ) {
      localMax = std::max(localMax, vertexdof::macrovertex::getMaxMagnitude< ValueType >( level, vertex, vertexDataID_ ));
    }
  }

  real_t globalMax = localMax;
  if( mpiReduce ) {
    globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
  }

  return globalMax;
}


template< typename ValueType >
inline real_t VertexDoFFunction< ValueType >::getMinValue( uint_t level, DoFType flag )
{
  real_t localMin = std::numeric_limits<real_t>::max();

  for( auto& it : this->getStorage()->getFaces() ) {
    Face& face = *it.second;
    const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if ( testFlag( faceBC, flag ) ) {
      localMin = std::min(localMin, vertexdof::macroface::getMinValue< ValueType >( level, face, faceDataID_ ));
    }
  }

  for( auto& it : this->getStorage()->getEdges() ) {
    Edge& edge = *it.second;
    const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if ( testFlag( edgeBC, flag ) ) {
      localMin = std::min(localMin, vertexdof::macroedge::getMinValue< ValueType >( level, edge, edgeDataID_ ));
    }
  }

  for( auto& it : this->getStorage()->getVertices() ) {
    Vertex& vertex = *it.second;
    const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if ( testFlag( vertexBC, flag ) ) {
      localMin = std::min(localMin, vertexdof::macrovertex::getMinValue< ValueType >( level, vertex, vertexDataID_ ));
    }
  }

  real_t globalMin = - walberla::mpi::allReduce( -localMin, walberla::mpi::MAX );

  return globalMin;
}


template< typename ValueType >
inline uint_t VertexDoFFunction< ValueType >::getNumLocalDoFs( const uint_t & level ) const
{
  const uint_t numDoFsOnSingleVertex = 1;
  const uint_t numDoFsOnSingleEdge   = levelinfo::num_microvertices_per_edge( level ) - 2;
  const uint_t numDoFsOnSingleFace   = levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microvertices_per_edge( level ) - 3 );
  const uint_t numDoFsOnSingleCell   = levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microvertices_per_edge( level ) - 4 );
  return   numDoFsOnSingleVertex * this->getStorage()->getNumberOfLocalVertices()
         + numDoFsOnSingleEdge   * this->getStorage()->getNumberOfLocalEdges()
         + numDoFsOnSingleFace   * this->getStorage()->getNumberOfLocalFaces()
         + numDoFsOnSingleCell   * this->getStorage()->getNumberOfLocalCells();
}

template< typename ValueType >
inline uint_t VertexDoFFunction< ValueType >::getNumGlobalDoFs( const uint_t & level ) const
{
  return walberla::mpi::allReduce( getNumLocalDoFs( level ), walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
}

} // namespace vertexdof
} // namespace hhg
