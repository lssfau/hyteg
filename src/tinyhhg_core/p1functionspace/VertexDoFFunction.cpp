#include "VertexDoFFunction.hpp"

#include <utility>

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/boundary/BoundaryConditions.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFunction.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFAdditivePackInfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFPackInfo.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"

namespace hhg {
namespace vertexdof {

template < typename ValueType >
VertexDoFFunction< ValueType >::VertexDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage )
: Function< VertexDoFFunction< ValueType > >( name, storage )
, vertexDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Vertex >, Vertex >() )
, edgeDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Edge >, Edge >() )
, faceDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Face >, Face >() )
, cellDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Cell >, Cell >() )
{}

template < typename ValueType >
VertexDoFFunction< ValueType >::VertexDoFFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel )
: VertexDoFFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create012BC() )
{}

template < typename ValueType >
VertexDoFFunction< ValueType >::VertexDoFFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel,
                                                   BoundaryCondition                          boundaryCondition )
: Function< VertexDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
, boundaryCondition_( std::move( boundaryCondition ) )
, boundaryTypeToSkipDuringAdditiveCommunication_( DoFType::DirichletBoundary )
{
   auto cellVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Cell > >(
       minLevel, maxLevel, vertexDoFMacroCellFunctionMemorySize );
   auto faceVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face > >(
       minLevel, maxLevel, vertexDoFMacroFaceFunctionMemorySize );
   auto edgeVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Edge > >(
       minLevel, maxLevel, vertexDoFMacroEdgeFunctionMemorySize );
   auto vertexVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Vertex > >(
       minLevel, maxLevel, vertexDoFMacroVertexFunctionMemorySize );

   storage->addCellData( cellDataID_, cellVertexDoFFunctionMemoryDataHandling, name );
   storage->addFaceData( faceDataID_, faceVertexDoFFunctionMemoryDataHandling, name );
   storage->addEdgeData( edgeDataID_, edgeVertexDoFFunctionMemoryDataHandling, name );
   storage->addVertexData( vertexDataID_, vertexVertexDoFFunctionMemoryDataHandling, name );

   for( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      communicators_[level]->addPackInfo( std::make_shared< VertexDoFPackInfo< ValueType > >(
          level, vertexDataID_, edgeDataID_, faceDataID_, cellDataID_, this->getStorage() ) );
      additiveCommunicators_[level]->addPackInfo(
          std::make_shared< VertexDoFAdditivePackInfo< ValueType > >( level,
                                                                      vertexDataID_,
                                                                      edgeDataID_,
                                                                      faceDataID_,
                                                                      cellDataID_,
                                                                      this->getStorage(),
                                                                      boundaryCondition_,
                                                                      boundaryTypeToSkipDuringAdditiveCommunication_ ) );
   }
}

template < typename ValueType >
BoundaryCondition VertexDoFFunction< ValueType >::getBoundaryCondition() const
{
   return boundaryCondition_;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const ValueType& constant, uint_t level, DoFType flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );

   interpolateByPrimitiveType< Vertex >( constant, level, flag );
   interpolateByPrimitiveType< Edge >( constant, level, flag );
   interpolateByPrimitiveType< Face >( constant, level, flag );
   interpolateByPrimitiveType< Cell >( constant, level, flag );

   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                  uint_t                                              level,
                                                  DoFType                                             flag ) const
{
   if( isDummy() )
   {
      return;
   }
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hhg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolateExtended( exprExtended, {}, level, flag );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                  uint_t                                              level,
                                                  BoundaryUID                                         boundaryUID ) const
{
   if( isDummy() )
   {
      return;
   }
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hhg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolateExtended( exprExtended, {}, level, boundaryUID );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolateExtended(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< VertexDoFFunction* >                                              srcFunctions,
    uint_t                                                                               level,
    DoFType                                                                              flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for( const auto& function : srcFunctions )
   {
      srcVertexIDs.push_back( function->vertexDataID_ );
      srcEdgeIDs.push_back( function->edgeDataID_ );
      srcFaceIDs.push_back( function->faceDataID_ );
      srcCellIDs.push_back( function->cellDataID_ );
   }

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::interpolate( vertex, vertexDataID_, srcVertexIDs, expr, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolateExtended(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< VertexDoFFunction* >                                              srcFunctions,
    uint_t                                                                               level,
    BoundaryUID                                                                          boundaryUID ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for( const auto& function : srcFunctions )
   {
      srcVertexIDs.push_back( function->vertexDataID_ );
      srcEdgeIDs.push_back( function->edgeDataID_ );
      srcFaceIDs.push_back( function->faceDataID_ );
      srcCellIDs.push_back( function->cellDataID_ );
   }

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( boundaryCondition_.getBoundaryUIDFromMeshFlag( vertex.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macrovertex::interpolate( vertex, vertexDataID_, srcVertexIDs, expr, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( boundaryCondition_.getBoundaryUIDFromMeshFlag( edge.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( boundaryCondition_.getBoundaryUIDFromMeshFlag( face.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;

      if( boundaryCondition_.getBoundaryUIDFromMeshFlag( cell.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::assign(
    const std::vector< ValueType >&                                                      scalars,
    const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
    size_t                                                                               level,
    DoFType                                                                              flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Assign" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for( const VertexDoFFunction< ValueType >& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }
   this->startTiming( "Vertex" );
   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::assign< ValueType >( vertex, scalars, srcVertexIDs, vertexDataID_, level );
      }
   }
   this->stopTiming( "Vertex" );
   this->startTiming( "Edge" );
   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::assign< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }
   this->stopTiming( "Edge" );
   this->startTiming( "Face" );
   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::assign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
      }
   }
   this->stopTiming( "Face" );
   this->startTiming( "Cell" );
   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::assign< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
      }
   }
   this->stopTiming( "Cell" );
   this->stopTiming( "Assign" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::add( const ValueType& scalar, const uint_t& level, DoFType flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Add" );

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::add< ValueType >( vertex, scalar, vertexDataID_, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::add< ValueType >( level, edge, scalar, edgeDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::add< ValueType >( level, face, scalar, faceDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::add< ValueType >( level, cell, scalar, cellDataID_ );
      }
   }

   this->stopTiming( "Add" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::add(
    const std::vector< ValueType >&                                                      scalars,
    const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
    size_t                                                                               level,
    DoFType                                                                              flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Add" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for( const VertexDoFFunction< ValueType >& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::add( vertex, scalars, srcVertexIDs, vertexDataID_, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::add< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::add< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::add< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
      }
   }
   this->stopTiming( "Add" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::multElementwise( const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
                                                      const uint_t                                               level,
                                                      const DoFType                                              flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Multiply elementwise" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for( const VertexDoFFunction& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::multElementwise< ValueType >( vertex, srcVertexIDs, vertexDataID_, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::multElementwise< ValueType >( level, edge, srcEdgeIDs, edgeDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::multElementwise< ValueType >( level, face, srcFaceIDs, faceDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::multElementwise< ValueType >( level, cell, srcCellIDs, cellDataID_ );
      }
   }
   this->stopTiming( "Multiply elementwise" );
}

template < typename ValueType >
real_t VertexDoFFunction< ValueType >::dotGlobal(const VertexDoFFunction< ValueType >& rhs, size_t level, DoFType flag ) const
{
   real_t scalarProduct = dotLocal( rhs, level, flag );
   this->startTiming( "Dot (reduce)" );
   walberla::mpi::allReduceInplace( scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   this->stopTiming( "Dot (reduce)" );
   return scalarProduct;
}

template < typename ValueType >
real_t VertexDoFFunction< ValueType >::dotLocal(const VertexDoFFunction< ValueType >& rhs, size_t level, DoFType flag ) const
{
   if( isDummy() )
   {
      return real_c( 0 );
   }
   this->startTiming( "Dot (local)" );
   real_t scalarProduct = 0.0;

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProduct += vertexdof::macrovertex::dot( vertex, vertexDataID_, rhs.vertexDataID_, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProduct += vertexdof::macroedge::dot< ValueType >( level, edge, edgeDataID_, rhs.edgeDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProduct += vertexdof::macroface::dot< ValueType >( level, face, faceDataID_, rhs.faceDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProduct += vertexdof::macrocell::dot< ValueType >( level, cell, cellDataID_, rhs.cellDataID_ );
      }
   }
   this->stopTiming( "Dot (local)" );
   return scalarProduct;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::enumerate( uint_t level ) const
{
   if( isDummy() )
   {
      return;
   }

   this->startTiming( "Enumerate" );

   uint_t counter = hhg::numberOfLocalDoFs< VertexDoFFunctionTag >( *( this->getStorage() ), level );

   std::vector< uint_t > dofs_per_rank = walberla::mpi::allGather( counter );

   ValueType startOnRank = 0;

   for( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
   {
      startOnRank += static_cast< ValueType >( dofs_per_rank[i] );
   }
   enumerate( level, startOnRank );
   this->stopTiming( "Enumerate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::enumerate( uint_t level, ValueType& offset ) const
{
   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;
      vertexdof::macrovertex::enumerate( level, vertex, vertexDataID_, offset );
   }

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;
      vertexdof::macroedge::enumerate< ValueType >( level, edge, edgeDataID_, offset );
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;
      vertexdof::macroface::enumerate< ValueType >( level, face, faceDataID_, offset );
   }

   for( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      vertexdof::macrocell::enumerate< ValueType >( level, cell, cellDataID_, offset );
   }

   /// in contrast to other methods in the function class enumerate needs to communicate due to its usage in the PETSc solvers
   communication::syncFunctionBetweenPrimitives( *this, level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::integrateDG( DGFunction< ValueType >&        rhs,
                                                  VertexDoFFunction< ValueType >& rhsP1,
                                                  uint_t                          level,
                                                  DoFType                         flag )
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "integrateDG" );

   rhsP1.startCommunication< Edge, Vertex >( level );
   rhsP1.startCommunication< Face, Edge >( level );

   rhs.template startCommunication< Face, Edge >( level );
   rhs.template endCommunication< Face, Edge >( level );

   rhs.template startCommunication< Edge, Vertex >( level );
   rhs.template endCommunication< Edge, Vertex >( level );

   rhsP1.endCommunication< Edge, Vertex >( level );

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::integrateDG< ValueType >(
             vertex, this->getStorage(), rhs.getVertexDataID(), rhsP1.getVertexDataID(), vertexDataID_, level );
      }
   }

   communicators_[level]->template startCommunication< Vertex, Edge >();
   rhsP1.endCommunication< Face, Edge >( level );

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::integrateDG< ValueType >(
             level, edge, this->getStorage(), rhs.getEdgeDataID(), rhsP1.getEdgeDataID(), edgeDataID_ );
      }
   }

   communicators_[level]->template endCommunication< Vertex, Edge >();
   communicators_[level]->template startCommunication< Edge, Face >();

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::integrateDG< ValueType >( level, face, rhs.getFaceDataID(), rhsP1.getFaceDataID(), faceDataID_ );
      }
   }

   communicators_[level]->template endCommunication< Edge, Face >();

   this->stopTiming( "integrateDG" );
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::getMaxValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMax = -std::numeric_limits< ValueType >::max();

   for( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMax = std::max( localMax, vertexdof::macrocell::getMaxValue< ValueType >( level, cell, cellDataID_ ) );
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroface::getMaxValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroedge::getMaxValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macrovertex::getMaxValue< ValueType >( level, vertex, vertexDataID_ ) );
      }
   }

   ValueType globalMax = localMax;
   if( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::getMaxMagnitude( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMax = ValueType( 0.0 );

   for( auto& it : this->getStorage()->getCells() )
   {
     Cell& cell = *it.second;
     localMax = std::max( localMax, vertexdof::macrocell::getMaxMagnitude< ValueType >( level, cell, cellDataID_ ) );
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroface::getMaxMagnitude< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroedge::getMaxMagnitude< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macrovertex::getMaxMagnitude< ValueType >( level, vertex, vertexDataID_ ) );
      }
   }

   ValueType globalMax = localMax;
   if( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::getMinValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMin = std::numeric_limits< ValueType >::max();

   for( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMin = std::min( localMin, vertexdof::macrocell::getMinValue< ValueType >( level, cell, cellDataID_ ) );
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macroface::getMinValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macroedge::getMinValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macrovertex::getMinValue< ValueType >( level, vertex, vertexDataID_ ) );
      }
   }

   ValueType globalMin = localMin;
   if( mpiReduce )
   {
      globalMin = -walberla::mpi::allReduce( -localMin, walberla::mpi::MAX );
   }

   return globalMin;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::setLocalCommunicationMode(
    const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode )
{
   if( isDummy() )
   {
      return;
   }
   for( auto& communicator : communicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
   for( auto& communicator : additiveCommunicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
}

template < typename ValueType >
template < typename PrimitiveType >
void VertexDoFFunction< ValueType >::interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );

   if( std::is_same< PrimitiveType, Vertex >::value )
   {
      for( const auto& it : this->getStorage()->getVertices() )
      {
         Vertex& vertex = *it.second;

         if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macrovertex::interpolate( level, vertex, vertexDataID_, constant );
         }
      }
   } else if( std::is_same< PrimitiveType, Edge >::value )
   {
      for( const auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroedge::interpolate( level, edge, edgeDataID_, constant );
         }
      }
   } else if( std::is_same< PrimitiveType, Face >::value )
   {
      for( const auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroface::interpolate( level, face, faceDataID_, constant );
         }
      }
   } else if( std::is_same< PrimitiveType, Cell >::value )
   {
      for( const auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macrocell::interpolate( level, cell, cellDataID_, constant );
         }
      }
   }

   this->stopTiming( "Interpolate" );
}

template class VertexDoFFunction< float >;
template class VertexDoFFunction< double >;
template class VertexDoFFunction< int >;

} // namespace vertexdof
} // namespace hhg
