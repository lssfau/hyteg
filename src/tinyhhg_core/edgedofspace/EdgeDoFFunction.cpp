#include "EdgeDoFFunction.hpp"

namespace hhg{

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::interpolate( const ValueType & constant, uint_t level, DoFType flag )
{
   if ( isDummy() ) { return; }
   this->startTiming( "Interpolate" );

   for ( auto & it : this->getStorage()->getEdges() )
   {
      Edge & edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, constant );
      }
   }

   for ( auto & it : this->getStorage()->getFaces() )
   {
      Face & face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroface::interpolate< ValueType >( level, face, faceDataID_, constant );
      }
   }

   for ( auto & it : this->getStorage()->getCells() )
   {
      Cell & cell = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, constant );
      }
   }
   this->stopTiming( "Interpolate" );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::interpolate(const std::function< ValueType( const Point3D& ) >& expr,
                                                      uint_t level, DoFType flag)
{
   if ( isDummy() ) { return; }
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
if ( isDummy() ) { return; }
this->startTiming( "Interpolate" );
// Collect all source IDs in a vector
std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>>   srcEdgeIDs;
std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>>   srcFaceIDs;
std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Cell>>   srcCellIDs;

for (auto& function : srcFunctions)
{
srcEdgeIDs.push_back(function->edgeDataID_);
srcFaceIDs.push_back(function->faceDataID_);
srcCellIDs.push_back(function->cellDataID_);
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

for ( auto & it : this->getStorage()->getCells() )
{
Cell & cell = *it.second;

if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
{
edgedof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
}
}
this->stopTiming( "Interpolate" );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::assign(const std::vector<ValueType> scalars, const std::vector<EdgeDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
if ( isDummy() ) { return; }
this->startTiming( "Assign" );
std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge >>     srcEdgeIDs;
std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face >>     srcFaceIDs;
std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Cell >>     srcCellIDs;

for ( auto& function : functions )
{
srcEdgeIDs.push_back( function->edgeDataID_);
srcFaceIDs.push_back( function->faceDataID_ );
srcCellIDs.push_back( function->cellDataID_ );
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

for ( auto & it : this->getStorage()->getCells() )
{
Cell & cell = *it.second;

if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
{
edgedof::macrocell::assign< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
}
}


this->stopTiming( "Assign" );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::add(const std::vector<ValueType> scalars, const std::vector<EdgeDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
if ( isDummy() ) { return; }
this->startTiming( "Add" );
std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge >>     srcEdgeIDs;
std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face >>     srcFaceIDs;
std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Cell >>     srcCellIDs;

for ( auto& function : functions )
{
srcEdgeIDs.push_back( function->edgeDataID_);
srcFaceIDs.push_back( function->faceDataID_ );
srcCellIDs.push_back( function->cellDataID_ );
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

for ( auto & it : this->getStorage()->getCells() )
{
Cell & cell = *it.second;

if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
{
edgedof::macrocell::add< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
}
}

this->stopTiming( "Add" );
}

template< typename ValueType >
inline real_t EdgeDoFFunction< ValueType >::dotLocal(EdgeDoFFunction< ValueType >& rhs, size_t level, DoFType flag)
{
   if ( isDummy() ) { return real_c(0); }
   this->startTiming( "Dot (local)" );
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

   for ( auto & it : this->getStorage()->getCells() )
   {
      Cell & cell = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProduct += edgedof::macrocell::dot< ValueType >( level, cell, cellDataID_, rhs.cellDataID_ );
      }
   }

   this->stopTiming( "Dot (local)" );

   return scalarProduct;
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::enumerate(uint_t level)
{
   if ( isDummy() ) { return; }
   this->startTiming( "Enumerate" );

   uint_t counter = hhg::numberOfLocalDoFs< EdgeDoFFunctionTag >( *( this->getStorage() ), level );

   std::vector< uint_t > dofs_per_rank = walberla::mpi::allGather( counter );

   ValueType startOnRank = 0;

   for( uint_t i = 0; i < uint_c(walberla::MPIManager::instance()->rank()); ++i )
   {
      startOnRank += dofs_per_rank[i];
   }

   enumerate( level, startOnRank);
   this->stopTiming( "Enumerate" );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::enumerate(uint_t level, ValueType& offset)
{

   for (auto& it : this->getStorage()->getEdges()) {
      Edge& edge = *it.second;
      edgedof::macroedge::enumerate< ValueType >(level, edge, edgeDataID_, offset);
   }

   communicators_[level]->template startCommunication<Edge, Face>();


   for (auto& it : this->getStorage()->getFaces()) {
      Face& face = *it.second;
      edgedof::macroface::enumerate< ValueType >(level, face, faceDataID_, offset);
   }

   communicators_[level]->template endCommunication<Edge, Face>();

   communicators_[level]->template startCommunication<Face, Edge>();
   communicators_[level]->template endCommunication<Face, Edge>();

   communicators_[level]->template startCommunication<Edge, Vertex>();
   communicators_[level]->template endCommunication<Edge, Vertex>();

}


template< typename ValueType >
inline ValueType EdgeDoFFunction< ValueType >::getMaxMagnitude( uint_t level, DoFType flag, bool mpiReduce )
{
   if ( isDummy() ) { return ValueType(0); }
   auto localMax = ValueType(0.0);

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

template class EdgeDoFFunction< float >;
template class EdgeDoFFunction< double >;
template class EdgeDoFFunction< int >;

}
