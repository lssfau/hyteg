
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

namespace hhg {

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;

static void testVertexDoFFunction( const communication::BufferedCommunicator::LocalCommunicationMode & localCommunicationMode,
                                   const std::string & meshFile )
{
  const uint_t level = 3;

  auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

  auto x = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "x", storage, level, level );
  auto y = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "y", storage, level, level );
  auto z = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "z", storage, level, level );

  x->getCommunicator( level )->setLocalCommunicationMode( localCommunicationMode );
  y->getCommunicator( level )->setLocalCommunicationMode( localCommunicationMode );

  std::function< real_t( const hhg::Point3D & ) > expr = []( const hhg::Point3D & xx ) -> real_t { return real_c( (1.0L/2.0L)*sin(2*xx[0])*sinh(xx[1]) ) * real_c( xx[2] ); };
  std::function< real_t( const hhg::Point3D & ) > ones = []( const hhg::Point3D &    ) -> real_t { return 1.0; };

  x->interpolate( expr, level );

  x->getCommunicator( level )->template communicate< Vertex, Edge >();
  x->getCommunicator( level )->template communicate< Edge, Face >();
  x->getCommunicator( level )->template communicate< Face, Cell >();

  y->assign( { 7.0 }, { x.get() }, level );
  y->add( { 6.0 }, { x.get() }, level );

  y->getCommunicator( level )->template communicate< Vertex, Edge >();
  y->getCommunicator( level )->template communicate< Edge, Face >();
  y->getCommunicator( level )->template communicate< Face, Cell >();

  for ( const auto & cellIt : storage->getCells() )
  {
    const Cell & cell = *cellIt.second;
    const auto xData = cell.getData( x->getCellDataID() )->getPointer( level );
    const auto yData = cell.getData( y->getCellDataID() )->getPointer( level );
    for ( const auto & it : vertexdof::macrocell::Iterator( level ) )
    {
      const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
      WALBERLA_CHECK_FLOAT_EQUAL( 13.0 * xData[ idx ], yData[ idx ] );
    }
  }

  z->interpolate( ones, level );
  const real_t zScalarProduct = z->dot( *z, level );
  WALBERLA_CHECK_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(), 1, "Test only works with 1 process currently." )
  WALBERLA_CHECK_FLOAT_EQUAL( zScalarProduct, real_c( z->getNumLocalDoFs( level ) ) );

  // *****************
  // Apply

  // Apply test I
  // Simple check if apply is working at all on macro-cells.
  // The apply test II below probably also covers this case but we want to keep old tests anyway right? :)
  {
    auto operatorHandling  = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Cell > >( level, level, vertexDoFMacroCellStencilMemorySize );
    PrimitiveDataID< StencilMemory< real_t >, Cell > cellOperatorID;
    storage->addCellData( cellOperatorID, operatorHandling, "cell operator" );

    auto src = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "src", storage, level, level );
    auto dst = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "dst", storage, level, level );

    src->getCommunicator( level )->setLocalCommunicationMode( localCommunicationMode );
    dst->getCommunicator( level )->setLocalCommunicationMode( localCommunicationMode );

    for ( const auto & cellIt : storage->getCells() )
    {
      auto operatorData = cellIt.second->getData( cellOperatorID )->getPointer( level );
      operatorData[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ]  = 1.0;
      operatorData[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BC ) ] = 2.0;
    }

    src->interpolate( ones, level );
    src->getCommunicator( level )->template communicate< Vertex, Edge >();
    src->getCommunicator( level )->template communicate< Edge, Face >();
    src->getCommunicator( level )->template communicate< Face, Cell >();

    for ( const auto & cellIt : storage->getCells() )
    {
      vertexdof::macrocell::apply< real_t >( level, *cellIt.second, cellOperatorID, src->getCellDataID(), dst->getCellDataID(), UpdateType::Replace );

      auto dstData = cellIt.second->getData( dst->getCellDataID() )->getPointer( level );
      for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
      {
        const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
        WALBERLA_CHECK_FLOAT_EQUAL( dstData[ idx ], 3.0 );
      }
    }
  }

  // Apply test II
  // 1. Set all stencil weights of all primitives to 1.0
  // 2. Interpolate 1.0 everywhere
  // 3. Apply the operator
  // Expected result: each point's value equals its stencil size
  // If either lower to upper or upper to lower primitive-dimension communication is defect,
  // the result should be wrong at the borders.
  {
    auto op = std::make_shared< P1ZeroOperator >( storage, level, level );

    for ( const auto & it : storage->getVertices() )
    {
      StencilMemory< real_t > * vertexStencil = it.second->getData( op->getVertexStencilID() );
      WALBERLA_CHECK_EQUAL( vertexStencil->getSize( level ), it.second->getNumNeighborEdges() + 1 );
      for ( uint_t i = 0; i < it.second->getNumNeighborEdges() + 1; i++ ) vertexStencil->getPointer( level )[ i ] = 1.0;
    }
    for ( const auto & it : storage->getEdges() )
    {
      StencilMemory< real_t > * edgeStencil = it.second->getData( op->getEdgeStencilID() );
      WALBERLA_CHECK_EQUAL( edgeStencil->getSize( level ), 3 + 2 * it.second->getNumNeighborFaces() );
      for ( uint_t i = 0; i < 3 + 2 * it.second->getNumNeighborFaces(); i++ ) edgeStencil->getPointer( level )[ i ] = 1.0;
    }
    for ( const auto & it : storage->getFaces() )
    {
      StencilMemory< real_t > * faceStencil = it.second->getData( op->getFaceStencilID() );
      WALBERLA_CHECK_EQUAL( faceStencil->getSize( level ), 7 + 4 * it.second->getNumNeighborCells() );
      for ( uint_t i = 0; i < 7 + 4 * it.second->getNumNeighborCells(); i++ ) faceStencil->getPointer( level )[ i ] = 1.0;
    }
    for ( const auto & it : storage->getCells() )
    {
      StencilMemory< real_t > * cellStencil = it.second->getData( op->getCellStencilID() );
      WALBERLA_CHECK_EQUAL( cellStencil->getSize( level ), 15 );
      for ( uint_t i = 0; i < 15; i++ ) cellStencil->getPointer( level )[ i ] = 1.0;
    }

    auto src = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "src", storage, level, level );
    auto dst = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "dst", storage, level, level );

    src->getCommunicator( level )->setLocalCommunicationMode( localCommunicationMode );
    dst->getCommunicator( level )->setLocalCommunicationMode( localCommunicationMode );

    src->interpolate( ones, level );
    op->apply( *src, *dst, level, DoFType::All );

    for ( const auto & it : storage->getVertices() )
    {
      auto vertexDst = it.second->getData( dst->getVertexDataID() )->getPointer( level );
      WALBERLA_CHECK_FLOAT_EQUAL( vertexDst[ 0 ], real_c( it.second->getNumNeighborEdges() + 1 ) );
    }
    for ( const auto & it : storage->getEdges() )
    {
      auto edgeDst = it.second->getData( dst->getEdgeDataID() )->getPointer( level );
      for ( const auto & idxIt : vertexdof::macroedge::Iterator( level, 1 ) )
      {
        WALBERLA_CHECK_FLOAT_EQUAL( edgeDst[vertexdof::macroedge::indexFromVertex( level, idxIt.x(), stencilDirection::VERTEX_C )], real_c( 3 + 2 * it.second->getNumNeighborFaces() ) );
      }
    }
    for ( const auto & it : storage->getFaces() )
    {
      auto faceDst = it.second->getData( dst->getFaceDataID() )->getPointer( level );
      for ( const auto & idxIt : vertexdof::macroface::Iterator( level, 1 ) )
      {
        WALBERLA_CHECK_FLOAT_EQUAL( faceDst[vertexdof::macroface::indexFromVertex( level, idxIt.x(),
                                                                                   idxIt.y(),
                                                                                   stencilDirection::VERTEX_C )], real_c( 7 + 4 * it.second->getNumNeighborCells() ) );
      }
    }
    for ( const auto & it : storage->getCells() )
    {
      auto cellDst = it.second->getData( dst->getCellDataID() )->getPointer( level );
      for ( const auto & idxIt : vertexdof::macrocell::Iterator( level, 1 ) )
      {
        WALBERLA_CHECK_FLOAT_EQUAL( cellDst[vertexdof::macrocell::indexFromVertex( level, idxIt.x(), idxIt.y(), idxIt.z(), stencilDirection::VERTEX_C )], 15.0 );
      }
    }
  }
}

} // namespace hhg


int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  hhg::testVertexDoFFunction( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/3D/tet_1el.msh" );
  hhg::testVertexDoFFunction( communication::BufferedCommunicator::LocalCommunicationMode::DIRECT      , "../../data/meshes/3D/tet_1el.msh" );

  return EXIT_SUCCESS;
}
