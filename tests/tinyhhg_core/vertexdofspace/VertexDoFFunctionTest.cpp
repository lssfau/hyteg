
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/forms/P1ZeroForm.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/FunctionTraits.hpp"
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

  x->setLocalCommunicationMode( localCommunicationMode );
  y->setLocalCommunicationMode( localCommunicationMode );

  std::function< real_t( const hhg::Point3D & ) > expr = []( const hhg::Point3D & xx ) -> real_t { return real_c( (1.0L/2.0L)*sin(2*xx[0])*sinh(xx[1]) ) * real_c( xx[2] ); };
  std::function< real_t( const hhg::Point3D & ) > ones = []( const hhg::Point3D &    ) -> real_t { return 1.0; };

  x->interpolate( expr, level );

  x->communicate< Vertex, Edge >( level );
  x->communicate< Edge, Face >( level );
  x->communicate< Face, Cell >( level );

  y->assign( { 7.0 }, { *x }, level );
  y->add( { 6.0 }, { *x }, level );

  y->communicate< Vertex, Edge >( level );
  y->communicate< Edge, Face >( level );
  y->communicate< Face, Cell >( level );

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
  real_t zScalarProduct = z->dotGlobal( *z, level );
  WALBERLA_CHECK_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(), 1, "Test only works with 1 process currently." )
  WALBERLA_CHECK_FLOAT_EQUAL( zScalarProduct, real_c( numberOfLocalDoFs< P1FunctionTag >( *storage, level ) ) );

  z->add( real_c( 1 ), level, All );
  zScalarProduct = z->dotGlobal( *z, level );
  WALBERLA_CHECK_FLOAT_EQUAL( zScalarProduct, real_c( 4 * numberOfLocalDoFs< P1FunctionTag >( *storage, level ) ) );

  // *****************
  // Apply

  // Apply test I
  // Simple check if apply is working at all on macro-cells.
  // The apply test II below probably also covers this case but we want to keep old tests anyway right? :)
  {
    auto operatorHandling  = std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell > >( level, level );
    PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell > cellOperatorID;
    storage->addCellData( cellOperatorID, operatorHandling, "cell operator" );

    auto src = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "src", storage, level, level );
    auto dst = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "dst", storage, level, level );

    src->setLocalCommunicationMode( localCommunicationMode );
    dst->setLocalCommunicationMode( localCommunicationMode );

    for ( const auto & cellIt : storage->getCells() )
    {
      auto & operatorData = cellIt.second->getData( cellOperatorID )->getData( level );
      for ( const auto & neighbor : vertexdof::macrocell::neighborsWithCenter )
      {
        operatorData[ vertexdof::logicalIndexOffsetFromVertex( neighbor ) ] = 0.0;
      }
      operatorData[ vertexdof::logicalIndexOffsetFromVertex( stencilDirection::VERTEX_C ) ]  = 1.0;
      operatorData[ vertexdof::logicalIndexOffsetFromVertex( stencilDirection::VERTEX_TC ) ] = 2.0;
    }

    src->interpolate( ones, level );
    src->communicate< Vertex, Edge >( level );
    src->communicate< Edge, Face >( level );
    src->communicate< Face, Cell >( level );

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
      WALBERLA_CHECK_EQUAL( edgeStencil->getSize( level ), 3 + 2 * it.second->getNumNeighborFaces() + it.second->getNumNeighborCells() );
      for ( uint_t i = 0; i < 3 + 2 * it.second->getNumNeighborFaces() + it.second->getNumNeighborCells(); i++ ) edgeStencil->getPointer( level )[ i ] = 1.0;
    }
    for ( const auto & it : storage->getFaces() )
    {
      auto & faceStencil = it.second->getData( op->getFaceStencil3DID() )->getData( level );
      for ( uint_t neighborCellID = 0; neighborCellID < it.second->getNumNeighborCells(); neighborCellID++ )
      {
        P1ZeroForm zeroForm;
        auto face = it.second;
        auto neighborCell = storage->getCell( face->neighborCells().at( neighborCellID ) );
        auto vertexAssemblyIndexInCell =
        vertexdof::macroface::getIndexInNeighboringMacroCell( {1, 1, 0}, *face, neighborCellID, *storage, level );
        faceStencil[neighborCellID] = P1Elements::P1Elements3D::assembleP1LocalStencilNew(
           storage, *neighborCell, vertexAssemblyIndexInCell, level, zeroForm );

        for ( auto & stencilIt : faceStencil[neighborCellID] )
        {
          auto leafIndexInMacroFace = vertexdof::macrocell::getIndexInNeighboringMacroFace(
          vertexAssemblyIndexInCell + stencilIt.first, *neighborCell, neighborCell->getLocalFaceID( face->getID() ), *storage, level );
          if ( leafIndexInMacroFace.z() == 0 )
            stencilIt.second = 0.5;
          else
            stencilIt.second = 1.0;
        }
      }
    }
    for ( const auto & it : storage->getCells() )
    {
      auto & cellStencil = it.second->getData( op->getCellStencilID() )->getData( level );
      for ( const auto & neighbor : vertexdof::macrocell::neighborsWithCenter )
        cellStencil[ vertexdof::logicalIndexOffsetFromVertex( neighbor ) ] = 1.0;
    }

    auto src = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "src", storage, level, level );
    auto dst = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "dst", storage, level, level );

    src->setLocalCommunicationMode( localCommunicationMode );
    dst->setLocalCommunicationMode( localCommunicationMode );

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
        WALBERLA_CHECK_FLOAT_EQUAL( edgeDst[vertexdof::macroedge::indexFromVertex( level, idxIt.x(), stencilDirection::VERTEX_C )], real_c( 3 + 2 * it.second->getNumNeighborFaces() + it.second->getNumNeighborCells() ) );
      }
    }
    for ( const auto & it : storage->getFaces() )
    {
      // face stencil size:
      // on face: 7
      // on top and bottom: each 6
      // => 13 for one, 19 for two neighboring cells
      auto faceDst = it.second->getData( dst->getFaceDataID() )->getPointer( level );
      for ( const auto & idxIt : vertexdof::macroface::Iterator( level, 1 ) )
      {
        WALBERLA_CHECK_FLOAT_EQUAL( faceDst[vertexdof::macroface::indexFromVertex( level, idxIt.x(),
                                                                                   idxIt.y(),
                                                                                   stencilDirection::VERTEX_C )], 7.5 * real_c( it.second->getNumNeighborCells() ) );
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
  hhg::testVertexDoFFunction( hhg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/3D/tet_1el.msh" );
  hhg::testVertexDoFFunction( hhg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT      , "../../data/meshes/3D/tet_1el.msh" );

  return EXIT_SUCCESS;
}
