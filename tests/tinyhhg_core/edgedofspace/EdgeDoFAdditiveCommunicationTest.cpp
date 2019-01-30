#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/FunctionTraits.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

namespace hhg {

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;


static void testEdgeDoFAdditiveCommunication2D( const communication::BufferedCommunicator::LocalCommunicationMode & localCommunicationMode,
                                                const std::string & meshFile )
{
  const uint_t level = 3;
  const real_t testValue = 1.0;
  const real_t someConstant = 6.345;

  auto meshInfo = MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupPrimitiveStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  auto storage = std::make_shared< PrimitiveStorage >( setupPrimitiveStorage );

  EdgeDoFFunction< real_t > x( "x", storage, level, level );
  x.setLocalCommunicationMode( localCommunicationMode );

  // 1. To avoid masking, interpolate any scalar != 0 everywhere
  x.interpolate( real_c( someConstant ), level, All );

  // 2. Set all entries (incl. ghost layers) on all cells to a test value.
  for ( const auto & it : storage->getFaces() )
  {
    auto facePtr = it.second->getData( x.getFaceDataID() )->getPointer( level );
    for ( const auto & idxIt : edgedof::macroface::Iterator( level, 0 ) )
    {
      for ( auto orientation : edgedof::faceLocalEdgeDoFOrientations )
      {
        const auto idx = edgedof::macroface::index( level, idxIt.x(), idxIt.y(), orientation );
        facePtr[ idx ] = testValue;
      }
    }
  }

  // 3. Communicate additively. Each unknown on each primitive should now be equal to the number of neighbor faces times the test value.
  x.communicateAdditively< Face, Edge >( level );

  for ( const auto & it : storage->getEdges() )
  {
    auto edgePtr = it.second->getData( x.getEdgeDataID() )->getPointer( level );
    for ( const auto & idxIt : edgedof::macroedge::Iterator( level, 0 ) )
    {
      const auto idx = edgedof::macroedge::index( level, idxIt.x() );
      if ( x.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
      {
        WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx], someConstant, "Face -> Edge additive comm failed (on Dirichlet boundary)" );
      }
      else
      {
        WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx], testValue * real_c( it.second->getNumNeighborFaces() ), "Face -> Edge additive comm failed (inner or Neumann)" );
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

  hhg::testEdgeDoFAdditiveCommunication2D( hhg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/tri_1el.msh" );
  hhg::testEdgeDoFAdditiveCommunication2D( hhg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,       "../../data/meshes/tri_1el.msh" );
  hhg::testEdgeDoFAdditiveCommunication2D( hhg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/tri_2el.msh" );
  hhg::testEdgeDoFAdditiveCommunication2D( hhg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,       "../../data/meshes/tri_2el.msh" );
  hhg::testEdgeDoFAdditiveCommunication2D( hhg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/quad_8el.msh" );
  hhg::testEdgeDoFAdditiveCommunication2D( hhg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,       "../../data/meshes/quad_8el.msh" );
  hhg::testEdgeDoFAdditiveCommunication2D( hhg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/annulus_coarse.msh" );
  hhg::testEdgeDoFAdditiveCommunication2D( hhg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,       "../../data/meshes/annulus_coarse.msh" );

  return EXIT_SUCCESS;
}
