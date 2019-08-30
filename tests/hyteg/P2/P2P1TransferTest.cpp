#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/VTKWriter.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

static void testP2P1Transfer()
{
  const uint_t level = 4;

  MeshInfo mesh  = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  auto p1Function = std::make_shared< P1Function< real_t > >( "p1Function", storage, level, level );
  auto p2Function = std::make_shared< P2Function< real_t > >( "p2Function", storage, level, level );

  VTKOutput vtkOutput("../../output", "P2P1TransferTest", storage);
  vtkOutput.add( *p1Function );
  vtkOutput.add( *p2Function );

  // To test the transfer we
  // 1. set one vertex unknown in the middle of a P1 macro-face to testValue (rest 0.0)
  // 2. we prolongate to P2 and should get testValue at the corresponding vertex unknown and 0.5 * testValue at all (six, directly) neighboring edge unknowns
  // 3. we restrict again and should get a 2.5 * testValue at the corresponding vertex unknown and 0.25 * testValue at all neighboring vertex unknowns

  // Step 1:

  std::function< real_t( const Point3D & )> zeros = []( const Point3D & ) { return 0; };
  p1Function->interpolate( zeros, level );
  p2Function->interpolate( zeros, level );

  const uint_t x = 4;
  const uint_t y = 4;
  const real_t testValue = 1.0;

  std::vector< PrimitiveID > faceIDs;
  storage->getFaceIDs( faceIDs );

  WALBERLA_CHECK_EQUAL( storage->getNumberOfLocalFaces(), 1 );
  WALBERLA_CHECK_EQUAL( faceIDs.size(),                   1 );

  const auto p1FaceDataID = p1Function->getFaceDataID();
        auto p1FaceData   = storage->getFace( faceIDs[ 0 ] )->getData( p1FaceDataID )->getPointer( level );

  const uint_t idx = vertexdof::macroface::indexFromVertex( level, x, y, stencilDirection::VERTEX_C );

  p1FaceData[ idx ] = testValue;

  vtkOutput.write( level, 1 );

  // Step 2:

  p2Function->prolongateP1ToP2( *p1Function, level );

  const auto p2VertexDoFFaceDataID = p2Function->getVertexDoFFunction().getFaceDataID();
  const auto p2EdgeDoFFaceDataID   = p2Function->getEdgeDoFFunction().getFaceDataID();

        auto p2VertexDoFFaceData   = storage->getFace( faceIDs[ 0 ] )->getData( p2VertexDoFFaceDataID )->getPointer( level );
        auto p2EdgeDoFFaceData     = storage->getFace( faceIDs[ 0 ] )->getData( p2EdgeDoFFaceDataID )->getPointer( level );

  WALBERLA_CHECK_FLOAT_EQUAL( p2VertexDoFFaceData[ idx ], testValue )

  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_HO_W )], 0.5 * testValue )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_HO_E )], 0.5 * testValue )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_DI_NW )], 0.5 * testValue )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_DI_SE )], 0.5 * testValue )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_VE_N )], 0.5 * testValue )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_VE_S )], 0.5 * testValue )

  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_HO_NW )], 0.0 )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_HO_SE )], 0.0 )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_DI_SW )], 0.0 )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_DI_NE )], 0.0 )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_VE_NW )], 0.0 )
  WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_VE_SE )], 0.0 )

  vtkOutput.write( level, 2 );

  // Step 3:

  p2Function->restrictP2ToP1( *p1Function, level );

  WALBERLA_CHECK_FLOAT_EQUAL( p1FaceData[ idx ], 2.5 * testValue );

  for ( const auto & neighbor : vertexdof::macroface::neighborsWithoutCenter )
  {
    WALBERLA_CHECK_FLOAT_EQUAL( p1FaceData[vertexdof::macroface::indexFromVertex( level, x, y,
                                                                                  neighbor )], 0.25 * testValue );
  }

  vtkOutput.write( level, 3 );

}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testP2P1Transfer();

   return EXIT_SUCCESS;
}
