// Test generation of MeshInfo object for inline mesh generators
//
// Simple test:
// - check that MeshInfo object can be generated
// - use it to setup a PrimitiveStorage
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::math::PI;
using walberla::uint_t;

typedef enum { RECTANGLE, ANNULUS, PARTIAL_ANNULUS, SPHERICAL_SHELL, FACE_CHAIN } testDomainType;

namespace hhg {

static void testMeshGenerator( testDomainType testDomain, MeshInfo::meshFlavour flavour )
{

  std::shared_ptr< MeshInfo > meshInfo;

  // perform mesh generation
  switch( testDomain )
    {
    case RECTANGLE:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshRectangle( Point2D( {-2.0, 1.0} ), Point2D( {0.0, 3.0} ),
                                                        flavour, 3, 2 ) );
      break;

    case PARTIAL_ANNULUS:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshAnnulus( 1.0, 2.0, 0.25*PI, 0.75*PI, flavour, 4, 2 ) );
      break;

    case ANNULUS:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshAnnulus( 1.0, 2.0, 4, 2 ) );
      break;

    case SPHERICAL_SHELL:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshSphericalShell( 5, 3, 1.0, 2.0 ) );
      break;

    case FACE_CHAIN:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshFaceChain( 3 ) );

    }

  // test primitive counts
  uint_t numVerts = meshInfo->getVertices().size();
  uint_t numFaces = meshInfo->getFaces().size();
  uint_t numEdges = meshInfo->getEdges().size();
  uint_t numCells = meshInfo->getCells().size();

  if( testDomain == RECTANGLE )
    {
      switch( flavour )
       {
       case hhg::MeshInfo::CRISS :
       case hhg::MeshInfo::CROSS :
         WALBERLA_CHECK_EQUAL( numVerts, 12 );
         WALBERLA_CHECK_EQUAL( numFaces, 12 );
         WALBERLA_CHECK_EQUAL( numEdges, 23 );
         WALBERLA_CHECK_EQUAL( numCells,  0 );
         break;

       case hhg::MeshInfo::CRISSCROSS :
       case hhg::MeshInfo::DIAMOND :
         WALBERLA_CHECK_EQUAL( numVerts, 18 );
         WALBERLA_CHECK_EQUAL( numFaces, 24 );
         WALBERLA_CHECK_EQUAL( numEdges, 41 );
         WALBERLA_CHECK_EQUAL( numCells,  0 );
         break;
       }
    }
  else if( testDomain == FACE_CHAIN )
    {
      WALBERLA_CHECK_EQUAL( numVerts, 5 );
      WALBERLA_CHECK_EQUAL( numEdges, 7 );
      WALBERLA_CHECK_EQUAL( numFaces, 3 );
      WALBERLA_CHECK_EQUAL( numCells, 0 );
    }
  else if( testDomain == SPHERICAL_SHELL )
    {
      WALBERLA_CHECK_EQUAL( numVerts,  486 );
      WALBERLA_CHECK_EQUAL( numFaces, 4160 );
      WALBERLA_CHECK_EQUAL( numCells, 1920 );
    }

  // use generated MeshInfo
  SetupPrimitiveStorage setupStorage( *meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage );
  WALBERLA_LOG_INFO_ON_ROOT( setupStorage );
  PrimitiveStorage storage( setupStorage );

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Check mesh generator for rectangles
   hhg::testMeshGenerator( RECTANGLE, hhg::MeshInfo::CRISS );
   hhg::testMeshGenerator( RECTANGLE, hhg::MeshInfo::CROSS );
   hhg::testMeshGenerator( RECTANGLE, hhg::MeshInfo::CRISSCROSS );
   hhg::testMeshGenerator( RECTANGLE, hhg::MeshInfo::DIAMOND );

   // Check mesh generator for annuli
   hhg::testMeshGenerator( ANNULUS, hhg::MeshInfo::CROSS );
   hhg::testMeshGenerator( PARTIAL_ANNULUS, hhg::MeshInfo::CRISS );
   hhg::testMeshGenerator( PARTIAL_ANNULUS, hhg::MeshInfo::CROSS );
   hhg::testMeshGenerator( PARTIAL_ANNULUS, hhg::MeshInfo::CRISSCROSS );
   hhg::testMeshGenerator( PARTIAL_ANNULUS, hhg::MeshInfo::DIAMOND );

   // Check mesh generator for face chains
   // (note that flavour argument is meaningless here)
   hhg::testMeshGenerator( FACE_CHAIN, hhg::MeshInfo::DIAMOND );

   // Check mesh generator for thick spherical shell
   // (note that flavour argument is meaningless here)
   hhg::testMeshGenerator( SPHERICAL_SHELL, hhg::MeshInfo::DIAMOND );

   return EXIT_SUCCESS;
}
