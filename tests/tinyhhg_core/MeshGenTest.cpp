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
#include "tinyhhg_core/tinyhhg.hpp"

using walberla::math::PI;

typedef enum { RECTANGLE, ANNULUS, PARTIAL_ANNULUS } testDomainType;

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

   return EXIT_SUCCESS;
}
