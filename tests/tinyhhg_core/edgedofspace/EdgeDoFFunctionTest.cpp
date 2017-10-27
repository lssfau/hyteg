
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"

namespace hhg {

static void testEdgeDoFFunction()
{
  const uint_t minLevel = 2;
  const uint_t maxLevel = 4;

  MeshInfo mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, walberla::mpi::MPIManager::instance()->numProcesses() );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  EdgeDoFFunction< real_t > edgeDoFFunction( "x", storage, minLevel, maxLevel );
  WALBERLA_UNUSED( edgeDoFFunction );

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testEdgeDoFFunction();

   return EXIT_SUCCESS;
}
