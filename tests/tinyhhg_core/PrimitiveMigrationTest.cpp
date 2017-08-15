#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

static void testPrimitiveMigration()
{
  uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  const std::string meshFileName = "../../data/meshes/tri_2el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  loadbalancing::roundRobin( setupStorage );

  WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

  WALBERLA_LOG_INFO_ON_ROOT( "Building PrimitiveStorage" );

  std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( setupStorage ) );

  writeDomainPartitioningVTK( storage, "../../output/", "domain_decomposition_before_migration" );

  WALBERLA_MPI_SECTION()
  {
    std::vector< PrimitiveID > primitiveIDs;
    storage->getPrimitiveIDsGenerically< Primitive >( primitiveIDs );

    std::map< PrimitiveID::IDType, uint_t > migrationInfo;
    for ( const auto & id : primitiveIDs )
    {
      migrationInfo[ id.getID() ] = 0;
    }

    storage->migratePrimitives( migrationInfo );
  }

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testPrimitiveMigration();

   return EXIT_SUCCESS;
}
