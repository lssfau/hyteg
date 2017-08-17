#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

struct Data
{
  uint_t primitiveIDInData = 345678;
};

class DataHandling : PrimitiveDataHandling< Data, Primitive >
{
public:

	std::shared_ptr< Data > initialize( const Primitive * const primitive ) const
	{
      auto data = std::make_shared< Data >();
      data->primitiveIDInData = uint_c( primitive->getID().getID() );
      return data;
	}

	virtual void serialize( const Primitive * const primitive, const PrimitiveDataID< Data, Primitive > & id, SendBuffer & buffer ) const
	{
      auto data = primitive->getData( id );
      buffer << data->primitiveIDInData;
	}

	virtual void deserialize( const Primitive * const primitive, const PrimitiveDataID< Data, Primitive > & id, RecvBuffer & buffer ) const
	{
	  auto data = primitive->getData( id );
	  buffer >> data->primitiveIDInData;
	}
};

static void testPrimitiveMigration()
{
  uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  const std::string meshFileName = "../../data/meshes/tri_2el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  loadbalancing::roundRobin( setupStorage );

  WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

  WALBERLA_LOG_INFO_ON_ROOT( "Building PrimitiveStorage" );

  std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( setupStorage ) );

  writeDomainPartitioningVTK( storage, "../../output/", "domain_decomposition_before_migration" );

  PrimitiveDataID< Data, Primitive > dataID;
  auto dataHandling = std::make_shared< DataHandling >();
  storage->addPrimitiveData( dataID, dataHandling, "test data" );

  WALBERLA_MPI_SECTION()
  {
    std::vector< PrimitiveID > primitiveIDs;
    storage->getPrimitiveIDsGenerically< Primitive >( primitiveIDs );

    std::map< PrimitiveID::IDType, uint_t > migrationInfo;
    uint_t lel = 0;
    for ( const auto & id : primitiveIDs )
    {
      uint_t targetRank = ++lel % numProcesses;
      WALBERLA_LOG_INFO( "Migrating " << id.getID() << " to rank " << targetRank );
      migrationInfo[ id.getID() ] = targetRank;
    }

    storage->migratePrimitives( migrationInfo );

    PrimitiveStorage::PrimitiveMap primitives;
    storage->getPrimitives( primitives );
    for ( const auto & it : primitives )
    {
      WALBERLA_CHECK_EQUAL( it.first, it.second->getData( dataID )->primitiveIDInData );
    }

  }

  writeDomainPartitioningVTK( storage, "../../output/", "domain_decomposition_after_migration" );

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();
   hhg::testPrimitiveMigration();

   return EXIT_SUCCESS;
}
