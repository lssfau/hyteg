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
  //uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  const std::string meshFileName = "../../data/meshes/tri_2el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  const uint_t globalNumPrimitives = setupStorage.getNumberOfPrimitives();

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
    for ( const auto & id : primitiveIDs )
    {
      migrationInfo[ id.getID() ] = 0;
    }

    storage->migratePrimitives( migrationInfo );

    WALBERLA_ROOT_SECTION()
    {
      WALBERLA_CHECK_EQUAL( globalNumPrimitives, storage->getNumberOfLocalPrimitives() );
      PrimitiveStorage::PrimitiveMap primitives;
      storage->getPrimitives( primitives );
      for ( const auto & it : primitives )
      {
        WALBERLA_CHECK_EQUAL( it.first, it.second->getData( dataID )->primitiveIDInData );
      }
    }
    WALBERLA_NON_ROOT_SECTION()
    {
      WALBERLA_CHECK_EQUAL( 0, storage->getNumberOfLocalPrimitives() );
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
  walberla::debug::enterTestMode();
   hhg::testPrimitiveMigration();

   return EXIT_SUCCESS;
}
