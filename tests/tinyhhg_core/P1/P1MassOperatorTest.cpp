// test that the product one^T*M*one with mass matrix M and vector of ones gives area of domain
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

using walberla::real_t;
using walberla::uint_t;
using namespace hhg;

void checkArea( std::shared_ptr<PrimitiveStorage> storage, real_t area )
{

  const size_t minLevel = 2;
  const size_t maxLevel = 4;

  P1ConstantMassOperator massOp( storage, minLevel, maxLevel );

  P1Function< real_t > aux( "aux", storage, minLevel, maxLevel );
  P1Function< real_t > vecOfOnes( "vecOfOnes", storage, minLevel, maxLevel );
  std::function< real_t( const Point3D& ) > ones = []( const Point3D& ) { return 1.0; };

  for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
    {
      vecOfOnes.interpolate( ones, lvl, All );
      massOp.apply( vecOfOnes, aux, lvl, All );
      real_t measure = vecOfOnes.dotGlobal( aux, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "measure = " << std::scientific << measure );
      WALBERLA_CHECK_FLOAT_EQUAL( measure, area );
    }
}

int main(int argc, char **argv)
{
  walberla::debug::enterTestMode();

  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  // Test with rectangle
  WALBERLA_LOG_INFO_ON_ROOT( "Testing with RECTANGLE" );
  MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {0.0, -1.0} ), Point2D( {2.0, 3.0} ),
                                               MeshInfo::CRISSCROSS, 1, 2 );
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);
  checkArea( storage, 8.0 );

  // Test with backward facing step
  WALBERLA_LOG_INFO_ON_ROOT( "Testing with BFS" );
  meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
  SetupPrimitiveStorage setupStorageBFS( meshInfo,
                                         uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr<PrimitiveStorage> storageBFS = std::make_shared<PrimitiveStorage>( setupStorageBFS );
  checkArea( storageBFS, 1.75 );

  return EXIT_SUCCESS;
}
