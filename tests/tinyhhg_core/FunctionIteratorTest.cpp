#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/DataTypes.h"

#include "tinyhhg_core/FunctionIterator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"

using walberla::uint_t;

namespace hhg {

static void testFunctionIterator()
{

  const std::string meshFileName = "../../data/meshes/3D/pyramid_2el.msh";

  const uint_t level = 3;

  auto storage = PrimitiveStorage::createFromGmshFile( meshFileName );

  P1Function< real_t > x( "x", storage, level, level );
  P1Function< int    > n( "n", storage, level, level );

  n.enumerate( level );

  auto storageInfo = storage->getGlobalInfo();
  WALBERLA_LOG_INFO_ON_ROOT( storageInfo );

  std::set< uint_t > testEnumerateSet;

  WALBERLA_LOG_INFO_ON_ROOT( "global dofs: " << numberOfGlobalDoFs< P1FunctionTag >( *storage, level ) );

  for ( const auto & dof : FunctionIterator< P1Function< int > >( n, level ) )
  {
    // WALBERLA_CHECK( dof.isVertexDoF() )
    auto value = getDoFValueFromFunction( n, dof );

    WALBERLA_LOG_DEVEL( dof );
    WALBERLA_LOG_DEVEL( value );

    WALBERLA_CHECK_EQUAL( testEnumerateSet.count( value ), 0 );
    testEnumerateSet.insert( value );
  }
  WALBERLA_LOG_INFO_ON_ROOT( testEnumerateSet.size() );

}

} // namespace hhg


int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  hhg::testFunctionIterator();

  return EXIT_SUCCESS;
}
