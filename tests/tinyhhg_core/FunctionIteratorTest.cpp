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
  const std::string meshFileName = "../../data/meshes/3D/cube_24el.msh";

  const uint_t level = 3;

  auto storage = PrimitiveStorage::createFromGmshFile( meshFileName );

  P1Function< real_t > x( "x", storage, level, level );
  P1Function< int    > n( "n", storage, level, level );

  n.enumerate( level );

  auto storageInfo = storage->getGlobalInfo();
  WALBERLA_LOG_INFO_ON_ROOT( storageInfo );

  std::set< uint_t > testEnumerateSet;

  const uint_t numGlobalDoFs = numberOfGlobalDoFs< P1FunctionTag >( *storage, level );
  const uint_t numLocalDoFs  = numberOfLocalDoFs< P1FunctionTag >( *storage, level );
  WALBERLA_LOG_INFO_ON_ROOT( "global dofs: " << numGlobalDoFs );
  WALBERLA_LOG_INFO( "local dofs: " << numLocalDoFs );

  for ( const auto & dof : FunctionIterator< P1Function< int > >( n, level ) )
  {
    WALBERLA_CHECK( dof.isVertexDoF() )
    WALBERLA_LOG_INFO( dof );
    WALBERLA_CHECK_EQUAL( testEnumerateSet.count( dof.value() ), 0 );
    testEnumerateSet.insert( dof.value() );
  }
  WALBERLA_CHECK_EQUAL( testEnumerateSet.size(), numLocalDoFs );

}

} // namespace hhg


int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  hhg::testFunctionIterator();

  return EXIT_SUCCESS;
}
