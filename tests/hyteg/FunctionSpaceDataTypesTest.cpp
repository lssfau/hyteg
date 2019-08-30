
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/DataTypes.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

static void testFunctionSpaceDataTypes()
{
  const std::string meshFileName = "../../data/meshes/tri_2el.msh";

  const uint_t minLevel = 2;
  const uint_t maxLevel = 5;

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( setupStorage ) );

  P1Function< real_t > realP1Function ( "real_t P1Function", storage, minLevel, maxLevel );
  P1Function< int >    intP1Function  ( "int    P1Function", storage, minLevel, maxLevel );

}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testFunctionSpaceDataTypes();

   return EXIT_SUCCESS;
}
