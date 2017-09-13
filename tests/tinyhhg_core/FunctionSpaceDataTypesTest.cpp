
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

static void testFunctionSpaceDataTypes()
{
  const std::string meshFileName = "../../data/meshes/tri_2el.msh";

  const uint_t minLevel = 2;
  const uint_t maxLevel = 5;

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( setupStorage ) );

  P1Function< real_t > realFunction ( "real_t P1Function", storage, minLevel, maxLevel );
  P1Function< uint_t > uintFunction ( "uint_t P1Function", storage, minLevel, maxLevel );
  P1Function< int >    intFunction  ( "int    P1Function", storage, minLevel, maxLevel );
  P1Function< float >  floatFunction( "float  P1Function", storage, minLevel, maxLevel );

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testFunctionSpaceDataTypes();

   return EXIT_SUCCESS;
}
