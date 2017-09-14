
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

  P1Function< real_t > realP1Function ( "real_t P1Function", storage, minLevel, maxLevel );
  P1Function< uint_t > uintP1Function ( "uint_t P1Function", storage, minLevel, maxLevel );
  P1Function< int >    intP1Function  ( "int    P1Function", storage, minLevel, maxLevel );
  P1Function< float >  floatP1Function( "float  P1Function", storage, minLevel, maxLevel );

  BubbleFunction< real_t > realBubbleFunction ( "real_t BubbleFunction", storage, minLevel, maxLevel );
  BubbleFunction< uint_t > uintBubbleFunction ( "uint_t BubbleFunction", storage, minLevel, maxLevel );
  BubbleFunction< int >    intBubbleFunction  ( "int    BubbleFunction", storage, minLevel, maxLevel );
  BubbleFunction< float >  floatBubbleFunction( "float  BubbleFunction", storage, minLevel, maxLevel );

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
