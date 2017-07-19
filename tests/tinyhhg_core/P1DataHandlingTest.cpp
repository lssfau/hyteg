
#include "tinyhhg_core/tinyhhg.hpp"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

namespace hhg {

static void testP1DataHandling()
{
  uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
  uint_t rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  std::string meshFileName = "../../data/meshes/tri_2el.msh";

  MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
  PrimitiveStorage      storage( rank, setupStorage );

  VertexP1FunctionMemoryDataHandling vertexP1FunctionMemoryDataHandling;
  EdgeP1FunctionMemoryDataHandling   edgeP1FunctionMemoryDataHandling;
  FaceP1FunctionMemoryDataHandling   faceP1FunctionMemoryDataHandling;

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testP1DataHandling();

   return EXIT_SUCCESS;
}
