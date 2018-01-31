#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"


namespace hhg {

static void testP2Smooth() {
  const uint_t sourceLevel = 3;

  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/quad_2el.msh");

  SetupPrimitiveStorage setupStorage(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  std::shared_ptr <PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  auto x = std::make_shared < P2Function < real_t > > ("x", storage, sourceLevel-1, sourceLevel);

  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1; };


  x->interpolate(ones,sourceLevel);


  x->restrict(sourceLevel,hhg::All);

  for (auto &edgeIT : storage->getEdges()) {
    auto edge = edgeIT.second;
    hhg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel-1>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
  }

  //hhg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel-1>(*storage->getEdge(PrimitiveID(6)), x->getVertexDoFFunction()->getEdgeDataID());


}

}/// namespace hhg

int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  hhg::testP2Smooth();

  return EXIT_SUCCESS;
}