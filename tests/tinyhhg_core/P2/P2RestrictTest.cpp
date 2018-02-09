#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"


namespace hhg {

///this test check restrict on specific points on the grid
static void testP2Smooth() {
  const uint_t sourceLevel = 3;

  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");

  SetupPrimitiveStorage setupStorage(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  std::shared_ptr <PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  auto x = std::make_shared < P2Function < real_t > > ("x", storage, sourceLevel-1, sourceLevel);

  uint_t num = 1;
  x->enumerate(sourceLevel,num);

  for (auto &faceIT : storage->getFaces()) {
    auto face = faceIT.second;
    hhg::edgedof::macroface::printFunctionMemory<real_t, sourceLevel>(*face, x->getEdgeDoFFunction()->getFaceDataID());
  }

    for (auto &faceIT : storage->getFaces()) {
    auto face = faceIT.second;
    hhg::edgedof::macroface::printFunctionMemory<real_t, sourceLevel>(*face, x->getEdgeDoFFunction()->getFaceDataID());
  }


//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hhg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }

  x->restrict(sourceLevel,hhg::All);



//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hhg::vertexdof::macroface::printFunctionMemory<real_t, sourceLevel - 1>(*face, x->getVertexDoFFunction()->getFaceDataID());
//  }
//
//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hhg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel-1>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }

  ///calculate expected entry for idx = 1 on edge 0
  real_t expected = 5;
  expected += (1./8.) * (-46 + 3 * 47 + 3 * 48 - 49);
  expected += (1./8.) * (-98 + 3 * 99 + 3 *127 - 128);
  expected += (1./8.) * (-70 - 72);
  expected += (1./8.) * (-105 - 133 - 106 - 134);

  WALBERLA_CHECK_FLOAT_EQUAL(storage->getEdge(PrimitiveID(3))->getData(x->getVertexDoFFunction()->getEdgeDataID())->getPointer(sourceLevel - 1)[1], expected);


  ///calculate expected entry for idx = 2 on edge 1
  expected = 21;
  expected += (1./8.) * (- 64 + 3 * 65 + 3 * 66 - 67);
  expected += (1./8.) * (- 143 + 3 * 147 + 3 * 91 - 94);
  expected += (1./8.) * (- 115 - 122);
  expected += (1./8.) * (- 142 - 146 - 86 - 90);

  WALBERLA_CHECK_FLOAT_EQUAL(storage->getEdge(PrimitiveID(5))->getData(x->getVertexDoFFunction()->getEdgeDataID())->getPointer(sourceLevel - 1)[2], expected);

  ///calculate expected entry for idx = 2 on edge 1
  expected = 16;
  expected += (1./8.) * (- 61 + 3 * 60 + 3 * 59 - 58);
  expected += (1./8.) * (- 120 + 3 * 123 + 3 * 95 - 97);
  expected += (1./8.) * (- 148 - 153);
  expected += (1./8.) * (- 124 - 121 - 96 - 93);

  WALBERLA_CHECK_FLOAT_EQUAL(storage->getEdge(PrimitiveID(4))->getData(x->getVertexDoFFunction()->getEdgeDataID())->getPointer(sourceLevel - 1)[3], expected);


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