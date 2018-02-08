#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"


namespace hhg {

static void testP2Smooth() {
  const uint_t sourceLevel = 2;

  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");

  SetupPrimitiveStorage setupStorage(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  auto x = std::make_shared<P2Function<real_t> >("x", storage, sourceLevel, sourceLevel + 1);

  std::function<real_t(const hhg::Point3D &)> zeros = [](const hhg::Point3D &) { return 0; };

  typedef stencilDirection sD;

  x->interpolate(zeros, sourceLevel);

  real_t* edgeDoFFineData = storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(sourceLevel + 1);
  real_t* vertexDoFFineData = storage->getFace(PrimitiveID(6))->getData(x->getVertexDoFFunction()->getFaceDataID())->getPointer(sourceLevel + 1);

  std::vector<std::pair<uint_t, uint_t > > vertical = { {1,0},{2,0},{1,1},{2,1},{1,2}};

  /// VERTICAL EDGE ///

  for(auto p : vertical) {

    storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(
      sourceLevel)[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(p.first, p.second, stencilDirection::EDGE_VE_N)] = 16.0;

    x->prolongate(sourceLevel, hhg::All);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::EDGE_VE_N)], 12.,p.first << " " << p.second);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::EDGE_VE_S)], 12.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::EDGE_HO_E)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::EDGE_HO_W)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::EDGE_DI_NW)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::EDGE_DI_SE)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::EDGE_VE_SE)], 4.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::EDGE_VE_NW)], 4.);

    WALBERLA_CHECK_FLOAT_EQUAL(vertexDoFFineData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::VERTEX_C)], 16.);

    storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(
      sourceLevel)[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(p.first, p.second, stencilDirection::EDGE_VE_N)] = 0.0;


  }
  //////////////////////////

  /// OTHER VERTICAL EDGE ///

  storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(sourceLevel)[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(1,1,stencilDirection::EDGE_VE_N)] = 0.0;
  storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(sourceLevel)[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(0,0,stencilDirection::EDGE_VE_N)] = 16.0;
  x->prolongate(sourceLevel, hhg::All);


//  WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1>(0,1,sD::EDGE_VE_N)],12.);
//  WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1>(0,1,sD::EDGE_VE_S)],12.);
//
//  WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1>(0,1,sD::EDGE_HO_E)],8.);
//  WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1>(0,1,sD::EDGE_HO_W)],8.);
//
//  WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1>(0,1,sD::EDGE_DI_NW)],8.);
//  WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1>(0,1,sD::EDGE_DI_SE)],8.);
//
//  WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1>(0,1,sD::EDGE_VE_SE)],4.);
//  WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1>(0,1,sD::EDGE_VE_NW)],4.);
//
//  WALBERLA_CHECK_FLOAT_EQUAL(vertexDoFFineData[hhg::vertexdof::macroface::indexFromVertex< sourceLevel + 1>(0,1,sD::VERTEX_C)],16.);

  //////////////////////////


  WALBERLA_LOG_INFO_ON_ROOT(storage->getFace(PrimitiveID(6))->getData(x->getVertexDoFFunction()->getFaceDataID())->getPointer(sourceLevel + 1)[hhg::vertexdof::macroface::indexFromVertex< sourceLevel + 1 >(2,3,stencilDirection::VERTEX_C)]);

  for (auto &faceIT : storage->getFaces()) {
    auto face = faceIT.second;
    hhg::vertexdof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getVertexDoFFunction()->getFaceDataID());
  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hhg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }

  for (auto &faceIT : storage->getFaces()) {
    auto face = faceIT.second;
    hhg::edgedof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getEdgeDoFFunction()->getFaceDataID());
  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hhg::edgedof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getEdgeDoFFunction()->getEdgeDataID());
//  }



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