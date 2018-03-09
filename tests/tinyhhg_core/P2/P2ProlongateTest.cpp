#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"


namespace hhg {

/// this test checks if a constant value is prolongated correctly to all points on the finer level
static void testP2Prolongate() {
  const uint_t sourceLevel = 2;

  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");
  std::shared_ptr<SetupPrimitiveStorage> setupStorage =
    std::make_shared<SetupPrimitiveStorage>(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(*setupStorage);
  auto x = std::make_shared<P2Function<real_t> >("x", storage, sourceLevel, sourceLevel + 2);
  typedef stencilDirection sD;
  std::function<real_t(const hhg::Point3D &)> values = [](const hhg::Point3D &) { return 13; };
  x->interpolate(values, sourceLevel, hhg::All);
  x->prolongate(sourceLevel, hhg::All);

  real_t* edgeDoFFineData = storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(sourceLevel + 1);
  real_t* vertexDoFFineData = storage->getFace(PrimitiveID(6))->getData(x->getVertexDoFFunction()->getFaceDataID())->getPointer(sourceLevel + 1);

  for( const auto & it : hhg::vertexdof::macroface::Iterator( sourceLevel + 1, 1)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      vertexDoFFineData[hhg::vertexdof::macroface::indexFromVertex< sourceLevel + 1 >(it.col(), it.row(), sD::VERTEX_C)],
      13.,
      it.col() << " " << it.row());
  }

  for( const auto & it : hhg::edgedof::macroface::Iterator( sourceLevel + 1, 0)) {
    if(it.row() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1 >(it.col(), it.row(), sD::EDGE_HO_E)],
        13.,
        it.col() << " " << it.row());
    }
    if(it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( sourceLevel + 1 ) - 1)) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1 >(it.col(), it.row(), sD::EDGE_DI_NE)],
        13.,
        it.col() << " " << it.row());
    }
    if(it.col() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 1 >(it.col(), it.row(), sD::EDGE_VE_N)],
        13.,
        it.col() << " " << it.row());
    }
  }


  for (auto &edgeIT : storage->getEdges()) {
    auto edge = edgeIT.second;
    edgeDoFFineData = edge->getData(x->getEdgeDoFFunction()->getEdgeDataID())->getPointer(
      sourceLevel + 1);
    vertexDoFFineData = edge->getData(x->getVertexDoFFunction()->getEdgeDataID())->getPointer(
      sourceLevel + 1);

    for (const auto &it : hhg::vertexdof::macroedge::Iterator(sourceLevel + 1, 1)) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        vertexDoFFineData[hhg::vertexdof::macroedge::indexFromVertex<sourceLevel + 1>(it.col(), sD::VERTEX_C)],
        13.,
        it.col() << " " << it.row());
    }

    for (const auto &it : hhg::edgedof::macroedge::Iterator(sourceLevel + 1, 0)) {
        WALBERLA_CHECK_FLOAT_EQUAL(
          edgeDoFFineData[hhg::edgedof::macroedge::indexFromVertex<sourceLevel + 1>(it.col(), sD::EDGE_HO_E)],
          13.,
          it.col() << " " << it.row());
    }
  }


  x->prolongate(sourceLevel + 1, hhg::All);

  edgeDoFFineData = storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(sourceLevel + 2);
  vertexDoFFineData = storage->getFace(PrimitiveID(6))->getData(x->getVertexDoFFunction()->getFaceDataID())->getPointer(sourceLevel + 2);

  for( const auto & it : hhg::vertexdof::macroface::Iterator( sourceLevel + 2, 1)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      vertexDoFFineData[hhg::vertexdof::macroface::indexFromVertex< sourceLevel + 2 >(it.col(), it.row(), sD::VERTEX_C)],
      13.,
      it.col() << " " << it.row());
  }

  for( const auto & it : hhg::edgedof::macroface::Iterator( sourceLevel + 2, 0)) {
    if(it.row() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 2 >(it.col(), it.row(), sD::EDGE_HO_E)],
        13.,
        it.col() << " " << it.row());
    }
    if(it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( sourceLevel + 2 ) - 1)) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 2 >(it.col(), it.row(), sD::EDGE_DI_NE)],
        13.,
        it.col() << " " << it.row());
    }
    if(it.col() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel + 2 >(it.col(), it.row(), sD::EDGE_VE_N)],
        13.,
        it.col() << " " << it.row());
    }
  }


//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hhg::vertexdof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getVertexDoFFunction()->getFaceDataID());
//  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hhg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }

//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hhg::edgedof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getEdgeDoFFunction()->getFaceDataID());
//  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hhg::edgedof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getEdgeDoFFunction()->getEdgeDataID());
//  }

}


/// this test writes specific values at certain points and
/// checks wether these values are propagated correctly after one prolongation step
static void testP2Prolongate2() {

  const uint_t sourceLevel = 2;
  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");
  std::shared_ptr<SetupPrimitiveStorage> setupStorage =
    std::make_shared<SetupPrimitiveStorage>(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(*setupStorage);
  auto x = std::make_shared<P2Function<real_t> >("x", storage, sourceLevel, sourceLevel + 1);
  typedef stencilDirection sD;

  /// this should not be necessary but just to be save
  std::function<real_t(const hhg::Point3D &)> zeros = [](const hhg::Point3D &) { return 0; };
  x->interpolate(zeros, sourceLevel);

  real_t* edgeDoFFineData = storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(sourceLevel + 1);
  real_t* vertexDoFFineData = storage->getFace(PrimitiveID(6))->getData(x->getVertexDoFFunction()->getFaceDataID())->getPointer(sourceLevel + 1);

  /// all possible vertical edge Dof locations that need to be updated by the face
  std::vector<std::pair<uint_t, uint_t > > vertical   = { {1,0},{2,0},{3,0},{1,1},{2,1},{1,2} };
  /// all possible horizontal edge Dof locations that need to be updated by the face
  std::vector<std::pair<uint_t, uint_t > > horizontal = { {0,1},{1,1},{2,1},{0,2},{1,2},{0,3} };
  /// all possible diagonal edge Dof locations that need to be updated by the face
  std::vector<std::pair<uint_t, uint_t > > diagonal   = { {0,0},{1,0},{2,0},{0,1},{1,1},{0,2} };
  /// all possible vertex dof locations that need to be updated by the face
  std::vector<std::pair<uint_t, uint_t > > vertex     = { {1,1},{2,1},{1,2} };

///////////////////////////
/// CHECH VERTICAL EDGE ///
///////////////////////////
  for(auto p : vertical) {

    storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(
      sourceLevel)[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(p.first, p.second, stencilDirection::EDGE_VE_N)] = 16.0;

    x->prolongate(sourceLevel, hhg::All);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2, p.second *2 + 1, sD::EDGE_VE_N)], 12.);
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
/////////////////////////////
/// CHECK HORIZONTAL EDGE ///
/////////////////////////////
  for(auto p : horizontal) {

    storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(
      sourceLevel)[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(p.first, p.second, stencilDirection::EDGE_HO_E)] = 16.0;

    x->prolongate(sourceLevel, hhg::All);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2, sD::EDGE_VE_N)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2, sD::EDGE_VE_S)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2, sD::EDGE_HO_E)], 12.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2, sD::EDGE_HO_W)], 12.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2, sD::EDGE_DI_NW)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2, sD::EDGE_DI_SE)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2, sD::EDGE_HO_SE)], 4.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2, sD::EDGE_HO_NW)], 4.);

    WALBERLA_CHECK_FLOAT_EQUAL(vertexDoFFineData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2, sD::VERTEX_C)], 16.);

    storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(
      sourceLevel)[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(p.first, p.second, stencilDirection::EDGE_HO_E)] = 0.0;

  }
///////////////////////////
/// CHECK DIAGONAL EDGE ///
///////////////////////////
  for(auto p : diagonal) {

    storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(
      sourceLevel)[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(p.first, p.second, stencilDirection::EDGE_DI_NE)] = 16.0;

    x->prolongate(sourceLevel, hhg::All);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_VE_N)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_VE_S)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_HO_E)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_HO_W)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_DI_NW)], 12.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_DI_SE)], 12.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_DI_SW)], 4.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hhg::edgedof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_DI_NE)], 4.);

    WALBERLA_CHECK_FLOAT_EQUAL(vertexDoFFineData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 + 1, p.second * 2 + 1, sD::VERTEX_C)], 16.);

    storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction()->getFaceDataID())->getPointer(
      sourceLevel)[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(p.first, p.second, stencilDirection::EDGE_DI_NE)] = 0.0;

  }
////////////////////
/// CHECK VERTEX ///
////////////////////
  for( auto p : vertex ){
    storage->getFace(PrimitiveID(6))->getData(x->getVertexDoFFunction()->getFaceDataID())->getPointer(
      sourceLevel)[hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(p.first, p.second, stencilDirection::VERTEX_C)] = 16.0;


    x->prolongate(sourceLevel, hhg::All);

    WALBERLA_CHECK_FLOAT_EQUAL(vertexDoFFineData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel + 1>(p.first * 2 , p.second * 2 , sD::VERTEX_C)],
                               16.,
                               p.first << " " << p.second);

    storage->getFace(PrimitiveID(6))->getData(x->getVertexDoFFunction()->getFaceDataID())->getPointer(
      sourceLevel)[hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(p.first, p.second, stencilDirection::VERTEX_C)] = 0.0;

  }
///////////////////////////



//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hhg::vertexdof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getVertexDoFFunction()->getFaceDataID());
//  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hhg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }

//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hhg::edgedof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getEdgeDoFFunction()->getFaceDataID());
//  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hhg::edgedof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getEdgeDoFFunction()->getEdgeDataID());
//  }

}

static void testP2InterpolateAndProlongate() {

  const uint_t sourceLevel = 2;
  const uint_t targetLevel = sourceLevel + 3;
  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
  std::shared_ptr<SetupPrimitiveStorage> setupStorage =
    std::make_shared<SetupPrimitiveStorage>(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(*setupStorage);
  auto x = P2Function<real_t> ("x", storage, sourceLevel, targetLevel);
  auto y = P2Function<real_t> ("y", storage, sourceLevel, targetLevel);
  auto error = P2Function<real_t> ("x", storage, sourceLevel, targetLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return 2. * xx[0] * xx[0] + 3. * xx[0] + 13. + 4. * xx[1] + 5. *  xx[1] * xx[1]; };
  x.interpolate(exact,sourceLevel    ,hhg::All);
  y.interpolate(exact,targetLevel,hhg::All);

  x.prolongate(sourceLevel    ,hhg::All);
  x.prolongate(sourceLevel + 1,hhg::All);
  x.prolongate(sourceLevel + 2,hhg::All);

  error.assign({1.0, -1.0}, {&x, &y}, targetLevel, hhg::All);

  WALBERLA_CHECK_FLOAT_EQUAL(error.dot(error,targetLevel,hhg::All),0.);


}



}/// namespace hhg

int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  hhg::testP2Prolongate();
  hhg::testP2Prolongate2();
  hhg::testP2InterpolateAndProlongate();

  return EXIT_SUCCESS;
}