#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hhg {

static void testOperator() {
  const uint_t level = 3;

  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");

  SetupPrimitiveStorage setupStorage(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  auto vertexDof = std::make_shared<P1Function<real_t> >("vertexDof", storage, level, level);
  auto edgeDof = std::make_shared<EdgeDoFFunction<real_t> >("edgeDof", storage, level, level);

  GenericVertexDoFToEdgeDoFOperator vertexToEdgeOperator(storage, level, level);

  std::function<real_t(const hhg::Point3D &)> ones = [](const hhg::Point3D &) { return 1; };
  std::function<real_t(const hhg::Point3D &,const std::vector<real_t> &)> onesVec = [](const hhg::Point3D &,const std::vector<real_t> &) { return 1; };

  for (auto &faceIT : storage->getFaces()) {
    auto face = faceIT.second;
    const uint_t faceStencilSize = face->getData(vertexToEdgeOperator.getFaceStencilID())->getSize(level);
    real_t *faceStencilData = face->getData(vertexToEdgeOperator.getFaceStencilID())->getPointer(level);
    for (uint_t i = 0; i < faceStencilSize; ++i) {
      faceStencilData[i] = 1;
    }
    hhg::vertexdof::macroface::interpolateTmpl< real_t, level >(*face,vertexDof->getFaceDataID(),{},onesVec);
  }

  for (auto &edgeIT : storage->getEdges()) {
    auto edge = edgeIT.second;
    const uint_t edgeStencilSize = edge->getData(vertexToEdgeOperator.getEdgeStencilID())->getSize(level);
    real_t *edgeStencilData = edge->getData(vertexToEdgeOperator.getEdgeStencilID())->getPointer(level);
    for (uint_t i = 0; i < edgeStencilSize; ++i) {
      edgeStencilData[i] = 1;
    }
    hhg::vertexdof::macroedge::interpolateTmpl< real_t, level >(*edge,vertexDof->getEdgeDataID(),{},onesVec);
  }

  for (auto &vertexIT : storage->getVertices()) {
    auto vertex = vertexIT.second;
    hhg::vertexdof::macrovertex::interpolate< real_t >(*vertex,vertexDof->getVertexDataID(),{},onesVec,level);
  }

  vertexToEdgeOperator.apply(*vertexDof, *edgeDof, level, hhg::All);

  auto edgeDoFcommunicator = edgeDof->getCommunicator( level );

  // Pull all halos
  edgeDoFcommunicator->communicate< Edge, Face >();
  edgeDoFcommunicator->communicate< Face, Edge >();
  edgeDoFcommunicator->communicate< Edge, Vertex >();
  edgeDoFcommunicator->communicate< Vertex, Edge >();



//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hhg::edgedof::macroface::printFunctionMemory<real_t, level>(*face, edgeDof->getFaceDataID());
//  }
//
//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hhg::edgedof::macroedge::printFunctionMemory<real_t, level>(*edge, edgeDof->getEdgeDataID());
//  }


  for (auto &faceIT : storage->getFaces()) {
    auto face = faceIT.second;
    const uint_t size = face->getData(edgeDof->getFaceDataID())->getSize(level);
    real_t *data = face->getData(edgeDof->getFaceDataID())->getPointer(level);
    for (uint_t i = 0; i < size; ++i) {
      ///the values on boundary can be 4 or 3 depending wether there is an adjacent face or not
      ///this check could be better but would be much more complicated
      WALBERLA_CHECK(walberla::floatIsEqual(data[i],4.0) || walberla::floatIsEqual(data[i],3.0));
    }
  }

  for (auto &edgeIT : storage->getEdges()) {
    auto edge = edgeIT.second;
    const uint_t size = edge->getData(edgeDof->getEdgeDataID())->getSize(level);
    real_t *data = edge->getData(edgeDof->getEdgeDataID())->getPointer(level);
    for (uint_t i = 0; i < size; ++i) {
      ///this check could also be imporoved
      WALBERLA_CHECK(walberla::floatIsEqual(data[i],4.0) || walberla::floatIsEqual(data[i],3.0));
    }
    hhg::vertexdof::macroedge::interpolateTmpl< real_t, level >(*edge,vertexDof->getEdgeDataID(),{},onesVec);
  }

}

}// namespace hhg

int main( int argc, char* argv[] ) {
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel(walberla::logging::Logging::PROGRESS);
  walberla::MPIManager::instance()->useWorldComm();
  hhg::testOperator();

  return EXIT_SUCCESS;
}