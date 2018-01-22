#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"


namespace hhg {

static void testP2Smooth() {
  const uint_t level = 3;

  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");

  SetupPrimitiveStorage setupStorage(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  std::shared_ptr <PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  auto x = std::make_shared < P2Function < real_t > > ("x", storage, level, level);
  auto rhs = std::make_shared < P2Function < real_t > > ("rhs", storage, level, level);

  P2ConstantLaplaceOperator p2operator( storage, level, level );

  auto face = storage->getFaces().begin()->second.get();

  real_t * vertexToVertexStencil = face->getData(p2operator.getVertexToVertexOpr().getFaceStencilID())->getPointer( level );
  real_t * edgeToVertexStencil = face->getData(p2operator.getEdgeToVertexOpr().getFaceStencilID())->getPointer( level );

  real_t * edgeToEdgeStencil = face->getData(p2operator.getEdgeToEdgeOpr().getFaceStencilID())->getPointer( level );
  real_t * vertexToEdgeStencil = face->getData(p2operator.getVertexToEdgeOpr().getFaceStencilID())->getPointer( level );
  //real_t edgeTovertexStencilSize = face->getData(p2operator.getEdgeToVertexOpr().getFaceStencilID())->getSize( level );

  /// vertex dofs
  for(uint_t k = 0; k < vertexdof::macroface::neighborsWithCenter.size(); ++k){
    vertexToVertexStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithCenter[k])] = 1;
  }
  vertexToVertexStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C )] = 19.   ;

  for(uint_t k = 0; k < edgedof::macroface::neighborsFromVertex.size(); ++k){
    edgeToVertexStencil[edgedof::stencilIndexFromVertex(edgedof::macroface::neighborsFromVertex[k])] = 1;
  }

  /// horizontal edges
  for(uint_t k = 0; k < edgedof::macroface::neighborsFromHorizontalEdge.size(); ++k){
    edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(edgedof::macroface::neighborsFromHorizontalEdge[k])] = 1;
  }
  edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)] = 9;

  for(uint_t k = 0; k < vertexdof::macroface::neighborsFromHorizontalEdge.size(); ++k){
    vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(vertexdof::macroface::neighborsFromHorizontalEdge[k])] = 1;
  }

  /// diagonal edges
  for(uint_t k = 0; k < edgedof::macroface::neighborsFromDiagonalEdge.size(); ++k){
    edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge(edgedof::macroface::neighborsFromDiagonalEdge[k])] = 1;
  }
  edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_DI_C)] = 9;

  for(uint_t k = 0; k < vertexdof::macroface::neighborsFromDiagonalEdge.size(); ++k){
    vertexToEdgeStencil[vertexdof::stencilIndexFromDiagonalEdge(vertexdof::macroface::neighborsFromDiagonalEdge[k])] = 1;
  }


  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1; };

  x->interpolate(ones,level);
  rhs->interpolate(ones,level);

  vertexdof::macroface::printFunctionMemory< real_t, level >(*(storage->getFaces().begin()->second),x->getVertexDoFFunction()->getFaceDataID());

  P2::face::smoothGSvertexDoFTmpl< level >(*face, p2operator.getVertexToVertexOpr().getFaceStencilID(),
                                           x->getVertexDoFFunction()->getFaceDataID(),
                                           p2operator.getEdgeToVertexOpr().getFaceStencilID(),
                            x->getEdgeDoFFunction()->getFaceDataID(),
                            rhs->getVertexDoFFunction()->getFaceDataID());

  vertexdof::macroface::printFunctionMemory< real_t, level >(*(storage->getFaces().begin()->second),x->getVertexDoFFunction()->getFaceDataID());

  edgedof::macroface::printFunctionMemory< real_t, level >(*(storage->getFaces().begin()->second),x->getEdgeDoFFunction()->getFaceDataID());

  P2::face::smoothGSedgeDoFTmpl< level >(*face, p2operator.getVertexToEdgeOpr().getFaceStencilID(),
                                           x->getVertexDoFFunction()->getFaceDataID(),
                                           p2operator.getEdgeToEdgeOpr().getFaceStencilID(),
                                           x->getEdgeDoFFunction()->getFaceDataID(),
                                           rhs->getVertexDoFFunction()->getFaceDataID());

  edgedof::macroface::printFunctionMemory< real_t, level >(*(storage->getFaces().begin()->second),x->getEdgeDoFFunction()->getFaceDataID());

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