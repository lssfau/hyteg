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

  real_t * vertexStencil = face->getData(p2operator.getVertexToVertexOpr().getFaceStencilID())->getPointer( level );
  real_t vertexStencilSize = face->getData(p2operator.getVertexToVertexOpr().getFaceStencilID())->getSize( level );

  vertexStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)] =  1 ;
  vertexStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N )] = -1 ;
  vertexStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W )] = -2;
  vertexStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C )] = 1  ;
  vertexStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E )] = 2;
  vertexStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S )] = -1 ;
  vertexStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)] =  1 ;


  std::function<real_t(const hhg::Point3D&)> rhsfunc = [](const hhg::Point3D&) { return 1; };


  real_t* faceMemory = face->getData(x->getVertexDoFFunction()->getFaceDataID())->getPointer( level );
  for ( const auto & it : vertexdof::macroface::Iterator( level, 0 ) ){
    faceMemory[hhg::vertexdof::macroface::indexFromVertex< level >(it.col(),it.row(), stencilDirection::VERTEX_C)] = 1;

  }

  rhs->interpolate(rhsfunc,level);

  vertexdof::macroface::printFunctionMemory< real_t, level >(*(storage->getFaces().begin()->second),x->getVertexDoFFunction()->getFaceDataID());

  P2::face::smoothGSvertexDoFTmpl< level >(*face, p2operator.getVertexToVertexOpr().getFaceStencilID(),
                                           x->getVertexDoFFunction()->getFaceDataID(),
                            p2operator.getEdgeToVertexOpr().getFaceStencilID(),
                            x->getEdgeDoFFunction()->getFaceDataID(),
                            rhs->getVertexDoFFunction()->getFaceDataID());

  vertexdof::macroface::printFunctionMemory< real_t, level >(*(storage->getFaces().begin()->second),x->getVertexDoFFunction()->getFaceDataID());
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