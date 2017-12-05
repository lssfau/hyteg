#include "core/mpi/all.h"
#include "core/debug/all.h"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"

using namespace hhg;

class dummyUFCOperator;

//template<typename OperatorType>
void checkOperator() {
  MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr <PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  const uint_t level = 2;

  EdgeDoFToVertexDoFOperator testOperator(storage,level,level);


}


int main(int argc, char** argv){
  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();

  checkOperator();
  //checkOperator<EdgeDoFToVertexDoFOperator>();
}