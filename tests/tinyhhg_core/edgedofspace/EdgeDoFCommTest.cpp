#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/indexing/VertexDoFIndexing.hpp"

#include "core/mpi/all.h"
#include "core/debug/all.h"

using namespace hhg;

using walberla::real_t;

template<uint_t level>
void checkComm(std::string meshfile, bool bufferComm = false){

  MeshInfo meshInfo = MeshInfo::fromGmshFile(meshfile);
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);


  //const uint_t level = 4;
  hhg::EdgeDoFFunction< uint_t > x("x", storage, level, level);
  if(bufferComm) {
    x.getCommunicator(level).get()->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
  }

  size_t num = 1;
  x.enumerate(level,num);

  uint_t numberOfChecks = 0;
  uint_t totalExpectedChecks = 0;

  for(auto &edgeIt : storage->getEdges()){
    if(edgeIt.second.get()->getNumHigherDimNeighbors() == 1){

      totalExpectedChecks += 4;
    } else if(edgeIt.second.get()->getNumHigherDimNeighbors() == 2){
      totalExpectedChecks += 8;
    } else {
      WALBERLA_CHECK(false);
    }
  }

  for (auto &faceIt : storage->getFaces()) {
    Face &face = *faceIt.second;
    uint_t *faceData = face.getData(x.getFaceDataID())->getPointer(level);
    std::vector<PrimitiveID> nbrEdges;
    face.getNeighborEdges(nbrEdges);
    for (uint_t i = 0; i < nbrEdges.size(); ++i) {
      Edge *edge = storage->getEdge(nbrEdges[0].getID());
      uint_t *edgeData = edge->getData(x.getEdgeDataID())->getPointer(level);
      uint_t idxCounter = 0;
      uint_t faceIdOnEdge = edge->face_index(face.getID());
      idxCounter = 0;
      auto it = hhg::indexing::edgedof::macroface::BorderIterator< level,0 >(indexing::getFaceBorderDirection(0,face.edge_orientation[0]));

    }
  }

  WALBERLA_CHECK_EQUAL(totalExpectedChecks,numberOfChecks);

}

int main (int argc, char ** argv ) {

  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();

  checkComm<4>("../../data/meshes/quad_4el.msh", true);

//  checkComm("../../data/meshes/quad_4el.msh", 5, true);
//
//  checkComm("../../data/meshes/quad_4el.msh", 4, false);
//
//  checkComm("../../data/meshes/quad_4el.msh", 5, false);
//
//  checkComm("../../data/meshes/bfs_12el.msh", 3, true);
//
//  checkComm("../../data/meshes/bfs_12el.msh", 3, false);

  return 0;


}