#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/tinyhhg.hpp"

#include "core/mpi/all.h"
#include "core/debug/all.h"

using namespace hhg;

using walberla::real_t;

template<uint_t Level>
void checkComm(std::string meshfile, bool bufferComm = false){

  MeshInfo meshInfo = MeshInfo::fromGmshFile(meshfile);
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  hhg::EdgeDoFFunction< uint_t > x("x", storage, Level, Level);
  if(bufferComm) {
    x.getCommunicator(Level).get()->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
  }

  size_t num = 0;
  uint_t check = 0;
  uint_t sum = 0;

  //x.enumerate(Level,num);

  uint_t totalDoFs = (hhg::levelinfo::num_microfaces_per_face(Level) * storage->getNumberOfLocalFaces());
  uint_t expectedSum = (totalDoFs * (totalDoFs + 1))/2;


  for ( auto &faceIt : storage->getFaces() ) {
    hhg::edgedof::macroface::enumerateTmpl< uint_t,Level >(*faceIt.second,x.getFaceDataID(),num);
    size_t idxCounter = 0;
    Face &face = *faceIt.second;
    uint_t *faceData = face.getData(x.getFaceDataID())->getPointer(Level);
    uint_t rowsize = levelinfo::num_microedges_per_edge(Level);
    uint_t inner_rowsize = rowsize;
    for(uint_t i = 0; i < rowsize ; ++i){
      for(uint_t j = 0; j <  inner_rowsize  ; ++j){
        if(i == 0) {
          WALBERLA_CHECK_EQUAL(faceData[idxCounter], 0, "idxCounter was: " << idxCounter);
          ++idxCounter;
        } else {
          WALBERLA_CHECK_EQUAL(faceData[idxCounter], check, "idxCounter was: " << idxCounter);
          sum += check;
          ++idxCounter;
          ++check;
        }
      }
      --inner_rowsize;
    }
    rowsize = levelinfo::num_microedges_per_edge(Level);
    inner_rowsize = rowsize;
    for(uint_t i = 0; i < rowsize ; ++i){
      for(uint_t j = 0; j <  inner_rowsize  ; ++j){
        if((j + i) == (levelinfo::num_microedges_per_edge( Level ) - 1)) {
          WALBERLA_CHECK_EQUAL(faceData[idxCounter], 0, "idxCounter was: " << idxCounter);
          ++idxCounter;
        } else {
          WALBERLA_CHECK_EQUAL(faceData[idxCounter], check, "idxCounter was: " << idxCounter);
          sum += check;
          ++idxCounter;
          ++check;
        }
      }
      --inner_rowsize;
    }
    rowsize = levelinfo::num_microedges_per_edge(Level);
    inner_rowsize = rowsize;
    for(uint_t i = 0; i < rowsize ; ++i){
      for(uint_t j = 0; j <  inner_rowsize  ; ++j){
        if(j == 0) {
          WALBERLA_CHECK_EQUAL(faceData[idxCounter], 0, "idxCounter was: " << idxCounter);
          ++idxCounter;
        } else {
          WALBERLA_CHECK_EQUAL(faceData[idxCounter], check, "idxCounter was: " << idxCounter);
          sum += check;
          ++idxCounter;
          ++check;
        }
      }
      --inner_rowsize;
    }
  }

  //  for (auto &edgeIt : storage->getEdges()) {
//    Edge &edge = *edgeIt.second;
//    uint_t *edgeData = edge.getData(x.getEdgeDataID())->getPointer(Level);
//    uint_t FaceDoFonFace = hhg::Levelinfo::num_microvertices_per_edge(Level) * 2 - 3 ;
//    //this only works with the linear default layout; can be changed to use index function
//    for(uint_t i = 0; i < edge.getNumHigherDimNeighbors(); ++i){
//      uint_t start = FaceDoFonFace * i;
//      for(uint_t j = 2; j < hhg::Levelinfo::num_microvertices_per_edge(Level) * 2 - 5; j += 2){
//        WALBERLA_CHECK_EQUAL(edgeData[start + j],check);
//        sum += check;
//        check++;
//      }
//    }
//  }

  --check;
  //WALBERLA_CHECK_EQUAL(check,totalDoFs);
  //WALBERLA_CHECK_EQUAL(sum,expectedSum);

}


int main (int argc, char ** argv ) {

  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();

  const uint_t Level = 3;

  checkComm< Level > ("../../data/meshes/tri_1el.msh");

}

