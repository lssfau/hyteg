#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/tinyhhg.hpp"

#include "core/mpi/all.h"
#include "core/debug/all.h"

using namespace hhg;

using walberla::real_t;

void checkComm(std::string meshfile,const uint_t level, bool bufferComm = false){

  MeshInfo meshInfo = MeshInfo::fromGmshFile(meshfile);
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  hhg::EdgeDoFFunction< uint_t > x("x", storage, level, level);
  if(bufferComm) {
    x.getCommunicator(level).get()->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
  }

  size_t num = 0;
  uint_t check = 0;
  uint_t sum = 0;
  x.enumerate(level,num);

  uint_t totalDoFs = (hhg::levelinfo::num_microfaces_per_face(level) * storage->getNumberOfLocalFaces());
  uint_t expectedSum = (totalDoFs * (totalDoFs + 1))/2;

//  for (auto &vertexIt : storage->getVertices()){
//    Vertex &vertex = *vertexIt.second;
//    uint_t *vertexData = vertex.getData(x.getVertexDataID())->getPointer(level);
//    for(uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i){
//      WALBERLA_CHECK_EQUAL(vertexData[i*2],check);
//      sum += check;
//      check++;
//    }
//  }
//
//
//
//  for (auto &edgeIt : storage->getEdges()) {
//    Edge &edge = *edgeIt.second;
//    uint_t *edgeData = edge.getData(x.getEdgeDataID())->getPointer(level);
//    uint_t FaceDoFonFace = hhg::levelinfo::num_microvertices_per_edge(level) * 2 - 3 ;
//    //this only works with the linear default layout; can be changed to use index function
//    for(uint_t i = 0; i < edge.getNumHigherDimNeighbors(); ++i){
//      uint_t start = FaceDoFonFace * i;
//      for(uint_t j = 2; j < hhg::levelinfo::num_microvertices_per_edge(level) * 2 - 5; j += 2){
//        WALBERLA_CHECK_EQUAL(edgeData[start + j],check);
//        sum += check;
//        check++;
//      }
//    }
//  }

  for ( auto &faceIt : storage->getFaces() ) {
    size_t idxCounter = 0;
    Face &face = *faceIt.second;
    uint_t *faceData = face.getData(x.getFaceDataID())->getPointer(level);
    uint_t rowsize = levelinfo::num_microedges_per_edge(level);
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
  }

  --check;
  WALBERLA_CHECK_EQUAL(check,totalDoFs);
  WALBERLA_CHECK_EQUAL(sum,expectedSum);

}


int main (int argc, char ** argv ) {

  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();

  checkComm("../../data/meshes/tri_1el.msh",3);

}

