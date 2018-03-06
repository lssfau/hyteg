#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"

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

  //x.enumerate(Level,num);

  uint_t totalDoFs = hhg::levelinfo::num_microedges_per_face( Level ) * storage->getNumberOfLocalFaces();

  for (auto &edgeIt : storage->getEdges()) {
    hhg::edgedof::macroedge::enumerateTmpl< uint_t,Level >(*edgeIt.second,x.getEdgeDataID(),num);
    Edge &edge = *edgeIt.second;
    uint_t *edgeData = edge.getData(x.getEdgeDataID())->getPointer(Level);

    if(edgeIt.second->getNumNeighborFaces() == 2) totalDoFs -= levelinfo::num_microedges_per_edge( Level );

    for(uint_t i = 0; i < levelinfo::num_microedges_per_edge( Level ); ++i){
      WALBERLA_CHECK_EQUAL(edgeData[i],check);
      WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_DI_S )], 0);
      WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_VE_SE )], 0);
      if( i != 0){
        WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_SE )], 0);
      }
      if(edgeIt.second->getNumNeighborFaces() == 2) {
        WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_DI_N )], 0);
        WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_VE_NW )], 0);
        if (i != 0) {
          WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_NW )], 0);
        }
      }
      check++;
    }
  }

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
          ++idxCounter;
          ++check;
        }
      }
      --inner_rowsize;
    }
  }



  //--check;
  WALBERLA_CHECK_EQUAL(check,totalDoFs);

}


int main (int argc, char ** argv ) {

  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();


  checkComm< 2 > ("../../data/meshes/quad_2el.msh");
  checkComm< 2 > ("../../data/meshes/quad_2el.msh");
  checkComm< 2 > ("../../data/meshes/bfs_12el.msh");

  checkComm< 3 > ("../../data/meshes/tri_1el.msh");
  checkComm< 3 > ("../../data/meshes/quad_2el.msh");
  checkComm< 3 > ("../../data/meshes/bfs_12el.msh");

  checkComm< 4 > ("../../data/meshes/tri_1el.msh");
  checkComm< 4 > ("../../data/meshes/quad_2el.msh");
  checkComm< 4 > ("../../data/meshes/bfs_12el.msh");


}

