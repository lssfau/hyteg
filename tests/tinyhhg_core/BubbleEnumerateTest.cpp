#include "tinyhhg_core/tinyhhg.hpp"

#include "core/mpi/all.h"
#include "core/debug/all.h"

using namespace hhg;
using walberla::real_t;

int main (int argc, char ** argv )
{

  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();

  MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);


  const uint_t minLevel = 2;
  const uint_t maxLevel = 4;

  hhg::BubbleFunction x("x", storage, minLevel, maxLevel);

  size_t num = 1;
  x.enumerate(maxLevel,num);

  //check face to edge comm; face inner vertex Dofs and Cells are communicated
  for (auto &faceIt : storage->getFaces()) {
    Face &face = *faceIt.second;

    using namespace BubbleEdge::EdgeCoordsVertex;
    real_t *faceData = face.getData(x.getFaceDataID())->data[maxLevel].get();
    std::vector<PrimitiveID> nbrEdges;
    face.getNeighborEdges(nbrEdges);
    for(uint_t i = 0; i < nbrEdges.size(); ++i){
      Edge* edge = storage->getEdge(nbrEdges[0].getID());
      real_t* edgeData = edge->getData(x.getEdgeDataID())->data[maxLevel].get();
      uint_t idxCounter = 0;
      uint_t faceIdOnEdge = edge->face_index(face.getID());
//////////////////// GRAY CELL //////////////////////
      idxCounter = 0;
      auto it = BubbleFace::indexIterator(face.edge_index(edge->getID()),
                                            face.edge_orientation[face.edge_index(edge->getID())],
                                            BubbleFace::CELL_GRAY,
                                            maxLevel);
      for(; it != BubbleFace::indexIterator(); ++it){
        if(faceIdOnEdge == 0) {
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, CELL_GRAY_SE)], faceData[*it]);
        } else if(faceIdOnEdge == 1){
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, CELL_GRAY_NE)], faceData[*it]);
        } else{
          WALBERLA_CHECK(false);
        }
        idxCounter++;
      }
//////////////////// BLUE CELL //////////////////////
      idxCounter = 0;
      it = BubbleFace::indexIterator(face.edge_index(edge->getID()),
                                            face.edge_orientation[face.edge_index(edge->getID())],
                                            BubbleFace::CELL_BLUE,
                                            maxLevel);
      for(; it != BubbleFace::indexIterator(); ++it){
        if(faceIdOnEdge == 0) {
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter + 1, CELL_BLUE_SE)], faceData[*it]);
        } else if(faceIdOnEdge == 1){
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter + 1, CELL_BLUE_NW)], faceData[*it]);
        } else{
          WALBERLA_CHECK(false);
        }
        idxCounter++;
      }
    }
  }

  return 0;
}
