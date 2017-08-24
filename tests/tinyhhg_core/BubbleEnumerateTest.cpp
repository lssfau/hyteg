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

  uint_t numberOfChecks = 0;
  uint_t totalExpectedChecks = (2 * hhg::levelinfo::num_microvertices_per_edge(maxLevel) - 3) * 3 * storage->getNumberOfLocalFaces();

  for (auto &faceIt : storage->getFaces()) {
    Face &face = *faceIt.second;
    //BubbleFace::printFunctionMemory<maxLevel>(face,x.getFaceDataID());
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
          numberOfChecks++;
        } else if(faceIdOnEdge == 1){
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, CELL_GRAY_NE)], faceData[*it]);
          numberOfChecks++;
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
          numberOfChecks++;
        } else if(faceIdOnEdge == 1){
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter + 1, CELL_BLUE_NW)], faceData[*it]);
          numberOfChecks++;
        } else{
          WALBERLA_CHECK(false);
        }
        idxCounter++;
      }
    }
  }

  WALBERLA_CHECK_EQUAL(totalExpectedChecks,numberOfChecks);

  numberOfChecks = 0;
  totalExpectedChecks = 0;
  for(auto &edgeIt : storage->getEdges()){
    if(edgeIt.second.get()->getNumHigherDimNeighbors() == 1){
      totalExpectedChecks += 2;
    } else if(edgeIt.second.get()->getNumHigherDimNeighbors() == 2){
      totalExpectedChecks += 4;
    } else {
      WALBERLA_CHECK(false);
    }
  }


  for (auto &edgeIt : storage->getEdges()) {
    using namespace BubbleEdge::EdgeCoordsVertex;
    Edge &edge = *edgeIt.second;
    //BubbleEdge::printFunctionMemory(edge,x.getEdgeDataID(),maxLevel);
    real_t *edgeData = edge.getData(x.getEdgeDataID())->data[maxLevel].get();
    std::vector<PrimitiveID> nbrVertices;
    edge.getNeighborVertices(nbrVertices);
    for(uint_t i = 0; i < nbrVertices.size(); ++i)
    {
      Vertex* vertex = storage->getVertex(nbrVertices[i].getID());
      real_t* vertexData = vertex->getData(x.getVertexDataID())->data[maxLevel].get();
      uint_t vertexIdOnEdge = edge.vertex_index(vertex->getID());
      uint_t vPerEdge = levelinfo::num_microvertices_per_edge(maxLevel);
      if(vertexIdOnEdge == 0){
        WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, 0, CELL_GRAY_SE )],
                             vertexData[vertex->face_index(edge.neighborFaces()[0])]);
        numberOfChecks++;
        if(edge.getNumHigherDimNeighbors() == 2)
        {
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, 0, CELL_GRAY_NE)],
                               vertexData[vertex->face_index(edge.neighborFaces()[1])]);
          numberOfChecks++;
        }
      } else if( vertexIdOnEdge == 1){
        WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, vPerEdge - 1, CELL_GRAY_SW )],
                             vertexData[vertex->face_index(edge.neighborFaces()[0])]);
        numberOfChecks++;
        if(edge.getNumHigherDimNeighbors() == 2)
        {
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, vPerEdge - 1, CELL_GRAY_NW)],
                               vertexData[vertex->face_index(edge.neighborFaces()[1])]);
          numberOfChecks++;
        }
      } else {
        WALBERLA_CHECK(false);
      }
    }
  }

  WALBERLA_CHECK_EQUAL(totalExpectedChecks,numberOfChecks);

  return 0;
}
