#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1PackInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/p1functionspace/P1EdgeIndex.hpp"
#include "tinyhhg_core/p1functionspace/P1FaceIndex.hpp"
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

  hhg::P1Function x("x", storage, minLevel, maxLevel);

  size_t num = 1;
  x.enumerate(maxLevel,num);

  uint_t numberOfChecks = 0;
  uint_t totalExpectedChecks = (3 * hhg::levelinfo::num_microvertices_per_edge(maxLevel)
                                + 3 * (hhg::levelinfo::num_microvertices_per_edge(maxLevel) - 1))
                               * storage->getNumberOfLocalFaces();

  for (auto &faceIt : storage->getFaces()) {
    Face &face = *faceIt.second;
    //BubbleFace::printFunctionMemory<maxLevel>(face,x.getFaceDataID());
    using namespace P1Edge::EdgeCoordsVertex;
    real_t *faceData = face.getData(x.getFaceDataID())->data[maxLevel].get();
    std::vector<PrimitiveID> nbrEdges;
    face.getNeighborEdges(nbrEdges);
    for(uint_t i = 0; i < nbrEdges.size(); ++i){
      Edge* edge = storage->getEdge(nbrEdges[0].getID());
      real_t* edgeData = edge->getData(x.getEdgeDataID())->data[maxLevel].get();
      uint_t idxCounter = 0;
      uint_t faceIdOnEdge = edge->face_index(face.getID());
//////////////////// INNER VERTEX //////////////////////
      idxCounter = 0;
      auto it = P1Face::indexIterator(face.edge_index(edge->getID()),
                                      face.edge_orientation[face.edge_index(edge->getID())],
                                      P1Face::VERTEX_INNER,
                                      maxLevel);
      for(; it != P1Face::indexIterator(); ++it){
        if(faceIdOnEdge == 0) {
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, VERTEX_SE)], faceData[*it]);
          numberOfChecks++;
        } else if(faceIdOnEdge == 1){
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, VERTEX_N)], faceData[*it]);
          numberOfChecks++;
        } else{
          WALBERLA_CHECK(false);
        }
        idxCounter++;
      }
//////////////////// VERTEX //////////////////////
      idxCounter = 0;
      it = P1Face::indexIterator(face.edge_index(edge->getID()),
                                 face.edge_orientation[face.edge_index(edge->getID())],
                                 P1Face::VERTEX,
                                 maxLevel);
      for(; it != P1Face::indexIterator(); ++it){
        WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, VERTEX_C)], faceData[*it]);
        numberOfChecks++;
        idxCounter++;
      }
    }
  }

  WALBERLA_CHECK_EQUAL(totalExpectedChecks,numberOfChecks);


  numberOfChecks = 0;
  totalExpectedChecks = 0;
  for(auto &vertexIt : storage->getVertices())
  {
    //totalExpectedChecks++;
    totalExpectedChecks += vertexIt.second->getNumHigherDimNeighbors() * 2;
  }

  for (auto &edgeIt : storage->getEdges()) {
    using namespace P1Edge::EdgeCoordsVertex;
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
        WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, 0, VERTEX_C )],
                             vertexData[0]);
        numberOfChecks++;

        WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, 1, VERTEX_C)],
                             vertexData[1 + vertex->edge_index(edge.getID())]);
        numberOfChecks++;

      } else if( vertexIdOnEdge == 1){
        WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, vPerEdge - 1, VERTEX_C )],
                             vertexData[0]);
        numberOfChecks++;

        WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, vPerEdge - 2, VERTEX_C)],
                             vertexData[1 + vertex->edge_index(edge.getID())]);
        numberOfChecks++;

      } else {
        WALBERLA_CHECK(false);
      }
    }
  }

  WALBERLA_CHECK_EQUAL(totalExpectedChecks,numberOfChecks);




  return 0;
}
