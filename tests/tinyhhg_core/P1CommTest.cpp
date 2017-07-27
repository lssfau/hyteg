#include "tinyhhg_core/p1functionspace/P1PackInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "core/mpi/all.h"
#include "core/debug/all.h"

using namespace hhg;
using walberla::real_t;
using namespace hhg::P1Edge;

int main (int argc, char ** argv )
{

  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();

  MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  PrimitiveStorage storage(uint_c(walberla::mpi::MPIManager::instance()->rank()), setupStorage);


  const uint_t minLevel = 2;
  const uint_t maxLevel = 3;

  //size_t v_perFace =  hhg::levelinfo::num_microvertices_per_face(maxLevel);
  size_t vPerEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
  //size_t nbr_v_perEdge = v_perEdge - 1;
  //size_t v_perVertex = hhg::levelinfo::num_microvertices_per_vertex(maxLevel);

  //hhg::P1Function x("x", mesh, minLevel, maxLevel);
//
//  const uint_t id = x.memory_id;
//  size_t num = 1;
//  x.enumerate(maxLevel,num);
//
//  using namespace hhg::P1BubbleEdge::EdgeCoordsVertex;
//  //check vertex to edge communication; vertex center and cell dofs are communicated
//  for(Vertex& vertex : mesh.vertices){
//    auto& vertexData = getVertexFunctionMemory(vertex, id)->data[maxLevel];
//    for(Edge* edge : vertex.edges){
//      auto& edgeData = getEdgeFunctionMemory(*edge, id)->data[maxLevel];
//      uint_t vOnEdgeID = edge->vertex_index(vertex);
//      //check vertex center
//      if(vOnEdgeID == 0) WALBERLA_CHECK_EQUAL(vertexData[0],edgeData[0]);
//      if(vOnEdgeID == 1) WALBERLA_CHECK_EQUAL(vertexData[0],edgeData[vPerEdge-1]);
//      //check cells
//      for(uint_t i = 0; i < vertex.faces.size(); ++i){
//        uint_t fOnEdgeID = edge->face_index(*vertex.faces[i]);
//        if(fOnEdgeID == 0 && vOnEdgeID == 0) {
//          WALBERLA_CHECK_EQUAL(vertexData[vertex.edges.size() + 1 + i],
//                               edgeData[edge_index(maxLevel, 0, CELL_GRAY_SE)])
//        } else if(fOnEdgeID == 0 && vOnEdgeID == 1){
//          WALBERLA_CHECK_EQUAL(vertexData[vertex.edges.size() + 1 + i],
//                               edgeData[edge_index(maxLevel, vPerEdge - 1, CELL_GRAY_SW)])
//        } else if(fOnEdgeID == 1 && vOnEdgeID == 0){
//          WALBERLA_CHECK_EQUAL(vertexData[vertex.edges.size() + 1 + i],
//                               edgeData[edge_index(maxLevel, 0, CELL_GRAY_NE)])
//        } else if(fOnEdgeID == 1 && vOnEdgeID == 1){
//          WALBERLA_CHECK_EQUAL(vertexData[vertex.edges.size() + 1 + i],
//                               edgeData[edge_index(maxLevel, vPerEdge - 1, CELL_GRAY_NW)])
//        }
//      }
//    }
//  }
//
//  //check edge to vertex comm; edge vertex Dofs are communicated
//  for(Edge& edge: mesh.edges){
//    auto& edgeData = getEdgeFunctionMemory(edge, id)->data[maxLevel];
//    WALBERLA_CHECK_EQUAL(edgeData[1],
//                         getVertexFunctionMemory(*edge.v0,id)->data[maxLevel][edge.v0->edge_index(edge) + 1]);
//    WALBERLA_CHECK_EQUAL(edgeData[vPerEdge-2],
//                         getVertexFunctionMemory(*edge.v1,id)->data[maxLevel][edge.v1->edge_index(edge) + 1]);
//  }
//
//  //check edge to face comm; edge vertex Dofs are communicated
//  for(Edge& edge: mesh.edges) {
//    auto &edgeData = getEdgeFunctionMemory(edge, id)->data[maxLevel];
//    for(Face* face: edge.faces){
//      auto &faceData = getFaceFunctionMemory(*face, id)->data[maxLevel];
//      uint_t idxCounter = 0;
//      for(auto it = P1BubbleFace::indexIterator(face->edge_index(edge),
//                                                face->edge_orientation[face->edge_index(edge)],
//                                                P1BubbleFace::VERTEX,
//                                                maxLevel); it != P1BubbleFace::indexIterator(); ++it){
//        WALBERLA_CHECK_EQUAL(edgeData[idxCounter],faceData[*it]);
//        idxCounter++;
//      }
//    }
//  }
//
//  //check face to edge comm; face inner vertex Dofs and Cells are communicated
//  for(Face& face: mesh.faces) {
//    using namespace P1BubbleEdge::EdgeCoordsVertex;
//    auto &faceData = getFaceFunctionMemory(face, id)->data[maxLevel];
//    for(Edge* edge: face.edges){
//      auto &edgeData = getEdgeFunctionMemory(*edge, id)->data[maxLevel];
//      uint_t idxCounter = 0;
//      uint_t faceIdOnEdge = edge->face_index(face);
///////////////////// VERTEX CELL /////////////////////
//      for(auto it = P1BubbleFace::indexIterator(face.edge_index(*edge),
//                                                face.edge_orientation[face.edge_index(*edge)],
//                                                P1BubbleFace::VERTEX_INNER,
//                                                maxLevel); it != P1BubbleFace::indexIterator(); ++it){
//        if(faceIdOnEdge == 0) {
//          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, VERTEX_SE)], faceData[*it]);
//        } else if(faceIdOnEdge == 1){
//          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, VERTEX_N)], faceData[*it]);
//        } else{
//          WALBERLA_CHECK(false);
//        }
//        idxCounter++;
//      }
////////////////////// GRAY CELL //////////////////////
//      idxCounter = 0;
//      for(auto it = P1BubbleFace::indexIterator(face.edge_index(*edge),
//                                                face.edge_orientation[face.edge_index(*edge)],
//                                                P1BubbleFace::CELL_GRAY,
//                                                maxLevel); it != P1BubbleFace::indexIterator(); ++it){
//        if(faceIdOnEdge == 0) {
//          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, CELL_GRAY_SE)], faceData[*it]);
//        } else if(faceIdOnEdge == 1){
//          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, CELL_GRAY_NE)], faceData[*it]);
//        } else{
//          WALBERLA_CHECK(false);
//        }
//        idxCounter++;
//      }
////////////////////// BLUE CELL //////////////////////
//      idxCounter = 0;
//      for(auto it = P1BubbleFace::indexIterator(face.edge_index(*edge),
//                                                face.edge_orientation[face.edge_index(*edge)],
//                                                P1BubbleFace::CELL_BLUE,
//                                                maxLevel); it != P1BubbleFace::indexIterator(); ++it){
//        if(faceIdOnEdge == 0) {
//          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter + 1, CELL_BLUE_SE)], faceData[*it]);
//        } else if(faceIdOnEdge == 1){
//          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter + 1, CELL_BLUE_NW)], faceData[*it]);
//        } else{
//          WALBERLA_CHECK(false);
//        }
//        idxCounter++;
//      }
//    }
//  }

  return 0;
}
