#include "tinyhhg_core/tinyhhg.hpp"

using namespace hhg;
using walberla::real_t;
using namespace hhg::P1Bubble;

int main (int argc, char ** argv )
{

  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();
  hhg::Mesh mesh("../../data/meshes/quad_4el.msh");


  const uint_t minLevel = 2;
  const uint_t maxLevel = 3;

  //size_t v_perFace =  hhg::levelinfo::num_microvertices_per_face(maxLevel);
  size_t vPerEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
  //size_t nbr_v_perEdge = v_perEdge - 1;
  //size_t v_perVertex = hhg::levelinfo::num_microvertices_per_vertex(maxLevel);

  hhg::P1BubbleFunction x("x", mesh, minLevel, maxLevel);

  const uint_t id = x.memory_id;
  size_t num = 1;
  x.enumerate(maxLevel,num);

  using namespace hhg::P1BubbleEdge::EdgeCoordsVertex;
  //check vertex to edge communication; vertex center and cell dofs are communicated
  for(Vertex& vertex : mesh.vertices){
    auto& vertexData = getVertexFunctionMemory(vertex, id)->data[maxLevel];
    for(Edge* edge : vertex.edges){
      auto& edgeData = getEdgeFunctionMemory(*edge, id)->data[maxLevel];
      uint_t vOnEdgeID = edge->vertex_index(vertex);
      //check vertex center
      if(vOnEdgeID == 0) WALBERLA_CHECK_EQUAL(vertexData[0],edgeData[0]);
      if(vOnEdgeID == 1) WALBERLA_CHECK_EQUAL(vertexData[0],edgeData[vPerEdge-1]);
      //check cells
      for(uint_t i = 0; i < vertex.faces.size(); ++i){
        uint_t fOnEdgeID = edge->face_index(*vertex.faces[i]);
        if(fOnEdgeID == 0 && vOnEdgeID == 0) {
          WALBERLA_CHECK_EQUAL(vertexData[vertex.edges.size() + 1 + i],
                               edgeData[edge_index(maxLevel, 0, CELL_GRAY_SE)])
        } else if(fOnEdgeID == 0 && vOnEdgeID == 1){
          WALBERLA_CHECK_EQUAL(vertexData[vertex.edges.size() + 1 + i],
                               edgeData[edge_index(maxLevel, vPerEdge - 1, CELL_GRAY_SW)])
        } else if(fOnEdgeID == 1 && vOnEdgeID == 0){
          WALBERLA_CHECK_EQUAL(vertexData[vertex.edges.size() + 1 + i],
                               edgeData[edge_index(maxLevel, 0, CELL_GRAY_NE)])
        } else if(fOnEdgeID == 1 && vOnEdgeID == 1){
          WALBERLA_CHECK_EQUAL(vertexData[vertex.edges.size() + 1 + i],
                               edgeData[edge_index(maxLevel, vPerEdge - 1, CELL_GRAY_NW)])
        }
      }
    }
  }

  //check edge to vertex comm; edge vertex Dofs are communicated
  for(Edge& edge: mesh.edges){
    auto& edgeData = getEdgeFunctionMemory(edge, id)->data[maxLevel];
    WALBERLA_CHECK_EQUAL(edgeData[1],
                         getVertexFunctionMemory(*edge.v0,id)->data[maxLevel][edge.v0->edge_index(edge) + 1]);
    WALBERLA_CHECK_EQUAL(edgeData[vPerEdge-2],
                         getVertexFunctionMemory(*edge.v1,id)->data[maxLevel][edge.v1->edge_index(edge) + 1]);
  }

  //check edge to face comm; edge vertex Dofs are communicated
  for(Edge& edge: mesh.edges) {
    auto &edgeData = getEdgeFunctionMemory(edge, id)->data[maxLevel];
    for(Face* face: edge.faces){
      auto &faceData = getFaceFunctionMemory(*face, id)->data[maxLevel];
      uint_t idxCounter = 0;
      for(auto it = P1BubbleFace::indexIterator(face->edge_index(edge),
                                                face->edge_orientation[face->edge_index(edge)],
                                                P1BubbleFace::VERTEX,
                                                maxLevel); it != P1BubbleFace::indexIterator(); ++it){
        WALBERLA_CHECK_EQUAL(edgeData[idxCounter],faceData[*it]);
        idxCounter++;
      }
    }
  }

  //check face to edge comm; face inner vertex Dofs and Cells are communicated
  for(Face& face: mesh.faces) {
    using namespace P1BubbleEdge::EdgeCoordsVertex;
    auto &faceData = getFaceFunctionMemory(face, id)->data[maxLevel];
    for(Edge* edge: face.edges){
      auto &edgeData = getEdgeFunctionMemory(*edge, id)->data[maxLevel];
      uint_t idxCounter = 0;
      uint_t faceIdOnEdge = edge->face_index(face);
      for(auto it = P1BubbleFace::indexIterator(face.edge_index(*edge),
                                                face.edge_orientation[face.edge_index(*edge)],
                                                P1BubbleFace::VERTEX_INNER,
                                                maxLevel); it != P1BubbleFace::indexIterator(); ++it){
        if(faceIdOnEdge == 0) {
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, VERTEX_SE)], faceData[*it]);
        } else if(faceIdOnEdge == 1){
          WALBERLA_CHECK_EQUAL(edgeData[edge_index(maxLevel, idxCounter, VERTEX_N)], faceData[*it]);
        } else{
          //TODO add walberla fail here
        }
        idxCounter++;
      }
    }
  }


//
//  hhg::P1BubbleFace::printFunctionMemory<maxLevel>(mesh.faces[0],id);
//  hhg::P1BubbleEdge::printFunctionMemory(mesh.edges[4],id,maxLevel);
//  hhg::P1BubbleEdge::printFunctionMemory(mesh.edges[5],id,maxLevel);
//  hhg::P1BubbleVertex::printFunctionMemory(mesh.vertices[4],id,maxLevel);

//  hhg::P1BubbleEdge::pull_halos(mesh.edges[4],0,maxLevel);
//
//  hhg::P1BubbleEdge::pull_vertices(mesh.edges[4],0,maxLevel);
//
//  hhg::P1BubbleVertex::pull_halos(mesh.vertices[3],0,maxLevel);
////  hhg::P1BubbleEdge::packDataforVertex(mesh.edges[4],0,sb,maxLevel,mesh.vertices[3]);
////  walberla::mpi::RecvBuffer rb3(sb);
////  hhg::P1BubbleVertex::unpackEdgeData(maxLevel,mesh.vertices[3],0,rb3,mesh.edges[4]);
//
//  hhg::P1BubbleVertex::print(mesh.vertices[3],0,maxLevel);
//
//  auto& face0mem = hhg::P1Bubble::getFaceFunctionMemory(face0, 0)->data[maxLevel];
//  auto& face1mem = hhg::P1Bubble::getFaceFunctionMemory(face1, 0)->data[maxLevel];
//  auto& edge5mem = hhg::P1Bubble::getEdgeFunctionMemory(mesh.edges[4], 0)->data[maxLevel];
//  real_t sumFace0 = 0;
//  real_t sumFace1 = 0;
//  std::cout << "Face 0: " << std::endl;
//  for(size_t i = 0; i < v_perEdge; ++i){
//    for(size_t j = 0; j < v_perEdge - i; ++j) {
//      sumFace0 += face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)];
//      //std::cout << face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
//      fmt::print("{0:.2f}  ",face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)]);
//    }
//    std::cout << std::endl;
//  }
//  std::cout << "=======================================" << std::endl;
//  std::cout << "Face 1: " << std::endl;
//  for(size_t i = 0; i < v_perEdge; ++i){
//    for(size_t j = 0; j < v_perEdge - i; ++j) {
//      sumFace1 += face1mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)];
//      //std::cout << face1mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
//      fmt::print("{0:.2f}  ",face1mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)]);
//    }
//    std::cout << std::endl;
//  }
//  std::cout << "=======================================" << std::endl;
//  std::cout << "Edge 5: " << std::endl;
//  for(size_t i = 0; i < v_perEdge - 1; ++i) {
//    fmt::print("{0:.2f}  ",edge5mem[hhg::P1BubbleEdge::EdgeCoordsVertex::index<maxLevel>(i, hhg::P1BubbleEdge::EdgeCoordsVertex::VERTEX_N)]);
//  }
//  std::cout << std::endl;
//  for(size_t i = 0; i < v_perEdge; ++i) {
//    fmt::print("{0:.2f}  ",edge5mem[hhg::P1BubbleEdge::EdgeCoordsVertex::index<maxLevel>(i, hhg::P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)]);
//  }
//  std::cout << std::endl;
//  for(size_t i = 0; i < v_perEdge -1 ; ++i) {
//    fmt::print("{0:.2f}  ",edge5mem[hhg::P1BubbleEdge::EdgeCoordsVertex::index<maxLevel>(i, hhg::P1BubbleEdge::EdgeCoordsVertex::VERTEX_SE)]);
//  }
//  std::cout << std::endl;
//  std::cout << "=======================================" << std::endl;

  //hhg::VTKWriter({ &x }, maxLevel, "output", "P1BubbleCommTest");

  return 0;
}
