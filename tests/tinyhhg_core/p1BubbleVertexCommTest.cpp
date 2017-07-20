//
// Created by thoennes on 13.04.17.
//

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubblefaceindex.hpp"
#include "core/mpi/SendBuffer.h"

using namespace hhg::P1BubbleFace;
using walberla::real_t;

int main (int argc, char ** argv )
{

  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  hhg::Mesh mesh("../../data/meshes/quad_4el.msh");


  size_t minLevel = 2;
  const size_t maxLevel = 3;

  //size_t v_perFace =  hhg::levelinfo::num_microvertices_per_face(maxLevel);
  //size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
  //size_t nbr_v_perEdge = v_perEdge - 1;
  //size_t v_perVertex = hhg::levelinfo::num_microvertices_per_vertex(maxLevel);

  hhg::P1BubbleFunction x("x", mesh, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> five = [](const hhg::Point3D&) { return 5; };
  std::function<real_t(const hhg::Point3D&)> six = [](const hhg::Point3D&) { return 6; };
  std::function<real_t(const hhg::Point3D&)> seven = [](const hhg::Point3D&) { return 7; };
  std::function<real_t(const hhg::Point3D&)> eight = [](const hhg::Point3D&) { return 8; };
  std::function<real_t(const hhg::Point3D&)> nine = [](const hhg::Point3D&) { return 9; };

  hhg::P1BubbleFace::interpolate(mesh.faces[0],0,six,maxLevel);
  hhg::P1BubbleFace::interpolate(mesh.faces[1],0,seven,maxLevel);
  hhg::P1BubbleFace::interpolate(mesh.faces[2],0,eight,maxLevel);
  hhg::P1BubbleFace::interpolate(mesh.faces[3],0,nine,maxLevel);

  for(auto edge : mesh.edges){
    if(edge.faces.size() == 2){
      hhg::P1BubbleEdge::interpolate(edge,0,five,maxLevel);
      hhg::P1BubbleEdge::pull_halos(edge,0,maxLevel);
    }
  }
  
  for(auto vertex : mesh.vertices){
    if(vertex.edges.size() == 4){
      std::cout << vertex << std::endl;
    }
  }

  hhg::P1BubbleVertex::pull_halos(mesh.vertices[4],0,maxLevel);
  
  auto& middleVertexMem = hhg::P1Bubble::getVertexFunctionMemory(mesh.vertices[4], 0)->data[maxLevel];

  std::cout << " Vertex: " << middleVertexMem[0] << std::endl;
  for(uint_t i = 1; i < 5; ++i){
    std::cout << " Edge " << i-1 << ": " << middleVertexMem[i] << std::endl;
  }
  for(uint_t i = 1; i < 5; ++i){
    std::cout << " Cell " << i-1 << ": " << middleVertexMem[i+4] << std::endl;
  }

//   int counter = 1;
//   for(auto edge : mesh.edges){
//     std::function<real_t(const hhg::Point3D&)> exact = [counter](const hhg::Point3D&) { return counter; };
//     hhg::P1BubbleEdge::interpolate(edge,0,exact,maxLevel);
//     counter++;
//   }

//   counter = 1;
//   for(auto vertex : mesh.vertices){
//     std::function<real_t(const hhg::Point3D&)> exact = [counter](const hhg::Point3D&) { return counter; };
//     hhg::P1BubbleVertex::interpolate(vertex,0,exact,maxLevel);
//     counter++;
//   }

//   walberla::mpi::SendBuffer sb;
//   auto& face0 = mesh.faces[0];
//   auto& face1 = mesh.faces[1];

//   for(hhg::Face& face : mesh.faces)
//   {
//     pull_edges(face,0,maxLevel);
//   }

//   hhg::P1BubbleEdge::pull_halos(mesh.edges[4],0,maxLevel);

//   hhg::P1BubbleEdge::pull_vertices(mesh.edges[4],0,maxLevel);

//   hhg::P1BubbleVertex::pull_halos(mesh.vertices[3],0,maxLevel);
// //  hhg::P1BubbleEdge::packDataforVertex(mesh.edges[4],0,sb,maxLevel,mesh.vertices[3]);
// //  walberla::mpi::RecvBuffer rb3(sb);
// //  hhg::P1BubbleVertex::unpackEdgeData(maxLevel,mesh.vertices[3],0,rb3,mesh.edges[4]);

//   hhg::P1BubbleVertex::print(mesh.vertices[3],0,maxLevel);

//   auto& face0mem = hhg::P1Bubble::getFaceFunctionMemory(face0, 0)->data[maxLevel];
//   auto& face1mem = hhg::P1Bubble::getFaceFunctionMemory(face1, 0)->data[maxLevel];
//   auto& edge5mem = hhg::P1Bubble::getEdgeFunctionMemory(mesh.edges[4], 0)->data[maxLevel];
//   real_t sumFace0 = 0;
//   real_t sumFace1 = 0;
//   std::cout << "Face 0: " << std::endl;
//   for(size_t i = 0; i < v_perEdge; ++i){
//     for(size_t j = 0; j < v_perEdge - i; ++j) {
//       sumFace0 += face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)];
//       //std::cout << face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
//       fmt::print("{0:.2f}  ",face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)]);
//     }
//     std::cout << std::endl;
//   }
//   std::cout << "=======================================" << std::endl;
//   std::cout << "Face 1: " << std::endl;
//   for(size_t i = 0; i < v_perEdge; ++i){
//     for(size_t j = 0; j < v_perEdge - i; ++j) {
//       sumFace1 += face1mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)];
//       //std::cout << face1mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
//       fmt::print("{0:.2f}  ",face1mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)]);
//     }
//     std::cout << std::endl;
//   }
//   std::cout << "=======================================" << std::endl;
//   std::cout << "Edge 5: " << std::endl;
//   for(size_t i = 0; i < v_perEdge - 1; ++i) {
//     fmt::print("{0:.2f}  ",edge5mem[hhg::P1BubbleEdge::EdgeCoordsVertex::index<maxLevel>(i, hhg::P1BubbleEdge::EdgeCoordsVertex::VERTEX_N)]);
//   }
//   std::cout << std::endl;
//   for(size_t i = 0; i < v_perEdge; ++i) {
//     fmt::print("{0:.2f}  ",edge5mem[hhg::P1BubbleEdge::EdgeCoordsVertex::index<maxLevel>(i, hhg::P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)]);
//   }
//   std::cout << std::endl;
//   for(size_t i = 0; i < v_perEdge -1 ; ++i) {
//     fmt::print("{0:.2f}  ",edge5mem[hhg::P1BubbleEdge::EdgeCoordsVertex::index<maxLevel>(i, hhg::P1BubbleEdge::EdgeCoordsVertex::VERTEX_SE)]);
//   }
//   std::cout << std::endl;
//   std::cout << "=======================================" << std::endl;

//   //hhg::VTKWriter({ &x }, maxLevel, "output", "P1BubbleCommTest");

  return 0;
}
