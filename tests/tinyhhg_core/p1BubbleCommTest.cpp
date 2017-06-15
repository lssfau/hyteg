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
  hhg::Mesh mesh("../../data/meshes/tri_2el.msh");


  size_t minLevel = 2;
  const size_t maxLevel = 3;

  //size_t v_perFace =  hhg::levelinfo::num_microvertices_per_face(maxLevel);
  size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
  //size_t nbr_v_perEdge = v_perEdge - 1;
  //size_t v_perVertex = hhg::levelinfo::num_microvertices_per_vertex(maxLevel);

  hhg::P1BubbleFunction x("x", mesh, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D& xx) { return 0; };

  for(auto face : mesh.faces){
    hhg::P1BubbleFace::interpolate(face,0,zero,maxLevel);
  }
  int counter = 1;
  for(auto edge : mesh.edges){
    std::function<real_t(const hhg::Point3D&)> exact = [counter](const hhg::Point3D& xx) { return counter; };
    hhg::P1BubbleEdge::interpolate(edge,0,exact,maxLevel);
    counter++;
  }

  walberla::mpi::SendBuffer sb;
  auto& face0 = mesh.faces[0];
  auto& face1 = mesh.faces[1];

  for(uint_t i = 0; i < mesh.edges.size(); ++i) {
    if(face0.edge_index(mesh.edges[i]) <= 2){
      hhg::P1BubbleEdge::packData(mesh.edges[i], 0, sb, maxLevel);
      walberla::mpi::RecvBuffer rb(sb);
      hhg::P1BubbleFace::unpackEdgeData(maxLevel, face0, 0, rb, mesh.edges[i]);
    }
  }
  for(uint_t i = 0; i < mesh.edges.size(); ++i) {
    if(face1.edge_index(mesh.edges[i]) <= 2){
      hhg::P1BubbleEdge::packData(mesh.edges[i], 0, sb, maxLevel);
      walberla::mpi::RecvBuffer rb(sb);
      hhg::P1BubbleFace::unpackEdgeData(maxLevel, face1, 0, rb, mesh.edges[i]);
    }
  }

  auto& face0mem = hhg::P1Bubble::getFaceFunctionMemory(face0, 0)->data[maxLevel];
  auto& face1mem = hhg::P1Bubble::getFaceFunctionMemory(face1, 0)->data[maxLevel];
  real_t sumFace0 = 0;
  real_t sumFace1 = 0;
  std::cout << "Face 0: " << std::endl;
  for(size_t i = 0; i < v_perEdge; ++i){
    for(size_t j = 0; j < v_perEdge - i; ++j) {
      sumFace0 += face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)];
      //std::cout << face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
      fmt::print("{0:.2f}  ",face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)]);
    }
    std::cout << std::endl;
  }
  std::cout << "=======================================" << std::endl;
  std::cout << "Face 1: " << std::endl;
  for(size_t i = 0; i < v_perEdge; ++i){
    for(size_t j = 0; j < v_perEdge - i; ++j) {
      sumFace1 += face1mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)];
      //std::cout << face1mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
      fmt::print("{0:.2f}  ",face1mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)]);
    }
    std::cout << std::endl;
  }
  std::cout << "=======================================" << std::endl;


  return 0;
}
