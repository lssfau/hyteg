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

  size_t v_perFace =  hhg::levelinfo::num_microvertices_per_face(maxLevel);
  size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
  //size_t nbr_v_perEdge = v_perEdge - 1;
  //size_t v_perVertex = hhg::levelinfo::num_microvertices_per_vertex(maxLevel);

  hhg::P1BubbleFunction x("x", mesh, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return xx[0]; };

  for(auto face : mesh.faces){
    hhg::P1BubbleFace::interpolate(face,0,exact,maxLevel);
  }
  for(auto edge : mesh.edges){
    hhg::P1BubbleEdge::interpolate(edge,0,exact,maxLevel);
  }

  //walberla::mpi::SendBuffer sb;

  for(uint_t i = 0; i < mesh.edges.size(); ++i) {
    walberla::mpi::SendBuffer sb;
    hhg::P1BubbleEdge::packData(mesh.edges[i], 0, sb, maxLevel);
    walberla::mpi::RecvBuffer rb(sb);
    hhg::P1BubbleFace::unpackData(maxLevel,mesh.faces[0], 0, rb, mesh.edges[i]);
  }
  for(uint_t i = 0; i < mesh.edges.size(); ++i) {
    walberla::mpi::SendBuffer sb;
    hhg::P1BubbleEdge::packData(mesh.edges[i], 0, sb, maxLevel);
    walberla::mpi::RecvBuffer rb(sb);
    hhg::P1BubbleFace::unpackData(maxLevel,mesh.faces[1], 0, rb, mesh.edges[i]);
  }


  real_t* face0 = getFaceP1BubbleFunctionMemory(mesh.faces[0], 0)->data[maxLevel];
  real_t* face1 = getFaceP1BubbleFunctionMemory(mesh.faces[1], 0)->data[maxLevel];
  real_t sumFace0 = 0;
  real_t sumFace1 = 0;
  std::cout << "Face 0: " << std::endl;
  for(size_t i = 0; i < v_perEdge; ++i){
    for(size_t j = 0; j < v_perEdge - i; ++j) {
      sumFace0 += face0[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)];
      //std::cout << face0[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
      fmt::print("{0:.2f}  ",face0[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)]);
    }
    std::cout << std::endl;
  }
  std::cout << "=======================================" << std::endl;
  std::cout << "Face 1: " << std::endl;
  for(size_t i = 0; i < v_perEdge; ++i){
    for(size_t j = 0; j < v_perEdge - i; ++j) {
      sumFace1 += face1[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)];
      //std::cout << face1[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
      fmt::print("{0:.2f}  ",face1[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)]);
    }
    std::cout << std::endl;
  }
  std::cout << "=======================================" << std::endl;
  WALBERLA_LOG_INFO_ON_ROOT(v_perFace);
  WALBERLA_LOG_INFO_ON_ROOT(v_perEdge);
  WALBERLA_LOG_INFO_ON_ROOT(sumFace0);
  WALBERLA_LOG_INFO_ON_ROOT(sumFace1);


  return 0;
}
