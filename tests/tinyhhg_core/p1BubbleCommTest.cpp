//
// Created by thoennes on 13.04.17.
//

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"
#include <core/debug/TestSubsystem.h>
#include <core/mpi/MPIManager.h>
#include "core/DataTypes.h"

using namespace hhg::P1BubbleFace;
using walberla::real_t;

int main (int argc, char ** argv )
{

  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  hhg::Mesh mesh("../../data/meshes/tri_2el.msh");


  size_t minLevel = 2;
  const size_t maxLevel = 5;

  size_t v_perFace =  hhg::levelinfo::num_microvertices_per_face(maxLevel);
  size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
  size_t nbr_v_perEdge = v_perEdge - 1;
  size_t v_perVertex = hhg::levelinfo::num_microvertices_per_vertex(maxLevel);

  hhg::P1BubbleFunction x("x", mesh, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D&) { return 1.0; };

  for(auto face : mesh.faces){
    hhg::P1BubbleFace::interpolate(face,0,exact,maxLevel);
  }
  real_t* face0 = getFaceP1BubbleFunctionMemory(mesh.faces[0], 0)->data[maxLevel];
  real_t* face1 = getFaceP1BubbleFunctionMemory(mesh.faces[0], 0)->data[maxLevel];
  real_t sumFace0 = 0;
  real_t sumFace1 = 0;
  for(size_t i = 0; i < v_perEdge; ++i){
    for(size_t j = 0; j < v_perEdge - i; ++j) {
      sumFace0 += face0[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)];
    }
  }
  WALBERLA_LOG_INFO_ON_ROOT(v_perFace);
  WALBERLA_LOG_INFO_ON_ROOT(v_perEdge);
  WALBERLA_LOG_INFO_ON_ROOT(sumFace0);

  return 0;
}
