//
// Created by thoennes on 13.04.17.
//

#include "tinyhhg_core/tinyhhg.hpp"
#include <core/debug/TestSubsystem.h>
#include <core/mpi/MPIManager.h>
#include "core/DataTypes.h"

using walberla::real_t;

int main (int argc, char ** argv )
{

  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  hhg::Mesh mesh("../../data/meshes/tri_2el.msh");


  size_t minLevel = 2;
  size_t maxLevel = 5;

  size_t v_perFace =  hhg::levelinfo::num_microvertices_per_face(maxLevel);
  size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
  size_t nbr_v_perEdge = v_perEdge - 1;
  size_t v_perVertex = hhg::levelinfo::num_microvertices_per_vertex(maxLevel);

  hhg::P1BubbleFunction x("x", mesh, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D&) { return 13.0; };

  for(auto face : mesh.faces){
    hhg::P1BubbleFace::interpolate(face,0,exact,maxLevel);
  }



  return 0;
}
