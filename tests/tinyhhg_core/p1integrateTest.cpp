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
  hhg::Mesh mesh("./tri_1el.msh");

  size_t minLevel = 2;
  const size_t maxLevel = 5;

  hhg::P1Function x("x", mesh, minLevel, maxLevel);

  std::function<double(const hhg::Point3D&)> exact = [](const hhg::Point3D&) { return 13.0; };

  auto faceZero = mesh.faces[0];
  hhg::P1Face::interpolate<maxLevel>(faceZero,0,exact);

  size_t totalPoints =  hhg::levelinfo::num_microvertices_per_face(maxLevel);
  size_t length = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
  for(size_t i = 0; i < totalPoints; ++i){
    if(hhg::P1Face::is_boundary(i,length)) {
      WALBERLA_CHECK_FLOAT_EQUAL(mesh.faces[0].data[0][maxLevel - 2][i],0.0)
    } else {
      WALBERLA_CHECK_FLOAT_EQUAL(mesh.faces[0].data[0][maxLevel - 2][i], 13.0)
    }
  }



//   auto edgeZero = mesh.edges[0];
//   hhg::P1Edge::interpolate(edgeZero,0,exact,maxLevel);
//
//   for(size_t i = 1; i < (micro_per_e-1); ++i){
//
//     WALBERLA_LOG_INFO()
//   }
//
//   for(size_t i = 1; i < (micro_per_e-1); ++i){
//     WALBERLA_CHECK_FLOAT_EQUAL(mesh.faces[0].data[0][maxLevel - 2][i], 13.0,
//                                "index was " << i)
//   }
//   WALBERLA_CHECK_FLOAT_EQUAL(mesh.faces[0].data[0][maxLevel - 2][0], 0.0)
//   WALBERLA_CHECK_FLOAT_EQUAL(mesh.faces[0].data[0][maxLevel - 2][micro_per_e-1], 0.0)


  return 0;
}
