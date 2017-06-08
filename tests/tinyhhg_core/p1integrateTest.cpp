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
  size_t maxLevel = 5;

  size_t v_perFace =  hhg::levelinfo::num_microvertices_per_face(maxLevel);
  size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
  size_t nbr_v_perEdge = v_perEdge - 1;
  size_t v_perVertex = hhg::levelinfo::num_microvertices_per_vertex(maxLevel);

  hhg::P1Function x("x", mesh, minLevel, maxLevel);

  for(auto face : mesh.faces){
    for(size_t i = 0; i < v_perFace; ++i){
      WALBERLA_CHECK_FLOAT_EQUAL(hhg::P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i],0.0)
    }
  }
  for(auto edge : mesh.edges){
    for(size_t i = 0; i < v_perEdge + edge.faces.size() * nbr_v_perEdge; ++i){
      WALBERLA_CHECK_FLOAT_EQUAL(hhg::P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i],0.0)
    }
  }
  for(auto vertex : mesh.vertices){
    //vertex have variable data sizes depending on the adjacent edges
    for(size_t i = 0; i < v_perVertex + vertex.edges.size();++i){
      WALBERLA_CHECK_FLOAT_EQUAL(hhg::P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i],0.0)
    }
  }


  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D&) { return 13.0; };

  auto faceZero = mesh.faces[0];
  hhg::P1Face::interpolate(faceZero,0,exact,maxLevel);


  for(size_t i = 0; i < v_perFace; ++i){
    if(hhg::P1Face::is_boundary(i,v_perEdge)) {
      WALBERLA_CHECK_FLOAT_EQUAL(hhg::P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i],0.0)
    } else {
      WALBERLA_CHECK_FLOAT_EQUAL(hhg::P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i], 13.0)
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
