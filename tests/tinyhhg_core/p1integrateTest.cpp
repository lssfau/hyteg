//
// Created by thoennes on 13.04.17.
//

#include "tinyhhg_core/tinyhhg.hpp"
#include <core/debug/CheckFunctions.h>
#include <core/debug/TestSubsystem.h>
#include <core/mpi/MPIManager.h>
#include "core/DataTypes.h"

using walberla::real_t;

namespace hhg {

static void testP1Integration()
{
  Mesh mesh("../../data/meshes/tri_1el.msh");

  MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  WALBERLA_LOG_INFO_ON_ROOT( setupStorage );
  PrimitiveStorage storage( uint_c( walberla::mpi::MPIManager::instance()->rank() ), setupStorage );

  size_t minLevel = 2;
  size_t maxLevel = 5;

  size_t v_perFace =  levelinfo::num_microvertices_per_face(maxLevel);
  size_t v_perEdge = levelinfo::num_microvertices_per_edge(maxLevel);
  size_t nbr_v_perEdge = v_perEdge - 1;
  size_t v_perVertex = levelinfo::num_microvertices_per_vertex(maxLevel);

  P1FunctionOld x("x", mesh, minLevel, maxLevel);

  for(auto face : mesh.faces){
    for(size_t i = 0; i < v_perFace; ++i){
      WALBERLA_CHECK_FLOAT_EQUAL( P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i], 0.0 );
    }
  }
  for(auto edge : mesh.edges){
    for(size_t i = 0; i < v_perEdge + edge.faces.size() * nbr_v_perEdge; ++i){
      WALBERLA_CHECK_FLOAT_EQUAL( P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i], 0.0 );
    }
  }
  for(auto vertex : mesh.vertices){
    //vertex have variable data sizes depending on the adjacent edges
    for(size_t i = 0; i < v_perVertex + vertex.edges.size();++i){
      WALBERLA_CHECK_FLOAT_EQUAL( P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i], 0.0 );
    }
  }


  std::function<real_t(const Point3D&)> exact = [](const Point3D&) { return 13.0; };

  auto faceZero = mesh.faces[0];
  P1Face::interpolate(faceZero,0,exact,maxLevel);


  for(size_t i = 0; i < v_perFace; ++i){
    if(P1Face::is_boundary(i,v_perEdge)) {
      WALBERLA_CHECK_FLOAT_EQUAL( P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i], 0.0);
    } else {
      WALBERLA_CHECK_FLOAT_EQUAL( P1::getFaceFunctionMemory(mesh.faces[0], x.memory_id)->data[maxLevel][i], 13.0);
    }
  }



  //   auto edgeZero = mesh.edges[0];
  //   P1Edge::interpolate(edgeZero,0,exact,maxLevel);
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


}
} // namespace hhg


int main ( int argc, char ** argv )
{
  walberla::debug::enterTestMode();
  walberla::mpi::Environment MPIenv( argc, argv );
  walberla::MPIManager::instance()->useWorldComm();

  hhg::testP1Integration();

  return EXIT_SUCCESS;

}


