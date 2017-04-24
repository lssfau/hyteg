//
// Created by thoennes on 13.04.17.
//

#include <tinyhhg_core/tinyhhg.hpp>
#include <core/debug/TestSubsystem.h>
#include <core/mpi/MPIManager.h>

int main (int argc, char ** argv )
{

  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  hhg::Mesh mesh("./tri_1el.msh");

  size_t minLevel = 2;
  size_t maxLevel = 6;

  hhg::P1Function x("x", mesh, minLevel, maxLevel);

  //std::function<double(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return x[0] + x[1]; };
  std::function<double(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return 13; };

  auto faceZero = mesh.faces[0];
  hhg::P1Face::interpolate(faceZero,0,exact,maxLevel);

  for(size_t i = 0)


  //x.interpolate(exact, maxLevel);
  WALBERLA_LOG_DEVEL(mesh.faces[0].data[0][maxLevel][113]);

  return 0;
}
