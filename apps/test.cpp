#include <tinyhhg_core/tinyhhg.hpp>

#include <fmt/format.h>

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG default test");

  hhg::Mesh mesh("../data/meshes/tri_1el.msh");

  size_t minLevel = 2;
  const size_t maxLevel = 5;

  std::function<walberla::real_t(const hhg::Point3D&)> expr = [](const hhg::Point3D&) { return 1.0; };
  hhg::P1Function u("u", mesh, minLevel, maxLevel);
  hhg::P1Function Lu("Lu", mesh, minLevel, maxLevel);

  u.interpolate<maxLevel>(expr);
  hhg::P1LaplaceOperator L(mesh, minLevel, maxLevel);
  hhg::P1MassOperator L_gen(mesh, minLevel, maxLevel);

  hhg::VTKWriter({&u}, maxLevel, "../output", "test");
  return 0;
}
