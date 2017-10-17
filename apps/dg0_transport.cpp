#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>
#include <core/Environment.h>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/quad_2el.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  const uint_t minLevel = 2;
  const uint_t maxLevel = 3;
  const uint_t timesteps = 100;
  real_t dt = 0.25 * std::pow(2.0, -walberla::real_c(maxLevel+1));
  WALBERLA_LOG_DEVEL("dt = " << dt)

  std::function<real_t(const hhg::Point3D&)> initialConcentration = [](const hhg::Point3D& x) {
    if ((x - Point3D{{1.0/3.0, 2.0/3.0, 0.0}}).norm() < 0.1) {
      return 1.0;
    } else {
      return 0.0;
    }
//    return x[0]+1;
  };

  std::function<real_t(const hhg::Point3D&)> vel_x = [](const hhg::Point3D& x) {
    return 1.0;
  };

  std::function<real_t(const hhg::Point3D&)> vel_y = [](const hhg::Point3D& x) {
    return 0.0;
  };

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::DGFunction<real_t> c_old("c_old", storage, minLevel, maxLevel);
  hhg::DGFunction<real_t> c("c", storage, minLevel, maxLevel);
  std::shared_ptr<hhg::P1Function<real_t>> u = std::make_shared<hhg::P1Function<real_t>>("u", storage, minLevel, maxLevel);
  std::shared_ptr<hhg::P1Function<real_t>> v = std::make_shared<hhg::P1Function<real_t>>("v", storage, minLevel, maxLevel);

  std::array<std::shared_ptr<hhg::P1Function<real_t>>, 2> velocity{{u,v}};

  hhg::DGUpwindOperator<hhg::P1Function<real_t>> N(storage, velocity, minLevel, maxLevel);

  u->interpolate(vel_x, maxLevel);
  v->interpolate(vel_y, maxLevel);
  c_old.interpolate(initialConcentration, maxLevel);

  hhg::VTKWriter<hhg::P1Function< real_t >, hhg::DGFunction< real_t >, maxLevel>({ u.get(), v.get() }, { &c_old, &c }, "../output", fmt::format("dg0_transport-{:0>4}", 0));

  for(uint_t i = 1; i <= timesteps; i++) {
    N.apply(c_old, c, maxLevel, hhg::Inner, Replace);
    c.assign({1.0, -dt}, {&c_old, &c}, maxLevel, hhg::Inner);

//    if (i % 10 == 0) {
      hhg::VTKWriter<hhg::P1Function<real_t>, hhg::DGFunction<real_t>, maxLevel>({u.get(), v.get()},
                                                                                 {&c_old, &c},
                                                                                 "../output",
                                                                                 fmt::format("dg0_transport-{:0>4}",
                                                                                             i));
//    }

    c_old.assign({1.0}, {&c}, maxLevel, hhg::Inner);
  }

  return EXIT_SUCCESS;
}
