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
  const uint_t maxLevel = 6;
  const uint_t solverMaxiter = 200;
  const uint_t timesteps = 50000;
  real_t dt = 0.75 * std::pow(2.0, -walberla::real_c(maxLevel));

  std::function<real_t(const hhg::Point3D&)> initialConcentration = [dt](const hhg::Point3D& x) {

    if (x[1] < 2.1/3.0 * dt/0.75) {
      return 100.0 * sin(M_PI*x[0]);
    } else {
      return 0.0;
    }

  };

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  // Setting up Functions
  auto c_old = std::make_shared<hhg::DGFunction<real_t>>("c", storage, minLevel, maxLevel);
  auto c = std::make_shared<hhg::DGFunction<real_t>>("c", storage, minLevel, maxLevel);
  auto r = std::make_shared<hhg::P1StokesFunction<real_t>>("r", storage, minLevel, maxLevel);
  auto f = std::make_shared<hhg::P1StokesFunction<real_t>>("f", storage, minLevel, maxLevel);
  auto u = std::make_shared<hhg::P1StokesFunction<real_t>>("u", storage, minLevel, maxLevel);

  // Setting up Operators


  std::array<std::shared_ptr<hhg::P1Function<real_t>>, 2> velocity{{std::shared_ptr<hhg::P1Function<real_t>>(&u->u), std::shared_ptr<hhg::P1Function<real_t>>(&u->v)}};
  hhg::DGUpwindOperator<hhg::P1Function<real_t>> N(storage, velocity, minLevel, maxLevel);
  hhg::P1StokesOperator L(storage, minLevel, maxLevel);

  // Interpolate initial functions
  c_old->interpolate(initialConcentration, maxLevel);
  c->assign({1.0}, {c_old.get()}, maxLevel);
  f->v.integrateDG(*c_old, maxLevel, hhg::All);

  auto solver = hhg::MinResSolver<hhg::P1StokesFunction<real_t>, hhg::P1StokesOperator>(storage, minLevel, maxLevel);

  for (uint_t t = 1; t <= timesteps; ++t) {

    if (t % 40 == 0) {
      f->v.integrateDG(*c_old, maxLevel, hhg::All);
      solver.solve(L, *u, *f, *r, maxLevel, 1e-3, solverMaxiter, hhg::Inner | hhg::NeumannBoundary, true);
      hhg::VTKWriter<hhg::P1Function< real_t >, hhg::DGFunction< real_t >, maxLevel>({ &u->u, &u->v, &u->p, &f->u, &f->v }, { c_old.get() }, "../output", fmt::format("plume-{:0>6}", t));
    }

    N.apply(*c_old, *c, maxLevel, hhg::Inner, Replace);
    c->assign({1.0, -dt}, {c_old.get(), c.get()}, maxLevel, hhg::Inner);

    c_old.swap(c);
  }

  return EXIT_SUCCESS;
}
