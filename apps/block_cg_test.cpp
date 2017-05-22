#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG CG Test\n");

  hhg::Mesh mesh("../data/meshes/quad_4el.msh");

  size_t minLevel = 2;
  const size_t maxLevel = 2;
  size_t maxiter = 10000;

  hhg::P1StokesFunction r("r", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction f("f", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction u("u", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction u_exact("u_exact", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction err("err", mesh, minLevel, maxLevel);
  hhg::P1Function npoints_helper("npoints_helper", mesh, minLevel, maxLevel);

  hhg::P1BlockLaplaceOperator L(mesh, minLevel, maxLevel);

  std::function<walberla::real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return xx[0]*xx[0] - xx[1]*xx[1]; };
  std::function<walberla::real_t(const hhg::Point3D&)> rhs   = [](const hhg::Point3D&) { return 0.0; };
  std::function<walberla::real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate<maxLevel>(exact, hhg::DirichletBoundary);
  u.v.interpolate<maxLevel>(exact, hhg::DirichletBoundary);
  u.p.interpolate<maxLevel>(exact, hhg::DirichletBoundary);

  u_exact.u.interpolate<maxLevel>(exact);
  u_exact.v.interpolate<maxLevel>(exact);
  u_exact.p.interpolate<maxLevel>(exact);

  auto solver = hhg::CGSolver<hhg::P1StokesFunction, hhg::P1BlockLaplaceOperator>(mesh, minLevel, maxLevel);
  walberla::WcTimer timer;
  solver.solve<maxLevel>(L, u, f, r, 1e-8, maxiter, hhg::Inner, true);
  timer.end();
  fmt::printf("time was: %e\n",timer.last());
  err.assign<maxLevel>({1.0, -1.0}, {&u, &u_exact});

  npoints_helper.interpolate<maxLevel>(ones);

  walberla::real_t npoints = npoints_helper.dot<maxLevel>(npoints_helper);

  walberla::real_t discr_l2_err = std::sqrt(err.dot<maxLevel>(err) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

//  hhg::VTKWriter({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "test");
  return 0;
}
