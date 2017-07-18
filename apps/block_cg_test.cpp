#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG CG Test\n");

  hhg::Mesh mesh("../data/meshes/quad_4el.msh");

  size_t minLevel = 2;
  size_t maxLevel = 2;
  size_t maxiter = 10000;

  hhg::P1StokesFunction r("r", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction f("f", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction u("u", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction u_exact("u_exact", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction err("err", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld npoints_helper("npoints_helper", mesh, minLevel, maxLevel);

  hhg::P1BlockLaplaceOperator L(mesh, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return xx[0]*xx[0] - xx[1]*xx[1]; };
  std::function<real_t(const hhg::Point3D&)> rhs   = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u.v.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u.p.interpolate(exact, maxLevel, hhg::DirichletBoundary);

  u_exact.u.interpolate(exact, maxLevel);
  u_exact.v.interpolate(exact, maxLevel);
  u_exact.p.interpolate(exact, maxLevel);

  auto solver = hhg::CGSolver<hhg::P1StokesFunction, hhg::P1BlockLaplaceOperator>(mesh, minLevel, maxLevel);
  walberla::WcTimer timer;
  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, true);
  timer.end();
  fmt::printf("time was: %e\n",timer.last());
  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);

  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

//  hhg::VTKWriter({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "test");
  return 0;
}
