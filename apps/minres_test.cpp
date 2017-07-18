#include <tinyhhg_core/tinyhhg.hpp>

#include <fmt/format.h>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  hhg::Mesh mesh("../data/meshes/quad_4el.msh");

  size_t minLevel = 2;
  size_t maxLevel = 5;
  size_t maxiter = 1000;

  hhg::P1FunctionOld r("r", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld f("f", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld u("u", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld u_exact("u_exact", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld err("err", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld npoints_helper("npoints_helper", mesh, minLevel, maxLevel);

  hhg::P1LaplaceOperator L(mesh, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) -> real_t { return x[0]*x[0] - x[1]*x[1]; };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u_exact.interpolate(exact, maxLevel);

  auto solver = hhg::MinResSolver<hhg::P1FunctionOld, hhg::P1LaplaceOperator>(mesh, minLevel, maxLevel);
  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, true);

  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("discrete L2 error = {:e}", discr_l2_err));

  hhg::VTKWriter({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "minres_test");
  return EXIT_SUCCESS;
}
