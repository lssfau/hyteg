#include <tinyhhg_core/tinyhhg.hpp>

#include <fmt/format.h>

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  
  hhg::Mesh mesh("../data/meshes/quad_4el.msh");

  size_t minLevel = 2;
  const size_t maxLevel = 5;
  size_t maxiter = 1000;

  hhg::P1Function r("r", mesh, minLevel, maxLevel);
  hhg::P1Function f("f", mesh, minLevel, maxLevel);
  hhg::P1Function u("u", mesh, minLevel, maxLevel);
  hhg::P1Function u_exact("u_exact", mesh, minLevel, maxLevel);
  hhg::P1Function err("err", mesh, minLevel, maxLevel);
  hhg::P1Function npoints_helper("npoints_helper", mesh, minLevel, maxLevel);

  hhg::P1LaplaceOperator L(mesh, minLevel, maxLevel);

  std::function<walberla::real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) -> walberla::real_t { return x[0]*x[0] - x[1]*x[1]; };
  std::function<walberla::real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<walberla::real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate<maxLevel>(exact, hhg::DirichletBoundary);
  u_exact.interpolate<maxLevel>(exact);

  auto solver = hhg::MinResSolver<hhg::P1Function, hhg::P1LaplaceOperator>(mesh, minLevel, maxLevel);
  solver.solve<maxLevel>(L, u, f, r, 1e-8, maxiter, hhg::Inner, true);

  err.assign<maxLevel>({1.0, -1.0}, {&u, &u_exact});

  npoints_helper.interpolate<maxLevel>(ones);
  walberla::real_t npoints = npoints_helper.dot<maxLevel>(npoints_helper);

  walberla::real_t discr_l2_err = std::sqrt(err.dot<maxLevel>(err) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("discrete L2 error = {:e}", discr_l2_err));

  hhg::VTKWriter({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "minres_test");
  return EXIT_SUCCESS;
}
