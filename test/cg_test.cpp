#include <tinyhhg.hpp>

#include <fmt/format.h>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int rk = hhg::Comm::get().rk;
  if(rk == 0)
  {
    fmt::printf("TinyHHG CG Test\n");
  }
  hhg::Mesh mesh("../data/meshes/bfs_126el.msh");

  size_t minLevel = 2;
  size_t maxLevel = 5;
  size_t maxiter = 10000;

  hhg::P1Function r("r", mesh, minLevel, maxLevel);
  hhg::P1Function f("f", mesh, minLevel, maxLevel);
  hhg::P1Function u("u", mesh, minLevel, maxLevel);
  hhg::P1Function u_exact("u_exact", mesh, minLevel, maxLevel);
  hhg::P1Function err("err", mesh, minLevel, maxLevel);
  hhg::P1Function npoints_helper("npoints_helper", mesh, minLevel, maxLevel);

  hhg::P1LaplaceOperator L(mesh, minLevel, maxLevel);

  std::function<double(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return x[0]*x[0] - x[1]*x[1]; };
  std::function<double(const hhg::Point3D&)> rhs = [](const hhg::Point3D& x) { return 0.0; };
  std::function<double(const hhg::Point3D&)> ones = [](const hhg::Point3D& x) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u_exact.interpolate(exact, maxLevel);

  auto solver = hhg::CGSolver<hhg::P1Function>(mesh, minLevel, maxLevel);
  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, true);

  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);
  double npoints = npoints_helper.dot(npoints_helper, maxLevel);

  double discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  if (rk == 0)
  {
    fmt::printf("discrete L2 error = %e\n", discr_l2_err);
  }
  hhg::VTKWriter({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "test");
  MPI_Finalize();
  return 0;
}