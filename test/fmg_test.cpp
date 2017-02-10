#include <tinyhhg.hpp>

#include <fmt/format.h>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv); 
  int rk = hhg::Comm::get().rk;
  int np = hhg::Comm::get().np;

  if (rk == 0)
  {
    fmt::printf("[%d] TinyHHG FMG Test\n", rk);
  }

  hhg::Mesh mesh("../data/meshes/quad_4el.msh");

  size_t minLevel = 2;
  size_t maxLevel = 3;
  size_t nu_pre = 2;
  size_t nu_post = 2;
  size_t outer = 50;

  size_t coarse_maxiter = 100;
  double coarse_tolerance = 1e-6;
  double mg_tolerance = 1e-8;

  hhg::P1FunctionSpace V = hhg::P1FunctionSpace(mesh);
  hhg::Function<hhg::P1FunctionSpace> r("r", V, minLevel, maxLevel);
  hhg::Function<hhg::P1FunctionSpace> b("b", V, minLevel, maxLevel);
  hhg::Function<hhg::P1FunctionSpace> x("x", V, minLevel, maxLevel);
  hhg::Function<hhg::P1FunctionSpace> x_exact("x_exact", V, minLevel, maxLevel);
  hhg::Function<hhg::P1FunctionSpace> ax("ax", V, minLevel, maxLevel);
  hhg::Function<hhg::P1FunctionSpace> tmp("tmp", V, minLevel, maxLevel);
  hhg::Function<hhg::P1FunctionSpace> err("err", V, minLevel, maxLevel);

  hhg::P1LaplaceOperator A(V, minLevel, maxLevel);

  std::function<double(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return x[0]*x[0] - x[1]*x[1]; };
  std::function<double(const hhg::Point3D&)> rhs = [](const hhg::Point3D& x) { return 0.0; };
  std::function<double(const hhg::Point3D&)> zero = [](const hhg::Point3D& x) { return 0.0; };
  std::function<double(const hhg::Point3D&)> ones = [](const hhg::Point3D& x) { return 1.0; };

  x.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  x_exact.interpolate(exact, maxLevel);

  tmp.interpolate(ones, maxLevel);
  double npoints = tmp.dot(tmp, maxLevel);

  auto solver = hhg::CGSolver(V, minLevel, minLevel);

  if (rk == 0)
  {
    fmt::print("[{}] Num dofs = {}\n", rk, (size_t)npoints);
    fmt::printf("[%d] Starting V cycles\n", rk);
    fmt::printf("[%d] iter  abs_res       rel_res       conv          L2-error\n", rk);
  }

  double rel_res = 1.0;

  x.apply(A, ax, maxLevel, hhg::Inner);
  r.assign<2>({1.0, -1.0}, {&b, &ax}, maxLevel, hhg::Inner);

  double begin_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner));
  double abs_res_old = begin_res;

  err.assign<2>({1.0, -1.0}, {&x, &x_exact}, maxLevel);
  double discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  if (rk == 0)
  {
    fmt::printf("[%d] %-3d   %e  %e  %e  %e\n", rk, 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err);
  }

  std::function<void(size_t)> cscycle;

  cscycle = [&](size_t level)
  {
    if (level == minLevel)
    {
      // fmt::printf("Coarse solve...\n");
      solver.solve(A, x, b, r, minLevel, coarse_tolerance, coarse_maxiter, hhg::Inner, false);
    }
    else
    {
      // fmt::printf("Level %d...\n", level);

      // pre-smooth
      for (size_t i = 0; i < nu_pre; ++i)
      {
        x.smooth_gs(A, b, level, hhg::Inner);
      }

      x.apply(A, ax, level, hhg::Inner);
      r.assign<2>({1.0, -1.0}, { &b, &ax }, level, hhg::Inner);

      // restrict
      r.restrict(level, hhg::Inner);
      b.assign<1>({1.0}, { &r }, level - 1, hhg::Inner);

      x.interpolate(zero, level-1);

      cscycle(level-1);

      // prolongate
      tmp.assign<1>({1.0}, { &x }, level, hhg::Inner);
      x.prolongate(level-1, hhg::Inner);
      x.add<1>({1.0}, { &tmp }, level, hhg::Inner);

      // post-smooth
      for (size_t i = 0; i < nu_post; ++i)
      {
        x.smooth_gs(A, b, level, hhg::Inner);
      }
    }
  };

  for (size_t i = 0; i < outer; ++i)
  {
    cscycle(maxLevel);
    x.apply(A, ax, maxLevel, hhg::Inner);
    r.assign<2>({1.0, -1.0}, { &b, &ax }, maxLevel, hhg::Inner);
    double abs_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner));
    rel_res = abs_res / begin_res;
    err.assign<2>({1.0, -1.0}, { &x, &x_exact }, maxLevel);
    discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

    if (rk == 0)
    {
      fmt::printf("[%d] %-3d   %e  %e  %e  %e\n", rk, i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err);
    }

    abs_res_old = abs_res;

    if (rel_res < mg_tolerance)
    {
      break;
    }
  }

  // hhg::VTKWriter({ &x }, maxLevel, "../output", "test");
  MPI_Finalize();

  return 0;
}