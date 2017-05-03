#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>
#include <tinyhhg_core/likwidwrapper.hpp>

int main(int argc, char* argv[])
{
  LIKWID_MARKER_INIT;
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  LIKWID_MARKER_THREADINIT;
  int rk = walberla::MPIManager::instance()->rank();
  int np = walberla::MPIManager::instance()->numProcesses();

  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG FMG Test");


  //hhg::Mesh mesh("C:/cygwin64/home/TUM/tinyhhg/build/apps/data/meshes/bfs_12el.msh");
  hhg::Mesh mesh("C:/cygwin64/home/TUM/tinyhhg/build/apps/data/meshes/quad_4el.msh");

  size_t minLevel = 2;
  size_t maxLevel = 8;
  size_t nu_pre = 2;
  size_t nu_post = 2;
  size_t outer = 50;

  size_t coarse_maxiter = 100;
  double coarse_tolerance = 1e-6;
  double mg_tolerance = 1e-8;

  hhg::P1Function r("r", mesh, minLevel, maxLevel);
  hhg::P1Function b("b", mesh, minLevel, maxLevel);
  hhg::P1Function x("x", mesh, minLevel, maxLevel);
  hhg::P1Function x_exact("x_exact", mesh, minLevel, maxLevel);
  hhg::P1Function ax("ax", mesh, minLevel, maxLevel);
  hhg::P1Function tmp("tmp", mesh, minLevel, maxLevel);
  hhg::P1Function err("err", mesh, minLevel, maxLevel);

  hhg::P1LaplaceOperator A(mesh, minLevel, maxLevel);
  
  WALBERLA_LOG_DEVEL("Creating functions")
  std::function<double(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return x[0]*x[0] - x[1]*x[1]; };
  std::function<double(const hhg::Point3D&)> rhs = [](const hhg::Point3D& x) { return 0.0; };
  std::function<double(const hhg::Point3D&)> zero = [](const hhg::Point3D& x) { return 0.0; };
  std::function<double(const hhg::Point3D&)> ones = [](const hhg::Point3D& x) { return 1.0; };

  WALBERLA_LOG_DEVEL("Interpolating x")
  x.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  WALBERLA_LOG_DEVEL("Interpolatin x_exact")
  x_exact.interpolate(exact, maxLevel);

  tmp.interpolate(ones, maxLevel);
  double npoints = tmp.dot(tmp, maxLevel);

  auto solver = hhg::CGSolver<hhg::P1Function>(mesh, minLevel, minLevel);

  if (rk == 0)
  {
    WALBERLA_LOG_INFO_ON_ROOT(fmt::format("Num dofs = {}", (size_t)npoints));
    WALBERLA_LOG_INFO_ON_ROOT("Starting V cycles");
    WALBERLA_LOG_INFO_ON_ROOT("iter  abs_res       rel_res       conv          L2-error");
  }

  double rel_res = 1.0;

  x.apply(A, ax, maxLevel, hhg::Inner);
  r.assign({1.0, -1.0}, {&b, &ax}, maxLevel, hhg::Inner);

  double begin_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner));
  double abs_res_old = begin_res;

  err.assign({1.0, -1.0}, {&x, &x_exact}, maxLevel);
  double discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  if (rk == 0)
  {
    WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err));
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
      r.assign({1.0, -1.0}, { &b, &ax }, level, hhg::Inner);

      // restrict
      r.restrict(level, hhg::Inner);
      b.assign({1.0}, { &r }, level - 1, hhg::Inner);

      x.interpolate(zero, level-1);

      cscycle(level-1);

      // prolongate
      tmp.assign({1.0}, { &x }, level, hhg::Inner);
      x.prolongate(level-1, hhg::Inner);
      x.add({1.0}, { &tmp }, level, hhg::Inner);

      // post-smooth
      for (size_t i = 0; i < nu_post; ++i)
      {
        x.smooth_gs(A, b, level, hhg::Inner);
      }
    }
  };

  LIKWID_MARKER_START("Compute");
  for (size_t i = 0; i < outer; ++i)
  {
    cscycle(maxLevel);
    x.apply(A, ax, maxLevel, hhg::Inner);
    r.assign({1.0, -1.0}, { &b, &ax }, maxLevel, hhg::Inner);
    double abs_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner));
    rel_res = abs_res / begin_res;
    err.assign({1.0, -1.0}, { &x, &x_exact }, maxLevel);
    discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

    if (rk == 0)
    {
      WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err));
    }

    abs_res_old = abs_res;

    if (rel_res < mg_tolerance)
    {
      break;
    }
  }
  LIKWID_MARKER_STOP("Compute");

  // hhg::VTKWriter({ &x }, maxLevel, "../output", "test");
  LIKWID_MARKER_CLOSE;
  return 0;
}
