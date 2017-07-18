#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>
#include <tinyhhg_core/likwidwrapper.hpp>

//using namespace walberla;
using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;
#define PI 3.14159265359

int main(int argc, char* argv[])
{
  LIKWID_MARKER_INIT;

  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  LIKWID_MARKER_THREADINIT;
  uint_t rk = uint_c(walberla::MPIManager::instance()->rank());

  if (walberlaEnv.config() == nullptr) {
    WALBERLA_ABORT("No parameter file was given");
  }

  auto parameters = walberlaEnv.config()->getOneBlock("Parameters");
  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG FMG Test");

  hhg::Mesh mesh(parameters.getParameter<std::string>("mesh"));

  size_t minLevel = parameters.getParameter<size_t>("minlevel");
  size_t maxLevel = parameters.getParameter<size_t>("maxlevel");
  size_t nu_pre = parameters.getParameter<size_t>("nu_pre");
  size_t nu_post = parameters.getParameter<size_t>("nu_post");
  size_t outer = parameters.getParameter<size_t>("outer_iter");

  size_t coarse_maxiter = 100;
  real_t coarse_tolerance = 1e-6;

  hhg::P1FunctionOld r("r", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld b("b", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld x("x", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld x_exact("x_exact", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld ax("ax", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld tmp("tmp", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld err("err", mesh, minLevel, maxLevel);

  hhg::P1LaplaceOperator A(mesh, minLevel, maxLevel);
  hhg::P1MassOperator M(mesh, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return sin(PI*xx[0])*sin(PI*xx[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs   = [](const hhg::Point3D& xx) { return 2*PI*PI*sin(PI*xx[0])*sin(PI*xx[1]); };
  std::function<real_t(const hhg::Point3D&)> zero  = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  for(size_t ll = minLevel; ll <= maxLevel;++ll)
  {
      x.interpolate(zero, ll, hhg::Inner);
      x.interpolate(exact, ll, hhg::DirichletBoundary);
      x_exact.interpolate(exact, ll);
      tmp.interpolate(rhs, ll, hhg::Inner);
      M.apply(tmp, b, ll, hhg::Inner);
  }

  auto solver = hhg::CGSolver<hhg::P1FunctionOld, hhg::P1LaplaceOperator>(mesh, minLevel, minLevel);

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
        A.smooth_gs(x, b, level, hhg::Inner);
      }

      A.apply(x, ax, level, hhg::Inner);
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
        A.smooth_gs(x, b, level, hhg::Inner);
      }
    }
  };
//    hhg::VTKWriter({ &x }, 0, "../output", "out");

  LIKWID_MARKER_START("Compute");
  for (size_t ll = minLevel; ll <= maxLevel; ++ll)
  {
      tmp.interpolate(ones, ll);
      real_t npoints = tmp.dot(tmp, ll);
      if (rk == 0)
      {
          WALBERLA_LOG_INFO_ON_ROOT(fmt::format("Level = {}", (size_t)ll));
          WALBERLA_LOG_INFO_ON_ROOT(fmt::format("Num dofs = {}", (size_t)npoints));
          WALBERLA_LOG_INFO_ON_ROOT("Starting V cycles");
          WALBERLA_LOG_INFO_ON_ROOT("iter  abs_res       rel_res       conv          L2-error           H1-semi");
      }
      real_t rel_res = 1.0;

      A.apply(x, ax, ll, hhg::Inner);
      r.assign({1.0, -1.0}, { &b, &ax }, ll, hhg::Inner);
      real_t abs_res_old = std::sqrt(r.dot(r, ll, hhg::Inner));
      real_t begin_res = abs_res_old;
      err.assign({1.0, -1.0}, {&x, &x_exact}, ll);
      real_t discr_l2_err = std::sqrt(err.dot(err, ll) / npoints);
      A.apply(err, tmp, ll, hhg::Inner);
      real_t discr_h1_err = std::sqrt(err.dot(tmp, ll));


      if (rk == 0)
      {
          WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  {:e}", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err, discr_h1_err));
      }

      for (size_t i = 0; i < outer; ++i)
      {
          cscycle(ll);
          A.apply(x, ax, ll, hhg::Inner);
          r.assign({1.0, -1.0}, { &b, &ax }, ll, hhg::Inner);
          real_t abs_res = std::sqrt(r.dot(r, ll, hhg::Inner));
          rel_res = abs_res / begin_res;
          err.assign({1.0, -1.0}, { &x, &x_exact }, ll);
          discr_l2_err = std::sqrt(err.dot(err, ll) / npoints);
          A.apply(err, tmp, ll, hhg::Inner);
          discr_h1_err = std::sqrt(err.dot(tmp, ll));

          if (rk == 0)
          {
              WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  {:e}", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err, discr_h1_err));
          }
          abs_res_old = abs_res;
      }
      if(ll < maxLevel)
      {
          x.prolongateQuadratic(ll, hhg::Inner);
      }

  }
  LIKWID_MARKER_STOP("Compute");

  LIKWID_MARKER_CLOSE;
  return EXIT_SUCCESS;
}
