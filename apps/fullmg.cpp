#include <tinyhhg_core/tinyhhg.hpp>
#include <tinyhhg_core/likwidwrapper.hpp>
#include <core/Environment.h>

#define PI 3.14159265359

int main(int argc, char* argv[])
{
  LIKWID_MARKER_INIT;

  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  LIKWID_MARKER_THREADINIT;

  walberla::shared_ptr<walberla::config::Config> cfg(new walberla::config::Config);
  if (walberlaEnv.config() == nullptr) {
    auto defaultFile = "../data/param/fullMGTest.prm";
    cfg->readParameterFile(defaultFile);
    if(!*cfg){
      WALBERLA_ABORT("could not open default file: " << defaultFile);
    }
  } else {
    cfg = walberlaEnv.config();
  }

  auto parameters = cfg->getOneBlock("Parameters");
  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG FMG Test");

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile(parameters.getParameter<std::string>("mesh"));
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  size_t minLevel = parameters.getParameter<size_t>("minlevel");
  size_t maxLevel = parameters.getParameter<size_t>("maxlevel");
  size_t nu_pre = parameters.getParameter<size_t>("nu_pre");
  size_t nu_post = parameters.getParameter<size_t>("nu_post");
  size_t outer = parameters.getParameter<size_t>("outer_iter");

  size_t coarse_maxiter = 100;
  real_t coarse_tolerance = 1e-6;

  hhg::P1Function< real_t > r("r", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > b("b", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > x("x", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > x_exact("x_exact", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > ax("ax", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > tmp("tmp", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > err("err", storage, minLevel, maxLevel);

  hhg::P1LaplaceOperator A(storage, minLevel, maxLevel);
  hhg::P1MassOperator M(storage, minLevel, maxLevel);

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  r.enableTiming( timingTree );
  b.enableTiming( timingTree );
  x.enableTiming( timingTree );
  x_exact.enableTiming( timingTree );
  ax.enableTiming( timingTree );
  tmp.enableTiming( timingTree );
  err.enableTiming( timingTree );

  A.enableTiming( timingTree );
  M.enableTiming( timingTree );

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

  auto solver = hhg::CGSolver<hhg::P1Function< real_t >, hhg::P1LaplaceOperator>(storage, minLevel, minLevel);

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
  //hhg::VTKWriter< hhg::P1Function< real_t > >({ &x, &x_exact, &b }, minLevel, "../output", "fullmg");

  LIKWID_MARKER_START("Compute");
  for (size_t ll = minLevel; ll <= maxLevel; ++ll)
  {
      tmp.interpolate(ones, ll);
      real_t npoints = tmp.dot(tmp, ll);
      WALBERLA_LOG_INFO_ON_ROOT(fmt::format("Level = {}", (size_t)ll));
      WALBERLA_LOG_INFO_ON_ROOT(fmt::format("Num dofs = {}", (size_t)npoints));
      WALBERLA_LOG_INFO_ON_ROOT("Starting V cycles");
      WALBERLA_LOG_INFO_ON_ROOT("iter  abs_res       rel_res       conv          L2-error           H1-semi");
      real_t rel_res = 1.0;

      A.apply(x, ax, ll, hhg::Inner);
      r.assign({1.0, -1.0}, { &b, &ax }, ll, hhg::Inner);
      real_t abs_res_old = std::sqrt(r.dot(r, ll, hhg::Inner));
      real_t begin_res = abs_res_old;
      err.assign({1.0, -1.0}, {&x, &x_exact}, ll);
      real_t discr_l2_err = std::sqrt(err.dot(err, ll) / npoints);
      A.apply(err, tmp, ll, hhg::Inner);
      real_t discr_h1_err = std::sqrt(err.dot(tmp, ll));

      WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  {:e}", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err, discr_h1_err));

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

          WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  {:e}", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err, discr_h1_err));
          abs_res_old = abs_res;
      }
      if(ll < maxLevel)
      {
          x.prolongateQuadratic(ll, hhg::Inner);
      }

  }
  LIKWID_MARKER_STOP("Compute");

  walberla::WcTimingTree tt = timingTree->getReduced();
  WALBERLA_LOG_INFO_ON_ROOT( tt );

  LIKWID_MARKER_CLOSE;
  return EXIT_SUCCESS;
}
