#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>
#include <tinyhhg_core/likwidwrapper.hpp>

#include <memory>

//using namespace walberla;
using walberla::uint_t;
using walberla::uint_c;

template<size_t Level, size_t MinLevel, class CSSolver, class O, class F>
struct CSCycle
{
  static void solve(CSSolver& solver, O& A, F& x, F& ax, F& b, F& r, F& tmp, walberla::real_t coarse_tolerance, size_t coarse_maxiter, size_t nu_pre, size_t nu_post)
  {
    std::function<double(const hhg::Point3D&)> zero  = [](const hhg::Point3D&) { return 0.0; };

    // fmt::printf("Level %d...\n", level);

    if (Level == MinLevel)
    {
      solver.template solve<MinLevel>(A, x, b, r, coarse_tolerance, coarse_maxiter, hhg::Inner, false);
    }
    else {
      // pre-smooth
      for (size_t i = 0; i < nu_pre; ++i) {
        A.template smooth_gs<Level>(x, b, hhg::Inner);
      }

      A.template apply<Level>(x, ax, hhg::Inner);
      r.template assign<Level>({1.0, -1.0}, {&b, &ax}, hhg::Inner);

      // restrict
      r.template restrict<Level>(hhg::Inner);
      b.template assign<Level - 1>({1.0}, {&r}, hhg::Inner);

      x.template interpolate<Level - 1>(zero);

      CSCycle<Level - 1, MinLevel, CSSolver, O, F>::solve(solver, A, x, ax, b, r, tmp, coarse_tolerance, coarse_maxiter, nu_pre,
                                                   nu_post);

      // prolongate
      tmp.template assign<Level>({1.0}, {&x}, hhg::Inner);
      x.template prolongate<Level - 1>(hhg::Inner);
      x.template add<Level>({1.0}, {&tmp}, hhg::Inner);

      // post-smooth
      for (size_t i = 0; i < nu_post; ++i) {
        A.template smooth_gs<Level>(x, b, hhg::Inner);
      }
    }
  }
};

template<size_t MinLevel, class CSSolver, class O, class F>
struct CSCycle<0, MinLevel, CSSolver, O, F>
{
  static void solve(CSSolver&, O&, F&, F&, F&, F&, F&, walberla::real_t, size_t, size_t, size_t)
  {
  };
};

int main(int argc, char* argv[])
{
  LIKWID_MARKER_INIT;

  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  LIKWID_MARKER_THREADINIT;

  walberla::shared_ptr<walberla::config::Config> cfg(new walberla::config::Config);
  if (walberlaEnv.config() == nullptr) {
    auto defaultFile = "../../data/param/fmg_test.prm";
    cfg->readParameterFile(defaultFile);
    if(!*cfg){
      WALBERLA_ABORT("could not open default file: " << defaultFile);
    }
  } else {
    cfg = walberlaEnv.config();
  }

  auto parameters = cfg->getOneBlock("Parameters");

  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG FMG Test");


  hhg::Mesh mesh(parameters.getParameter<std::string>("mesh"));

  const size_t minLevel = 2;
  const size_t maxLevel = 11;
  size_t nu_pre = parameters.getParameter<size_t>("nu_pre");
  size_t nu_post = parameters.getParameter<size_t>("nu_post");
  size_t outer = parameters.getParameter<size_t>("outer_iter");

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

  std::function<double(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return xx[0]*xx[0] - xx[1]*xx[1]; };
  std::function<double(const hhg::Point3D&)> rhs   = [](const hhg::Point3D&) { return 0.0; };
  std::function<double(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  x.interpolate<maxLevel>(exact, hhg::DirichletBoundary);
  x_exact.interpolate<maxLevel>(exact);

  tmp.interpolate<maxLevel>(ones);
  double npoints = tmp.dot<maxLevel>(tmp);

  auto solver = hhg::CGSolver<hhg::P1Function, hhg::P1LaplaceOperator>(mesh, minLevel, minLevel);

  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("Num dofs = {}", (size_t)npoints));
  WALBERLA_LOG_INFO_ON_ROOT("Starting V cycles");
  WALBERLA_LOG_INFO_ON_ROOT("iter  abs_res       rel_res       conv          L2-error");

  double rel_res = 1.0;

  A.apply<maxLevel>(x, ax, hhg::Inner);
  r.assign<maxLevel>({1.0, -1.0}, {&b, &ax}, hhg::Inner);

  double begin_res = std::sqrt(r.dot<maxLevel>(r, hhg::Inner));
  double abs_res_old = begin_res;

  err.assign<maxLevel>({1.0, -1.0}, {&x, &x_exact});
  double discr_l2_err = std::sqrt(err.dot<maxLevel>(err) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err));

  LIKWID_MARKER_START("Compute");
  size_t i = 0;
  for (; i < outer; ++i)
  {
    CSCycle<maxLevel, minLevel, decltype(solver), decltype(A), decltype(x)>::solve(solver, A, x, ax, b, r, tmp, coarse_tolerance, coarse_maxiter, nu_pre, nu_post);
    A.apply<maxLevel>(x, ax, hhg::Inner);
    r.assign<maxLevel>({1.0, -1.0}, { &b, &ax }, hhg::Inner);
    double abs_res = std::sqrt(r.dot<maxLevel>(r, hhg::Inner));
    rel_res = abs_res / begin_res;
    err.assign<maxLevel>({1.0, -1.0}, { &x, &x_exact });
    discr_l2_err = std::sqrt(err.dot<maxLevel>(err) / npoints);

    WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err));
    abs_res_old = abs_res;

    if (rel_res < mg_tolerance)
    {
      break;
    }
  }
  LIKWID_MARKER_STOP("Compute");

  WALBERLA_CHECK_LESS( i, outer );

  // hhg::VTKWriter({ &x }, maxLevel, "../output", "test");
  LIKWID_MARKER_CLOSE;
  return EXIT_SUCCESS;
}
