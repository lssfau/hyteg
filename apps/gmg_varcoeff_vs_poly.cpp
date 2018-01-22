#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>
#include <core/Environment.h>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  walberla::shared_ptr<walberla::config::Config> cfg(new walberla::config::Config);
  cfg->readParameterFile("../data/param/gmg_varcoeff_vs_poly.prm");
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

  uint_t level_H = parameters.getParameter<uint_t>("level_h_coarse");
  uint_t level_h = parameters.getParameter<uint_t>("level_h_fine");
  const uint_t minLevel = 2;
  const uint_t maxLevel = level_h - level_H;
  const uint_t maxPolyDegree = 1;
  const uint_t interpolationLevel = maxLevel-1;
  const uint_t max_outer_iter =  parameters.getParameter<uint_t>("max_outer_iter");
  const uint_t max_cg_iter =  parameters.getParameter<uint_t>("max_cg_iter");
  const real_t mg_tolerance = parameters.getParameter<real_t>("mg_tolerance");
  const real_t coarse_tolerance = parameters.getParameter<real_t>("coarse_tolerance");
  const bool polynomialOperator = parameters.getParameter<bool>("polynomialOperator");

  MeshInfo meshInfo = MeshInfo::unitSquareMesh(level_H);
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  if (polynomialOperator) {
    WALBERLA_LOG_INFO_ON_ROOT("Polynomial Operator enabled");
  } else {
    WALBERLA_LOG_INFO_ON_ROOT("Polynomial Operator disabled");
  }

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  hhg::P1Function< real_t > r("r", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > f("f", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > u("u", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > Lu("Lu", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > u_exact("u_exact", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > err("err", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > npoints_helper("npoints_helper", storage, minLevel, maxLevel);
  std::shared_ptr<P1Function< real_t >> coefficient = std::make_shared<P1Function< real_t >>("coeff", storage, minLevel, maxLevel);

  hhg::P1MassOperator M(storage, minLevel, maxLevel);

  typedef hhg::P1VariableCoefficientLaplaceOperator SolveOperatorNodal;
  typedef hhg::P1PolynomialLaplaceOperator<maxPolyDegree> SolveOperatorPoly;
  typedef Operator< P1Function< real_t >, P1Function< real_t > > GeneralOperator;
  typedef std::shared_ptr<GeneralOperator> SolveOperator;

  std::function<real_t(const hhg::Point3D&)> coeff = [](const hhg::Point3D& x) { return 0.78546668391367*pow(x[0], 5) + 0.510366500099077*pow(x[0], 4)*x[1] + 0.391730092723877*pow(x[0], 4) + 0.207487088173805*pow(x[0], 3)*pow(x[1], 2) + 0.0805663025864225*pow(x[0], 3)*x[1] + 0.536962603494453*pow(x[0], 3) + 0.531408246732952*pow(x[0], 2)*pow(x[1], 3) + 0.169273512762915*pow(x[0], 2)*pow(x[1], 2) + 0.650313969995998*pow(x[0], 2)*x[1] + 0.936768966999531*pow(x[0], 2) + 0.33272576173186*x[0]*pow(x[1], 4) + 0.467034357005117*x[0]*pow(x[1], 3) + 0.710869731775775*x[0]*pow(x[1], 2) + 0.0716850897166866*x[0]*x[1] + 0.952980913121537*x[0] + 0.80834990475715*pow(x[1], 5) + 0.736229430842718*pow(x[1], 4) + 0.908328209520154*pow(x[1], 3) + 0.864072771153344*pow(x[1], 2) + 0.475144899525835*x[1] + 1; };
  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return pow(x[0], 3)*pow(x[1], 2)/(x[0]*x[1] + 1); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D& x) { return x[0]*(-x[0]*x[1]*(x[0]*x[1] + 1)*(x[0]*(x[0]*x[1] + 2)*(0.510366500099077*pow(x[0], 4) + 0.41497417634761*pow(x[0], 3)*x[1] + 0.0805663025864225*pow(x[0], 3) + 1.59422474019885*pow(x[0], 2)*pow(x[1], 2) + 0.33854702552583*pow(x[0], 2)*x[1] + 0.650313969995998*pow(x[0], 2) + 1.33090304692744*x[0]*pow(x[1], 3) + 1.40110307101535*x[0]*pow(x[1], 2) + 1.42173946355155*x[0]*x[1] + 0.0716850897166866*x[0] + 4.04174952378575*pow(x[1], 4) + 2.94491772337087*pow(x[1], 3) + 2.72498462856046*pow(x[1], 2) + 1.72814554230669*x[1] + 0.475144899525835) + x[1]*(2*x[0]*x[1] + 3)*(3.92733341956835*pow(x[0], 4) + 2.04146600039631*pow(x[0], 3)*x[1] + 1.56692037089551*pow(x[0], 3) + 0.622461264521414*pow(x[0], 2)*pow(x[1], 2) + 0.241698907759268*pow(x[0], 2)*x[1] + 1.61088781048336*pow(x[0], 2) + 1.0628164934659*x[0]*pow(x[1], 3) + 0.33854702552583*x[0]*pow(x[1], 2) + 1.300627939992*x[0]*x[1] + 1.87353793399906*x[0] + 0.33272576173186*pow(x[1], 4) + 0.467034357005117*pow(x[1], 3) + 0.710869731775775*pow(x[1], 2) + 0.0716850897166866*x[1] + 0.952980913121537)) + 2*(pow(x[0], 2)*(-pow(x[0], 2)*pow(x[1], 2) + 2*x[0]*x[1]*(x[0]*x[1] + 1) - pow(x[0]*x[1] + 1, 2)) + pow(x[1], 2)*(-pow(x[0], 2)*pow(x[1], 2) + 3*x[0]*x[1]*(x[0]*x[1] + 1) - 3*pow(x[0]*x[1] + 1, 2)))*(0.78546668391367*pow(x[0], 5) + 0.510366500099077*pow(x[0], 4)*x[1] + 0.391730092723877*pow(x[0], 4) + 0.207487088173805*pow(x[0], 3)*pow(x[1], 2) + 0.0805663025864225*pow(x[0], 3)*x[1] + 0.536962603494453*pow(x[0], 3) + 0.531408246732952*pow(x[0], 2)*pow(x[1], 3) + 0.169273512762915*pow(x[0], 2)*pow(x[1], 2) + 0.650313969995998*pow(x[0], 2)*x[1] + 0.936768966999531*pow(x[0], 2) + 0.33272576173186*x[0]*pow(x[1], 4) + 0.467034357005117*x[0]*pow(x[1], 3) + 0.710869731775775*x[0]*pow(x[1], 2) + 0.0716850897166866*x[0]*x[1] + 0.952980913121537*x[0] + 0.80834990475715*pow(x[1], 5) + 0.736229430842718*pow(x[1], 4) + 0.908328209520154*pow(x[1], 3) + 0.864072771153344*pow(x[1], 2) + 0.475144899525835*x[1] + 1))/pow(x[0]*x[1] + 1, 3); };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);

  for (uint_t level = minLevel; level <= maxLevel; ++level) {
    coefficient->interpolate(coeff, level);
  }

  u_exact.interpolate(exact, maxLevel);
  npoints_helper.interpolate(rhs, maxLevel);
  M.apply(npoints_helper, f, maxLevel, hhg::All);

  SolveOperator L;

  auto start = walberla::timing::getWcTime();
  if (polynomialOperator) {
    L = std::make_shared<SolveOperatorPoly>(storage, coefficient, coeff, minLevel, maxLevel, interpolationLevel);
  } else {
    L = std::make_shared<SolveOperatorNodal>(storage, coefficient, coeff, minLevel, maxLevel);
  }
  auto end = walberla::timing::getWcTime();
  real_t setupTime = end - start;

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  typedef hhg::CGSolver<hhg::P1Function<real_t>, GeneralOperator> CoarseSolver;
  auto coarseLaplaceSolver = std::make_shared<CoarseSolver>(storage, minLevel, minLevel);
  typedef GMultigridSolver<hhg::P1Function<real_t>, GeneralOperator, CoarseSolver> LaplaceSover;
  LaplaceSover laplaceSolver(storage, coarseLaplaceSolver, minLevel, maxLevel);

  WALBERLA_LOG_INFO_ON_ROOT("Starting V cycles");
  WALBERLA_LOG_INFO_ON_ROOT("iter  abs_res       rel_res       conv          L2-error      Time");

  real_t rel_res = 1.0;

  L->apply(u, Lu, maxLevel, hhg::Inner);
  r.assign({1.0, -1.0}, {&f, &Lu}, maxLevel, hhg::Inner);

  real_t begin_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner));
  real_t abs_res_old = begin_res;

  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);
  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  -", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err));

  real_t totalTime = real_c(0.0);
  real_t averageConvergenceRate = real_c(0.0);
  const uint_t convergenceStartIter = 3;

  uint_t i = 0;
  for (; i < max_outer_iter; ++i)
  {
    start = walberla::timing::getWcTime();
    laplaceSolver.solve(*L, u, f, r, maxLevel, coarse_tolerance, max_cg_iter, hhg::Inner, LaplaceSover::CycleType::VCYCLE, false);
    end = walberla::timing::getWcTime();
    L->apply(u, Lu, maxLevel, hhg::Inner);
    r.assign({1.0, -1.0}, { &f, &Lu }, maxLevel, hhg::Inner);
    real_t abs_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner));
    rel_res = abs_res / begin_res;
    err.assign({1.0, -1.0}, { &u, &u_exact }, maxLevel);
    discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

    WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  {:e}", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err, end-start));
    totalTime += end-start;

    if (i >= convergenceStartIter) {
      averageConvergenceRate += abs_res/abs_res_old;
    }

    abs_res_old = abs_res;

    if (rel_res < mg_tolerance)
    {
      break;
    }
  }

  WALBERLA_LOG_INFO_ON_ROOT("L(H): " << level_H);
  WALBERLA_LOG_INFO_ON_ROOT("L(h): " << level_h);
  if (polynomialOperator) {
    WALBERLA_LOG_INFO_ON_ROOT("Polynomial degree: " << maxPolyDegree);
    WALBERLA_LOG_INFO_ON_ROOT("LSQP level: " << interpolationLevel);
  }
  WALBERLA_LOG_INFO_ON_ROOT("Setup time: " << std::defaultfloat << setupTime);
  WALBERLA_LOG_INFO_ON_ROOT("Time to solution: " << std::defaultfloat << totalTime);
  WALBERLA_LOG_INFO_ON_ROOT("Avg. convergence rate: " << std::scientific << averageConvergenceRate / real_c(i-convergenceStartIter));
  WALBERLA_LOG_INFO_ON_ROOT("L^2 error: " << std::scientific << discr_l2_err);
  WALBERLA_LOG_INFO_ON_ROOT("DoFs: " << (uint_t) npoints);

  return 0;
}
