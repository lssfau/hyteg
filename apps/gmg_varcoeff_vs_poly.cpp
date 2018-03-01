#define FP_FAST_FMA
#define FP_FAST_FMAF
#define FP_FAST_FMAL

#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <core/Environment.h>
#include <core/config/Create.h>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  walberla::shared_ptr<walberla::config::Config> cfg;

  if (argc == 1) {
    walberla::shared_ptr<walberla::config::Config> cfg_(new walberla::config::Config);
    cfg_->readParameterFile("../data/param/gmg_varcoeff_vs_poly.prm");
    cfg = cfg_;
  } else {
    cfg = walberla::config::create(argc, argv);
  }
  WALBERLA_LOG_INFO("config = " << *cfg);
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

  uint_t level_H = parameters.getParameter<uint_t>("level_h_coarse");
  uint_t level_h = parameters.getParameter<uint_t>("level_h_fine");
  const uint_t minLevel = 2;
  const uint_t maxLevel = level_h - level_H;
  const uint_t maxPolyDegree = parameters.getParameter<uint_t>("maxPolyDegree");
  const uint_t interpolationLevel = parameters.getParameter<uint_t>("interpolationLevel");
  const uint_t max_outer_iter =  parameters.getParameter<uint_t>("max_outer_iter");
  const uint_t max_cg_iter =  parameters.getParameter<uint_t>("max_cg_iter");
  const real_t mg_tolerance = parameters.getParameter<real_t>("mg_tolerance");
  const real_t coarse_tolerance = parameters.getParameter<real_t>("coarse_tolerance");
  const bool polynomialOperator = parameters.getParameter<bool>("polynomialOperator");

  uint_t numMacroEdgesPerRectangleEdge = uint_c(std::pow(2, level_H));
  MeshInfo meshInfo = MeshInfo::meshRectangle({{0.0, 0.0}}, {{1.0, 1.0}}, MeshInfo::meshFlavour::CRISS,
                                              numMacroEdgesPerRectangleEdge, numMacroEdgesPerRectangleEdge);

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
  std::shared_ptr<P1Function< real_t >> coefficient_11 = std::make_shared<P1Function< real_t >>("coeff_11", storage, minLevel, maxLevel);
  std::shared_ptr<P1Function< real_t >> coefficient_12 = std::make_shared<P1Function< real_t >>("coeff_12", storage, minLevel, maxLevel);
  std::shared_ptr<P1Function< real_t >> coefficient_22 = std::make_shared<P1Function< real_t >>("coeff_22", storage, minLevel, maxLevel);
  std::vector<std::shared_ptr<P1Function< real_t >>> coefficients;
  coefficients.push_back(coefficient_11);
  coefficients.push_back(coefficient_12);
  coefficients.push_back(coefficient_22);
  auto coordX = std::make_shared<hhg::P1Function<real_t>>("x", storage, minLevel, maxLevel);
  auto coordY = std::make_shared<hhg::P1Function<real_t>>("y", storage, minLevel, maxLevel);

  hhg::P1MassOperator M(storage, minLevel, maxLevel);

  typedef hhg::P1TensorCoefficientLaplaceOperator SolveOperatorNodal;
  typedef hhg::P1PolynomialTensorCoefficientLaplaceOperator SolveOperatorPoly;
  typedef Operator< P1Function< real_t >, P1Function< real_t > > GeneralOperator;
  typedef std::shared_ptr<GeneralOperator> SolveOperator;

  std::function<real_t(const hhg::Point3D&)> map_x = [](const hhg::Point3D& x) { return (1.0*x[1] + 1.0)*cos(1.5707963267949*x[0] - 2.35619449019234); };
  std::function<real_t(const hhg::Point3D&)> map_y = [](const hhg::Point3D& x) { return -(1.0*x[1] + 1.0)*sin(1.5707963267949*x[0] - 2.35619449019234); };
  std::function<real_t(const hhg::Point3D&)> coeff_11 = [](const hhg::Point3D& x) { return 0.636619772367581/(pow(x[1] + 1, 2)*sqrt(pow(x[1] + 1, -2))); };
  std::function<real_t(const hhg::Point3D&)> coeff_12 = [](const hhg::Point3D& x) { return 0; };
  std::function<real_t(const hhg::Point3D&)> coeff_22 = [](const hhg::Point3D& x) { return 1.5707963267949/sqrt(pow(x[1] + 1, -2)); };
  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D& x) { return -(1.5707963267949*pow(x[1] + 1, 2)*sinh(x[1]) + 1.5707963267949*(x[1] + 1)*cosh(x[1]) - 0.636619772367581*sinh(x[1]))*sin(x[0])/(pow(x[1] + 1, 4)*pow(pow(x[1] + 1, -2), 3.0L/2.0L)); };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  std::vector<std::function<real_t(const hhg::Point3D&)>> analyticCoefficients;
  analyticCoefficients.push_back(coeff_11);
  analyticCoefficients.push_back(coeff_12);
  analyticCoefficients.push_back(coeff_22);

  coordX->interpolate(map_x, maxLevel);
  coordY->interpolate(map_y, maxLevel);

  WALBERLA_LOG_INFO_ON_ROOT("Interpolating u");
  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);

  for (uint_t level = minLevel; level <= maxLevel; ++level) {
    coefficient_11->interpolate(coeff_11, level);
    coefficient_12->interpolate(coeff_12, level);
    coefficient_22->interpolate(coeff_22, level);
  }

  WALBERLA_LOG_INFO_ON_ROOT("Interpolating exact function");
  u_exact.interpolate(exact, maxLevel);
  WALBERLA_LOG_INFO_ON_ROOT("Interpolating and integrating rhs");
  npoints_helper.interpolate(rhs, maxLevel);
  M.apply(npoints_helper, f, maxLevel, hhg::All);

  WALBERLA_LOG_INFO_ON_ROOT("Setting up stiffness operator");
  SolveOperator L;

  auto start = walberla::timing::getWcTime();
  if (polynomialOperator) {
    L = std::make_shared<SolveOperatorPoly>(storage, coefficients, analyticCoefficients, minLevel, maxLevel, maxPolyDegree, interpolationLevel);
  } else {
    L = std::make_shared<SolveOperatorNodal>(storage, coefficients, minLevel, maxLevel);
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
  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6s|%10s|%10s|%10s|%10s|%10s","iter","abs_res","rel_res","conv","L2-error","Time"));

  real_t rel_res = 1.0;

  L->apply(u, Lu, maxLevel, hhg::Inner);
  r.assign({1.0, -1.0}, {&f, &Lu}, maxLevel, hhg::Inner);

  real_t begin_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner));
  real_t abs_res_old = begin_res;

  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);
  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err,0));

  real_t solveTime = real_c(0.0);
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

    WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err,end - start));
    solveTime += end-start;

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
  WALBERLA_LOG_INFO_ON_ROOT("Solve time " << std::defaultfloat << solveTime);
  WALBERLA_LOG_INFO_ON_ROOT("Time to solution: " << std::defaultfloat << setupTime + solveTime);
  WALBERLA_LOG_INFO_ON_ROOT("Avg. convergence rate: " << std::scientific << averageConvergenceRate / real_c(i-convergenceStartIter));
  WALBERLA_LOG_INFO_ON_ROOT("L^2 error: " << std::scientific << discr_l2_err);
  WALBERLA_LOG_INFO_ON_ROOT("DoFs: " << (uint_t) npoints);

  if (parameters.getParameter<bool>("vtkOutput")) {
    hhg::VTKOutput vtkOutput("../output", "gmg_varcoeff_vs_poly");
    vtkOutput.add(coordX.get());
    vtkOutput.add(coordY.get());
    vtkOutput.write(maxLevel, 0);
  }

  return 0;
}
