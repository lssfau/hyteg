#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <core/Environment.h>
#include <core/config/Config.h>

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
  cfg->readParameterFile("../data/param/gmg_P2.prm");
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

  const uint_t minLevel = parameters.getParameter<uint_t>("minLevel");
  const uint_t maxLevel = parameters.getParameter<uint_t>("maxLevel");
  const uint_t max_outer_iter =  parameters.getParameter<uint_t>("max_outer_iter");
  const uint_t max_cg_iter =  parameters.getParameter<uint_t>("max_cg_iter");
  const real_t mg_tolerance = parameters.getParameter<real_t>("mg_tolerance");
  const real_t coarse_tolerance = parameters.getParameter<real_t>("coarse_tolerance");
  const uint_t nuPre = 3;
  const uint_t nuPost = 3;


  MeshInfo meshInfo = MeshInfo::fromGmshFile( parameters.getParameter<std::string>("mesh") );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage, timingTree);

  hhg::P1Function< real_t > r_p1("r_p1", storage, minLevel, maxLevel);
  std::shared_ptr<hhg::P1Function< real_t > > f_p1 = std::make_shared<hhg::P1Function< real_t >>("f_p1", storage, minLevel, maxLevel);
  std::shared_ptr<hhg::P1Function< real_t > > u_p1 = std::make_shared<hhg::P1Function< real_t >>("u_p1", storage, minLevel, maxLevel);

  hhg::P2Function< real_t > r_p2("r_p2", storage, maxLevel, maxLevel);
  hhg::P2Function< real_t > f_p2("f_p2", storage, maxLevel, maxLevel);
  hhg::P2Function< real_t > u_p2("u_p2", storage, maxLevel, maxLevel);
  hhg::P2Function< real_t > Lu_p2("Lu_p2", storage, maxLevel, maxLevel);
  hhg::P2Function< real_t > tmp_p2("tmp_p2", storage, maxLevel, maxLevel);
  hhg::P2Function< real_t > u_exact_p2("u_exact_p2", storage, maxLevel, maxLevel);
  hhg::P2Function< real_t > err_p2("err_p2", storage, maxLevel, maxLevel);
  hhg::P2Function< real_t > npoints_helper_p2("npoints_helper_p2", storage, maxLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0; };
  std::function<real_t(const hhg::Point3D&)> zero  = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  WALBERLA_LOG_INFO_ON_ROOT("Interpolating u");
  u_p2.interpolate(exact, maxLevel, hhg::DirichletBoundary);

  WALBERLA_LOG_INFO_ON_ROOT("Interpolating exact function");
  u_exact_p2.interpolate(exact, maxLevel);
//  WALBERLA_LOG_INFO_ON_ROOT("Interpolating and integrating rhs");
//  npoints_helper.interpolate(rhs, maxLevel);
//  M.apply(npoints_helper, f, maxLevel, hhg::All);

  WALBERLA_LOG_INFO_ON_ROOT("Setting up stiffness operator");
  auto start = walberla::timing::getWcTime();
  hhg::P1LaplaceOperator L_p1(storage, minLevel, maxLevel);
  hhg::P2ConstantLaplaceOperator L_p2(storage, maxLevel, maxLevel);
  auto end = walberla::timing::getWcTime();
  real_t setupTime = end - start;

  npoints_helper_p2.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper_p2.dot(npoints_helper_p2, maxLevel);

  typedef hhg::CGSolver<hhg::P1Function<real_t>, hhg::P1LaplaceOperator> CoarseSolver;
  auto coarseLaplaceSolver = std::make_shared<CoarseSolver>(storage, minLevel, minLevel);
  typedef GMultigridSolver<hhg::P1Function<real_t>, hhg::P1LaplaceOperator, CoarseSolver> LaplaceSover;
  LaplaceSover laplaceSolver(storage, coarseLaplaceSolver, minLevel, maxLevel);

  WALBERLA_LOG_INFO_ON_ROOT("Starting V cycles");
  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6s|%10s|%10s|%10s|%10s|%10s","iter","abs_res","rel_res","conv","L2-error","Time"));

  real_t rel_res = 1.0;

  L_p2.apply(u_p2, Lu_p2, maxLevel, hhg::Inner);
  r_p2.assign({1.0, -1.0}, {&f_p2, &Lu_p2}, maxLevel, hhg::Inner);

  real_t begin_res = std::sqrt(r_p2.dot(r_p2, maxLevel, hhg::Inner));
  real_t abs_res_old = begin_res;

  err_p2.assign({1.0, -1.0}, {&u_p2, &u_exact_p2}, maxLevel);
  real_t discr_l2_err = std::sqrt(err_p2.dot(err_p2, maxLevel) / npoints);

  //WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  -", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err));
  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err,0))

  real_t solveTime = real_c(0.0);
  real_t averageConvergenceRate = real_c(0.0);
  const uint_t convergenceStartIter = 3;

  uint_t i = 0;
  for (; i < max_outer_iter; ++i)
  {
    start = walberla::timing::getWcTime();

    // pre-smooth
    for (size_t nu = 0; nu < nuPre; ++nu)
    {
      L_p2.smooth_gs(u_p2, f_p2, maxLevel, hhg::Inner);
    }

    // compute residuum
    L_p2.apply(u_p2, Lu_p2, maxLevel, hhg::Inner);
    r_p2.assign({1.0, -1.0}, { &f_p2, &Lu_p2 }, maxLevel, hhg::Inner);

    // restrict
    r_p2.restrictP2ToP1(f_p1, maxLevel, hhg::Inner);

    u_p1->interpolate(zero, maxLevel);

    // Apply P1 geometric multigrid solver
    laplaceSolver.solve(L_p1, *u_p1, *f_p1, r_p1, maxLevel, coarse_tolerance, max_cg_iter, hhg::Inner, LaplaceSover::CycleType::VCYCLE, false);

    // prolongate
    tmp_p2.assign({1.0}, { &u_p2 }, maxLevel, hhg::Inner);
    u_p2.prolongateP1ToP2(u_p1, maxLevel, hhg::Inner);
    u_p2.add({1.0}, { &tmp_p2 }, maxLevel, hhg::Inner);

    // post-smooth
    for (size_t nu = 0; nu < nuPost; ++nu)
    {
      L_p2.smooth_gs(u_p2, f_p2, maxLevel, hhg::Inner);
    }

    end = walberla::timing::getWcTime();


    L_p2.apply(u_p2, Lu_p2, maxLevel, hhg::Inner);
    r_p2.assign({1.0, -1.0}, { &f_p2, &Lu_p2 }, maxLevel, hhg::Inner);
    real_t abs_res = std::sqrt(r_p2.dot(r_p2, maxLevel, hhg::Inner));
    rel_res = abs_res / begin_res;
    err_p2.assign({1.0, -1.0}, { &u_p2, &u_exact_p2 }, maxLevel);
    discr_l2_err = std::sqrt(err_p2.dot(err_p2, maxLevel) / npoints);

    //WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  {:e}", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err, end-start));
    WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", i+1, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err,end - start))
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

  WALBERLA_LOG_INFO_ON_ROOT("Setup time: " << std::defaultfloat << setupTime);
  WALBERLA_LOG_INFO_ON_ROOT("Solve time " << std::defaultfloat << solveTime);
  WALBERLA_LOG_INFO_ON_ROOT("Time to solution: " << std::defaultfloat << setupTime + solveTime);
  WALBERLA_LOG_INFO_ON_ROOT("Avg. convergence rate: " << std::scientific << averageConvergenceRate / real_c(i-convergenceStartIter));
  WALBERLA_LOG_INFO_ON_ROOT("L^2 error: " << std::scientific << discr_l2_err);
  WALBERLA_LOG_INFO_ON_ROOT("DoFs: " << (uint_t) npoints);

  if (parameters.getParameter<bool>("vtkOutput")) {
    VTKOutput vtkOutput("../output", "gmg_P2");
    vtkOutput.add(&u_p2);
    vtkOutput.add(&u_exact_p2);
    vtkOutput.add(&f_p2);
    vtkOutput.add(&r_p2);
    vtkOutput.add(&err_p2);
    vtkOutput.add(&npoints_helper_p2);
    vtkOutput.write(maxLevel);
  }

  walberla::WcTimingTree tt = timingTree->getReduced();
  WALBERLA_LOG_INFO_ON_ROOT( tt );

  return 0;
}
