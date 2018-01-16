#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>
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

  walberla::Config::BlockHandle parameters;
  if(walberlaEnv.config()) {
    parameters = walberlaEnv.config()->getOneBlock("Parameters");
  } else {
    walberla::shared_ptr<walberla::config::Config> cfg(new walberla::config::Config);
    cfg->readParameterFile("../data/param/cg_P2.prm");
    parameters = walberlaEnv.config()->getOneBlock("Parameters");
  }

  size_t minLevel = parameters.getParameter<size_t>("minlevel");
  size_t maxLevel = parameters.getParameter<size_t>("maxlevel");
  size_t maxiter = parameters.getParameter<size_t>("maxiter");
  real_t tolerance = parameters.getParameter<real_t>("tolerance");

  MeshInfo meshInfo = MeshInfo::fromGmshFile( parameters.getParameter<std::string>("mesh") );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage, timingTree);

  hhg::P2ConstantLaplaceOperator L(storage, minLevel, maxLevel);

  hhg::P2Function< real_t > r("r", storage, minLevel, maxLevel);
  hhg::P2Function< real_t > f("f", storage, minLevel, maxLevel);
  hhg::P2Function< real_t > u("u", storage, minLevel, maxLevel);
  hhg::P2Function< real_t > u_exact("u_exact", storage, minLevel, maxLevel);
  hhg::P2Function< real_t > err("err", storage, minLevel, maxLevel);
  hhg::P2Function< real_t > npoints_helper("npoints_helper", storage, minLevel, maxLevel);



  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u_exact.interpolate(exact, maxLevel);

  auto solver = hhg::CGSolver<hhg::P2Function< real_t >, hhg::P2ConstantLaplaceOperator>(storage, minLevel, maxLevel, 30);
  walberla::WcTimer timer;
  solver.solve(L, u, f, r, maxLevel, tolerance, maxiter, hhg::Inner, true);
  timer.end();

  WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

  VTKOutput vtkOutput( "../output", "cg_P2" );
  vtkOutput.add( &u );
  vtkOutput.add( &u_exact );
  vtkOutput.add( &f );
  vtkOutput.add( &r );
  vtkOutput.add( &err );
  vtkOutput.add( &npoints_helper );
  vtkOutput.write( maxLevel );

  walberla::WcTimingTree tt = timingTree->getReduced();
  WALBERLA_LOG_INFO_ON_ROOT( tt );

  return 0;
}
