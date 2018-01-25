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
  cfg->readParameterFile("../data/param/cg_P2.prm");
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

  size_t level = parameters.getParameter<size_t>("level");
  size_t maxiter = parameters.getParameter<size_t>("maxiter");
  real_t tolerance = parameters.getParameter<real_t>("tolerance");

  MeshInfo meshInfo = MeshInfo::fromGmshFile( parameters.getParameter<std::string>("mesh") );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage, timingTree);

  hhg::P2ConstantLaplaceOperator L(storage, level, level);

  hhg::P2Function< real_t > r("r", storage, level, level);
  hhg::P2Function< real_t > f("f", storage, level, level);
  hhg::P2Function< real_t > u("u", storage, level, level);
  hhg::P2Function< real_t > u_exact("u_exact", storage, level, level);
  hhg::P2Function< real_t > err("err", storage, level, level);
  hhg::P2Function< real_t > npoints_helper("npoints_helper", storage, level, level);



  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, level, hhg::DirichletBoundary);
  u_exact.interpolate(exact, level);

  auto solver = hhg::CGSolver<hhg::P2Function< real_t >, hhg::P2ConstantLaplaceOperator>(storage, level, level);
  walberla::WcTimer timer;
  solver.solve(L, u, f, r, level, tolerance, maxiter, hhg::Inner, true);
  timer.end();

  WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
  err.assign({1.0, -1.0}, {&u, &u_exact}, level);

  npoints_helper.interpolate(ones, level);
  real_t npoints = npoints_helper.dot(npoints_helper, level);

  real_t discr_l2_err = std::sqrt(err.dot(err, level) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

  VTKOutput vtkOutput( "../output", "cg_P2" );
  vtkOutput.add( &u );
  vtkOutput.add( &u_exact );
  vtkOutput.add( &f );
  vtkOutput.add( &r );
  vtkOutput.add( &err );
  vtkOutput.add( &npoints_helper );
  vtkOutput.write( level );

  walberla::WcTimingTree tt = timingTree->getReduced();
  WALBERLA_LOG_INFO_ON_ROOT( tt );

  return 0;
}
