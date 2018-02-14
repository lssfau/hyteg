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
  cfg->readParameterFile("../../data/param/gs_P2.prm");
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
  hhg::P2Function< real_t > Lu("u", storage, level, level);
  hhg::P2Function< real_t > u_exact("u_exact", storage, level, level);
  hhg::P2Function< real_t > err("err", storage, level, level);
  hhg::P2Function< real_t > npoints_helper("npoints_helper", storage, level, level);



  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, level, hhg::DirichletBoundary);
  u_exact.interpolate(exact, level);

  real_t begin_res, abs_res_old, rel_res, abs_res = 0;

  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6s|%10s|%10s|%10s","iter","abs_res","rel_res","conv"));

  L.apply(u, Lu, level, hhg::Inner);
  r.assign({1.0, -1.0}, { &f, &Lu }, level, hhg::Inner);
  begin_res = std::sqrt(r.dot(r, level, hhg::Inner));
  abs_res_old = begin_res;

  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res/abs_res_old))
  walberla::WcTimer timer;
  for(uint_t i = 0; i < maxiter; ++i) {
    L.smooth_gs(u, f, level, hhg::Inner);
    L.apply(u, Lu, level, hhg::Inner);
    r.assign({1.0, -1.0}, { &f, &Lu }, level, hhg::Inner);
    abs_res = std::sqrt(r.dot(r, level, hhg::Inner));
    rel_res = abs_res / begin_res;
    WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e", i+1, abs_res, rel_res, abs_res/abs_res_old))
    abs_res_old = abs_res;
  }
  timer.end();

  WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
  err.assign({1.0, -1.0}, {&u, &u_exact}, level);

  npoints_helper.interpolate(ones, level);
  real_t npoints = npoints_helper.dot(npoints_helper, level);

  real_t discr_l2_err = std::sqrt(err.dot(err, level) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

  if (parameters.getParameter<bool>("vtkOutput")) {
    VTKOutput vtkOutput( "../../output", "gs_P2" );
    vtkOutput.add( &u );
    vtkOutput.add( &u_exact );
    vtkOutput.add( &f );
    vtkOutput.add( &r );
    vtkOutput.add( &err );
    vtkOutput.add( &npoints_helper );
    vtkOutput.write( level );
  }

  if (parameters.getParameter<bool>("printTiming")) {
    walberla::WcTimingTree tt = timingTree->getReduced();
    WALBERLA_LOG_INFO_ON_ROOT(tt);
  }

  return 0;
}
