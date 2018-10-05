#include "core/timing/Timer.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Random.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/Format.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

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
  MeshInfo meshInfo = MeshInfo::fromGmshFile( parameters.getParameter<std::string>("mesh") );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  hhg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage, timingTree);

  hhg::P2ConstantLaplaceOperator L(storage, level, level);

  hhg::P2Function< real_t > residuum  ("residuum",  storage, level, level);
  hhg::P2Function< real_t > rhs       ("rhs",       storage, level, level);
  hhg::P2Function< real_t > p2function("p2Function",storage, level, level);
  hhg::P2Function< real_t > Lu        ("Lu",        storage, level, level);
  hhg::P2Function< real_t > p2Exact   ("p2Exact",   storage, level, level);
  hhg::P2Function< real_t > error     ("error",     storage, level, level);
  hhg::P2Function< real_t > helperFun ("helperFun", storage, level, level);


  std::function<real_t(const hhg::Point3D&)> exactFunction = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  //std::function<real_t(const hhg::Point3D&)> zeros = [](const hhg::Point3D&) { return 0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };
  walberla::math::seedRandomGenerator(0);
  std::function<real_t(const Point3D &)> rand = [](const Point3D &) { return walberla::math::realRandom(0.0, 1.0); };

  p2function.interpolate(exactFunction, level, hhg::DirichletBoundary);
  p2function.interpolate(rand, level, hhg::Inner);
  p2Exact.interpolate(exactFunction, level);

  real_t begin_res, abs_res_old, rel_res, abs_res = 0;

  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6s|%10s|%10s|%10s","iter","abs_res","rel_res","conv"));

  L.apply(p2function, Lu, level, hhg::Inner);
  residuum.assign({1.0, -1.0}, { &rhs, &Lu }, level, hhg::Inner);
  begin_res = std::sqrt(residuum.dotGlobal(residuum, level, hhg::Inner));
  abs_res_old = begin_res;

  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res/abs_res_old))
  walberla::WcTimer timer;
  for(uint_t i = 0; i < maxiter; ++i) {
    L.smooth_gs(p2function, rhs, level, hhg::Inner);
    L.apply(p2function, Lu, level, hhg::Inner);
    residuum.assign({1.0, -1.0}, { &rhs, &Lu }, level, hhg::Inner);
    abs_res = std::sqrt(residuum.dotGlobal(residuum, level, hhg::Inner));
    rel_res = abs_res / begin_res;
    WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e", i+1, abs_res, rel_res, abs_res/abs_res_old))
    WALBERLA_CHECK_LESS(abs_res,abs_res_old);
    abs_res_old = abs_res;
  }
  timer.end();

  WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
  error.assign({1.0, -1.0}, {&p2function, &p2Exact}, level);

  helperFun.interpolate(ones, level);
  real_t npoints = helperFun.dotGlobal(helperFun, level);

  real_t discr_l2_err = std::sqrt(error.dotGlobal(error, level) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

  if (parameters.getParameter<bool>("vtkOutput")) {
    VTKOutput vtkOutput("../../output", "gs_P2", storage);
    vtkOutput.add( &p2function );
    vtkOutput.add( &p2Exact );
    vtkOutput.add( &rhs );
    vtkOutput.add( &residuum );
    vtkOutput.add( &error );
    vtkOutput.add( &helperFun );
    vtkOutput.write( level );
  }

  if (parameters.getParameter<bool>("printTiming")) {
    walberla::WcTimingTree tt = timingTree->getReduced();
    WALBERLA_LOG_INFO_ON_ROOT(tt);
  }

  return 0;
}
