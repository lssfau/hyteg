#include "core/timing/Timer.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/Format.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hyteg;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  walberla::shared_ptr<walberla::config::Config> cfg(new walberla::config::Config);
  cfg->readParameterFile("../../data/param/jacobi_P2.prm");
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");
  size_t level = parameters.getParameter<size_t>("level");
  size_t maxiter = parameters.getParameter<size_t>("maxiter");
  MeshInfo meshInfo = MeshInfo::fromGmshFile( parameters.getParameter<std::string>("mesh") );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  hyteg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage, timingTree);

  hyteg::P2ConstantLaplaceOperator L(storage, level, level);

  hyteg::P2Function< real_t > residuum  ("residuum",  storage, level, level);
  hyteg::P2Function< real_t > rhs       ("rhs",       storage, level, level);
  hyteg::P2Function< real_t > p2function("p2Function",storage, level, level);
  hyteg::P2Function< real_t > Lu        ("Lu",        storage, level, level);
  hyteg::P2Function< real_t > p2Exact   ("p2Exact",   storage, level, level);
  hyteg::P2Function< real_t > error     ("error",     storage, level, level);
  hyteg::P2Function< real_t > helperFun ("helperFun", storage, level, level);



  std::function<real_t(const hyteg::Point3D&)> exactFunction = [](const hyteg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  //std::function<real_t(const hyteg::Point3D&)> zeros = [](const hyteg::Point3D&) { return 0.0; };
  std::function<real_t(const hyteg::Point3D&)> ones  = [](const hyteg::Point3D&) { return 1.0; };
  std::function<real_t(const hyteg::Point3D&)> random  = [](const hyteg::Point3D&) { return real_c(rand()); };
  p2function.interpolate(random, level, hyteg::Inner);
  helperFun.interpolate(random, level, hyteg::Inner);

  p2function.interpolate(exactFunction, level, hyteg::DirichletBoundary);
  helperFun.interpolate(exactFunction, level, hyteg::DirichletBoundary);
  p2Exact.interpolate(exactFunction, level);

  real_t begin_res, abs_res_old, rel_res, abs_res = 0;

  WALBERLA_LOG_INFO_ON_ROOT( hyteg::format("%6s|%10s|%10s|%10s","iter","abs_res","rel_res","conv"));

  L.apply(p2function, Lu, level, hyteg::Inner);
  residuum.assign({1.0, -1.0}, { rhs, Lu }, level, hyteg::Inner);
  begin_res = std::sqrt(residuum.dotGlobal(residuum, level, hyteg::Inner));
  abs_res_old = begin_res;

  WALBERLA_LOG_INFO_ON_ROOT( hyteg::format("%6d|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res/abs_res_old))
  walberla::WcTimer timer;
  for(uint_t i = 0; i < maxiter; ++i) {
    helperFun.assign({1.0},{p2function},level, hyteg::Inner);
    L.smooth_jac(p2function, rhs, helperFun, level, hyteg::Inner);
    L.apply(p2function, Lu, level, hyteg::Inner);
    residuum.assign({1.0, -1.0}, { rhs, Lu }, level, hyteg::Inner);
    abs_res = std::sqrt(residuum.dotGlobal(residuum, level, hyteg::Inner));
    rel_res = abs_res / begin_res;
    WALBERLA_LOG_INFO_ON_ROOT( hyteg::format("%6d|%10.3e|%10.3e|%10.3e", i+1, abs_res, rel_res, abs_res/abs_res_old))
    WALBERLA_CHECK_LESS(abs_res,abs_res_old);
    abs_res_old = abs_res;
  }
  timer.end();

  if (parameters.getParameter<bool>("vtkOutput")) {
    VTKOutput vtkOutput("../../output", "gs_P2", storage);
    vtkOutput.add( p2function );
    vtkOutput.add( p2Exact );
    vtkOutput.add( rhs );
    vtkOutput.add( residuum );
    vtkOutput.add( error );
    vtkOutput.add( helperFun );
    vtkOutput.write( level );
  }

  WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
  error.assign({1.0, -1.0}, {p2function, p2Exact}, level);

  helperFun.interpolate(ones, level);
  real_t npoints = helperFun.dotGlobal(helperFun, level);

  real_t discr_l2_err = std::sqrt(error.dotGlobal(error, level) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);


  if (parameters.getParameter<bool>("printTiming")) {
    walberla::WcTimingTree tt = timingTree->getReduced();
    WALBERLA_LOG_INFO_ON_ROOT(tt);
  }

  return 0;
}
