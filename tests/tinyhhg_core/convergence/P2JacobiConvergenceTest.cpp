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
  cfg->readParameterFile("../../data/param/jacobi_P2.prm");
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



  std::function<real_t(const hhg::Point3D&)> exactFunction = [](const hhg::Point3D& x) { return 50.0; }; //sin(x[0])*sinh(x[1]); };
  //std::function<real_t(const hhg::Point3D&)> zeros = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };
  std::function<real_t(const hhg::Point3D&)> random  = [](const hhg::Point3D&) { return real_c(rand()); };
//  std::function<real_t(const hhg::Point3D&)> exactFunction2 = [](const hhg::Point3D& x) { return x[0]; };
//
//  p2function.interpolate(random, level, hhg::Inner);
//  helperFun.interpolate(random, level, hhg::Inner);

  //p2function.getVertexDoFFunction()->interpolate(exactFunction, level, hhg::Inner);
  //helperFun.getVertexDoFFunction()->interpolate(exactFunction, level, hhg::Inner);
  p2function.interpolate(exactFunction, level, hhg::DirichletBoundary);
  helperFun.interpolate(exactFunction, level, hhg::DirichletBoundary);
  p2Exact.interpolate(exactFunction, level);

  real_t begin_res, abs_res_old, rel_res, abs_res = 0;

  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6s|%10s|%10s|%10s","iter","abs_res","rel_res","conv"));

  L.apply(p2function, Lu, level, hhg::Inner);
  residuum.assign({1.0, -1.0}, { &rhs, &Lu }, level, hhg::Inner);
  begin_res = std::sqrt(residuum.dot(residuum, level, hhg::Inner));
  abs_res_old = begin_res;


  //if (parameters.getParameter<bool>("vtkOutput")) {
    VTKOutput vtkOutput( "../../output", "gs_P2" );
    vtkOutput.add( &p2function );
    vtkOutput.add( &p2Exact );
    vtkOutput.add( &rhs );
    vtkOutput.add( &residuum );
    vtkOutput.add( &error );
    vtkOutput.add( &helperFun );
  //}

  //auto faceDataPtr = storage->getFace(PrimitiveID(5))->getData(p2function.getVertexDoFFunction()->getFaceDataID())->getPointer(level );

  auto face0 = storage->getFace(PrimitiveID(6));
  hhg::vertexdof::macroface::printFunctionMemory< real_t >(level,*face0, p2function.getVertexDoFFunction()->getFaceDataID());
  hhg::edgedof::macroface::printFunctionMemory< real_t >(level,*face0, p2function.getEdgeDoFFunction()->getFaceDataID());
  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res/abs_res_old))
  walberla::WcTimer timer;
  for(uint_t i = 0; i < maxiter; ++i) {

    helperFun.assign({1.0},{&p2function},level,hhg::Inner);
    //vtkOutput.write( level, i );
    L.smooth_jac(p2function, rhs, helperFun, level, hhg::Inner);
    hhg::vertexdof::macroface::printFunctionMemory< real_t > (level,*face0, p2function.getVertexDoFFunction()->getFaceDataID());
    hhg::edgedof::macroface::printFunctionMemory< real_t >(level,*face0, p2function.getEdgeDoFFunction()->getFaceDataID());
    L.apply(p2function, Lu, level, hhg::Inner);
    residuum.assign({1.0, -1.0}, { &rhs, &Lu }, level, hhg::Inner);
    abs_res = std::sqrt(residuum.dot(residuum, level, hhg::Inner));
    rel_res = abs_res / begin_res;
    WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e", i+1, abs_res, rel_res, abs_res/abs_res_old))
    //WALBERLA_CHECK_LESS(abs_res,abs_res_old);
    abs_res_old = abs_res;
  }
  timer.end();

  WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
  error.assign({1.0, -1.0}, {&p2function, &p2Exact}, level);

  helperFun.interpolate(ones, level);
  real_t npoints = helperFun.dot(helperFun, level);

  real_t discr_l2_err = std::sqrt(error.dot(error, level) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);


  if (parameters.getParameter<bool>("printTiming")) {
    walberla::WcTimingTree tt = timingTree->getReduced();
    WALBERLA_LOG_INFO_ON_ROOT(tt);
  }

  return 0;
}
