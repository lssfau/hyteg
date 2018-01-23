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

  size_t level = 5;

  MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_2el.msh" );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage, timingTree);

  hhg::P2ConstantLaplaceOperator L(storage, level, level);

  hhg::P2Function< real_t > res("r", storage, level, level);
  hhg::P2Function< real_t > r("r", storage, level, level);
  hhg::P2Function< real_t > u("u", storage, level, level);
  hhg::P2Function< real_t > u_exact("u_exact", storage, level, level);
  hhg::P2Function< real_t > err("err", storage, level, level);
  hhg::P2Function< real_t > npoints_helper("npoints_helper", storage, level, level);



  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, level, hhg::DirichletBoundary);
  u_exact.interpolate(exact, level);


  WALBERLA_LOG_INFO_ON_ROOT("residual: = " << std::sqrt(res.dot(res, level, hhg::Inner)));
  for(uint i = 0; i < 100; ++i) {

    L.smooth_gs(u, r, level, hhg::Inner);
    L.apply(u,res,level,hhg::Inner);
    res.add({-1},{&r},level,hhg::Inner);
    WALBERLA_LOG_INFO_ON_ROOT("residual: = " << std::sqrt(res.dot(res, level, hhg::Inner)));


  }

  VTKOutput vtkOutput( "../../output", "cg_P2" );
  vtkOutput.add( &u );
  vtkOutput.write( level );

  err.assign({1.0, -1.0}, {&u, &u_exact}, level);

  npoints_helper.interpolate(ones, level);
  real_t npoints = npoints_helper.dot(npoints_helper, level);

  real_t discr_l2_err = std::sqrt(err.dot(err, level) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);



  //walberla::WcTimingTree tt = timingTree->getReduced();
  //WALBERLA_LOG_INFO_ON_ROOT( tt );

  return 0;
}
