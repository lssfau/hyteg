#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  std::string meshFileName = "../data/meshes/quad_4el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  RoundRobin loadbalancer;
  setupStorage.balanceLoad( loadbalancer, 0.0 );

  size_t minLevel = 2;
  size_t maxLevel = 2;
  size_t maxiter = 10000;

  PrimitiveStorage storage( rank, setupStorage );
  //////OLD STUFF

  hhg::P1Function r("r", storage, minLevel, maxLevel);
  hhg::P1Function f("f", storage, minLevel, maxLevel);
  hhg::P1Function u("u", storage, minLevel, maxLevel);
  hhg::P1Function u_exact("u_exact", storage, minLevel, maxLevel);
  hhg::P1Function err("err", storage, minLevel, maxLevel);
  hhg::P1Function npoints_helper("npoints_helper", storage, minLevel, maxLevel);

  //hhg::P1LaplaceOperator L(storage, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return xx[0]*xx[0] - xx[1]*xx[1]; };
  std::function<real_t(const hhg::Point3D&)> rhs   = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u_exact.interpolate(exact, maxLevel);

//  auto solver = hhg::CGSolver<hhg::P1FunctionOld, hhg::P1LaplaceOperator>(mesh, minLevel, maxLevel);
//  walberla::WcTimer timer;
//  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, true);
//  timer.end();
//  fmt::printf("time was: %e\n",timer.last());
  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

  hhg::VTKWriter({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "test");
  return 0;
}
