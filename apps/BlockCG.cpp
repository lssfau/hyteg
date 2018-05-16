#include "core/timing/Timer.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"

#include "tinyhhg_core/composites/P1BlockLaplaceOperator.hpp"

#include "tinyhhg_core/solvers/CGSolver.hpp"

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG CG Test\n");

  std::string meshFileName = "../data/meshes/quad_4el.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t minLevel = 2;
  size_t maxLevel = 2;
  size_t maxiter = 10000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P1StokesFunction<real_t> r("r", storage, minLevel, maxLevel);
  hhg::P1StokesFunction<real_t> f("f", storage, minLevel, maxLevel);
  hhg::P1StokesFunction<real_t> u("u", storage, minLevel, maxLevel);
  hhg::P1StokesFunction<real_t> u_exact("u_exact", storage, minLevel, maxLevel);
  hhg::P1StokesFunction<real_t> err("err", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > npoints_helper("npoints_helper", storage, minLevel, maxLevel);

  hhg::P1BlockLaplaceOperator L(storage, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return xx[0]*xx[0] - xx[1]*xx[1]; };
  std::function<real_t(const hhg::Point3D&)> rhs   = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u.v.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u.p.interpolate(exact, maxLevel, hhg::DirichletBoundary);

  u_exact.u.interpolate(exact, maxLevel);
  u_exact.v.interpolate(exact, maxLevel);
  u_exact.p.interpolate(exact, maxLevel);

  auto solver = hhg::CGSolver<hhg::P1StokesFunction<real_t>, hhg::P1BlockLaplaceOperator>(storage, minLevel, maxLevel);
  walberla::WcTimer timer;
  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, true);
  timer.end();
  WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);

  real_t npoints = npoints_helper.dotGlobal(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dotGlobal(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

//  hhg::VTKWriter({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "test");
  return 0;
}
