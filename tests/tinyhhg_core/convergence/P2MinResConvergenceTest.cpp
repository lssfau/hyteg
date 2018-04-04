#include <string>

#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/minressolver.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"

using walberla::real_t;

int main(int argc, char* argv[]) {
  walberla::MPIManager::instance()->initializeMPI(&argc, &argv);
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../../data/meshes/quad_4el.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile(meshFileName);
  hhg::SetupPrimitiveStorage setupStorage(meshInfo, walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  hhg::loadbalancing::roundRobin(setupStorage);

  size_t minLevel = 2;
  size_t maxLevel = 3;
  size_t maxiter = 1000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P2Function<real_t> r("r", storage, minLevel, maxLevel);
  hhg::P2Function<real_t> f("f", storage, minLevel, maxLevel);
  hhg::P2Function<real_t> u("u", storage, minLevel, maxLevel);
  hhg::P2Function<real_t> u_exact("u_exact", storage, minLevel, maxLevel);
  hhg::P2Function<real_t> err("err", storage, minLevel, maxLevel);
  hhg::P2Function<real_t> npoints_helper("npoints_helper", storage, minLevel, maxLevel);

  hhg::P2ConstantLaplaceOperator L(storage, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D &)> exact = [](const hhg::Point3D &x) -> real_t { return x[0] * x[0] - x[1] * x[1]; };
  std::function<real_t(const hhg::Point3D &)> rhs = [](const hhg::Point3D &) { return 0.0; };
  std::function<real_t(const hhg::Point3D &)> ones = [](const hhg::Point3D &) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u_exact.interpolate(exact, maxLevel);

  auto solver = hhg::MinResSolver<hhg::P2Function<real_t>, hhg::P2ConstantLaplaceOperator>(storage, minLevel, maxLevel);
  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, true);

  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << std::scientific << discr_l2_err);
  WALBERLA_CHECK_LESS( discr_l2_err, 2e-8 );
  return EXIT_SUCCESS;
}