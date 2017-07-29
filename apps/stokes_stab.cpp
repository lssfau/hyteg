#include <tinyhhg_core/tinyhhg.hpp>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/bfs_12el_neumann.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::RoundRobin loadbalancer;
  setupStorage.balanceLoad( loadbalancer, 0.0 );

  size_t minLevel = 2;
  size_t maxLevel = 4;
  size_t maxiter = 10000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P1StokesFunction r("r", storage, minLevel, maxLevel);
  hhg::P1StokesFunction f("f", storage, minLevel, maxLevel);
  hhg::P1StokesFunction u("u", storage, minLevel, maxLevel);

  hhg::P1StokesOperator L(storage, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> bc_x = [](const hhg::Point3D& x) {
    if (x[0] < 1e-8)
    {
      return 16.0 * (x[1]-0.5) * (1.0 - x[1]);
    }
    else
    {
      return 0.0;
    }
  };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate(bc_x, maxLevel, hhg::DirichletBoundary);
  u.v.interpolate(zero, maxLevel, hhg::DirichletBoundary);

  auto solver = hhg::MinResSolver<hhg::P1StokesFunction, hhg::P1StokesOperator>(storage, minLevel, maxLevel);
  solver.solve(L, u, f, r, maxLevel, 1e-12, maxiter, hhg::Inner | hhg::NeumannBoundary, true);

  hhg::VTKWriter({ &u.u, &u.v, &u.p }, maxLevel, "../output", "stokes_stab_test");
  return EXIT_SUCCESS;
}
