#include <tinyhhg_core/tinyhhg.hpp>

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/quad_4el.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::RoundRobin loadbalancer;
  setupStorage.balanceLoad( loadbalancer, 0.0 );

  size_t minLevel = 2;
  size_t maxLevel = 5;
  size_t maxiter = 1000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P1Function r("r", storage, minLevel, maxLevel);
  hhg::P1Function f("f", storage, minLevel, maxLevel);
  hhg::P1Function u("u", storage, minLevel, maxLevel);
  hhg::P1Function u_exact("u_exact", storage, minLevel, maxLevel);
  hhg::P1Function err("err", storage, minLevel, maxLevel);
  hhg::P1Function npoints_helper("npoints_helper", storage, minLevel, maxLevel);

  hhg::P1LaplaceOperator L(storage, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) -> real_t { return x[0]*x[0] - x[1]*x[1]; };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u_exact.interpolate(exact, maxLevel);

  auto solver = hhg::MinResSolver<hhg::P1Function, hhg::P1LaplaceOperator>(storage, minLevel, maxLevel);
  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, true);

  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("discrete L2 error = {:e}", discr_l2_err));

  hhg::VTKWriter<hhg::P1Function>({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "minres");
  return EXIT_SUCCESS;
}
