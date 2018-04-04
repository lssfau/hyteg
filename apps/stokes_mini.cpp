

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/quad_4el_neumann.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile(meshFileName);
  hhg::SetupPrimitiveStorage setupStorage(meshInfo, walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t minLevel = 2;
  size_t maxLevel = 4;
  size_t maxiter = 1000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::MiniStokesFunction<real_t> r("r", storage, minLevel, maxLevel);
  hhg::MiniStokesFunction<real_t> f("f", storage, minLevel, maxLevel);
  hhg::MiniStokesFunction<real_t> u("u", storage, minLevel, maxLevel);

//  hhg::MiniStokesFunction numerator("numerator", storage, minLevel, maxLevel);

  hhg::MiniStokesOperator L(storage, minLevel, maxLevel);

//  size_t num = 1;
//  numerator.enumerate(maxLevel, num);
//
//  std::ofstream fileL("../output/L.txt");
//  L.save(numerator, numerator, fileL, maxLevel, hhg::All);
//  return 0;

  std::function<real_t(const hhg::Point3D&)> bc_x = [](const hhg::Point3D& x) {
    return 4.0 * (1.0-x[1]) * x[1];
  };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate(zero, maxLevel);
  u.u.interpolate(bc_x, maxLevel, hhg::DirichletBoundary);
  u.v.interpolate(zero, maxLevel, hhg::DirichletBoundary);


  auto solver = hhg::MinResSolver<hhg::MiniStokesFunction<real_t>, hhg::MiniStokesOperator>(storage, minLevel, maxLevel);
  solver.solve(L, u, f, r, maxLevel, 1e-12, maxiter, hhg::Inner | hhg::NeumannBoundary, true);

  return EXIT_SUCCESS;
}
