#include <tinyhhg_core/tinyhhg.hpp>

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  hhg::Mesh mesh("../data/meshes/bfs_12el_neumann.msh");

  size_t minLevel = 2;
  const size_t maxLevel = 4;
  size_t maxiter = 10000;

  hhg::P1StokesFunction r("r", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction f("f", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction u("u", mesh, minLevel, maxLevel);

  hhg::P1StokesOperator L(mesh, minLevel, maxLevel);

  std::function<walberla::real_t(const hhg::Point3D&)> bc_x = [](const hhg::Point3D& x) {
    if (x[0] < 1e-8)
    {
      return 16.0 * (x[1]-0.5) * (1.0 - x[1]);
    }
    else
    {
      return 0.0;
    }
  };
  std::function<walberla::real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<walberla::real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
  std::function<walberla::real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate<maxLevel>(bc_x, hhg::DirichletBoundary);
  u.v.interpolate<maxLevel>(zero, hhg::DirichletBoundary);

  auto solver = hhg::MinResSolver<hhg::P1StokesFunction, hhg::P1StokesOperator>(mesh, minLevel, maxLevel);
  solver.solve<maxLevel>(L, u, f, r, 1e-12, maxiter, hhg::Inner | hhg::NeumannBoundary, true);

  hhg::VTKWriter({ &u.u, &u.v, &u.p }, maxLevel, "../output", "stokes_stab_test");
  return EXIT_SUCCESS;
}
