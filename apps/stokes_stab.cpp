#include <tinyhhg_core/tinyhhg.hpp>

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  
  hhg::Mesh mesh("../data/meshes/quad_2el_neumann.msh");

  size_t minLevel = 2;
  size_t maxLevel = 5;
  size_t maxiter = 10000;

  hhg::P1StokesFunction r("r", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction f("f", mesh, minLevel, maxLevel);
  hhg::P1StokesFunction u("u", mesh, minLevel, maxLevel);

  hhg::P1StokesOperator L(mesh, minLevel, maxLevel);

  std::function<double(const hhg::Point3D&)> u_exact = [](const hhg::Point3D& x) -> double { return 4.0 * x[1] * (1.0 - x[1]); };
  std::function<double(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<double(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
  std::function<double(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate(u_exact, maxLevel, hhg::DirichletBoundary);
  u.v.interpolate(zero, maxLevel, hhg::DirichletBoundary);

  auto solver = hhg::MinResSolver<hhg::P1StokesFunction, hhg::P1StokesOperator>(mesh, minLevel, maxLevel);
  solver.solve(L, u, f, r, maxLevel, 1e-12, maxiter, hhg::Inner | hhg::NeumannBoundary, true);

  hhg::VTKWriter({ &u.u, &u.v, &u.p }, maxLevel, "../output", "stokes_stab_test");
  return EXIT_SUCCESS;
}
