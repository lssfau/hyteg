#include <tinyhhg_core/tinyhhg.hpp>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/quad_4el.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t minLevel = 2;
  size_t maxLevel = 5;
  size_t maxiter = 10000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P1Function<real_t> tmp("tmp", storage, minLevel, maxLevel);
  hhg::P1StokesFunction<real_t> r("r", storage, minLevel, maxLevel);
  hhg::P1StokesFunction<real_t> f("f", storage, minLevel, maxLevel);
  hhg::P1StokesFunction<real_t> u("u", storage, minLevel, maxLevel);
  std::shared_ptr<hhg::P1Function< real_t >> coefficient = std::make_shared<hhg::P1Function< real_t >>("coeff", storage, minLevel, maxLevel);

  hhg::P1MassOperator M(storage, minLevel, maxLevel);
  hhg::P1CoefficientStokesOperator L(storage, coefficient, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> coeff = [](const hhg::Point3D& x) { return ((0.000112225535684453*exp(17*x[1]) + 1)*(0.89*exp(-10*pow(x[0] - 0.5, 2) - 30*pow(x[1] - (-sqrt(pow(x[0] - 0.5, 2)) + 0.5)*(1.5*sqrt(pow(x[0] - 0.5, 2)) + 0.75) - 0.3, 2)) + 1)*exp(100*pow(-x[0] + 0.5, 2)) + 0.98)*exp(-100*pow(-x[0] + 0.5, 2))/(0.000112225535684453*exp(17*x[1]) + 1); };
  std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate(zero, maxLevel, hhg::DirichletBoundary);
  u.v.interpolate(zero, maxLevel, hhg::DirichletBoundary);
  coefficient->interpolate(coeff, maxLevel);
  tmp.interpolate(coeff, maxLevel);
  M.apply(tmp, f.v, maxLevel, hhg::All);

  auto solver = hhg::MinResSolver<hhg::P1StokesFunction<real_t>, hhg::P1CoefficientStokesOperator>(storage, minLevel, maxLevel);
  solver.solve(L, u, f, r, maxLevel, 1e-12, maxiter, hhg::Inner | hhg::NeumannBoundary, true);

  // u_u*iHat + u_v*jHat
  hhg::VTKWriter<hhg::P1Function< real_t >>({ &u.u, &u.v, &u.p, coefficient.get() }, maxLevel, "../output", "stokes_stab_varcoeff");
  return EXIT_SUCCESS;
}
