#include <tinyhhg_core/tinyhhg.hpp>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/quad_4el_neumann_freeslip.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t minLevel = 2;
  size_t maxLevel = 5;
//  size_t maxiter = 10000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P1Function<real_t> p("p", storage, minLevel, maxLevel);

  hhg::P1Function<real_t> tmp_x("tmp_x", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> tmp_y("tmp_y", storage, minLevel, maxLevel);

  hhg::P1Function<real_t> n_x("n_x", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> n_y("n_y", storage, minLevel, maxLevel);

//  hhg::P1StokesFunction<real_t> r("r", storage, minLevel, maxLevel);
//  hhg::P1StokesFunction<real_t> f("f", storage, minLevel, maxLevel);
//  hhg::P1StokesFunction<real_t> u("u", storage, minLevel, maxLevel);

  hhg::P1StokesOperator L(storage, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> one = [](const hhg::Point3D&) {
    return 1.0;
  };

  std::function<real_t(const hhg::Point3D&, const std::vector<real_t>&)> normalizer = [](const hhg::Point3D&, const std::vector<real_t>& val) {
    return -val[0] / std::sqrt(val[0]*val[0] + val[1]*val[1]);
  };

  p.interpolate(one, maxLevel, hhg::All);
  L.divT_x.apply(p, tmp_x, maxLevel, hhg::All);
  L.divT_y.apply(p, tmp_y, maxLevel, hhg::All);

  n_x.interpolateExtended(normalizer, { &tmp_x, &tmp_y }, maxLevel, hhg::Boundary);
  n_y.interpolateExtended(normalizer, { &tmp_y, &tmp_x }, maxLevel, hhg::Boundary);

//  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
//  std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
//  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };
//
//  u.u.interpolate(bc_x, maxLevel, hhg::DirichletBoundary);
//  u.v.interpolate(zero, maxLevel, hhg::DirichletBoundary);
//
//  auto solver = hhg::MinResSolver<hhg::P1StokesFunction<real_t>, hhg::P1StokesOperator>(storage, minLevel, maxLevel);
//  solver.solve(L, u, f, r, maxLevel, 1e-12, maxiter, hhg::Inner | hhg::NeumannBoundary, true);

  hhg::VTKWriter<hhg::P1Function<real_t>, hhg::DGFunction<real_t >>({&n_x, &n_y}, {}, maxLevel, "../output", "stokes_freeslip");
  return EXIT_SUCCESS;
}
