#include <tinyhhg_core/tinyhhg.hpp>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/quad_2el.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t minLevel = 2;
  size_t maxLevel = 5;
  size_t coarseMaxiter = 20;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  auto r = std::make_shared<hhg::P1StokesFunction<real_t>>("r", storage, minLevel, maxLevel);
  auto f = std::make_shared<hhg::P1StokesFunction<real_t>>("f", storage, minLevel, maxLevel);
  auto u = std::make_shared<hhg::P1StokesFunction<real_t>>("u", storage, minLevel, maxLevel);

  auto tmp = std::make_shared<hhg::P1Function<real_t>>("tmp", storage, minLevel, maxLevel);

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
  std::function<real_t(const hhg::Point3D&)> rand = [](const hhg::Point3D&) { return static_cast <real_t> (std::rand()) / static_cast <real_t> (RAND_MAX); };

  u->u.interpolate(rand, maxLevel, hhg::Inner);
  u->v.interpolate(rand, maxLevel, hhg::Inner);
  u->p.interpolate(rand, maxLevel, hhg::All);

  u->u.interpolate(zero, maxLevel, hhg::DirichletBoundary);
  u->v.interpolate(zero, maxLevel, hhg::DirichletBoundary);

  L.apply(*u, *r, maxLevel, hhg::Inner | hhg::NeumannBoundary);
  r->assign({1.0, -1.0}, { f.get(), r.get() }, maxLevel, hhg::Inner | hhg::NeumannBoundary);
  WALBERLA_LOG_DEVEL("[Uzawa] residuum: " << std::scientific << std::sqrt(r->dot(*r, maxLevel, hhg::Inner | hhg::NeumannBoundary)));

  hhg::VTKWriter<hhg::P1Function< real_t >, hhg::DGFunction< real_t >>({ &u->u, &u->v, &u->p, &f->u, &f->v, &f->p }, { }, "../output", "before", maxLevel);

  auto solver = hhg::UzawaSolver<hhg::P1StokesFunction<real_t>, hhg::P1StokesOperator>(storage, minLevel, maxLevel);

  for (uint_t outer = 0; outer < 10; ++outer) {
    solver.solve(L, *u, *f, *r, maxLevel, 1e-15, coarseMaxiter, hhg::Inner | hhg::NeumannBoundary, true);
    hhg::projectMean(u->p, *tmp, maxLevel);


    L.apply(*u, *r, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    r->assign({1.0, -1.0}, { f.get(), r.get() }, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    WALBERLA_LOG_DEVEL("[Uzawa] residuum: " << std::scientific << std::sqrt(r->dot(*r, maxLevel, hhg::Inner | hhg::NeumannBoundary)));
  }

  hhg::VTKWriter<hhg::P1Function< real_t >, hhg::DGFunction< real_t >>({ &u->u, &u->v, &u->p, &f->u, &f->v, &f->p }, { }, "../output", "after", maxLevel);
  return EXIT_SUCCESS;
}
