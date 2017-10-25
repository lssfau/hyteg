#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>
#include <core/Environment.h>

#include <boost/core/null_deleter.hpp>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

  timingTree->start("Global");

  std::string meshFileName = "../data/meshes/annulus.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  const uint_t minLevel = 2;
  const uint_t maxLevel = 3;
  const uint_t solverMaxiter = 100;
  const uint_t timesteps = 5000;
  real_t dt = 0.75 * std::pow(2.0, -walberla::real_c(maxLevel));

  std::function<real_t(const hhg::Point3D&,const std::vector<real_t>&)> initialConcentration = [dt](const hhg::Point3D& x,const std::vector<real_t>&) {
    if (sqrt(x[0] * x[0] + x[1] * x[1]) < 1.1){
      return 1.0;
    } else {
      return 0.0;
    }
  };

  std::function<real_t(const hhg::Point3D&,const std::vector<real_t>&)> expr_f_x = [](const hhg::Point3D& x,const std::vector<real_t>& val) {
    return val[0] * std::cos(std::atan2 (x[1], x[0]));
  };

  std::function<real_t(const hhg::Point3D&,const std::vector<real_t>&)> expr_f_y = [](const hhg::Point3D& x,const std::vector<real_t>& val) {
    return val[0] * std::sin(std::atan2 (x[1], x[0]));
  };

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  #ifdef WALBERLA_BUILD_WITH_PARMETIS
    loadbalancing::distributed::parmetis( *storage );
  #endif

  storage->enableGlobalTiming(timingTree);

  // Setting up Functions
  auto c_old = std::make_shared<hhg::DGFunction<real_t>>("c", storage, minLevel, maxLevel);
  auto c = std::make_shared<hhg::DGFunction<real_t>>("c", storage, minLevel, maxLevel);

  auto f_dg_x = std::make_shared<hhg::DGFunction<real_t>>("f_dg_x", storage, minLevel, maxLevel);
  auto f_dg_y = std::make_shared<hhg::DGFunction<real_t>>("f_dg_y", storage, minLevel, maxLevel);

  auto r = std::make_shared<hhg::P1StokesFunction<real_t>>("r", storage, minLevel, maxLevel);
  auto f = std::make_shared<hhg::P1StokesFunction<real_t>>("f", storage, minLevel, maxLevel);
  auto u = std::make_shared<hhg::P1StokesFunction<real_t>>("u", storage, minLevel, maxLevel);

  auto tmp = std::make_shared<hhg::P1Function<real_t>>("tmp", storage, minLevel, maxLevel);

  // Setting up Operators
  std::array<std::shared_ptr<hhg::P1Function<real_t>>, 2> velocity{{std::shared_ptr<hhg::P1Function<real_t>>(&u->u, boost::null_deleter()), std::shared_ptr<hhg::P1Function<real_t>>(&u->v, boost::null_deleter())}};
  hhg::DGUpwindOperator<hhg::P1Function<real_t>> N(storage, velocity, minLevel, maxLevel);
  hhg::P1StokesOperator L(storage, minLevel, maxLevel);
  hhg::P1MassOperator M(storage, minLevel, maxLevel);

  // Interpolate initial functions
  c_old->interpolate(initialConcentration,{}, maxLevel);

  auto solver = hhg::UzawaSolver<hhg::P1StokesFunction<real_t>, hhg::P1StokesOperator>(storage, minLevel, maxLevel);

  for (uint_t t = 0; t <= timesteps; ++t) {
    WALBERLA_LOG_PROGRESS_ON_ROOT("Current timestep: " << t);

    if (t % 3 == 0) {
      WALBERLA_LOG_PROGRESS_ON_ROOT("Solving Stokes system...")

      f_dg_x->interpolate(expr_f_x, { c_old.get() }, maxLevel);
      f_dg_y->interpolate(expr_f_y, { c_old.get() }, maxLevel);

      f->u.integrateDG(*f_dg_x, maxLevel, hhg::All);
      f->v.integrateDG(*f_dg_y, maxLevel, hhg::All);

      for (uint_t outer = 0; outer < 2; ++outer) {
        solver.solve(L, *u, *f, *r, maxLevel, 1e-4, solverMaxiter, hhg::Inner | hhg::NeumannBoundary, true);
        hhg::projectMean(u->p, *tmp, maxLevel);

        L.apply(*u, *r, maxLevel, hhg::Inner | hhg::NeumannBoundary);
        projectMean(u->p, *tmp, maxLevel);

        r->assign({1.0, -1.0}, { f.get(), r.get() }, maxLevel, hhg::Inner | hhg::NeumannBoundary);
        real_t residuum = std::sqrt(r->dot(*r, maxLevel, hhg::Inner | hhg::NeumannBoundary));
        WALBERLA_LOG_PROGRESS_ON_ROOT("[Uzawa] residuum: " << std::scientific << residuum);
      }
    }

    if (t % 25 == 0) {
      timingTree->start("VTK");
      hhg::VTKWriter<hhg::P1Function<real_t>, hhg::DGFunction<real_t >>({&u->u, &u->v, &u->p, &f->u, &f->v}, {c_old.get()}, maxLevel,
                                                                        "../output", fmt::format("plume-{:0>6}", t));
      timingTree->stop("VTK");
    }

    WALBERLA_LOG_PROGRESS_ON_ROOT("Advecting temperature...")
    N.apply(*c_old, *c, maxLevel, hhg::Inner, Replace);
    c->assign({1.0, -dt}, {c_old.get(), c.get()}, maxLevel, hhg::Inner);

    c_old.swap(c);
  }

  timingTree->stop("Global");
  auto reduced_tt = timingTree->getReduced();
  WALBERLA_LOG_INFO_ON_ROOT(reduced_tt);

  return EXIT_SUCCESS;
}
