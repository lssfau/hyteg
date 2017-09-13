#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>
#include <core/Environment.h>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

#ifdef HHG_BUILD_WITH_PETSC
  PETScManager petscManager;
#endif

  std::string meshFileName = "../data/meshes/bfs_12el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t minLevel = 2;
  size_t maxLevel = 4;
  size_t maxiter = 10000;

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  hhg::P1Function< real_t > r("r", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > f("f", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > u("u", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > u_exact("u_exact", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > err("err", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > npoints_helper("npoints_helper", storage, minLevel, maxLevel);

  hhg::P1LaplaceOperator L(storage, minLevel, maxLevel);

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  r.enableTiming( timingTree );
  f.enableTiming( timingTree );
  u.enableTiming( timingTree );
  u_exact.enableTiming( timingTree );
  err.enableTiming( timingTree );
  npoints_helper.enableTiming( timingTree );

  L.enableTiming( timingTree );

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return xx[0]; };
  std::function<real_t(const hhg::Point3D&)> rhs   = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u_exact.interpolate(exact, maxLevel);

#ifdef HHG_BUILD_WITH_PETSC
  typedef hhg::PETScPreconditioner<hhg::P1Function< real_t >, hhg::P1LaplaceOperator> PreconditionerType;
  std::shared_ptr< hhg::P1Function< real_t > > numerator = std::make_shared< hhg::P1Function< real_t > >("numerator", storage, maxLevel, maxLevel);
  uint_t globalSize = 0;
  uint_t localSize = numerator->enumerate(maxLevel, globalSize);
  auto prec = std::make_shared<PreconditionerType>(L, numerator, localSize, globalSize);
#else
  typedef hhg::GaussSeidelPreconditioner<hhg::P1Function< real_t >, hhg::P1LaplaceOperator> PreconditionerType;
  auto prec = std::make_shared<PreconditionerType>(L, 30);
#endif
  auto solver = hhg::CGSolver<hhg::P1Function< real_t >, hhg::P1LaplaceOperator, PreconditionerType>(storage, minLevel, maxLevel, prec);
  walberla::WcTimer timer;
  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, true);
  timer.end();
  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("time was: {}",timer.last()));
  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

//  auto face0data = *storage->beginFaces().operator*().second->getData(u.getFaceDataID());

//  hhg::P1Face::printFunctionMemory(face0data,maxLevel);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

  hhg::VTKWriter< P1Function< real_t > >({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "cg_P1");

  walberla::WcTimingTree tt = timingTree->getReduced();
  WALBERLA_LOG_INFO_ON_ROOT( tt );

  return 0;
}
