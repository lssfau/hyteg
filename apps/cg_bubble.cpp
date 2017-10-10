#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <fmt/format.h>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG CG Test\n");

  std::string meshFileName = "../data/meshes/tri_1el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t minLevel = 2;
  size_t maxLevel = 2;
  size_t maxiter = 10000;

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  hhg::P1BubbleFunction<real_t> r("r", storage, minLevel, maxLevel);
  hhg::P1BubbleFunction<real_t> f("f", storage, minLevel, maxLevel);
  hhg::P1BubbleFunction<real_t> u("u", storage, minLevel, maxLevel);
  hhg::P1BubbleFunction<real_t> u_exact("u_exact", storage, minLevel, maxLevel);
  hhg::P1BubbleFunction<real_t> err("err", storage, minLevel, maxLevel);
  hhg::P1BubbleFunction<real_t> npoints_helper("npoints_helper", storage, minLevel, maxLevel);

  hhg::P1BubbleLaplaceOperator L(storage, minLevel, maxLevel);

  hhg::BubbleToP1DivTxOperator divtx(storage, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return xx[0]*xx[0] - xx[1]*xx[1]; };
  std::function<real_t(const hhg::Point3D&)> rhs   = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);
  u_exact.interpolate(exact, maxLevel);

  auto solver = hhg::CGSolver<hhg::P1BubbleFunction<real_t>, hhg::P1BubbleLaplaceOperator>(storage, minLevel, maxLevel);
  walberla::WcTimer timer;
  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, true);
  timer.end();
  fmt::printf("time was: %e\n",timer.last());
  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

  hhg::VTKWriter<hhg::P1Function< real_t >>({ &u.p1, &u_exact.p1, &f.p1, &r.p1, &err.p1 }, maxLevel, "../output", "cg_bubble");
  return 0;
}
