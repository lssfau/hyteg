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

  std::string meshFileName = "../data/meshes/tri_1el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  const uint_t minLevel = 2;
  const uint_t maxLevel = 8;
  const uint_t maxPolyDegree = 2;
  const uint_t interpolationLevel = 6;
  const uint_t maxiter = 10000;

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  hhg::P1Function< real_t > r("r", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > f("f", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > u("u", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > u_exact("u_exact", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > err("err", storage, minLevel, maxLevel);
  hhg::P1Function< real_t > npoints_helper("npoints_helper", storage, minLevel, maxLevel);
  std::shared_ptr<P1Function< real_t >> coefficient = std::make_shared<P1Function< real_t >>("coeff", storage, minLevel, maxLevel);

  hhg::P1MassOperator M(storage, minLevel, maxLevel);

  typedef hhg::P1PolynomialLaplaceOperator<maxPolyDegree, interpolationLevel> SolveOperator;
//  typedef hhg::P1VariableCoefficientLaplaceOperator SolveOperator;

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  r.enableTiming( timingTree );
  f.enableTiming( timingTree );
  u.enableTiming( timingTree );
  u_exact.enableTiming( timingTree );
  err.enableTiming( timingTree );
  npoints_helper.enableTiming( timingTree );

  std::function<real_t(const hhg::Point3D&)> coeff = [](const hhg::Point3D& x) { return 7*pow(x[0], 2) + 3*x[0] + 4*x[1] + 1; };
  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D& x) { return -(14*x[0] + 3)*cos(x[0])*sinh(x[1]) - 4*sin(x[0])*cosh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, maxLevel, hhg::DirichletBoundary);

  for (uint_t level = minLevel; level <= maxLevel; ++level) {
    coefficient->interpolate(coeff, level);
  }

  u_exact.interpolate(exact, maxLevel);
  npoints_helper.interpolate(rhs, maxLevel);
  M.apply(npoints_helper, f, maxLevel, hhg::All);

  SolveOperator L(storage, coefficient, minLevel, maxLevel);
  L.enableTiming( timingTree );

//  typedef hhg::GaussSeidelPreconditioner<hhg::P1Function< real_t >, hhg::P1LaplaceOperator> PreconditionerType;
//  auto prec = std::make_shared<PreconditionerType>(L, 30);

  auto solver = hhg::CGSolver<hhg::P1Function< real_t >, SolveOperator>(storage, minLevel, maxLevel);
  walberla::WcTimer timer;
  solver.solve(L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner, false);
  timer.end();
  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("time was: {}",timer.last()));
  err.assign({1.0, -1.0}, {&u, &u_exact}, maxLevel);

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dot(npoints_helper, maxLevel);

  real_t discr_l2_err = std::sqrt(err.dot(err, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

  //hhg::VTKWriter< P1Function< real_t > >({ &u, &u_exact, &f, &r, &err, coefficient.get() }, maxLevel, "../output", "cg_P1_varcoeff");

  walberla::WcTimingTree tt = timingTree->getReduced();
  WALBERLA_LOG_INFO_ON_ROOT( tt );

  return 0;
}
