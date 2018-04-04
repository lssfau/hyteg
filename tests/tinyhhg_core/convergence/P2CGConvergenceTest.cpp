#include <core/timing/Timer.h>

#include <core/Environment.h>
#include <core/config/Config.h>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{
  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  const uint_t level         = 3;
  const std::string meshFile = "../../data/meshes/quad_8el.msh";
  const real_t tolerance     = 1e-15;
  const uint_t maxIter       = 1000;

  auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

  hhg::P2ConstantLaplaceOperator L(storage, level, level);

  hhg::P2Function< real_t > r("r", storage, level, level);
  hhg::P2Function< real_t > f("f", storage, level, level);
  hhg::P2Function< real_t > u("u", storage, level, level);
  hhg::P2Function< real_t > u_exact("u_exact", storage, level, level);
  hhg::P2Function< real_t > err("err", storage, level, level);
  hhg::P2Function< real_t > npoints_helper("npoints_helper", storage, level, level);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  u.interpolate(exact, level, hhg::DirichletBoundary);
  u_exact.interpolate(exact, level);

  auto solver = hhg::CGSolver<hhg::P2Function< real_t >, hhg::P2ConstantLaplaceOperator>(storage, level, level);
  solver.solve(L, u, f, r, level, tolerance, maxIter, hhg::Inner, true);

  err.assign({1.0, -1.0}, {&u, &u_exact}, level);
  npoints_helper.interpolate(ones, level);

  const real_t npoints = npoints_helper.dot(npoints_helper, level);
  const real_t discr_l2_err = std::sqrt(err.dot(err, level) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);
  WALBERLA_CHECK_LESS( discr_l2_err, 2e-8 );
  return 0;
}
