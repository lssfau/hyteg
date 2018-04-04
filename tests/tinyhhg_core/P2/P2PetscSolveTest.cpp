#include "core/timing/Timer.h"

#include <core/Environment.h>
#include "core/math/Random.h"
#include "tinyhhg_core/misc/ExactStencilWeights.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"

#ifndef HHG_BUILD_WITH_PETSC
#error "This test only works with PETSc enabled. Please enable it via -DHHG_BUILD_WITH_PETSC=ON"
#endif

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{
  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  PETScManager petscManager;

  std::string meshFileName = "../../data/meshes/quad_4el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  const size_t level = 2;

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  hhg::P2Function< real_t > x("x", storage, level, level + 1);
  hhg::P2Function< real_t > x_exact("x_exact", storage, level, level + 1);
  hhg::P2Function< real_t > b("x", storage, level, level + 1);
  hhg::P2Function< real_t > err("err", storage, level, level + 1);
  hhg::P2Function< real_t > residuum("err", storage, level, level + 1);
  std::shared_ptr<hhg::P2Function< PetscInt >> numerator = std::make_shared<hhg::P2Function< PetscInt >>("numerator", storage, level, level + 1);

  hhg::P2ConstantLaplaceOperator A(storage, level, level + 1);

  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return sin(xx[0])*sinh(xx[1]); };
  walberla::math::seedRandomGenerator(0);
  std::function<real_t(const Point3D &)> rand = [](const Point3D &) { return walberla::math::realRandom(0.0, 1.0); };

  x.interpolate(exact, level, hhg::DirichletBoundary);
  x.interpolate(rand, level, hhg::Inner);
  b.interpolate(exact, level, hhg::DirichletBoundary);
  x_exact.interpolate(exact, level);

  x.interpolate(exact, level + 1, hhg::DirichletBoundary);
  x.interpolate(rand, level + 1, hhg::Inner);
  b.interpolate(exact, level + 1, hhg::DirichletBoundary);
  x_exact.interpolate(exact, level + 1);

  uint_t num_1 = 0,num_2 = 0;
  uint_t dofsOnRank_1 = numerator->enumerate(level,num_1);
  uint_t dofsOnRank_2 = numerator->enumerate(level + 1,num_2);

  PETScLUSolver<real_t, hhg::P2Function, hhg::P2ConstantLaplaceOperator> solver_1(numerator, dofsOnRank_1, num_1);
  PETScLUSolver<real_t, hhg::P2Function, hhg::P2ConstantLaplaceOperator> solver_2(numerator, dofsOnRank_2, num_2);

  walberla::WcTimer timer;
  solver_1.solve(A,x,b,x,level,0,0);
  solver_2.solve(A,x,b,x,level + 1,0,0);
  timer.end();

  WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
  A.apply(x,residuum,level,hhg::Inner);
  A.apply(x,residuum,level + 1,hhg::Inner);

  err.assign({1.0, -1.0}, {&x, &x_exact}, level);
  err.assign({1.0, -1.0}, {&x, &x_exact}, level + 1);

  real_t discr_l2_err_1 = std::sqrt(err.dot(err, level) / (real_t)num_1);
  real_t discr_l2_err_2 = std::sqrt(err.dot(err, level + 1) / (real_t)num_2);
  real_t residuum_l2_1 = std::sqrt(residuum.dot(residuum, level) / (real_t)num_1);
  real_t residuum_l2_2 = std::sqrt(residuum.dot(residuum, level + 1) / (real_t)num_2);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error 1 = " << discr_l2_err_1);
  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error 2 = " << discr_l2_err_2);
  WALBERLA_LOG_INFO_ON_ROOT("error ratio = " << (discr_l2_err_1 / discr_l2_err_2));
  WALBERLA_LOG_INFO_ON_ROOT("residuum 1 = " << residuum_l2_1);
  WALBERLA_LOG_INFO_ON_ROOT("residuum 1 = " << residuum_l2_2);

  WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(residuum_l2_1,0.0,1e-15);
  WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(residuum_l2_2,0.0,1e-15);
  WALBERLA_CHECK_LESS(8.0, (discr_l2_err_1 / discr_l2_err_2));

  VTKOutput vtkOutput( "../../output", "P2PetscSolve" );
  vtkOutput.add( &x );
  vtkOutput.add( &x_exact );
  vtkOutput.add( &err );
  vtkOutput.add( &residuum );
  vtkOutput.write( level + 1 );


  return EXIT_SUCCESS;
}
