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

  PETScManager petscManager;

  std::string meshFileName = "../data/meshes/bfs_126el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t Level = 3;

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  hhg::P1Function x("x", storage, Level, Level);
  hhg::P1Function x_exact("x_exact", storage, Level, Level);
  hhg::P1Function numerator("numerator", storage, Level, Level);

  hhg::P1LaplaceOperator A(storage, Level, Level);


  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& xx) { return xx[0]; };
  std::function<real_t(const hhg::Point3D&)> rhs   = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };

  x.interpolate(exact, Level, hhg::DirichletBoundary);
  x_exact.interpolate(exact, Level);

  uint_t num = 0;
  numerator.enumerate(Level,num);
  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("Num dofs = {}", (size_t)num));

  WALBERLA_LOG_INFO_ON_ROOT("Creating Matrix")
  hhg::PETScSparseMatrix<hhg::P1LaplaceOperator,hhg::P1Function> Amat("A",A,Level,numerator,num);
  Amat.applyDirichletBC(numerator,Level);
  WALBERLA_LOG_INFO_ON_ROOT("Applying DBC")
  Amat.print("CreateMatrixDirichlet.m");


  //hhg::PETScVector<hhg::P1Function> x_exact_vector("x_exact",num);
  WALBERLA_LOG_INFO_ON_ROOT("Creating RHS")
  hhg::PETScVector<hhg::P1Function> rhs_vector("b",num);

  //x_exact_vector.createVectorFromFunction(x_exact,numerator,maxLevel);
  rhs_vector.createVectorFromFunction(x_exact,numerator,Level,hhg::DirichletBoundary);
  rhs_vector.print("CreateRHS.m");












  WALBERLA_LOG_INFO_ON_ROOT("Printing Solution")
  hhg::VTKWriter< P1Function >({ &x, &x_exact }, Level, "../output", "exact_solver");


  return 0;
}
