#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/quad_4el_neumann.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile(meshFileName);
  hhg::SetupPrimitiveStorage setupStorage(meshInfo, walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t level = 4;
  size_t maxiter = 10000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P2P1TaylorHoodFunction<real_t> r("r", storage, level, level);
  hhg::P2P1TaylorHoodFunction<real_t> f("f", storage, level, level);
  hhg::P2P1TaylorHoodFunction<real_t> u("u", storage, level, level);

  hhg::P2P1TaylorHoodStokesOperator L(storage, level, level);

  std::function<real_t(const hhg::Point3D&)> bc_x = [](const hhg::Point3D& x) {
    return 4.0 * (1.0-x[1]) * x[1];
  };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate(zero, level);
  u.u.interpolate(bc_x, level, hhg::DirichletBoundary);
  u.v.interpolate(zero, level, hhg::DirichletBoundary);


  auto solver = hhg::MinResSolver<hhg::P2P1TaylorHoodFunction<real_t>, hhg::P2P1TaylorHoodStokesOperator>(storage, level, level);
  solver.solve(L, u, f, r, level, 1e-12, maxiter, hhg::Inner | hhg::NeumannBoundary, true);

//  PETScManager petscManager;
//  f.u.interpolate(bc_x, level, hhg::DirichletBoundary);
//  auto numerator = std::make_shared<hhg::P2P1TaylorHoodFunction<PetscInt>>("numerator", storage, level, level);
//  uint_t num = 0;
//  uint_t localSize = numerator->enumerate(level, num);
//  PETScLUSolver<real_t, hhg::P2P1TaylorHoodFunction, hhg::P2P1TaylorHoodStokesOperator> solver(numerator, localSize, num);
//  solver.solve(L, u, f, r, level, 1e-12, maxiter, hhg::Inner | hhg::NeumannBoundary, true);

  VTKOutput vtkOutput( "../output", "StokesP2P1" );
  vtkOutput.add( &u.u );
  vtkOutput.add( &u.v );
  vtkOutput.add( &u.p );
  vtkOutput.write( level );

  return EXIT_SUCCESS;
}
