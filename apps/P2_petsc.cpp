#include <tinyhhg_core/tinyhhg.hpp>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  PETScManager petscManager;

  std::string meshFileName = "../data/meshes/tri_1el.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile(meshFileName);
  hhg::SetupPrimitiveStorage setupStorage(meshInfo, walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  hhg::loadbalancing::roundRobin( setupStorage );

  uint_t level = 2;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P2Function<PetscInt> numerator("numerator", storage, level, level);
  hhg::P2Function<real_t> ones("ones", storage, level, level);
  hhg::P2Function<real_t> dst("dst", storage, level, level);

  std::function<real_t(const hhg::Point3D&)> one  = [](const hhg::Point3D&) { return 1.0; };

  ones.interpolate(one, level);
  dst.interpolate(one, level);

  hhg::P2ConstantLaplaceOperator L(storage, level, level);
  L.apply(ones, dst, level, hhg::All, hhg::Replace);

  VTKOutput vtkOutput( "../output", "P2_petsc" );
  vtkOutput.add( &ones );
  vtkOutput.add( &dst );
  vtkOutput.write( level );

  uint_t num = 0;
  uint_t localSize = numerator.enumerate(level, num);

  fmt::print("num = {}\n", num);
  hhg::PETScSparseMatrix<hhg::P2ConstantLaplaceOperator, hhg::P2Function> Lpetsc(localSize, num);
  Lpetsc.createMatrixFromFunction(L, level, numerator, hhg::All);
  Lpetsc.print("../output/matrix.m");

  return 0;
}
