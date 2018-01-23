#include <tinyhhg_core/tinyhhg.hpp>

#include <core/math/Random.h>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  PETScManager petscManager;

  std::string meshFileName = "../../data/meshes/quad_16el.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile(meshFileName);
  hhg::SetupPrimitiveStorage setupStorage(meshInfo, walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  hhg::loadbalancing::roundRobin( setupStorage );

  uint_t level = 4;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P2Function<PetscInt> numerator("numerator", storage, level, level);
  hhg::P2Function<real_t> ones("ones", storage, level, level);
  hhg::P2Function<real_t> dst("dst", storage, level, level);

  std::function<real_t(const hhg::Point3D&)> one  = [](const hhg::Point3D&) { return 1.0; };
  std::function<real_t(const hhg::Point3D&)> rand  = [](const hhg::Point3D&) { return walberla::math::realRandom<real_t>(); };

  ones.interpolate(one, level);
  dst.interpolate(rand, level);

  hhg::P2ConstantLaplaceOperator L(storage, level, level);
  L.apply(ones, dst, level, hhg::All, hhg::Replace);

  real_t sqSum = dst.dot(dst, level, hhg::All);

  // Check if row sum is zero
  WALBERLA_CHECK_LESS( sqSum, 1e-14 );

  uint_t num = 0;
  uint_t localSize = numerator.enumerate(level, num);

  hhg::PETScSparseMatrix<hhg::P2ConstantLaplaceOperator, hhg::P2Function> Lpetsc(localSize, num);
  Lpetsc.createMatrixFromFunction(L, level, numerator, hhg::All);

  WALBERLA_CHECK_EQUAL( Lpetsc.isSymmetric(), true );

  return 0;
}
