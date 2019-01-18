#include "core/math/Random.h"
#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/petsc/PETScSparseMatrix.hpp"
#include "tinyhhg_core/FunctionTraits.hpp"


using walberla::real_t;
using walberla::uint_t;

namespace hhg {

static void test( const std::string & meshFile, const uint_t & level )
{
  PETScManager petscManager;

  auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

  P2Function< PetscInt > numerator( "numerator", storage, level, level );
  P2ConstantLaplaceOperator L( storage, level, level );

  const uint_t globalDoFs = numberOfGlobalDoFs< hhg::P2FunctionTag >( *storage, level );
  const uint_t localDoFs  = numberOfLocalDoFs< hhg::P2FunctionTag >( *storage, level );

  numerator.enumerate( level );

  hhg::PETScSparseMatrix< hhg::P2ConstantLaplaceOperator, hhg::P2Function > Lpetsc( localDoFs, globalDoFs );
  Lpetsc.createMatrixFromFunction( L, level, numerator, hhg::All );

  WALBERLA_CHECK( Lpetsc.isSymmetric(), "P2 Laplacian _NOT_ symmetric for: level = " << level << ", mesh: " << meshFile );
  WALBERLA_LOG_INFO_ON_ROOT( "P2 Laplacian symmetric for: level = " << level << ", mesh: " << meshFile );
}

}

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  for ( uint_t level = 2; level <= 3; level++ )
  {
    hhg::test( "../../data/meshes/annulus_coarse.msh",            level );

    hhg::test( "../../data/meshes/3D/tet_1el.msh",                level );
    hhg::test( "../../data/meshes/3D/pyramid_2el.msh",            level );
    hhg::test( "../../data/meshes/3D/pyramid_4el.msh",            level );
    hhg::test( "../../data/meshes/3D/regular_octahedron_8el.msh", level );
    hhg::test( "../../data/meshes/3D/cube_24el.msh",              level );
  }


  return 0;
}
