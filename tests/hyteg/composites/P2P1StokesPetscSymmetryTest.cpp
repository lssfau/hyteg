#include "core/math/Random.h"
#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/FunctionTraits.hpp"


using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

static void test( const std::string & meshFile, const uint_t & level )
{
  PETScManager petscManager;

  auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

  P2P1TaylorHoodFunction< PetscInt > numerator( "numerator", storage, level, level );
  P2P1TaylorHoodStokesOperator L( storage, level, level );

  const uint_t globalDoFs = numberOfGlobalDoFs< hyteg::P2P1TaylorHoodFunctionTag >( *storage, level );
  const uint_t localDoFs  = numberOfLocalDoFs< hyteg::P2P1TaylorHoodFunctionTag >( *storage, level );

  numerator.enumerate( level );

  hyteg::PETScSparseMatrix< P2P1TaylorHoodStokesOperator, P2P1TaylorHoodFunction > Lpetsc( localDoFs, globalDoFs );
  Lpetsc.createMatrixFromFunction( L, level, numerator, hyteg::All );

  WALBERLA_CHECK( Lpetsc.isSymmetric(), "P2P1 Stokes operator _NOT_ symmetric for: level = " << level << ", mesh: " << meshFile );
  WALBERLA_LOG_INFO_ON_ROOT( "P2P1 Stokes operator symmetric for: level = " << level << ", mesh: " << meshFile );
}

}

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  for ( uint_t level = 2; level <= 3; level++ )
  {
     hyteg::test( "../../data/meshes/annulus_coarse.msh",            level );

     hyteg::test( "../../data/meshes/3D/tet_1el.msh",                level );
     hyteg::test( "../../data/meshes/3D/pyramid_2el.msh",            level );
     hyteg::test( "../../data/meshes/3D/pyramid_4el.msh",            level );
     hyteg::test( "../../data/meshes/3D/regular_octahedron_8el.msh", level );
     hyteg::test( "../../data/meshes/3D/cube_24el.msh",              level );
  }


  return 0;
}
