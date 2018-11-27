
#include "core/math/Random.h"
#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/petsc/PETScSparseMatrix.hpp"
#include "tinyhhg_core/petsc/PETScVector.hpp"
#include "tinyhhg_core/FunctionTraits.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"

#include <numeric>

using walberla::real_t;

namespace hhg {

void p2PetscApplyTest( const uint_t & level, const std::string & meshFile, const DoFType & location, const real_t & eps )
{
  WALBERLA_LOG_INFO_ON_ROOT( "level: " << level << ", mesh file: " << meshFile );

  PETScManager petscManager;

  MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupStorage(meshInfo, walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  loadbalancing::roundRobin( setupStorage );
  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  writeDomainPartitioningVTK( storage, "../../output", "P2PetscApplyTestDomain" );

  P2Function< real_t >   src      ( "src",       storage, level, level );
  P2Function< real_t >   hhgDst   ( "hhgDst",    storage, level, level );
  P2Function< real_t >   petscDst ( "petscDst",  storage, level, level );
  P2Function< real_t >   err      ( "error",     storage, level, level );
  P2Function< real_t >   ones     ( "ones",      storage, level, level );
  P2Function< PetscInt > numerator( "numerator", storage, level, level );

  std::function<real_t(const hhg::Point3D&)> zero  = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> one   = [](const hhg::Point3D&) { return 1.0; };
  std::function<real_t(const hhg::Point3D&)> rand         = []( const hhg::Point3D &   ) { return walberla::math::realRandom<real_t>(); };
  std::function<real_t(const hhg::Point3D&)> srcFunction  = []( const hhg::Point3D & x ) { return x[0] * x[0] * x[0] * x[0] * std::sinh( x[1] ) * std::cos( x[2] ); };

  src.interpolate( srcFunction, level, hhg::All );
  hhgDst.interpolate( rand, level, location );
  petscDst.interpolate( rand, level, location );
  ones.interpolate( one, level, location );

  P2ConstantLaplaceOperator L( storage, level, level );

  numerator.enumerate( level );

  const uint_t globalDoFs = hhg::numberOfGlobalDoFs< hhg::P2FunctionTag >( *storage, level );
  const uint_t localDoFs  = hhg::numberOfLocalDoFs< hhg::P2FunctionTag >( *storage, level );

  WALBERLA_LOG_INFO_ON_ROOT( "Global DoFs: " << globalDoFs );

  // HyTeG apply
  L.apply( src, hhgDst, level, location );

  // PETSc apply
  PETScVector< real_t, P2Function > srcPetscVec( localDoFs );
  PETScVector< real_t, P2Function > dstPetscVec( localDoFs );
  PETScSparseMatrix< P2ConstantLaplaceOperator, P2Function > petscMatrix( localDoFs, globalDoFs );

  srcPetscVec.createVectorFromFunction( src, numerator, level, All );
  dstPetscVec.createVectorFromFunction( petscDst, numerator, level, All );
  petscMatrix.createMatrixFromFunction( L, level, numerator, All );

  WALBERLA_CHECK( petscMatrix.isSymmetric() );

  MatMult( petscMatrix.get(), srcPetscVec.get(), dstPetscVec.get() );

  dstPetscVec.createFunctionFromVector( petscDst, numerator, level, location );

  // compare
  err.assign( {1.0, -1.0}, {&hhgDst, &petscDst}, level, location );
  const auto absScalarProd = std::abs( err.dotGlobal( ones, level, location ) );

  WALBERLA_LOG_INFO_ON_ROOT( "Error sum = " << absScalarProd );

  // VTK
  VTKOutput vtkOutput( "../../output", "P2PetscApplyTest", storage );
  vtkOutput.add( &src );
  vtkOutput.add( &hhgDst );
  vtkOutput.add( &petscDst );
  vtkOutput.add( &err );
  vtkOutput.write( level, 0 );

  WALBERLA_CHECK_LESS( absScalarProd, eps );

}

}

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  hhg::p2PetscApplyTest( 3, "../../data/meshes/quad_4el.msh",       hhg::All,   1.6e-15 );
  hhg::p2PetscApplyTest( 3, "../../data/meshes/annulus_coarse.msh", hhg::All,   1.2e-13 );
  hhg::p2PetscApplyTest( 3, "../../data/meshes/3D/tet_1el.msh",     hhg::Inner, 1.1e-17 );
  hhg::p2PetscApplyTest( 3, "../../data/meshes/3D/pyramid_2el.msh", hhg::Inner, 9.1e-16 );
  hhg::p2PetscApplyTest( 3, "../../data/meshes/3D/pyramid_4el.msh", hhg::Inner, 1.5e-15 );
  hhg::p2PetscApplyTest( 3, "../../data/meshes/3D/regular_octahedron_8el.msh", hhg::Inner, 3.1e-15 );

  return EXIT_SUCCESS;
}
