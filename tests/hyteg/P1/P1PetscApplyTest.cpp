
#include "core/math/Random.h"
#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/FunctionTraits.hpp"
#include "hyteg/VTKWriter.hpp"
#include "hyteg/communication/Syncing.hpp"

#include <numeric>

using walberla::real_t;

namespace hyteg {

void p1PetscApplyTest( const uint_t & level, const std::string & meshFile, const DoFType & location, const real_t & eps )
{
  WALBERLA_LOG_INFO_ON_ROOT( "level: " << level << ", mesh file: " << meshFile );

  PETScManager petscManager;

  MeshInfo meshInfo = hyteg::MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupStorage(meshInfo, walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< hyteg::PrimitiveStorage> storage = std::make_shared< hyteg::PrimitiveStorage>(setupStorage);

  writeDomainPartitioningVTK( storage, "../../output", "P1PetscApplyTestDomain" );

  P1Function< real_t >   src      ( "src",       storage, level, level );
  P1Function< real_t >   hhgDst   ( "hhgDst",    storage, level, level );
  P1Function< real_t >   petscDst ( "petscDst",  storage, level, level );
  P1Function< real_t >   err      ( "error",     storage, level, level );
  P1Function< real_t >   ones     ( "ones",      storage, level, level );
  P1Function< PetscInt > numerator( "numerator", storage, level, level );

  std::function<real_t(const hyteg::Point3D&)> zero  = [](const hyteg::Point3D&) { return 0.0; };
  std::function<real_t(const hyteg::Point3D&)> one   = [](const hyteg::Point3D&) { return 1.0; };
  std::function<real_t(const hyteg::Point3D&)> rand         = []( const hyteg::Point3D &   ) { return walberla::math::realRandom<real_t>(); };
  std::function<real_t(const hyteg::Point3D&)> srcFunction  = []( const hyteg::Point3D & x ) { return x[0] * x[0] * x[0] * x[0] * std::sinh( x[1] ) * std::cos( x[2] ); };

  src.interpolate( srcFunction, level, hyteg::All );
  hhgDst.interpolate( rand, level, location );
  petscDst.interpolate( rand, level, location );
  ones.interpolate( one, level, location );

  P1ConstantLaplaceOperator L( storage, level, level );

  numerator.enumerate( level );

  const uint_t globalDoFs = hyteg::numberOfGlobalDoFs< hyteg::P1FunctionTag >( *storage, level );
  const uint_t localDoFs  = hyteg::numberOfLocalDoFs< hyteg::P1FunctionTag >( *storage, level );

  WALBERLA_LOG_INFO_ON_ROOT( "Global DoFs: " << globalDoFs );

  // HyTeG apply
  L.apply( src, hhgDst, level, location );

  // PETSc apply
  PETScVector< real_t, P1Function > srcPetscVec( localDoFs );
  PETScVector< real_t, P1Function > dstPetscVec( localDoFs );
  PETScSparseMatrix< P1ConstantLaplaceOperator, P1Function > petscMatrix( localDoFs, globalDoFs );

  srcPetscVec.createVectorFromFunction( src, numerator, level, All );
  dstPetscVec.createVectorFromFunction( petscDst, numerator, level, All );
  petscMatrix.createMatrixFromFunction( L, level, numerator, All );

  WALBERLA_CHECK( petscMatrix.isSymmetric() );

  MatMult( petscMatrix.get(), srcPetscVec.get(), dstPetscVec.get() );

  dstPetscVec.createFunctionFromVector( petscDst, numerator, level, location );

  // compare
  err.assign( {1.0, -1.0}, {hhgDst, petscDst}, level, location );
  const auto absScalarProd = std::abs( err.dotGlobal( ones, level, location ) );

  WALBERLA_LOG_INFO_ON_ROOT( "Error sum = " << absScalarProd );

  // VTK
  VTKOutput vtkOutput( "../../output", "P1PetscApplyTest", storage );
  vtkOutput.add( src );
  vtkOutput.add( hhgDst );
  vtkOutput.add( petscDst );
  vtkOutput.add( err );
  vtkOutput.write( level, 0 );

  WALBERLA_CHECK_LESS( absScalarProd, eps );

}

}

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  hyteg::p1PetscApplyTest( 3, "../../data/meshes/quad_4el.msh", hyteg::All,   1.0e-15 );
  hyteg::p1PetscApplyTest( 3, "../../data/meshes/annulus_coarse.msh", hyteg::All,   1.7e-13 );
  hyteg::p1PetscApplyTest( 3, "../../data/meshes/3D/tet_1el.msh", hyteg::Inner, 7.6e-18 );
  hyteg::p1PetscApplyTest( 3, "../../data/meshes/3D/pyramid_2el.msh", hyteg::Inner, 2.0e-16 );
  hyteg::p1PetscApplyTest( 3, "../../data/meshes/3D/pyramid_4el.msh", hyteg::Inner, 1.0e-15 );
  hyteg::p1PetscApplyTest( 3, "../../data/meshes/3D/regular_octahedron_8el.msh", hyteg::Inner, 1.0e-15 );


  return EXIT_SUCCESS;
}
