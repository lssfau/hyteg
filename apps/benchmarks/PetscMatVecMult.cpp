
#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Petsc.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/petsc/PETScSparseMatrix.hpp"
#include "tinyhhg_core/petsc/PETScVector.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"

using walberla::real_t;

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;

   PETScManager petscManager;

   std::string meshFileName = "../data/meshes/quad_4el_neumann.msh";

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hhg::loadbalancing::roundRobin( setupStorage );

   size_t level = 11;

   std::function< real_t( const hhg::Point3D& ) > ones  = []( const hhg::Point3D& ) { return 1.0; };
   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& xx ) { return xx[0] + xx[1] + 1; };

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   hhg::P1Function< double > oneFunc( "x", storage, level, level );
   oneFunc.interpolate( ones, level );
   hhg::P1Function< double > x( "x", storage, level, level );
   hhg::P1Function< double > y( "x", storage, level, level );
   hhg::P1Function< double > z( "x", storage, level, level );
   x.interpolate( exact, level );
   hhg::P1Function< PetscInt > numerator( "numerator", storage, level, level );
   hhg::P1MassOperator         mass( storage, level, level );

   uint_t num       = 0;
   uint_t localSize = numerator.enumerate( level, num );

   hhg::PETScSparseMatrix< hhg::P1MassOperator, hhg::P1Function > matPetsc( localSize, num );
   matPetsc.createMatrixFromFunction( mass, level, numerator );
   hhg::PETScVector< real_t, hhg::P1Function > vecPetsc( localSize );
   vecPetsc.createVectorFromFunction( x, numerator, level );
   hhg::PETScVector< real_t, hhg::P1Function > dstvecPetsc( localSize );

   LIKWID_MARKER_START( "HyTeG-apply" );
   mass.apply( x, y, level, hhg::All );
   LIKWID_MARKER_STOP( "HyTeG-apply" );

   //vecPetsc.print("../output/vector0.vec");
   PetscErrorCode ierr;

   LIKWID_MARKER_START( "Petsc-MatMult" );
   ierr = MatMult( matPetsc.get(), vecPetsc.get(), dstvecPetsc.get() );
   LIKWID_MARKER_STOP( "Petsc-MatMult" );
   CHKERRQ( ierr );

   dstvecPetsc.createFunctionFromVector( z, numerator, level );

   WALBERLA_LOG_INFO_ON_ROOT( y.dotGlobal( oneFunc, level ) );
   WALBERLA_LOG_INFO_ON_ROOT( z.dotGlobal( oneFunc, level ) );

   //dstvecPetsc.print("../output/vector1.vec");

   LIKWID_MARKER_CLOSE;
}
