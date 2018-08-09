
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Petsc.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScSparseMatrix.hpp"
#include "tinyhhg_core/petsc/PETScVector.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;

   PETScManager petscManager;

   std::string meshFileName = "../../../data/meshes/annulus_coarse.msh";

   hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );

   //hhg::MeshInfo::meshUnitSquare( 2 );

   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hhg::loadbalancing::roundRobin( setupStorage );

   auto cfg = std::make_shared< walberla::config::Config >();
   if( env.config() == nullptr )
   {
      auto defaultFile = "./PetscCompare.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   } else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t level = mainConf.getParameter< uint_t >( "level" );

   std::function< real_t( const hhg::Point3D& ) > ones  = []( const hhg::Point3D& ) { return 1.0; };
   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& xx ) {
      //return 5.0;
      return std::sin( walberla::math::PI * xx[0] ) + std::cos( walberla::math::PI * xx[1] );
      //return ( real_c(std::rand()) / real_c(RAND_MAX));
   };

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   storage->setTimingTree( timingTree );

   hhg::P1Function< double > oneFunc( "x", storage, level, level );
   oneFunc.interpolate( ones, level );
   hhg::P1Function< double > x( "x", storage, level, level );
   hhg::P1Function< double > y( "y", storage, level, level );
   hhg::P1Function< double > z( "z", storage, level, level );
   hhg::P1Function< double > diff( "diff", storage, level, level );
   x.interpolate( exact, level, hhg::Inner );
   hhg::P1Function< PetscInt >    numerator( "numerator", storage, level, level );
   hhg::P1ConstantLaplaceOperator mass( storage, level, level );

   uint_t num       = 0;
   uint_t localSize = numerator.enumerate( level, num );

   hhg::PETScSparseMatrix< hhg::P1ConstantLaplaceOperator, hhg::P1Function > matPetsc( localSize, num );
   matPetsc.createMatrixFromFunction( mass, level, numerator, hhg::Inner );
   hhg::PETScVector< real_t, hhg::P1Function > vecPetsc( localSize );
   vecPetsc.createVectorFromFunction( x, numerator, level, hhg::Inner );
   hhg::PETScVector< real_t, hhg::P1Function > dstvecPetsc( localSize );

   LIKWID_MARKER_START( "HyTeG-apply" );
   mass.apply( x, y, level, hhg::Inner );
   LIKWID_MARKER_STOP( "HyTeG-apply" );

   //vecPetsc.print("../output/vector0.vec");
   PetscErrorCode ierr;

   LIKWID_MARKER_START( "Petsc-MatMult" );
   ierr = MatMult( matPetsc.get(), vecPetsc.get(), dstvecPetsc.get() );
   LIKWID_MARKER_STOP( "Petsc-MatMult" );
   CHKERRQ( ierr );

   dstvecPetsc.createFunctionFromVector( z, numerator, level, hhg::Inner );

   WALBERLA_LOG_INFO_ON_ROOT( y.dotGlobal( oneFunc, level, hhg::Inner ) );
   WALBERLA_LOG_INFO_ON_ROOT( z.dotGlobal( oneFunc, level, hhg::Inner ) );

   //dstvecPetsc.print("../output/vector1.vec");

   walberla::WcTimingTree tt  = timingTree->getReduced();
   auto                   tt2 = tt.getCopyWithRemainder();

   //WALBERLA_LOG_INFO_ON_ROOT(tt2);

   diff.assign( {1.0, -1.0}, {&z, &y}, level, hhg::All );

   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT("writing VTK output");
      hhg::VTKOutput vtkOutput( "./output", "petscCompare" );
      vtkOutput.add( &x );
      vtkOutput.add( &z );
      vtkOutput.add( &y );
      vtkOutput.add( &diff );
      vtkOutput.write( level );
   }

   LIKWID_MARKER_CLOSE;
}
