#include <complex>
#include <fstream>
#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/ExactStencilWeights.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScExportLinearSystem.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScVersion.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

/* Multi sinker benchmark described in section 2.1 of [1]. Specify solver type, 
   number of elements in each spatial direction, number of sinkers, viscosity 
   contrast, decay rate of the sinkers and radius of the sinkers.

   [1]: "WEIGHTED BFBT PRECONDITIONER FOR STOKES FLOW PROBLEMS WITH HIGHLY HETEROGENEOUS VISCOSITY" 
   by JOHANN RUDI, GEORG STADLER, AND OMAR GHATTAS 
*/

void MultiSinker2D( const uint_t& level,
                    const uint_t& solver,
                    const uint_t& nxy,
                    const uint_t& nSinkers,
                    const real_t& visc_min,
                    const real_t& visc_max,
                    const real_t& delta,
                    const real_t& omega )
{
   // storage and domain
   auto meshInfo = MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), MeshInfo::CRISSCROSS, nxy, nxy );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "MultiSinker2DBenchmark_Domain" );

   WALBERLA_LOG_INFO_ON_ROOT( "#Sinkers: " << nSinkers << ", visc_max: " << visc_max );

   // function setup
   hyteg::P2P1TaylorHoodFunction< real_t >          x( "x", storage, level, level );
   hyteg::P2P1TaylorHoodFunction< real_t >          btmp( "btmp", storage, level, level );
   hyteg::P2P1TaylorHoodFunction< real_t >          b( "b", storage, level, level );
   hyteg::P2P1TaylorHoodFunction< real_t >          residuum( "res", storage, level, level );
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return real_c( 0 ); };

   // generate sinker centers randomly (without them exiting the domain)
   std::uniform_real_distribution< real_t > unif( -1 + omega, 1 - omega );
   std::default_random_engine               re;
   re.seed( 1312412415 );
   //re.seed( 213512512 );
   
   std::vector< Point3D > centers;
   for ( uint_t c = 0; c < nSinkers; c++ )
   {
      centers.push_back( Point3D( { unif( re ), unif( re ) } ) );
   }

   // characteristic function for sinkers
   std::function< real_t( const hyteg::Point3D& ) > Xi = [centers, visc_max, visc_min, delta, omega]( const hyteg::Point3D& xx ) {
      real_t val = 1;
      for ( auto& c : centers )
      {
         auto distance = c - xx;
         val *= 1 - exp( -delta * std::pow( std::max( 0.0, distance.norm() - omega / 2 ), 2 ) );
      }
      return val;
   };

   // viscosity function
   std::function< real_t( const hyteg::Point3D& ) > viscosity = [visc_max, visc_min, Xi]( const hyteg::Point3D& xx ) {
      return ( visc_max - visc_min ) * ( 1 - Xi( xx ) ) + visc_min;
   };

   // right hand side: "force sinkers downward": negative v velocity at sinker locations
   std::function< real_t( const hyteg::Point3D& ) > rhsV = [Xi]( const hyteg::Point3D& xx ) { return 10 * ( Xi( xx ) - 1 ); };
   b.uvw().interpolate( { zero, rhsV }, level );

   // operator setup
   hyteg::P2P1ElementwiseAffineEpsilonStokesOperator A( storage, level, level, viscosity );

   // Visualization
   VTKOutput vtkOutput( "../../output", "MultiSinker2DBenchmark", storage );
   vtkOutput.add( x.uvw() );
   vtkOutput.add( x.p() );
   vtkOutput.add( b.uvw() );
   vtkOutput.add( b.p() );
   vtkOutput.write( level, 0 );

   // DoFs
   uint_t localDoFs1  = hyteg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );

   // Initial errors and residual
   x.interpolate( zero, level, hyteg::DirichletBoundary );
   A.apply( x, btmp, level, hyteg::Inner | hyteg::NeumannBoundary );
   residuum.assign( { 1.0, -1.0 }, { b, btmp }, level, hyteg::Inner | hyteg::NeumannBoundary );
   real_t residuum_l2_1 =
       std::sqrt( residuum.dotGlobal( residuum, level, hyteg::Inner | hyteg::NeumannBoundary ) / (real_t) globalDoFs1 );

   WALBERLA_LOG_INFO_ON_ROOT( "initial residual = " << residuum_l2_1 );

   // Solve
   PETScLUSolver< P2P1ElementwiseAffineEpsilonStokesOperator >                        LU( storage, level );
   PETScBlockPreconditionedStokesSolver< P2P1ElementwiseAffineEpsilonStokesOperator > StdBlkdiagPMINRES_PETSC(
       storage, level, 1e-6, std::numeric_limits< PetscInt >::max(), 0 );
   PETScBlockPreconditionedStokesSolver< P2P1ElementwiseAffineEpsilonStokesOperator > GKB(
       storage, level, 1e-6, std::numeric_limits< PetscInt >::max(), 5, 1, 2 );

   auto ViscWeightedPMINRES = solvertemplates::varViscStokesMinResSolver( storage, level, viscosity, 1, 1e-15, 1e-9, 10000, true );
   auto OnlyPressurePMINRES =
       solvertemplates::stokesMinResSolver< P2P1ElementwiseAffineEpsilonStokesOperator >( storage, level, 1e-14, 10000, true );
   auto StdBlkdiagPMINRES = solvertemplates::blkdiagPrecStokesMinResSolver( storage, 2, level, 1e-15, 1e-11, 10000, true );
   auto BFBT_PMINRES     = solvertemplates::BFBTStokesMinResSolver( storage, level, viscosity, 1e-15, 10000, true, x.uvw()[0].getBoundaryCondition() );

   //auto bfbtop = std::make_shared<BFBT_P2P1>( storage, level, level, viscosity );
   //bfbtop->printComponentMatrices(level, storage);

   walberla::WcTimer timer;
   switch ( solver )
   {
   case 0:
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: StdBlkdiagPMINRES_PETSC" );
      StdBlkdiagPMINRES_PETSC.solve( A, x, b, level );
      break;
   case 1:
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: GKB" );
      GKB.solve( A, x, b, level );
      break;
   case 2:
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: LU" );
      LU.solve( A, x, b, level );
      break;
   case 3:
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: ViscWeightedPMINRES" );
      ViscWeightedPMINRES->solve( A, x, b, level );
      break;
   case 4:
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: StdBlkdiagPMINRES" );
      StdBlkdiagPMINRES->solve( A, x, b, level );
      break;
   case 5:
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: BFBT_PMINRES" );
      BFBT_PMINRES->solve( A, x, b, level );
      break;
   case 6:
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: OnlyPressurePMINRES" );
      OnlyPressurePMINRES->solve( A, x, b, level );
      break;

   default:
      WALBERLA_ABORT( "No solver chosen! aborting..." );
   }
   timer.end();

   hyteg::vertexdof::projectMean( x.p(), level );

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   A.apply( x, btmp, level, hyteg::Inner | hyteg::NeumannBoundary );
   residuum.assign( { 1.0, -1.0 }, { b, btmp }, level );
   residuum_l2_1 =
       std::sqrt( residuum.dotGlobal( residuum, level, hyteg::Inner | hyteg::NeumannBoundary ) / (real_t) globalDoFs1 );
   WALBERLA_LOG_INFO_ON_ROOT( "final residual = " << residuum_l2_1 );

   vtkOutput.write( level, 1 );
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );
   /*
  const uint_t & level, 
  const uint_t & solver, 
  const uint_t & nxy,
  const uint_t & nSinkers, 
  const uint_t & visc_min, 
  const uint_t & visc_max, 
  const real_t & delta, 
  const real_t & omega)
  */

   MultiSinker2D( 4, atoi( argv[2] ), 1, atoi( argv[3] ), 1, atoi( argv[4] ), 200, 0.1 );

   return EXIT_SUCCESS;
}
