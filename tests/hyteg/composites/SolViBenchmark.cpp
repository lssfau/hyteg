#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

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
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/ExactStencilWeights.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVersion.hpp"
#include "hyteg/petsc/PETScExportLinearSystem.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

// SolVi benchmark (Circular inclusion)
/*  Described in section 5.3 of [1] with analytical solution from [2] (equations 26, 34).
The viscosity contains a jump from visc_matrix to visc_inclusion in a circle located at the center of the domain. 
Difficulty for FE methods: the mesh can not align with the viscosity jump.

[1]: "On the choice of finite element for applications in geodynamics" 
by Cedric Thieulot and Wolfgang Bangerth

[2]: "Analytical solutions for deformable elliptical inclusions in general shear" 
by Daniel W. Schmid and Yuri Yu. Podladchikov
*/

void SolViBenchmark( const uint_t& level,
                     const uint_t& solver,
                     const uint_t& nxy,
                     const real_t  r_inclusion,
                     const real_t& visc_inclusion,
                     const real_t& visc_matrix )
{
   // storage and domain
   auto meshInfo = MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), MeshInfo::CRISSCROSS, nxy, nxy );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "SolViBenchmark_Domain" );

   // function setup
   hyteg::P2P1TaylorHoodFunction< real_t > x( "x", storage, level, level );
   hyteg::P2P1TaylorHoodFunction< real_t > x_exact( "x_exact", storage, level, level );
   hyteg::P2P1TaylorHoodFunction< real_t > btmp( "btmp", storage, level, level );
   hyteg::P2P1TaylorHoodFunction< real_t > b( "b", storage, level, level );
   hyteg::P2P1TaylorHoodFunction< real_t > err( "err", storage, level, level );
   hyteg::P2P1TaylorHoodFunction< real_t > residuum( "res", storage, level, level );
  std::function< real_t( const hyteg::Point3D& ) > zero =   []( const hyteg::Point3D&    ) { return real_c(0); };
  std::function< real_t( const hyteg::Point3D& ) > ones =   []( const hyteg::Point3D&    ) { return real_c(1); };

   // radius helper function
   std::function< real_t( const hyteg::Point3D& ) > rad = []( const hyteg::Point3D& xx ) {
      return  sqrt(std::pow( xx[0], 2.0 ) + std::pow( xx[1], 2.0 ) );
   };
   
   // viscosity function and operator setup
   std::function< real_t( const hyteg::Point3D& ) > viscosity =
       [r_inclusion, visc_inclusion, visc_matrix, rad]( const hyteg::Point3D& xx ) {
          if ( rad( xx ) < r_inclusion )
             return visc_inclusion;
          else
             return visc_matrix;
       };

  
   hyteg::P2P1ElementwiseAffineEpsilonStokesOperator A( storage, level, level, viscosity );
   // exact solution
   const real_t                                     C_visc = visc_matrix / ( visc_inclusion + visc_matrix );
   std::function< real_t( const hyteg::Point3D& ) > exactU = [C_visc, r_inclusion, rad]( const hyteg::Point3D& xx ) {
        if ( rad( xx ) < r_inclusion )
            return real_c( 2 ) * C_visc * xx[0];
        else
            return real_c( 1 ) * xx[0];
   };
   std::function< real_t( const hyteg::Point3D& ) > exactV = [C_visc, r_inclusion, rad]( const hyteg::Point3D& xx ) {
        if ( rad( xx ) < r_inclusion )
            return real_c( -2 ) * C_visc * xx[1];
        else
            return real_c( -1 )  * xx[1];
   };

   const real_t C_p = real_c(-4) * C_visc * (visc_inclusion - visc_matrix) * std::pow(r_inclusion,2.0);
   std::function< real_t( const hyteg::Point3D& ) > exactP =
       [C_p, r_inclusion, rad]( const hyteg::Point3D& xx ) {
          if ( rad( xx ) < r_inclusion )
             return real_c( 0 );
          else
             return  C_p * cos( 2.0 * atan( xx[1]/ xx[0] ) ) /
                    ( std::pow( xx[0], 2.0 ) + std::pow( xx[1], 2.0 ) ) ;
       };
   x_exact.uvw().interpolate( { exactU, exactV }, level );
   x_exact.p().interpolate( exactP, level );
   x.uvw().interpolate( { exactU, exactV }, level, hyteg::DirichletBoundary );
   //x.p().interpolate( { exactP }, level, hyteg::DirichletBoundary );

   // Right-hand-side
   std::function< real_t( const hyteg::Point3D& ) > rhsU =
    [r_inclusion, C_p, rad]( const hyteg::Point3D& xx ) {
        if ( rad( xx ) < r_inclusion )
            return real_c( 0 );
        else
            return real_c(-2) * C_p * (xx[0]*cos(2.0*atan(xx[1]/xx[0])) - xx[1]*sin(2.0*atan(xx[1]/xx[0])))/std::pow(xx[0]*xx[0] + xx[1]*xx[1],2.0);
    };

    std::function< real_t( const hyteg::Point3D& ) > rhsV =
    [r_inclusion, C_p, rad]( const hyteg::Point3D& xx ) {
        if ( rad( xx ) < r_inclusion )
            return real_c( 0 );
        else
            return real_c(-2) * C_p * (xx[0]*sin(2.0*atan(xx[1]/xx[0])) + xx[1]*cos(2.0*atan(xx[1]/xx[0])))/std::pow(xx[0]*xx[0] + xx[1]*xx[1],2.0);
    };
    
    btmp.uvw().interpolate({rhsU, rhsV},level, hyteg::Inner);
    P2ConstantMassOperator VelMassOp(storage, level, level);
    VelMassOp.apply(btmp.uvw()[0], b.uvw()[0], level,All);
    VelMassOp.apply(btmp.uvw()[1], b.uvw()[1], level,All);  
    b.uvw().interpolate( { exactU, exactV}, level, DirichletBoundary );
    b.p().interpolate( zero, level, All );

    /*
     hyteg::P2P1TaylorHoodFunction< real_t > mask( "mask", storage, level, level );
         std::function< real_t( const hyteg::Point3D& ) > mask_func =
       [r_inclusion, visc_inclusion, visc_matrix, rad]( const hyteg::Point3D& xx ) {
          if ( rad( xx ) < r_inclusion )
             return 0.0;
          else
             return 1.0;
       };
    mask.interpolate( mask_func, level );
    b.multElementwise({mask}, level, All);
*/
    
 // nullspace.p().interpolate( ones, level, All );
    //b.uvw().interpolate({rhsU, rhsV},level, hyteg::Inner);
   //A.apply( x_exact, b, level, hyteg::Inner | hyteg::NeumannBoundary );



   // Visualization
   VTKOutput vtkOutput( "../../output", "SolViBenchmark", storage );
   vtkOutput.add( x.uvw() );
   vtkOutput.add( x.p() );
   vtkOutput.add( x_exact.uvw() );
   vtkOutput.add( x_exact.p() );
   vtkOutput.add( err.uvw() );
   vtkOutput.add( err.p() );
   vtkOutput.add( b.uvw() );
   vtkOutput.add( b.p() );
//   vtkOutput.add( mask.uvw());
   
   vtkOutput.write( level, 0 );

   // DoFs
   uint_t localDoFs1  = hyteg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );

   // Initial errors and residual
   //A.apply( x, btmp, level, hyteg::Inner | hyteg::NeumannBoundary );
   residuum.assign( { 1.0, -1.0 }, { b, btmp }, level, hyteg::Inner  );
   err.assign( { 1.0, -1.0 }, { x, x_exact }, level, hyteg::Inner  );
   real_t discr_l2_err_u =
       std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level, hyteg::Inner ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_v =
       std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level, hyteg::Inner) / (real_t) globalDoFs1 );
   real_t discr_l2_err_p =
       std::sqrt( err.p().dotGlobal( err.p(), level, hyteg::Inner ) / (real_t) globalDoFs1 );
   real_t residuum_l2 =
       std::sqrt( residuum.dotGlobal( residuum, level, hyteg::Inner  ) / (real_t) globalDoFs1 );
   WALBERLA_LOG_INFO_ON_ROOT( "initial errors and residual:" );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residual = " << residuum_l2 );

   // Solve
   PETScLUSolver< P2P1ElementwiseAffineEpsilonStokesOperator >                        LU( storage, level );
   PETScBlockPreconditionedStokesSolver< P2P1ElementwiseAffineEpsilonStokesOperator > BlockprecMINRES(
       storage, level, 1e-12, std::numeric_limits< PetscInt >::max(), 0 );
   PETScBlockPreconditionedStokesSolver< P2P1ElementwiseAffineEpsilonStokesOperator > GKB(
       storage, level, 1e-12, std::numeric_limits< PetscInt >::max(), 5, 1, 2 );


   walberla::WcTimer timer;
   switch ( solver )
   {
   case 0:
      BlockprecMINRES.solve( A, x, b, level );
      break;
   case 1:
      GKB.solve( A, x, b, level );
      break;
      //   case 2: GMG.solve( A, x, b, level ); break;
   default:
      LU.solve( A, x, b, level );
   }
   timer.end();

   hyteg::vertexdof::projectMean( x.p(), level );
   hyteg::vertexdof::projectMean( x_exact.p(), level );

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   A.apply( x, btmp, level, hyteg::Inner  );
   residuum.assign( { 1.0, -1.0 }, { b, btmp }, level, hyteg::Inner );
   err.assign( { 1.0, -1.0 }, { x, x_exact }, level, hyteg::Inner  );
   discr_l2_err_u =
       std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level, hyteg::Inner  )/ (real_t) globalDoFs1  );
   discr_l2_err_v =
       std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level, hyteg::Inner)/ (real_t) globalDoFs1 );
   discr_l2_err_p =
       std::sqrt( err.p().dotGlobal( err.p(), level, hyteg::Inner ) / (real_t) globalDoFs1 );
   residuum_l2 =
       std::sqrt( residuum.dotGlobal( residuum, level, hyteg::Inner )/ (real_t) globalDoFs1 );
   WALBERLA_LOG_INFO_ON_ROOT( "final errors and residual:" );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residual = " << residuum_l2 );

   vtkOutput.write( level, 1 );
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   // const uint_t & level, const uint_t & solver, const uint_t & nxy, const real_t r_inclusion, const real_t & visc_inclusion, const real_t & visc_matrix

   SolViBenchmark( 7, 2, 1, 0.2, 1000.0, 1.0 );

   return EXIT_SUCCESS;
}
