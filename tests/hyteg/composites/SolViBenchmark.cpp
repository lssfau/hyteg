#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"
#include <complex>
#include <iostream>
#include <fstream>

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
#include "hyteg/MeshQuality.hpp"

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

std::tuple<real_t, real_t, real_t> SolViBenchmark( const uint_t& level,
                     const uint_t& solver,
                     const uint_t& nxy,
                     const real_t  r_inclusion,
                     const real_t& visc_inclusion,
                     const real_t& visc_matrix )
{
   using namespace std::complex_literals;

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
   hyteg::P2P1ElementwiseAffineEpsilonStokesOperator OP( storage, level, level, viscosity );


   // exact solution
   const real_t                                     C_visc = visc_matrix / ( visc_inclusion + visc_matrix );
   const real_t A = C_visc* (visc_inclusion - visc_matrix);
   std::function<  hyteg::Point3D( const hyteg::Point3D& ) > analytic_uvp = [A,r_inclusion,visc_matrix,visc_inclusion]( const hyteg::Point3D& xx ) {
    
        std::complex<real_t> phi, psi, dphi;
        real_t r2_inclusion = r_inclusion * r_inclusion;

        real_t x = xx[0];
        real_t y = xx[1];
        real_t r2 = x*x+y*y;

        std::complex<real_t> z(x,y);
        if (r2<r2_inclusion)
          {
            //inside the inclusion
            phi=0;
            dphi=0;
            psi=-4.0*(visc_inclusion*visc_matrix/(visc_matrix+visc_inclusion))*z;
          }
        else
          {
            //outside the inclusion
            phi=-2*A*r2_inclusion/z;
            dphi=-phi/z;
            psi=-2.0*(visc_matrix*z+A*r2_inclusion*r2_inclusion/(z*z*z));
          }
        real_t visc = (r2<r2_inclusion)? visc_inclusion : 1.0;
        std::complex<real_t> v = (phi - z*conj(dphi) - conj(psi))/(2.0*visc); 
        return Point3D({v.real(), v.imag(), -2*dphi.real()});
    };
   std::function< real_t( const hyteg::Point3D& ) > exactU = [analytic_uvp]( const hyteg::Point3D& xx ) {
        auto uvp =  analytic_uvp(xx);
        return uvp[0];
   };
   std::function< real_t( const hyteg::Point3D& ) > exactV = [analytic_uvp]( const hyteg::Point3D& xx ) {
        auto uvp =  analytic_uvp(xx);
        return uvp[1];
   };
   std::function< real_t( const hyteg::Point3D& ) > exactP = [analytic_uvp]( const hyteg::Point3D& xx ) {
        auto uvp =  analytic_uvp(xx);
        return uvp[2];
    };
   x_exact.uvw().interpolate( { exactU, exactV }, level );
   x_exact.p().interpolate( exactP, level );
   x.uvw().interpolate( { exactU, exactV }, level, hyteg::DirichletBoundary );
   //x.p().interpolate( { exactP }, level, hyteg::DirichletBoundary );


    // Right-hand-side
    std::function< real_t( const hyteg::Point3D& ) > ddx_u = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
       return A*std::pow(r_inclusion,2.0)*xx[0]*(
                std::pow(r_inclusion,2.0)
                *(
                    12.0*std::pow(xx[0],4.0)-120.0*std::pow(xx[0]*xx[1],2.0)+60.0*std::pow(xx[1],4.0)
                 )
                - 4.0*std::pow(xx[0],6.0) 
                + 52.0*std::pow(xx[0],4.0)*std::pow(xx[1],2.0)
                + 20.0*std::pow(xx[0],2.0)*std::pow(xx[1],4.0) 
                - 36.0*std::pow(xx[1],6.0)
            )/(std::pow(xx[0]*xx[0] + xx[1]*xx[1],5.0));
    };
     std::function< real_t( const hyteg::Point3D& ) > ddy_u = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
            return A*std::pow(r_inclusion,2.0)*xx[0]*(
                std::pow(r_inclusion,2.0)
                *(
                    -12.0*std::pow(xx[0],4.0)+120.0*std::pow(xx[0]*xx[1],2.0)-60.0*std::pow(xx[1],4.0)
                 )
                + 12.0*std::pow(xx[0],6.0) 
                - 60.0*std::pow(xx[0],4.0)*std::pow(xx[1],2.0)
                - 60.0*std::pow(xx[0],2.0)*std::pow(xx[1],4.0) 
                +12.0* std::pow(xx[1],6.0)
            )/(std::pow(xx[0]*xx[0] + xx[1]*xx[1],5.0));
    };
    std::function< real_t( const hyteg::Point3D& ) > ddx_v = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
          return A*std::pow(r_inclusion,2.0)*xx[1]*(
                std::pow(r_inclusion,2.0)
                *(
                   60.0*std::pow(xx[0],4.0) - 120.0*std::pow(xx[0],2.0)*std::pow(xx[1],2.0)+12.0*std::pow(xx[1],4.0)
                 )
                -12.0* std::pow(xx[0],6.0) 
                +60.0*std::pow(xx[0],4.0)*std::pow(xx[1],2.0)
                +60.0*std::pow(xx[0],2.0)*std::pow(xx[1],4.0)
                -12.0* std::pow(xx[1],6.0)
            )/(std::pow(xx[0]*xx[0] + xx[1]*xx[1],5.0));
    };
    std::function< real_t( const hyteg::Point3D& ) > ddy_v = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
       return A*std::pow(r_inclusion,2.0)*xx[1]*(
                std::pow(r_inclusion,2.0)
                *(
                    -60.0*std::pow(xx[0],4.0)+120.0*std::pow(xx[0]*xx[1],2.0)-12.0*std::pow(xx[1],4.0)
                 )
                + 36*std::pow(xx[0],6.0) 
                - 20.0*std::pow(xx[0],4.0)*std::pow(xx[1],2.0)
                - 52.0*std::pow(xx[0],2.0)*std::pow(xx[1],4.0) 
                + 4.0*std::pow(xx[1],6.0)
            )/(std::pow(xx[0]*xx[0] + xx[1]*xx[1],5.0));
    }; 
    std::function< real_t( const hyteg::Point3D& ) > dydx_v = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
          return A*std::pow(r_inclusion,2.0)*xx[0]*(
                std::pow(r_inclusion,2.0)
                *(
                    -12.0*std::pow(xx[0],4.0)+120.0*std::pow(xx[0]*xx[1],2.0)-60.0*std::pow(xx[1],4.0)
                 )
                +4.0*std::pow(xx[0],6.0) 
                -52.0*std::pow(xx[0],4.0)*std::pow(xx[1],2.0)
                -20.0*std::pow(xx[0],2.0)*std::pow(xx[1],4.0) 
                +36.0* std::pow(xx[1],6.0)
            )/(std::pow(xx[0]*xx[0] + xx[1]*xx[1],5.0));
    }; 
  std::function< real_t( const hyteg::Point3D& ) > dxdy_u = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
          return A*std::pow(r_inclusion,2.0)*xx[1]*(
                std::pow(r_inclusion,2.0)
                *(
                    60.0*std::pow(xx[0],4.0)-120.0*std::pow(xx[0]*xx[1],2.0)+12.0*std::pow(xx[1],4.0)
                 )
                - 36.0*std::pow(xx[0],6.0) 
                + 20.0*std::pow(xx[0],4.0)*std::pow(xx[1],2.0)
                + 52.0*std::pow(xx[0],2.0)*std::pow(xx[1],4.0) 
                - 4.0*std::pow(xx[1],6.0)
            )/(std::pow(xx[0]*xx[0] + xx[1]*xx[1],5.0));
    }; 
     std::function< real_t( const hyteg::Point3D& ) > ddx_p = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
          return (12.0*A*std::pow(r_inclusion,2.0)*(std::pow(xx[0],4.0) - 6.0*std::pow(xx[0],2.0)*std::pow(xx[1],2.0)+std::pow(xx[1],4.0))
            )/(std::pow(xx[0]*xx[0] + xx[1]*xx[1],4.0));
    }; 
       std::function< real_t( const hyteg::Point3D& ) > ddy_p = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
          return (-12.0*A*std::pow(r_inclusion,2.0)*(std::pow(xx[0],4.0) - 6.0*std::pow(xx[0],2.0)*std::pow(xx[1],2.0)+std::pow(xx[1],4.0))
            )/(std::pow(xx[0]*xx[0] + xx[1]*xx[1],4.0));
    }; 
   std::function< real_t( const hyteg::Point3D& ) > rhsU =
    [r_inclusion, A, rad, ddx_u, ddy_u, dydx_v, ddx_p, visc_matrix]( const hyteg::Point3D& xx ) {
        if ( rad( xx ) < r_inclusion )
            return 0.0;
        else 
            return  -(2*ddx_u(xx) + ddy_u(xx) + dydx_v(xx)) + ddx_p(xx);
    
    };
    std::function< real_t( const hyteg::Point3D& ) > rhsV =
    [r_inclusion, A, rad, ddx_v, ddy_v,dxdy_u, ddy_p, visc_matrix]( const hyteg::Point3D& xx ) {
        if (   rad( xx ) < r_inclusion )
            return 0.0;
        else 
            return  -(ddx_v(xx) + 2*ddy_v(xx) + dxdy_u(xx)) + ddy_p(xx);

    };
    btmp.uvw().interpolate({rhsU, rhsV},level, hyteg::Inner);
    /*P2ConstantMassOperator VelMassOp(storage, level, level);
    VelMassOp.apply(btmp.uvw()[0], b.uvw()[0], level,All);
    VelMassOp.apply(btmp.uvw()[1], b.uvw()[1], level,All);  
   */// b.uvw().interpolate( { exactU, exactV}, level, DirichletBoundary );


    /*
    hyteg::P2Function< real_t > mask( "mask", storage, level, level );
         std::function< real_t( const hyteg::Point3D& ) > mask_func =
       [r_inclusion, visc_inclusion, visc_matrix, rad]( const hyteg::Point3D& xx ) {
          if ( rad( xx ) <= r_inclusion )
             return 0.0;
          else
             return 1.0;
       };
    mask.interpolate( {mask_func}, level );
    b.uvw()[0].multElementwise({mask,b.uvw()[0]}, level, All);
    b.uvw()[1].multElementwise({mask, b.uvw()[1]}, level, All);
*/
    
 // nullspace.p().interpolate( ones, level, All );
    //b.uvw().interpolate({rhsU, rhsV},level, hyteg::Inner);
   //A.apply( x_exact, b, level, hyteg::Inner | hyteg::NeumannBoundary );



   // Visualization
   /*
   VTKOutput vtkOutput( "../../output", "SolViBenchmark", storage );
   vtkOutput.add( x.uvw() );
   vtkOutput.add( x.p() );
   vtkOutput.add( x_exact.uvw() );
   vtkOutput.add( x_exact.p() );
   vtkOutput.add( err.uvw() );
   vtkOutput.add( err.p() );
   vtkOutput.add( b.uvw() );
   vtkOutput.add( b.p() );
   //vtkOutput.add( mask);   
   vtkOutput.write( level, 0 );
   */

   // DoFs
   uint_t localDoFs1  = hyteg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );


   // Initial errors and residual
   OP.apply( x, btmp, level, hyteg::Inner | hyteg::NeumannBoundary );
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
      BlockprecMINRES.solve( OP, x, b, level );
      break;
   case 1:
      GKB.solve( OP, x, b, level );
      break;
      //   case 2: GMG.solve( A, x, b, level ); break;
   default:
      LU.solve( OP, x, b, level );
   }
   timer.end();

   hyteg::vertexdof::projectMean( x.p(), level );
   hyteg::vertexdof::projectMean( x_exact.p(), level );

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   OP.apply( x, btmp, level, hyteg::Inner  );
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

   //vtkOutput.write( level, 1 );
   
   
   real_t h = MeshQuality::getMaximalEdgeLength(storage, level);
   real_t err_vel = std::sqrt( err.uvw().dotGlobal( err.uvw(), level, hyteg::Inner)/ (real_t) globalDoFs1 );
   return std::make_tuple<real_t&, real_t&, real_t&>(h, err_vel, discr_l2_err_p);
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

    // compute errors
    std::list<uint_t> levels = {2,3,4,5,6,7};
    std::list<std::tuple<real_t, real_t, real_t>> errors_per_h;
   for(auto lvl : levels) {
        auto error_per_h = SolViBenchmark( lvl, 2, 1, 0.2, 1000.0, 1.0 );
        errors_per_h.push_back(error_per_h);
   }

    
    // write to plot file
    
    std::ofstream err_file;
    err_file.open ("err_file.txt");
    for(auto err : errors_per_h) {
        err_file << std::get<0>(err) << ", " << std::get<1>(err) << ", " << std::get<2>(err) << "\n";
    }
    err_file.close();
    
    return EXIT_SUCCESS;
}
