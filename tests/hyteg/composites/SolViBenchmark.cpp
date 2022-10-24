
#include <complex>
#include <fstream>
#include <iostream>
#include <type_traits>
#include <typeinfo>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
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
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"
#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0EpsilonStokesOperator;

using hyteg::P2P1ElementwiseAffineEpsilonStokesOperator;
using hyteg::Point3D;
using hyteg::dg::eg::EGLaplaceOperator;
using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0ConstEpsilonStokesOperator;
using hyteg::dg::eg::EGP0EpsilonStokesOperator;
using hyteg::dg::eg::EGP0StokesOperator;
namespace hyteg {

// SolVi benchmark (Circular inclusion)
/*  Described in section 5.3 of [1] with analytical solution from [2] (equations 26, 34).
The viscosity contains a jump from visc_matrix to visc_inclusion in a circle located at the center of the domain. 
Difficulty for FE methods: the mesh can not align with the viscosity jump well.

[1]: "On the choice of finite element for applications in geodynamics" 
by Cedric Thieulot and Wolfgang Bangerth

[2]: "Analytical solutions for deformable elliptical inclusions in general shear" 
by Daniel W. Schmid and Yuri Yu. Podladchikov
*/

auto copyBdry = []( EGP0StokesFunction< real_t > fun ) { fun.p().setBoundaryCondition( fun.uvw().getBoundaryCondition() ); };

template < typename StokesOperatorType >
constexpr bool isEGP0Discr()
{
   return std::is_same< StokesOperatorType, EGP0EpsilonStokesOperator >::value ||
          std::is_same< StokesOperatorType, EGP0StokesOperator >::value ||
          std::is_same< StokesOperatorType, EGP0ConstEpsilonStokesOperator >::value;
}

template < typename StokesOperatorType >
constexpr bool isP2P1Discr()
{
   return std::is_same< StokesOperatorType, hyteg::P2P1ElementwiseAffineEpsilonStokesOperator >::value ||
          std::is_same< StokesOperatorType, hyteg::P2P1TaylorHoodStokesOperator >::value;
}

template < typename StokesOperatorType >
std::tuple< real_t, real_t, real_t > RunSolVi( const std::string& name,
                                               const uint_t&      level,
                                               const uint_t&      nxy,
                                               const real_t       r_inclusion,
                                               const real_t&      visc_inclusion,
                                               const real_t&      visc_matrix,
                                               bool               writeVTK = true )
{
   using namespace std::complex_literals;
   using StokesFunctionType          = typename StokesOperatorType::srcType;
   using StokesFunctionNumeratorType = typename StokesFunctionType::template FunctionType< idx_t >;

   // storage and domain
   auto meshInfo = MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), MeshInfo::CRISSCROSS, nxy, nxy );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );
   //writeDomainPartitioningVTK( storage, "../../output", "SolViBenchmark_Domain" );

   // function setup
   StokesFunctionType   x( "x", storage, level, level );
   StokesFunctionType   x_exact( "x_exact", storage, level, level );
   StokesFunctionType   btmp( "btmp", storage, level, level );
   StokesFunctionType   b( "b", storage, level, level );
   StokesFunctionType   err( "err", storage, level, level );
   StokesFunctionType   Merr( "err", storage, level, level );
   P2Function< real_t > Visc_func( "visc_func", storage, level, level );
   //StokesFunctionType nullspace( "nullspace", storage, level, level );
   if constexpr ( isEGP0Discr< StokesOperatorType >() )
   {
      copyBdry( x );
      copyBdry( x_exact );
      copyBdry( err );
      copyBdry( Merr );
      copyBdry( b );
      copyBdry( btmp );
   }
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return real_c( 0 ); };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return real_c( 1 ); };

   // radius helper function
   std::function< real_t( const hyteg::Point3D& ) > rad = []( const hyteg::Point3D& xx ) {
      return sqrt( std::pow( xx[0], 2.0 ) + std::pow( xx[1], 2.0 ) );
   };

   // viscosity function and operator setup
   std::function< real_t( const hyteg::Point3D& ) > viscosity =
       [r_inclusion, visc_inclusion, visc_matrix, rad]( const hyteg::Point3D& xx ) {
          if ( rad( xx ) < r_inclusion )
             return visc_inclusion;
          else
             return visc_matrix;
       };
   StokesOperatorType Op( storage, level, level, viscosity );

   // analytic solution for u,v,p
   const real_t                                             C_visc = visc_matrix / ( visc_inclusion + visc_matrix );
   const real_t                                             A      = C_visc * ( visc_inclusion - visc_matrix );
   std::function< hyteg::Point3D( const hyteg::Point3D& ) > analytic_uvp =
       [A, r_inclusion, visc_matrix, visc_inclusion]( const hyteg::Point3D& xx ) {
          std::complex< real_t > phi, psi, dphi;
          real_t                 r2_inclusion = r_inclusion * r_inclusion;
          real_t                 x            = xx[0];
          real_t                 y            = xx[1];
          real_t                 r2           = x * x + y * y;

          std::complex< real_t > z( x, y );
          if ( r2 < r2_inclusion )
          {
             //inside the inclusion
             phi  = 0;
             dphi = 0;
             psi  = -4.0 * ( visc_inclusion * visc_matrix / ( visc_matrix + visc_inclusion ) ) * z;
          }
          else
          {
             //outside the inclusion
             phi  = -2 * A * r2_inclusion / z;
             dphi = -phi / z;
             psi  = -2.0 * ( visc_matrix * z + A * r2_inclusion * r2_inclusion / ( z * z * z ) );
          }
          real_t                 visc = ( r2 < r2_inclusion ) ? visc_inclusion : 1.0;
          std::complex< real_t > v    = ( phi - z * conj( dphi ) - conj( psi ) ) / ( 2.0 * visc );
          return Point3D( { v.real(), v.imag(), -2 * dphi.real() } );
       };
   std::function< real_t( const hyteg::Point3D& ) > analyticU = [analytic_uvp]( const hyteg::Point3D& xx ) {
      auto uvp = analytic_uvp( xx );
      return uvp[0];
   };
   std::function< real_t( const hyteg::Point3D& ) > analyticV = [analytic_uvp]( const hyteg::Point3D& xx ) {
      auto uvp = analytic_uvp( xx );
      return uvp[1];
   };
   std::function< real_t( const hyteg::Point3D& ) > analyticP = [analytic_uvp]( const hyteg::Point3D& xx ) {
      auto uvp = analytic_uvp( xx );
      return uvp[2];
   };
   //hyteg::vertexdof::projectMean( x_exact.p(), level );
   //nullspace.p().interpolate( ones, level, All );
   //x.p().interpolate( { analyticP }, level, hyteg::DirichletBoundary );

   // Right-hand-side: derivatives of u, v, p for x and y
   std::function< real_t( const hyteg::Point3D& ) > ddx_u = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[0] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( 12.0 * std::pow( xx[0], 4.0 ) - 120.0 * std::pow( xx[0] * xx[1], 2.0 ) + 60.0 * std::pow( xx[1], 4.0 ) ) -
               4.0 * std::pow( xx[0], 6.0 ) + 52.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) +
               20.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) - 36.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > ddy_u = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[0] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( -12.0 * std::pow( xx[0], 4.0 ) + 120.0 * std::pow( xx[0] * xx[1], 2.0 ) - 60.0 * std::pow( xx[1], 4.0 ) ) +
               12.0 * std::pow( xx[0], 6.0 ) - 60.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) -
               60.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) + 12.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > ddx_v = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[1] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( 60.0 * std::pow( xx[0], 4.0 ) - 120.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 2.0 ) +
                     12.0 * std::pow( xx[1], 4.0 ) ) -
               12.0 * std::pow( xx[0], 6.0 ) + 60.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) +
               60.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) - 12.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > ddy_v = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[1] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( -60.0 * std::pow( xx[0], 4.0 ) + 120.0 * std::pow( xx[0] * xx[1], 2.0 ) - 12.0 * std::pow( xx[1], 4.0 ) ) +
               36.0 * std::pow( xx[0], 6.0 ) - 20.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) -
               52.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) + 4.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > dydx_v = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[0] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( -12.0 * std::pow( xx[0], 4.0 ) + 120.0 * std::pow( xx[0] * xx[1], 2.0 ) - 60.0 * std::pow( xx[1], 4.0 ) ) +
               4.0 * std::pow( xx[0], 6.0 ) - 52.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) -
               20.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) + 36.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > dxdy_u = [r_inclusion, A, visc_matrix]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[1] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( 60.0 * std::pow( xx[0], 4.0 ) - 120.0 * std::pow( xx[0] * xx[1], 2.0 ) + 12.0 * std::pow( xx[1], 4.0 ) ) -
               36.0 * std::pow( xx[0], 6.0 ) + 20.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) +
               52.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) - 4.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };

   // right hand side: setup by epsilon operator on u,v and gradient of p
   std::function< real_t( const hyteg::Point3D& ) > rhsU =
       [r_inclusion, A, rad, ddx_u, ddy_u, dydx_v, visc_matrix]( const hyteg::Point3D& xx ) {
          if ( rad( xx ) < r_inclusion )
             return 0.0;
          else
             return -0.5 * ( 2 * ddx_u( xx ) + ddy_u( xx ) + dydx_v( xx ) ) +
                    real_c( 2 ) * A * std::pow( r_inclusion, 2.0 ) *
                        ( xx[0] * cos( 2.0 * atan( xx[1] / xx[0] ) ) - xx[1] * sin( 2.0 * atan( xx[1] / xx[0] ) ) ) /
                        std::pow( xx[0] * xx[0] + xx[1] * xx[1], 2.0 ); //+ ddx_p(xx);
       };
   std::function< real_t( const hyteg::Point3D& ) > rhsV =
       [r_inclusion, A, rad, ddx_v, ddy_v, dxdy_u, visc_matrix]( const hyteg::Point3D& xx ) {
          if ( rad( xx ) < r_inclusion )
             return 0.0;
          else
             return -0.5 * ( ddx_v( xx ) + 2 * ddy_v( xx ) + dxdy_u( xx ) ) +
                    real_c( 2 ) * A * std::pow( r_inclusion, 2.0 ) *
                        ( xx[0] * sin( 2.0 * atan( xx[1] / xx[0] ) ) + xx[1] * cos( 2.0 * atan( xx[1] / xx[0] ) ) ) /
                        std::pow( xx[0] * xx[0] + xx[1] * xx[1], 2.0 ); //+ ddy_p(xx);
       };

   // interpolate solution, "Integrate" rhs
   // btmp.uvw().interpolate( { rhsU, rhsV }, level, Inner );
   x_exact.uvw().interpolate( { analyticU, analyticV }, level, All );
   x_exact.p().interpolate( analyticP, level, All );
   btmp.uvw().interpolate( { rhsU, rhsV }, level, All );
   //b.uvw().interpolate( { analyticU, analyticV }, level, DirichletBoundary );

   if constexpr ( isP2P1Discr< StokesOperatorType >() )
   {
      P2ConstantMassOperator M_vel( storage, level, level );
      M_vel.apply( btmp.uvw()[0], b.uvw()[0], level, All );
      M_vel.apply( btmp.uvw()[1], b.uvw()[1], level, All );

      x.uvw().interpolate( { analyticU, analyticV }, level, hyteg::DirichletBoundary );
   }
   else
   {
      if constexpr ( isEGP0Discr< StokesOperatorType >() )
      {
         EGMassOperator M_vel( storage, level, level );
         M_vel.apply( btmp.uvw(), b.uvw(), level, All, Replace );
         x.uvw().getConformingPart()->interpolate( { analyticU, analyticV }, level, DirichletBoundary );
      }
      else
      {
         WALBERLA_ABORT( "SolVi Benchmark not implemented for other discretizations!" );
      }
   }

   StokesFunctionNumeratorType Numerator( "Num", storage, level, level );
   Numerator.enumerate( level );

   //auto Solver = solvertemplates::varViscStokesMinResSolver<StokesOperatorType>( storage, level, viscosity, 1, 1e-8, 100, true );
   PETScLUSolver< StokesOperatorType >     LU( storage, level, Numerator );
   PETScMinResSolver< StokesOperatorType > MINRES( storage, level, Numerator );
   StokesFunctionType                      nullSpace( "ns", storage, level, level );
   nullSpace.uvw().interpolate( 0, level, All );
   nullSpace.p().interpolate( 1, level, All );
   LU.setNullSpace( nullSpace );
   MINRES.setNullSpace( nullSpace );
   if constexpr ( isP2P1Discr< StokesOperatorType >() )
   {
      LU.solve( Op, x, b, level );
   }
   else
   {
      if constexpr ( isEGP0Discr< StokesOperatorType >() )
      {
         LU.solve( Op, x, b, level );
      }
   }

   // calculate the error in the L2 norm
   if constexpr ( isEGP0Discr< StokesOperatorType >() )
   {
      hyteg::dg::projectMean( x.p(), level );
      hyteg::dg::projectMean( x_exact.p(), level );
   }
   else if constexpr ( isP2P1Discr< StokesOperatorType >() )
   {
      hyteg::vertexdof::projectMean( x.p(), level );
      hyteg::vertexdof::projectMean( x_exact.p(), level );
   }

   Visc_func.interpolate( viscosity, level, All );
   if ( writeVTK )
   {
      // Visualization
      VTKOutput vtkOutput( "../../output", name, storage );
      vtkOutput.add( Visc_func );
      vtkOutput.add( x.uvw() );
      vtkOutput.add( x.p() );
      vtkOutput.add( x_exact.uvw() );
      vtkOutput.add( x_exact.p() );
      vtkOutput.add( err.uvw() );
      vtkOutput.add( err.p() );
      vtkOutput.add( b.uvw() );
      vtkOutput.add( b.p() );
      vtkOutput.write( level, 0 );
   }

   err.assign( { 1.0, -1.0 }, { x, x_exact }, level, Inner );
   real_t discrL2_velocity_err = 0.0;
   real_t discrL2_pressure_err = 0.0;

   if constexpr ( isP2P1Discr< StokesOperatorType >() )
   {
      P2ConstantMassOperator M_vel( storage, level, level );
      M_vel.apply( err.uvw()[0], Merr.uvw()[0], level, Inner, Replace );
      M_vel.apply( err.uvw()[1], Merr.uvw()[1], level, Inner, Replace );
      P1ConstantMassOperator M_pressure( storage, level, level );
      M_pressure.apply( err.p(), Merr.p(), level, All, Replace );
   }
   else if constexpr ( isEGP0Discr< StokesOperatorType >() )
   {
      EGMassOperator M_vel( storage, level, level );
      M_vel.apply( err.uvw(), Merr.uvw(), level, Inner, Replace );
      auto           mass_form = std::make_shared< dg::DGMassFormP0P0 >();
      dg::DGOperator M_pressure( storage, level, level, mass_form );
      M_pressure.apply( *err.p().getDGFunction(), *Merr.p().getDGFunction(), level, All, Replace );
   }
   discrL2_velocity_err = sqrt( err.uvw().dotGlobal( Merr.uvw(), level, Inner ) );
   discrL2_pressure_err = sqrt( err.p().dotGlobal( Merr.p(), level, Inner ) );
   real_t h = MeshQuality::getMaximalEdgeLength( storage, level );
   return std::make_tuple< real_t&, real_t&, real_t& >( h, discrL2_velocity_err, discrL2_pressure_err );
}

} // namespace hyteg

using namespace hyteg;

template < typename StokesOperatorType >
void SolViConvergenceTest( const std::string& name, const uint_t minLevel, const uint_t maxLevel )
{
   std::vector< std::tuple< real_t, real_t, real_t > > errors_per_h;
   WALBERLA_LOG_INFO_ON_ROOT( "Running " << name );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6s|%15s|%15s|%15s|%15s", "level", "error_v", "error_p", "rate_v", "rate_p" ) );
   real_t lastError_v    = std::nan( "" );
   real_t lastError_p    = std::nan( "" );
   real_t currentError_v = std::nan( "" );
   real_t currentError_p = std::nan( "" );
   real_t currentRate_v  = std::nan( "" );
   real_t currentRate_p  = std::nan( "" );
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      lastError_v      = currentError_v;
      lastError_p      = currentError_p;
      auto error_per_h = RunSolVi< StokesOperatorType >( name, level, 1, 0.2, 100.0, 1.0 );
      currentError_v   = std::get< 1 >( error_per_h );
      currentError_p   = std::get< 2 >( error_per_h );
      errors_per_h.push_back( error_per_h );
      currentRate_v = lastError_v / currentError_v;
      currentRate_p = lastError_p / currentError_p;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
          "%6d|%15.2e|%15.2e|%15.2e|%15.2e", level, currentError_v, currentError_p, currentRate_v, currentRate_p ) );
   }

   //const real_t expectedRate = 4.;
   //WALBERLA_CHECK_LESS( 0.9 * expectedRate, currentRate, "unexpected rate!" );
   //WALBERLA_CHECK_GREATER( 1.1 * expectedRate, currentRate, "unexpected rate!" );
   //WALBERLA_LOG_INFO_ON_ROOT( "Test " << testName << " converged correctly." );

   // write to plot file
   std::ofstream err_file;
   auto          err_file_name = "../../../hyteg-plots/EG_ConvOrders/" + name;
   err_file.open( err_file_name );
   for ( auto err : errors_per_h )
   {
      err_file << std::get< 0 >( err ) << ", " << std::get< 1 >( err ) << ", " << std::get< 2 >( err ) << "\n";
   }
   err_file.close();
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   /* commandline arguments for petsc solver:
   -ksp_monitor -ksp_rtol 1e-7 -ksp_type minres  -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type diag  -fieldsplit_0_ksp_type cg -fieldsplit_1_ksp_type cg -pc_fieldsplit_detect_saddle_point -fieldsplit_1_ksp_constant_null_space
   */

   SolViConvergenceTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator >( "P2P1_SolVi", 4, 8 );
   SolViConvergenceTest< hyteg::dg::eg::EGP0EpsilonStokesOperator >( "EGP0_SolVi", 4, 8 );

   return EXIT_SUCCESS;
}