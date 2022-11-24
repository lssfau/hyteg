
#include <complex>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/egfunctionspace/EGConvTestUtils.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.cpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using hyteg::P2P1ElementwiseAffineEpsilonStokesOperator;
using hyteg::Point3D;
using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0EpsilonStokesOperator;
using hyteg::dg::eg::EGP0StokesOperator;
using hyteg::dg::eg::isEGP0Discr;
using hyteg::dg::eg::isP2P1Discr;
using hyteg::dg::eg::LambdaTuple;
using hyteg::dg::eg::ScalarLambda;
using hyteg::dg::eg::StokesConvergenceOrderTest;

namespace hyteg {

// SolVi benchmark (Circular inclusion)
/*  Described in section 5.3 of [1] with analytical solution from [2] (equations 26, 34).
The viscosity contains a jump from visc_matrix to visc_inclusion in a circle located at the center of the domain. 
Difficulty for FE methods: the mesh can not align with the viscosity jump well.

Implementation in parts from [3]

[1]: "On the choice of finite element for applications in geodynamics" 
by Cedric Thieulot and Wolfgang Bangerth

[2]: "Analytical solutions for deformable elliptical inclusions in general shear" 
by Daniel W. Schmid and Yuri Yu. Podladchikov

[3]: https://github.com/geodynamics/aspect/blob/main/benchmarks/inclusion/inclusion.h
*/

// returns (u,p) analytical solution, (f,g) rhs and mu viscosity function of the SolVi problem
std::tuple< LambdaTuple, LambdaTuple, ScalarLambda >
    SetupSolViSolution( const real_t r_inclusion, const real_t& visc_inclusion, const real_t& visc_matrix )
{
   using namespace std::complex_literals;

   // radius helper function
   std::function< real_t( const hyteg::Point3D& ) > rad = []( const hyteg::Point3D& xx ) {
      return sqrt( std::pow( xx[0], 2.0 ) + std::pow( xx[1], 2.0 ) );
   };

   // viscosity function and operator setup
   std::function< real_t( const hyteg::Point3D& ) > viscosity =
       [r_inclusion, visc_inclusion, visc_matrix, rad]( const hyteg::Point3D& xx ) {
          Point3D offset( { 1.0, 1.0, 0.0 } );
          if ( rad( xx - offset ) < r_inclusion )
             return visc_inclusion;
          else
             return visc_matrix;
       };

   // analytic solution for u,v,p
   const real_t                                             C_visc = visc_matrix / ( visc_inclusion + visc_matrix );
   const real_t                                             A      = C_visc * ( visc_inclusion - visc_matrix );
   std::function< hyteg::Point3D( const hyteg::Point3D& ) > analytic_uvp =
       [A, r_inclusion, visc_matrix, visc_inclusion, rad]( const hyteg::Point3D& xx ) {
          std::complex< real_t > phi, psi, dphi;
          real_t                 r2_inclusion = r_inclusion * r_inclusion;
          Point3D                offset( { 1.0, 1.0, 0.0 } );
          Point3D                tmp = xx - offset;
          real_t                 x   = tmp[0];
          real_t                 y   = tmp[1];
          real_t                 r2  = std::pow( rad( tmp ), 2.0 ); //x * x + y * y;

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
          real_t                 visc = ( r2 < r2_inclusion ) ? visc_inclusion : visc_matrix;
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

   // Right-hand-side: derivatives of u, v, p for x and y
   std::function< real_t( const hyteg::Point3D& ) > ddx_u = [r_inclusion, A]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[0] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( 12.0 * std::pow( xx[0], 4.0 ) - 120.0 * std::pow( xx[0] * xx[1], 2.0 ) + 60.0 * std::pow( xx[1], 4.0 ) ) -
               4.0 * std::pow( xx[0], 6.0 ) + 52.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) +
               20.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) - 36.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > ddy_u = [r_inclusion, A]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[0] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( -12.0 * std::pow( xx[0], 4.0 ) + 120.0 * std::pow( xx[0] * xx[1], 2.0 ) - 60.0 * std::pow( xx[1], 4.0 ) ) +
               12.0 * std::pow( xx[0], 6.0 ) - 60.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) -
               60.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) + 12.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > ddx_v = [r_inclusion, A]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[1] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( 60.0 * std::pow( xx[0], 4.0 ) - 120.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 2.0 ) +
                     12.0 * std::pow( xx[1], 4.0 ) ) -
               12.0 * std::pow( xx[0], 6.0 ) + 60.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) +
               60.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) - 12.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > ddy_v = [r_inclusion, A]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[1] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( -60.0 * std::pow( xx[0], 4.0 ) + 120.0 * std::pow( xx[0] * xx[1], 2.0 ) - 12.0 * std::pow( xx[1], 4.0 ) ) +
               36.0 * std::pow( xx[0], 6.0 ) - 20.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) -
               52.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) + 4.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > dydx_v = [r_inclusion, A]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[0] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( -12.0 * std::pow( xx[0], 4.0 ) + 120.0 * std::pow( xx[0] * xx[1], 2.0 ) - 60.0 * std::pow( xx[1], 4.0 ) ) +
               4.0 * std::pow( xx[0], 6.0 ) - 52.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) -
               20.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) + 36.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > dxdy_u = [r_inclusion, A]( const hyteg::Point3D& xx ) {
      return A * std::pow( r_inclusion, 2.0 ) * xx[1] *
             ( std::pow( r_inclusion, 2.0 ) *
                   ( 60.0 * std::pow( xx[0], 4.0 ) - 120.0 * std::pow( xx[0] * xx[1], 2.0 ) + 12.0 * std::pow( xx[1], 4.0 ) ) -
               36.0 * std::pow( xx[0], 6.0 ) + 20.0 * std::pow( xx[0], 4.0 ) * std::pow( xx[1], 2.0 ) +
               52.0 * std::pow( xx[0], 2.0 ) * std::pow( xx[1], 4.0 ) - 4.0 * std::pow( xx[1], 6.0 ) ) /
             ( std::pow( xx[0] * xx[0] + xx[1] * xx[1], 5.0 ) );
   };

   // right hand side: setup by epsilon operator on u,v and gradient of p
   std::function< real_t( const hyteg::Point3D& ) > rhsU =
       [r_inclusion, A, rad, ddx_u, ddy_u, dydx_v]( const hyteg::Point3D& xx ) {
          if ( rad( xx ) < r_inclusion )
             return 0.0;
          else
             return -0.5 * ( 2 * ddx_u( xx ) + ddy_u( xx ) + dydx_v( xx ) ) +
                    real_c( 2 ) * A * std::pow( r_inclusion, 2.0 ) *
                        ( xx[0] * cos( 2.0 * atan( xx[1] / xx[0] ) ) - xx[1] * sin( 2.0 * atan( xx[1] / xx[0] ) ) ) /
                        std::pow( xx[0] * xx[0] + xx[1] * xx[1], 2.0 ); //+ ddx_p(xx);
       };
   std::function< real_t( const hyteg::Point3D& ) > rhsV =
       [r_inclusion, A, rad, ddx_v, ddy_v, dxdy_u]( const hyteg::Point3D& xx ) {
          if ( rad( xx ) < r_inclusion )
             return 0.0;
          else
             return -0.5 * ( ddx_v( xx ) + 2 * ddy_v( xx ) + dxdy_u( xx ) ) +
                    real_c( 2 ) * A * std::pow( r_inclusion, 2.0 ) *
                        ( xx[0] * sin( 2.0 * atan( xx[1] / xx[0] ) ) + xx[1] * cos( 2.0 * atan( xx[1] / xx[0] ) ) ) /
                        std::pow( xx[0] * xx[0] + xx[1] * xx[1], 2.0 ); //+ ddy_p(xx);
       };

   /*
   PETScLUSolver< StokesOperatorType > solver( storage, level );
   //PETScMinResSolver< StokesOperatorType > solver( storage, level );

   StokesFunctionType nullSpace( "ns", storage, level, level );

   nullSpace.uvw().interpolate( 0, level, All );
   nullSpace.p().interpolate( 1, level, All );
   solver.setNullSpace( nullSpace, level );

   solver.solve( Op, x_num, b, level );
*/

   // pack
   auto        zero     = []( const hyteg::Point3D& ) { return real_c( 0 ); };
   LambdaTuple solTuple = std::make_tuple( analyticU, analyticV, zero, analyticP );
   LambdaTuple rhsTuple = std::make_tuple( rhsU, rhsV, zero, zero );

   return std::make_tuple( solTuple, rhsTuple, viscosity );
}

} // namespace hyteg

using namespace hyteg;
/*
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
   real_t mean_rate_v    = 0.;
   real_t mean_rate_p    = 0.;
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
      if ( level > minLevel )
      {
         mean_rate_v += currentRate_v;
         mean_rate_p += currentRate_p;
      }

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
          "%6d|%15.2e|%15.2e|%15.2e|%15.2e", level, currentError_v, currentError_p, currentRate_v, currentRate_p ) );
   }
   mean_rate_v = mean_rate_v / ( real_c( maxLevel ) - real_c( minLevel ) );
   mean_rate_p = mean_rate_p / ( real_c( maxLevel ) - real_c( minLevel ) );
   WALBERLA_LOG_INFO_ON_ROOT( name << ": mean conv. rate v = " << mean_rate_v << ", mean conv. rate p = " << mean_rate_p );

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
*/
int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );
   /* commandline arguments for petsc solver:
   -ksp_monitor -ksp_rtol 1e-7 -ksp_type minres  -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type diag  -fieldsplit_0_ksp_type cg -fieldsplit_1_ksp_type cg -pc_fieldsplit_detect_saddle_point -fieldsplit_1_ksp_constant_null_space
   */
   uint_t minLevel = 3;
   uint_t maxLevel = 7;

   // storage setup
   auto meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 2, 2 } ), MeshInfo::CRISS, 2, 2 );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

   // SolVi solution setup
   auto [solTuple, rhsTuple, viscosity] = SetupSolViSolution( 0.2, 100.0, 1.0 );

   if ( false )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Running SolVi with EGP0 ###" );

      EGP0EpsilonStokesOperator EGP0EpsilonOp( storage, minLevel, maxLevel, viscosity );

      StokesConvergenceOrderTest< EGP0EpsilonStokesOperator >(
          "SolVi_EGP0", solTuple, rhsTuple, EGP0EpsilonOp, storage, minLevel, maxLevel, 1, true );
   }

   if ( true )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Running SolVi with P2P1 ###" );

      P2P1ElementwiseAffineEpsilonStokesOperator P2P1EpsilonOp( storage, minLevel, maxLevel, viscosity );

      StokesConvergenceOrderTest< P2P1ElementwiseAffineEpsilonStokesOperator >(
          "SolVi_P2P1", solTuple, rhsTuple, P2P1EpsilonOp, storage, minLevel, maxLevel-1, 1, true );
   }
   return EXIT_SUCCESS;
}