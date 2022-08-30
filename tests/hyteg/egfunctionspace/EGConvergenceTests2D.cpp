/*
* Copyright (c) 2017-2022 Nils Kohl.
*
* This file is part of HyTeG
* (see https://i10git.cs.fau.de/hyteg/hyteg).
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.cpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_t;

// scalar lambda for one component of analytical solution and rhs
typedef std::function< real_t( const hyteg::PointND< real_t, 3 >& p ) > ScalarLambda;

// tuple of function for solution (u,p) and rhs of vector values stokes equation
typedef std::tuple< ScalarLambda, ScalarLambda, ScalarLambda > LambdaTuple;

using hyteg::Point3D;
using hyteg::dg::eg::EGLaplaceOperator;
using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0ConstEpsilonStokesOperator;
using hyteg::dg::eg::EGP0EpsilonStokesOperator;
using hyteg::dg::eg::EGP0StokesOperator;

namespace hyteg {
auto copyBdry = []( EGP0StokesFunction< real_t > fun ) { fun.p().setBoundaryCondition( fun.uvw().getBoundaryCondition() ); };

template < typename StokesOperatorType >
class EGStokesConvergenceTest
{
 public:
   using StokesFunctionType          = typename StokesOperatorType::srcType;
   using StokesFunctionNumeratorType = typename StokesFunctionType::template FunctionType< idx_t >;

   EGStokesConvergenceTest( const std::string&                         testName,
                            LambdaTuple                                sol_tuple,
                            LambdaTuple                                rhs_tuple,
                            StokesOperatorType&                        Op,
                            const std::shared_ptr< PrimitiveStorage >& storage,
                            const uint_t                               minLevel,
                            const uint_t                               maxLevel,
                            bool                                       writeVTK = false )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Running " << testName );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6s|%15s|%15s", "level", "error", "rate" ) );
      real_t lastError    = std::nan( "" );
      real_t currentError = std::nan( "" );
      real_t currentRate  = std::nan( "" );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         lastError    = currentError;
         currentError = RunStokesTestOnLevel( testName, level, sol_tuple, rhs_tuple, Op, storage, writeVTK );
         currentRate  = lastError / currentError;
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%15.2e|%15.2e", level, currentError, currentRate ) );
      }
      const real_t expectedRate = 4.;
      //WALBERLA_CHECK_LESS( 0.9 * expectedRate, currentRate, "unexpected rate!" );
      //WALBERLA_CHECK_GREATER( 1.1 * expectedRate, currentRate, "unexpected rate!" );
      //WALBERLA_LOG_INFO_ON_ROOT( "Test " << testName << " converged correctly." );
   }

   real_t RunStokesTestOnLevel( const std::string&                         testName,
                                const uint_t&                              level,
                                LambdaTuple                                sol_tuple,
                                LambdaTuple                                rhs_tuple,
                                StokesOperatorType&                        Op,
                                const std::shared_ptr< PrimitiveStorage >& storage,
                                bool                                       writeVTK = false )
   {
      StokesFunctionNumeratorType numerator( "numerator", storage, level, level );
      numerator.enumerate( level );

      // solution, rhs as a lambda function
      auto [u_x_expr, u_y_expr, p_expr] = sol_tuple;
      auto [f_x_expr, f_y_expr, g_expr] = rhs_tuple;

      StokesFunctionType u( "u", storage, level, level );
      StokesFunctionType f( "f", storage, level, level );
      StokesFunctionType rhs( "rhs", storage, level, level );
      StokesFunctionType sol( "sol", storage, level, level );
      StokesFunctionType err( "err", storage, level, level );
      StokesFunctionType Merr( "Merr", storage, level, level );
      if constexpr ( std::is_same< StokesOperatorType, EGP0EpsilonStokesOperator >::value ||
                     std::is_same< StokesOperatorType, EGP0StokesOperator >::value ||
                     std::is_same< StokesOperatorType, EGP0ConstEpsilonStokesOperator >::value )
      {
         copyBdry( u );
         copyBdry( f );
         copyBdry( rhs );
         copyBdry( sol );
         copyBdry( err );
         copyBdry( Merr );
      }

      // interpolate analytical solution and rhs
      sol.uvw().interpolate( { u_x_expr, u_y_expr }, level, All );
      sol.p().interpolate( p_expr, level, All );
      f.uvw().interpolate( { f_x_expr, f_y_expr }, level, All );
      f.p().interpolate( g_expr, level, All );

      if constexpr ( std::is_same< StokesOperatorType, hyteg::P2P1ElementwiseAffineEpsilonStokesOperator >::value ||
                     std::is_same< StokesOperatorType, hyteg::P2P1TaylorHoodStokesOperator >::value )
      {
         P2ConstantMassOperator M_vel( storage, level, level );
         M_vel.apply( f.uvw()[0], rhs.uvw()[0], level, All );
         M_vel.apply( f.uvw()[1], rhs.uvw()[1], level, All );
         u.uvw().interpolate( { u_x_expr, u_y_expr }, level, hyteg::DirichletBoundary );

         P1ConstantMassOperator M_pressure( storage, level, level );
         M_pressure.apply( f.p(), rhs.p(), level, All, Replace );
      }
      else
      {
         if constexpr ( std::is_same< StokesOperatorType, EGP0EpsilonStokesOperator >::value ||
                        std::is_same< StokesOperatorType, EGP0StokesOperator >::value ||
                        std::is_same< StokesOperatorType, EGP0ConstEpsilonStokesOperator >::value )
         {
            EGMassOperator M_vel( storage, level, level );
            M_vel.apply( f.uvw(), rhs.uvw(), level, All, Replace );
            u.uvw().getConformingPart()->interpolate( { u_x_expr, u_y_expr }, level, DirichletBoundary );

            auto           mass_form = std::make_shared< dg::DGMassFormP0P0 >();
            dg::DGOperator M_pressure( storage, level, level, mass_form );
            M_pressure.apply( *f.p().getDGFunction(), *rhs.p().getDGFunction(), level, All, Replace );
         }
         else
         {
            WALBERLA_ABORT( "Benchmark not implemented for other discretizations!" );
         }
      }

      // solve
      //PETScMinResSolver< StokesOperatorType > solver( storage, level, numerator );
      PETScLUSolver< StokesOperatorType > solver( storage, level, numerator );
      StokesFunctionType                  nullSpace( "ns", storage, level, level );
      nullSpace.uvw().interpolate( 0, level, All );
      nullSpace.p().interpolate( 1, level, All );
      solver.setNullSpace( nullSpace );
      solver.solve( Op, u, rhs, level );
      // calculate the error in the L2 norm
      err.assign( { 1.0, -1.0 }, { u, sol }, level, Inner );
      if constexpr ( std::is_same< StokesOperatorType, hyteg::P2P1ElementwiseAffineEpsilonStokesOperator >::value ||
                     std::is_same< StokesOperatorType, hyteg::P2P1TaylorHoodStokesOperator >::value )
      {
         P2ConstantMassOperator M_vel( storage, level, level );
         M_vel.apply( err.uvw()[0], Merr.uvw()[0], level, Inner, Replace );
         M_vel.apply( err.uvw()[1], Merr.uvw()[1], level, Inner, Replace );
      }
      else
      {
         if constexpr ( std::is_same< StokesOperatorType, EGP0EpsilonStokesOperator >::value ||
                        std::is_same< StokesOperatorType, EGP0StokesOperator >::value ||
                        std::is_same< StokesOperatorType, EGP0ConstEpsilonStokesOperator >::value )
         {
            EGMassOperator M_vel( storage, level, level );
            M_vel.apply( err.uvw(), Merr.uvw(), level, Inner, Replace );
         }
      }
      auto discrL2_velocity = sqrt( err.uvw().dotGlobal( Merr.uvw(), level, Inner ) );

      if ( writeVTK )
      {
         VTKOutput vtk( "../../output", testName, storage );
         vtk.add( u );
         vtk.add( sol );
         vtk.add( err );
         vtk.add( f );
         vtk.write( level );
      }

      return discrL2_velocity;
   }
};

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   auto meshInfo = hyteg::MeshInfo::meshRectangle(
       hyteg::Point2D( { -1, -1 } ), hyteg::Point2D( { 1, 1 } ), hyteg::MeshInfo::CRISSCROSS, 1, 1 );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

   EGP0StokesOperator             EGP0StokesOp( storage, 3, 5 );
   EGP0ConstEpsilonStokesOperator EGP0ConstantEpsilonOp( storage, 3, 5 );
   EGP0EpsilonStokesOperator      EGP0EpsilonOp_mu_1( storage, 3, 5, []( const hyteg::Point3D& p ) { return 1; } );
   EGP0EpsilonStokesOperator      EGP0EpsilonOp_mu_smooth( storage, 3, 5, []( const hyteg::Point3D& p ) {
      const real_t x = p[0];
      const real_t y = p[1];
      return std::exp( x ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                 std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) +
             1;
   } );
   hyteg::P2P1ElementwiseAffineEpsilonStokesOperator P2P1ElementwiseEpsilonOp_mu_smooth(
       storage, 3, 5, []( const hyteg::Point3D& p ) {
          const real_t x = p[0];
          const real_t y = p[1];
          return std::exp( x ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                     std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) +
                 1;
       } );
   hyteg::P2P1ElementwiseAffineEpsilonStokesOperator P2P1ElementwiseEpsilonOp_mu_1(
       storage, 3, 5, []( const hyteg::Point3D& p ) { return 1; } );
   hyteg::P2P1TaylorHoodStokesOperator P2P1StokesOp( storage, 3, 5 );

   /* commandline arguments for petsc solver:
   -ksp_monitor -ksp_rtol 1e-7 -ksp_type minres  -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type diag  -fieldsplit_0_ksp_type cg -fieldsplit_1_ksp_type cg -pc_fieldsplit_detect_saddle_point -fieldsplit_1_ksp_constant_null_space
   */

   hyteg::EGStokesConvergenceTest< EGP0EpsilonStokesOperator >(
       "EGP0EpsilonOp_muSmooth_divFreeVelocity",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return y + 2 * std::sin( M_PI * ( x + y ) ) + 4;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return -x - 2 * std::sin( M_PI * ( x + y ) ) + 3;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return 2*x - y + 1;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( x + y );
              const real_t x1 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2 = std::exp( x ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
              const real_t x3 = x2 * std::sin( x1 );
              return 4.0 * std::pow( M_PI, 2 ) * ( x3 + 1 ) * std::sin( x0 ) -
                     4.0 * M_PI * ( ( 1.0 / 2.0 ) * M_PI * x2 * std::cos( x1 ) + x3 ) * std::cos( x0 ) + 2;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::exp( x );
              const real_t x1 = std::pow( M_PI, 2 );
              const real_t x2 = M_PI * ( x + y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              return 2.0 * x0 * x1 * x3 * std::cos( x2 ) * std::cos( x4 ) -
                     4.0 * x1 * ( x0 * x3 * std::sin( x4 ) + 1 ) * std::sin( x2 ) - 1;
           },
           []( const Point3D& p ) -> real_t { return 0; } ),
       EGP0EpsilonOp_mu_smooth,
       storage,
       3,
       5 );

       hyteg::EGStokesConvergenceTest<  hyteg::P2P1ElementwiseAffineEpsilonStokesOperator >(
      "P2P1EpsilonOp_muSmooth_divFreeVelocity",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return y + 2 * std::sin( M_PI * ( x + y ) ) + 4;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return -x - 2 * std::sin( M_PI * ( x + y ) ) + 3;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return 2*x - y + 1;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( x + y );
              const real_t x1 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2 = std::exp( x ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
              const real_t x3 = x2 * std::sin( x1 );
              return 4.0 * std::pow( M_PI, 2 ) * ( x3 + 1 ) * std::sin( x0 ) -
                     4.0 * M_PI * ( ( 1.0 / 2.0 ) * M_PI * x2 * std::cos( x1 ) + x3 ) * std::cos( x0 ) + 2;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::exp( x );
              const real_t x1 = std::pow( M_PI, 2 );
              const real_t x2 = M_PI * ( x + y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              return 2.0 * x0 * x1 * x3 * std::cos( x2 ) * std::cos( x4 ) -
                     4.0 * x1 * ( x0 * x3 * std::sin( x4 ) + 1 ) * std::sin( x2 ) - 1;
           },
           []( const Point3D& p ) -> real_t { return 0; } ),
       P2P1ElementwiseEpsilonOp_mu_smooth,
       storage,
       3,
       5 );
   
   
   hyteg::EGStokesConvergenceTest< EGP0StokesOperator >(
       "DefaultStokesHomogeneousDirichlet_asymmetric_u",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4 = x2 * std::sin( x3 );
              const real_t x5 = M_PI_2;
              return -x1 * x2 * x5 * std::cos( x3 ) - x1 * x4 - x4 * x5 * std::cos( x0 );
           } ),
       EGP0StokesOp,
       storage,
       3,
       5,
       true );

   hyteg::EGStokesConvergenceTest< EGP0StokesOperator >(
       "DefaultStokesHomogeneousDirichlet:pointsymmetric_u",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                         std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) +
                     1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                         std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) +
                     1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2 = M_PI_2;
              return -x2 * std::sin( x0 ) * std::cos( x1 ) - x2 * std::sin( x1 ) * std::cos( x0 );
           } ),
       EGP0StokesOp,
       storage,
       3,
       5 );

   hyteg::EGStokesConvergenceTest< EGP0StokesOperator >(
       "DefaultStokesHomogeneousDirichlet:asymmetric_p",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * x ) * std::sin( M_PI * y ) * std::sin( M_PI * ( x + y ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * x ) * std::sin( M_PI * y ) * std::sin( M_PI * ( x + y ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return 2 * x - y + 4;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * x;
              const real_t x1 = std::sin( x0 );
              const real_t x2 = M_PI * ( x + y );
              const real_t x3 = std::pow( M_PI, 2 );
              const real_t x4 = M_PI * y;
              const real_t x5 = x3 * std::sin( x4 );
              const real_t x6 = 2 * std::cos( x2 );
              return -x1 * x3 * x6 * std::cos( x4 ) + 4 * x1 * x5 * std::sin( x2 ) - x5 * x6 * std::cos( x0 ) + 2;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::pow( M_PI, 2 );
              const real_t x1 = M_PI * x;
              const real_t x2 = std::sin( x1 );
              const real_t x3 = M_PI * y;
              const real_t x4 = std::sin( x3 );
              const real_t x5 = M_PI * ( x + y );
              const real_t x6 = 2 * std::cos( x5 );
              return 4 * x0 * x2 * x4 * std::sin( x5 ) - x0 * x2 * x6 * std::cos( x3 ) - x0 * x4 * x6 * std::cos( x1 ) - 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * y;
              const real_t x1 = M_PI * x;
              const real_t x2 = std::sin( x1 );
              const real_t x3 = M_PI * ( x + y );
              const real_t x4 = M_PI * std::sin( x3 );
              const real_t x5 = std::sin( x0 );
              return -x2 * x4 * std::cos( x0 ) - 2 * M_PI * x2 * x5 * std::cos( x3 ) - x4 * x5 * std::cos( x1 );
           } ),
       EGP0StokesOp,
       storage,
       3,
       5 );

   hyteg::EGStokesConvergenceTest< EGP0StokesOperator >( "DefaultStokesInhomogeneousDirichlet:asymmetric_u",
                                                         std::make_tuple(
                                                             []( const Point3D& p ) -> real_t {
                                                                const real_t x = p[0];
                                                                const real_t y = p[1];
                                                                return y + 2 * std::sin( M_PI * ( x + y ) ) + 4;
                                                             },
                                                             []( const Point3D& p ) -> real_t {
                                                                const real_t x = p[0];
                                                                const real_t y = p[1];
                                                                return -x - 2 * std::sin( M_PI * ( x + y ) ) + 3;
                                                             },
                                                             []( const Point3D& p ) -> real_t {
                                                                const real_t x = p[0];
                                                                const real_t y = p[1];
                                                                return 2 * x - y + 1;
                                                             } ),
                                                         std::make_tuple(
                                                             []( const Point3D& p ) -> real_t {
                                                                const real_t x = p[0];
                                                                const real_t y = p[1];
                                                                return 4 * std::pow( M_PI, 2 ) * std::sin( M_PI * ( x + y ) ) + 2;
                                                             },
                                                             []( const Point3D& p ) -> real_t {
                                                                const real_t x = p[0];
                                                                const real_t y = p[1];
                                                                return -4 * std::pow( M_PI, 2 ) * std::sin( M_PI * ( x + y ) ) -
                                                                       1;
                                                             },
                                                             []( const Point3D& p ) -> real_t { return 0; } ),
                                                         EGP0StokesOp,
                                                         storage,
                                                         3,
                                                         5 );

   hyteg::EGStokesConvergenceTest< EGP0ConstEpsilonStokesOperator >(
       "ConstEpsilonHomogeneousDirichlet:asymmetric_u",
       std::make_tuple(
           []( const Point3D& xx ) -> real_t {
              return std::sin( M_PI * ( xx[0] + 1.0 ) / 2.0 ) * std::sin( M_PI * ( xx[1] + 1.0 ) / 2.0 ) * std::exp( xx[1] );
           },
           []( const Point3D& xx ) -> real_t {
              return std::sin( M_PI * ( xx[0] + 1.0 ) / 2.0 ) * std::sin( M_PI * ( xx[1] + 1.0 ) / 2.0 ) * std::exp( xx[1] );
           },
           []( const Point3D& xx ) -> real_t {
              return std::sin( M_PI * ( xx[0] + 1.0 ) / 2.0 ) * std::sin( M_PI * ( xx[1] + 1.0 ) / 2.0 ) * std::exp( xx[1] );
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4 = std::sin( x3 );
              const real_t x5 = x2 * x4;
              const real_t x6 = x2 * std::cos( x3 );
              const real_t x7 = std::cos( x0 );
              const real_t x8 = std::pow( M_PI, 2 );
              return 0.75 * x1 * x2 * x4 * x8 - 1.0 * x1 * x5 - 1.0 * M_PI * x1 * x6 - 0.5 * M_PI * x5 * x7 -
                     0.25 * x6 * x7 * x8 + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4 = std::sin( x3 );
              const real_t x5 = x2 * x4;
              const real_t x6 = x2 * std::cos( x3 );
              const real_t x7 = std::cos( x0 );
              const real_t x8 = std::pow( M_PI, 2 );
              return 0.75 * x1 * x2 * x4 * x8 - 2.0 * x1 * x5 - 2.0 * M_PI * x1 * x6 - 0.5 * M_PI * x5 * x7 -
                     0.25 * x6 * x7 * x8 + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4 = x2 * std::sin( x3 );
              const real_t x5 = M_PI_2;
              return -x1 * x2 * x5 * std::cos( x3 ) - x1 * x4 - x4 * x5 * std::cos( x0 );
              ;
           } ),
       EGP0ConstantEpsilonOp,
       storage,
       3,
       5 );

   hyteg::EGStokesConvergenceTest< EGP0ConstEpsilonStokesOperator >(
       "ConstEpsilonHomogeneousDirichlet:symmetric_u",
       std::make_tuple(
           []( const Point3D& xx ) -> real_t {
              return std::sin( M_PI * ( xx[0] + 1.0 ) / 2.0 ) * std::sin( M_PI * ( xx[1] + 1.0 ) / 2.0 );
           },
           []( const Point3D& xx ) -> real_t {
              return std::sin( M_PI * ( xx[0] + 1.0 ) / 2.0 ) * std::sin( M_PI * ( xx[1] + 1.0 ) / 2.0 );
           },
           []( const Point3D& xx ) -> real_t { return xx[0] + xx[1]; } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::pow( M_PI, 2 );
              const real_t x1 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              return 0.75 * x0 * std::sin( x1 ) * std::sin( x2 ) - 0.25 * x0 * std::cos( x1 ) * std::cos( x2 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::pow( M_PI, 2 );
              const real_t x1 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              return 0.75 * x0 * std::sin( x1 ) * std::sin( x2 ) - 0.25 * x0 * std::cos( x1 ) * std::cos( x2 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2 = M_PI_2;
              return -( x2 * std::sin( x0 ) * std::cos( x1 ) + x2 * std::sin( x1 ) * std::cos( x0 ) );
           } ),
       EGP0ConstantEpsilonOp,
       storage,
       3,
       5 );

   hyteg::EGStokesConvergenceTest< EGP0ConstEpsilonStokesOperator >(
       "ConstEpsilonInhomogeneousDirichlet",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) + std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) + std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& xx ) -> real_t { return xx[0] + xx[1]; } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::pow( M_PI, 2 );
              return 0.5 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     0.25 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::pow( M_PI, 2 );
              return 0.25 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     0.5 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI_2;
              return -x0 * std::cos( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) -
                     x0 * std::cos( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           } ),
       EGP0ConstantEpsilonOp,
       storage,
       3,
       5 );

   hyteg::EGStokesConvergenceTest< EGP0EpsilonStokesOperator >(
       "EpsilonHomogeneousDirichlet:ConstantViscosity",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) + std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) + std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& xx ) -> real_t { return xx[0] + xx[1]; } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::pow( M_PI, 2 );
              return 0.5 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     0.25 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::pow( M_PI, 2 );
              return 0.25 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     0.5 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI_2;
              return -x0 * std::cos( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) -
                     x0 * std::cos( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           } ),
       EGP0EpsilonOp_mu_1,
       storage,
       3,
       5 );

   hyteg::EGStokesConvergenceTest< EGP0EpsilonStokesOperator >(
       "EpsilonInhomogeneousDirichlet:constantViscosity",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) + std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) + std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& xx ) -> real_t { return xx[0] + xx[1]; } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::pow( M_PI, 2 );
              return 0.5 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     0.25 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::pow( M_PI, 2 );
              return 0.25 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     0.5 * x0 * std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI_2;
              return -x0 * std::cos( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) -
                     x0 * std::cos( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           } ),
       EGP0EpsilonOp_mu_1,
       storage,
       3,
       5 );

   hyteg::EGStokesConvergenceTest< EGP0EpsilonStokesOperator >(
       "EpsilonInhomogeneousDirichlet:smoothViscosity",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::exp( y ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::exp( y ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              const real_t x0  = std::exp( y );
              const real_t x1  = std::pow( M_PI, 2 );
              const real_t x2  = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x3  = std::sin( x2 );
              const real_t x4  = std::exp( x );
              const real_t x5  = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x6  = std::sin( x5 );
              const real_t x7  = x4 * x6;
              const real_t x8  = x3 * x7;
              const real_t x9  = x8 + 1;
              const real_t x10 = std::cos( x2 );
              const real_t x11 = x0 * x10;
              const real_t x12 = std::cos( x5 );
              const real_t x13 = 0.25 * M_PI;
              const real_t x14 = 0.5 * x0 * x3 + x11 * x13;
              return 0.5 * x0 * x1 * x3 * x9 - 1.0 * M_PI * x11 * ( ( 1.0 / 2.0 ) * M_PI * x10 * x7 + x8 ) -
                     M_PI * x12 * x3 * x4 * ( x12 * x13 + x14 ) - 2 * x9 * ( -0.125 * x1 * x6 + x14 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x   = p[0];
              const real_t y   = p[1];
              const real_t x0  = std::exp( x );
              const real_t x1  = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2  = std::sin( x1 );
              const real_t x3  = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4  = M_PI * std::cos( x3 );
              const real_t x5  = std::exp( y );
              const real_t x6  = x2 * x5;
              const real_t x7  = 1.0 * x6;
              const real_t x8  = std::pow( M_PI, 2 );
              const real_t x9  = std::sin( x3 );
              const real_t x10 = x0 * x9;
              const real_t x11 = x10 * x2;
              const real_t x12 = 2 * x11 + 2;
              const real_t x13 = std::cos( x1 );
              const real_t x14 = M_PI * x13;
              return -x0 * x2 * x4 * ( 0.5 * x4 + x7 ) - x12 * ( x7 - 0.25 * x8 * x9 ) -
                     x12 * ( 0.25 * M_PI * x13 * x5 - 0.125 * x6 * x8 ) -
                     2 * ( ( 1.0 / 2.0 ) * x10 * x14 + x11 ) * ( 0.25 * x14 * x5 + 0.25 * x4 + 0.5 * x6 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI_2;
              const real_t x1 = std::exp( y );
              const real_t x2 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              return -x0 * x1 * std::cos( x2 ) - x0 * std::cos( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) - x1 * std::sin( x2 );
           } ),
       EGP0EpsilonOp_mu_smooth,
       storage,
       3,
       5,
       false );

   hyteg::EGStokesConvergenceTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator >(
       "THEpsilonInhomogeneousDirichlet:smoothViscosity",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::exp( y ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::exp( y ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) +
                     std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              const real_t x0  = std::exp( y );
              const real_t x1  = std::pow( M_PI, 2 );
              const real_t x2  = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x3  = std::sin( x2 );
              const real_t x4  = std::exp( x );
              const real_t x5  = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x6  = std::sin( x5 );
              const real_t x7  = x4 * x6;
              const real_t x8  = x3 * x7;
              const real_t x9  = x8 + 1;
              const real_t x10 = std::cos( x2 );
              const real_t x11 = x0 * x10;
              const real_t x12 = std::cos( x5 );
              const real_t x13 = 0.25 * M_PI;
              const real_t x14 = 0.5 * x0 * x3 + x11 * x13;
              return 0.5 * x0 * x1 * x3 * x9 - 1.0 * M_PI * x11 * ( ( 1.0 / 2.0 ) * M_PI * x10 * x7 + x8 ) -
                     M_PI * x12 * x3 * x4 * ( x12 * x13 + x14 ) - 2 * x9 * ( -0.125 * x1 * x6 + x14 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x   = p[0];
              const real_t y   = p[1];
              const real_t x0  = std::exp( x );
              const real_t x1  = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2  = std::sin( x1 );
              const real_t x3  = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4  = M_PI * std::cos( x3 );
              const real_t x5  = std::exp( y );
              const real_t x6  = x2 * x5;
              const real_t x7  = 1.0 * x6;
              const real_t x8  = std::pow( M_PI, 2 );
              const real_t x9  = std::sin( x3 );
              const real_t x10 = x0 * x9;
              const real_t x11 = x10 * x2;
              const real_t x12 = 2 * x11 + 2;
              const real_t x13 = std::cos( x1 );
              const real_t x14 = M_PI * x13;
              return -x0 * x2 * x4 * ( 0.5 * x4 + x7 ) - x12 * ( x7 - 0.25 * x8 * x9 ) -
                     x12 * ( 0.25 * M_PI * x13 * x5 - 0.125 * x6 * x8 ) -
                     2 * ( ( 1.0 / 2.0 ) * x10 * x14 + x11 ) * ( 0.25 * x14 * x5 + 0.25 * x4 + 0.5 * x6 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI_2;
              const real_t x1 = std::exp( y );
              const real_t x2 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              return -x0 * x1 * std::cos( x2 ) - x0 * std::cos( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) - x1 * std::sin( x2 );
           } ),
       P2P1ElementwiseEpsilonOp_mu_smooth,
       storage,
       3,
       5,
       false );

   hyteg::EGStokesConvergenceTest< hyteg::P2P1TaylorHoodStokesOperator >(
       "P2P1StokesHomogeneousDirichlet_asymmetric_u",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4 = x2 * std::sin( x3 );
              const real_t x5 = M_PI_2;
              return -x1 * x2 * x5 * std::cos( x3 ) - x1 * x4 - x4 * x5 * std::cos( x0 );
           } ),
       P2P1StokesOp,
       storage,
       3,
       5,
       false );

   hyteg::EGStokesConvergenceTest< hyteg::P2P1TaylorHoodStokesOperator >(
       "P2P1StokesHomogeneousDirichlet:simple_u",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * x );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return std::sin( M_PI * y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::pow( M_PI, 2 ) * std::sin( M_PI * x ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              return std::pow( M_PI, 2 ) * std::sin( M_PI * y ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return -M_PI * std::cos( M_PI * x ) - M_PI * std::cos( M_PI * y );
           } ),
       P2P1StokesOp,
       storage,
       3,
       5,
       false );

   hyteg::EGStokesConvergenceTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator >(
       "P2P1EpsilonHomogeneousDirichlet_asymmetric_u",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::sin( M_PI * ( x + 1 ) / 2 ) * std::sin( M_PI * ( y + 1 ) / 2 ) * std::exp( y );
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return x + y;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];

              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = x2 * x3;
              return ( 1.0 / 2.0 ) * std::pow( M_PI, 2 ) * x1 * x2 * x3 - x1 * x4 - M_PI * x4 * std::cos( x0 ) + 1;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x1 = std::sin( x0 );
              const real_t x2 = std::exp( y );
              const real_t x3 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              const real_t x4 = x2 * std::sin( x3 );
              const real_t x5 = M_PI_2;
              return -x1 * x2 * x5 * std::cos( x3 ) - x1 * x4 - x4 * x5 * std::cos( x0 );
           } ),
       P2P1ElementwiseEpsilonOp_mu_1,
       storage,
       3,
       5,
       false );
       
   return 0;
}
