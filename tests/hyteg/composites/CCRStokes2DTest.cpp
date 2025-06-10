/*
 * Copyright (c) 2017-2025 Marcus Mohr.
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

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"

#include "hyteg/composites/CCRStokesFunction.hpp"
#include "hyteg/composites/CCRStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/eigen/EigenSparseDirectSolver.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "mixed_operator/VectorMassOperator.hpp"

// #ifndef HYTEG_BUILD_WITH_PETSC
// WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
// #endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

namespace hyteg {

void objectInstantionTest()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running \"objectInstantionTest\"" );

   uint_t level = 2;

   MeshInfo meshInfo = MeshInfo::meshRectangle(
       Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), hyteg::MeshInfo::CRISSCROSS, 2, 2 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   CCRStokesFunction< real_t > func( "CCR", storage, level, level );
   CCRStokesOperator           oper( storage, level, level );
}

enum class ExpressionType
{
   UX,
   UY,
   P,
   RHS_X,
   RHS_Y
};

std::function< real_t( const hyteg::Point3D& ) > getAnalyticalExpression( ExpressionType exType )
{
   // The test uses example D.3 for Steady-State flow problems
   // from John, "Finite Element Methods for Incompressible Flow Problems"
   //             _
   //            |  x^2 (1-x)^4 y^2 (1-y) (3-5y)
   // u = 1000 * |
   //            |_ -2x (1-x)^3 (1-3x) y^3 (1-y)^2
   //
   //
   // p = pi^2 ( x y^3 cos( 2 pi x^2 y ) - x^2 y sin( 2 pi x y ) + 1/8
   //
   // The problem is posed in the unit square with no-slip boundary conditions
   // for velocity.
   std::function< real_t( const hyteg::Point3D& ) > expression;
   switch ( exType )
   {
   case ( ExpressionType::UX ):
      expression = []( const hyteg::Point3D& coords ) {
         real_t x = coords[0];
         real_t y = coords[1];

         real_t tmp1 = ( real_c( 1 ) - x );
         real_t tmp2 = tmp1 * tmp1;
         real_t tmp3 = tmp2 * tmp2;

         return x * x * tmp3 * y * y * ( real_c( 1 ) - y ) * ( real_c( 3 ) - real_c( 5 ) * y );
      };
      break;

   case ( ExpressionType::UY ):
      expression = []( const hyteg::Point3D& coords ) {
         real_t x = coords[0];
         real_t y = coords[1];

         real_t tmp1 = ( real_c( 1 ) - x );
         real_t tmp2 = tmp1 * tmp1 * tmp1;
         real_t tmp3 = ( real_c( 1 ) - y ) * ( real_c( 1 ) - y );
         real_t y3   = y * y * y;
         return -real_c( 2 ) * x * tmp2 * ( real_c( 1 ) - real_c( 3 ) * x ) * y3 * tmp3;
      };
      break;

   case ( ExpressionType::P ):
      expression = []( const hyteg::Point3D& coords ) {
         real_t x = coords[0];
         real_t y = coords[1];

         real_t x2 = x * x;
         real_t y3 = y * y * y;

         return pi * pi * ( x * y3 * std::cos( real_c( 2 ) * pi * x2 * y ) - x2 * y * std::sin( real_c( 2 ) * pi * x * y ) ) +
                real_c( 1.0 / 8.0 );
      };
      break;

   case ( ExpressionType::RHS_X ):
      expression = []( const hyteg::Point3D& coords ) {
         real_t x = coords[0];
         real_t y = coords[1];

         real_t t1  = x * x;
         real_t t4  = 0.2e1 * y * t1 * pi;
         real_t t5  = std::cos( t4 );
         real_t t6  = pi * pi;
         real_t t8  = y * y;
         real_t t9  = y * t8;
         real_t t11 = std::sin( t4 );
         real_t t12 = pi * t6;
         real_t t14 = t8 * t8;
         real_t t20 = 0.2e1 * pi * x * y;
         real_t t21 = std::cos( t20 );
         real_t t26 = std::sin( t20 );
         real_t t32 = std::pow( -0.1e1 + x, 0.2e1 );
         real_t t33 = 0.4e1 / 0.5e1 * y;
         real_t t35 = t1 * t1;
         real_t retVal =
             t9 * t6 * t5 - 0.4e1 * t14 * t1 * t12 * t11 - 0.2e1 * t8 * t1 * t12 * t21 - 0.2e1 * x * y * t6 * t26 -
             0.60000e5 *
                 ( t35 * ( t8 - t33 + 0.1e1 / 0.10e2 ) + x * t1 * ( -0.2e1 * t8 + 0.8e1 / 0.5e1 * y - 0.1e1 / 0.5e1 ) +
                   t1 * ( 0.5e1 / 0.2e1 * t14 - 0.4e1 * t9 + 0.5e1 / 0.2e1 * t8 - t33 + 0.1e1 / 0.10e2 ) +
                   x * ( -0.5e1 / 0.3e1 * t14 + 0.8e1 / 0.3e1 * t9 - t8 ) +
                   ( -0.3e1 / 0.5e1 + y ) * ( -0.1e1 + y ) * t8 / 0.6e1 ) *
                 t32;

         return retVal;
      };
      break;

   case ( ExpressionType::RHS_Y ):
      expression = []( const hyteg::Point3D& coords ) {
         real_t x = coords[0];
         real_t y = coords[1];

         real_t t1     = x * x;
         real_t t4     = 0.2e1 * y * t1 * pi;
         real_t t5     = std::cos( t4 );
         real_t t6     = pi * pi;
         real_t t8     = y * y;
         real_t t12    = std::sin( t4 );
         real_t t13    = pi * t6;
         real_t t15    = x * t1;
         real_t t16    = y * t8;
         real_t t22    = 0.2e1 * pi * x * y;
         real_t t23    = std::cos( t22 );
         real_t t28    = std::sin( t22 );
         real_t t31    = -0.1e1 + x;
         real_t t34    = t8 * t8;
         real_t t40    = t1 * t1;
         real_t t46    = t31 * t31;
         real_t t48    = ( -0.1e1 / 0.3e1 + x ) * t46;
         real_t retVal = 0.3e1 * t8 * x * t6 * t5 - 0.2e1 * t16 * t15 * t13 * t12 - 0.2e1 * y * t15 * t13 * t23 - t1 * t6 * t28 +
                         0.120000e6 *
                             ( t34 * ( t1 - x + 0.1e1 / 0.5e1 ) + t16 * ( -0.2e1 * t1 + 0.2e1 * x - 0.2e1 / 0.5e1 ) +
                               t8 * ( t40 - 0.7e1 / 0.3e1 * t15 + 0.8e1 / 0.3e1 * t1 - 0.4e1 / 0.3e1 * x + 0.1e1 / 0.5e1 ) -
                               0.6e1 / 0.5e1 * x * y * t48 + 0.3e1 / 0.10e2 * x * t48 ) *
                             y * t31;

         return retVal;
      };
      break;
   }

   return expression;
}

void manufacturedSolutionTest( bool doVTKOutput = false )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running \"manufacturedSolutionTest\"" );

   // The test uses example D.3 for Steady-State flow problems
   // from John, "Finite Element Methods for Incompressible Flow Problems"
   //             _
   //            |  x^2 (1-x)^4 y^2 (1-y) (3-5y)
   // u = 1000 * |
   //            |_ -2x (1-x)^3 (1-3x) y^3 (1-y)^2
   //
   //
   // p = pi^2 ( x y^3 cos( 2 pi x^2 y ) - x^2 y sin( 2 pi x y ) + 1/8
   //
   // The problem is posed in the unit square with no-slip boundary conditions
   // for velocity.

   std::function< real_t( const hyteg::Point3D& ) > ux_analytic = []( const hyteg::Point3D& coords ) {
      real_t x = coords[0];
      real_t y = coords[1];

      real_t tmp1 = ( real_c( 1 ) - x );
      real_t tmp2 = tmp1 * tmp1;
      real_t tmp3 = tmp2 * tmp2;

      return x * x * tmp3 * y * y * ( real_c( 1 ) - y ) * ( real_c( 3 ) - real_c( 5 ) * y );
   };

   std::function< real_t( const hyteg::Point3D& ) > uy_analytic = []( const hyteg::Point3D& coords ) {
      real_t x = coords[0];
      real_t y = coords[1];

      real_t tmp1 = ( real_c( 1 ) - x );
      real_t tmp2 = tmp1 * tmp1 * tmp1;
      real_t tmp3 = ( real_c( 1 ) - y ) * ( real_c( 1 ) - y );
      real_t y3   = y * y * y;
      return -real_c( 2 ) * x * tmp2 * ( real_c( 1 ) - real_c( 3 ) * x ) * y3 * tmp3;
   };

   std::function< real_t( const hyteg::Point3D& ) > p_analytic = []( const hyteg::Point3D& coords ) {
      real_t x = coords[0];
      real_t y = coords[1];

      real_t x2 = x * x;
      real_t y3 = y * y * y;

      return pi * pi * ( x * y3 * std::cos( real_c( 2 ) * pi * x2 * y ) - x2 * y * std::sin( real_c( 2 ) * pi * x * y ) ) +
             real_c( 1.0 / 8.0 );
   };

   std::function< real_t( const hyteg::Point3D& ) > rhsMomentumX = []( const hyteg::Point3D& coords ) {
      real_t x = coords[0];
      real_t y = coords[1];

      real_t t1  = x * x;
      real_t t4  = 0.2e1 * y * t1 * pi;
      real_t t5  = std::cos( t4 );
      real_t t6  = pi * pi;
      real_t t8  = y * y;
      real_t t9  = y * t8;
      real_t t11 = std::sin( t4 );
      real_t t12 = pi * t6;
      real_t t14 = t8 * t8;
      real_t t20 = 0.2e1 * pi * x * y;
      real_t t21 = std::cos( t20 );
      real_t t26 = std::sin( t20 );
      real_t t32 = std::pow( -0.1e1 + x, 0.2e1 );
      real_t t33 = 0.4e1 / 0.5e1 * y;
      real_t t35 = t1 * t1;
      real_t retVal =
          t9 * t6 * t5 - 0.4e1 * t14 * t1 * t12 * t11 - 0.2e1 * t8 * t1 * t12 * t21 - 0.2e1 * x * y * t6 * t26 -
          0.60000e5 *
              ( t35 * ( t8 - t33 + 0.1e1 / 0.10e2 ) + x * t1 * ( -0.2e1 * t8 + 0.8e1 / 0.5e1 * y - 0.1e1 / 0.5e1 ) +
                t1 * ( 0.5e1 / 0.2e1 * t14 - 0.4e1 * t9 + 0.5e1 / 0.2e1 * t8 - t33 + 0.1e1 / 0.10e2 ) +
                x * ( -0.5e1 / 0.3e1 * t14 + 0.8e1 / 0.3e1 * t9 - t8 ) + ( -0.3e1 / 0.5e1 + y ) * ( -0.1e1 + y ) * t8 / 0.6e1 ) *
              t32;

      return retVal;
   };

   std::function< real_t( const hyteg::Point3D& ) > rhsMomentumY = []( const hyteg::Point3D& coords ) {
      real_t x = coords[0];
      real_t y = coords[1];

      real_t t1     = x * x;
      real_t t4     = 0.2e1 * y * t1 * pi;
      real_t t5     = std::cos( t4 );
      real_t t6     = pi * pi;
      real_t t8     = y * y;
      real_t t12    = std::sin( t4 );
      real_t t13    = pi * t6;
      real_t t15    = x * t1;
      real_t t16    = y * t8;
      real_t t22    = 0.2e1 * pi * x * y;
      real_t t23    = std::cos( t22 );
      real_t t28    = std::sin( t22 );
      real_t t31    = -0.1e1 + x;
      real_t t34    = t8 * t8;
      real_t t40    = t1 * t1;
      real_t t46    = t31 * t31;
      real_t t48    = ( -0.1e1 / 0.3e1 + x ) * t46;
      real_t retVal = 0.3e1 * t8 * x * t6 * t5 - 0.2e1 * t16 * t15 * t13 * t12 - 0.2e1 * y * t15 * t13 * t23 - t1 * t6 * t28 +
                      0.120000e6 *
                          ( t34 * ( t1 - x + 0.1e1 / 0.5e1 ) + t16 * ( -0.2e1 * t1 + 0.2e1 * x - 0.2e1 / 0.5e1 ) +
                            t8 * ( t40 - 0.7e1 / 0.3e1 * t15 + 0.8e1 / 0.3e1 * t1 - 0.4e1 / 0.3e1 * x + 0.1e1 / 0.5e1 ) -
                            0.6e1 / 0.5e1 * x * y * t48 + 0.3e1 / 0.10e2 * x * t48 ) *
                          y * t31;

      return retVal;
   };

   // prepare domain and mesh
   MeshInfo meshInfo = MeshInfo::meshRectangle(
       Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), hyteg::MeshInfo::CRISSCROSS, 2, 2 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   uint_t                      level = 4;
   CCRStokesFunction< real_t > sol_analytic( "Analytic Solution", storage, level, level );

   sol_analytic.uvw().interpolate( { ux_analytic, uy_analytic }, level );
   sol_analytic.p().interpolate( p_analytic, level );

   CCRStokesFunction< real_t > rhs_analytic( "Analytic RHS", storage, level, level );
   rhs_analytic.uvw().interpolate( { rhsMomentumX, rhsMomentumY }, level );

   // setup FE problem
   CCRStokesOperator           stokesOperator( storage, level, level );
   CCRStokesFunction< real_t > discrete_rhs( "Discrete RHS", storage, level, level );
   CCRStokesFunction< real_t > discrete_sol( "Discrete Solution", storage, level, level );

   // note: problem has no-slip boundary conditions, so nothing to do for this

   P2PlusBubbleElementwiseVectorMassOperator massOperator( storage, level, level );
   massOperator.apply( rhs_analytic.uvw(), discrete_rhs.uvw(), level, All );

   EigenSparseDirectSolver< CCRStokesOperator > EigenLU( storage, level );
   // EigenLU.setReassembleMatrix( true );
   EigenLU.solve( stokesOperator, discrete_sol, discrete_rhs, level );

   // PETScLUSolver< CCRStokesOperator > petscLU( storage, level );
   // petscLU.solve( stokesOperator, discrete_sol, discrete_rhs, level );

   // export data for post-processing
   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "ccr", storage );
      vtkOutput.add( sol_analytic );
      vtkOutput.add( rhs_analytic );
      vtkOutput.add( discrete_sol );
      vtkOutput.add( discrete_rhs );
      vtkOutput.write( level );
   }
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   objectInstantionTest();
   manufacturedSolutionTest( true );

   return EXIT_SUCCESS;
}
