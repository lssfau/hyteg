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

void manufacturedSolutionTest( bool doVTKOutput = false )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running \"manufacturedSolutionTest\"" );

   // The test uses example D.3 for Steady-State flow problems
   // from John, "Finite Element Methods for Incompressible Flow Problems"
   //             _
   //            |  x^2 (1-x)^4 y^2 (1-y) (3-5y)
   // u = 1000 * |
   //            |_ -2x (1-x)^3 (^-3x) y^3 (1-y)^2
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

   // setup FE problem
   CCRStokesOperator           stokesOperator( storage, level, level );
   CCRStokesFunction< real_t > discrete_rhs( "Discrete RHS", storage, level, level );
   CCRStokesFunction< real_t > discrete_sol( "Discrete Solution", storage, level, level );

   // note: problem has no-slip boundary conditions, so nothing to do for this

   EigenSparseDirectSolver< CCRStokesOperator > EigenLU( storage, level );
   // EigenLU.setReassembleMatrix( true );
   EigenLU.solve( stokesOperator, discrete_sol, discrete_rhs, level );

   // export data for post-processing
   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "ccr", storage );
      vtkOutput.add( sol_analytic );
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
