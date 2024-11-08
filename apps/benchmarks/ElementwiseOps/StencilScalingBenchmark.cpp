/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

// ----------------------------------------------------------------------
// Numeric benchmark testing using the 2D scenario from Sec. 7.1 of the
// Bauer et al. paper "A stencil scaling approach for accelerating
// matrix-free finite element implementations" (doi:10.1137/17M1148384)
// ----------------------------------------------------------------------

#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/timing/Timer.h>

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

typedef struct
{
   uint_t level;
   uint_t dofs;
   real_t maxNorm;
   real_t errNorm;
} caseResult;

// =============
//  analyseCase
// =============
template < typename funcType >
caseResult analyseCase( std::shared_ptr< PrimitiveStorage > storage,
                        uint_t                              level,
                        const funcType&                     u_exact,
                        const funcType&                     u,
                        funcType&                           error )
{
   typedef typename FunctionTrait< typename funcType::template FunctionType< idx_t > >::Tag funcTag;

   // embed numeric solution in next finer space
   P2toP2QuadraticProlongation embeddor;
   embeddor.prolongate( u, level, All );

   error.assign( {1.0, -1.0}, {u_exact, u}, level + 1, All );
   real_t     npoints = static_cast< real_t >( numberOfGlobalDoFs< funcTag >( *storage, level + 1 ) );
   caseResult testCase;
   testCase.level   = level;
   testCase.dofs    = numberOfGlobalInnerDoFs< funcTag >( *storage, level );
   testCase.maxNorm = error.getMaxDoFMagnitude( level + 1 );
   testCase.errNorm = std::sqrt( error.dotGlobal( error, level + 1, All ) / npoints );

   return testCase;
}

// ============
//  addToTable
// ============
void addToTable( const std::vector< caseResult >& results, uint_t curLevel, uint_t minLevel, uint_t maxLevel )
{
   std::stringstream tab;
   uint_t            idx = curLevel - minLevel;
   if ( curLevel == minLevel )
   {
      tab << "-------------------------------------------------------------------------------\n"
          << "| level |  #DoFs  |        disc. L2 norm        |         maximum norm        |\n"
          << "|       |         |     error    |   reduction  |     error    |   reduction  |\n"
          << "-------------------------------------------------------------------------------\n";
   }

   static real_t     errOld = 0.0;
   static real_t     maxOld = 0.0;
   const caseResult& entry  = results[idx];
   tab << "| " << std::scientific << std::setw( 3 ) << entry.level << "   | " << std::setw( 7 ) << entry.dofs << " | "
       << entry.errNorm << " | " << errOld / entry.errNorm << " | " << entry.maxNorm << " | " << maxOld / entry.maxNorm << " |";
   errOld = entry.errNorm;
   maxOld = entry.maxNorm;
   if ( curLevel == maxLevel )
   {
      tab << "\n-------------------------------------------------------------------------------\n";
   }
   WALBERLA_LOG_INFO_ON_ROOT( "" << tab.str() );
}

int main( int argc, char* argv[] )
{
   // Setup enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Going to use PETSc's LU solver
   PETScManager petscManager( &argc, &argv );

   // Set steering parameters
   uint_t minLevel  = 2;
   uint_t maxLevel  = 8;
   uint_t freqM     = 2;
   bool   outputVTK = true;

   std::function< real_t( const Point3D& ) > manufacturedSolution = []( const Point3D& x ) {
      return x[0] * x[0] * x[0] * x[0] * x[1] / ( real_c( 1.0 ) + x[0] * x[1] );
   };

   std::function< real_t( const Point3D& ) > parameterFunc = [freqM]( const Point3D& x ) {
      real_t m = real_c( freqM );
      return real_c( 2.0 ) + std::sin( m * pi * x[0] ) * std::sin( m * pi * x[1] );
   };

   std::function< real_t( const Point3D& ) > manufacturedRHS = [freqM]( const Point3D& x ) {
      real_t m     = real_c( freqM );
      real_t t1    = m * pi;
      real_t t2    = t1 * x[0];
      real_t t3    = std::cos( t2 );
      real_t t4    = t1 * x[1];
      real_t t5    = std::sin( t4 );
      real_t t7    = x[0] * x[0];
      real_t t8    = t7 * x[0];
      real_t t11   = x[0] * x[1] + 0.1e1;
      real_t t12   = 0.1e1 / t11;
      real_t t15   = t7 * t7;
      real_t t16   = x[1] * x[1];
      real_t t18   = t11 * t11;
      real_t t19   = 0.1e1 / t18;
      real_t t24   = std::sin( t2 );
      real_t t26   = t24 * t5 + 0.2e1;
      real_t t36   = 0.1e1 / t18 / t11;
      real_t t42   = std::cos( t4 );
      real_t t45   = t15 * x[0];
      real_t value = t1 * t3 * t5 * ( 0.4e1 * t8 * x[1] * t12 - t15 * t16 * t19 ) +
                     t26 * ( 0.2e1 * t15 * t16 * x[1] * t36 + 0.12e2 * t7 * x[1] * t12 - 0.8e1 * t8 * t16 * t19 ) +
                     t24 * m * pi * t42 * ( -t45 * x[1] * t19 + t15 * t12 ) +
                     t26 * ( 0.2e1 * t15 * t7 * x[1] * t36 - 0.2e1 * t45 * t19 );
      return ( -value );
   };

   // setup mesh and primitives
   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISSCROSS, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // make the form include the parameter function in the integral
   P2Form_divKgrad::callback = parameterFunc;

   typedef P2Function< real_t > funcType;
   funcType                     u_exact( "u_analytic", storage, minLevel, maxLevel + 1 );
   funcType                     u( "u_numeric", storage, minLevel, maxLevel + 1 );
   funcType                     rhs( "rhs", storage, minLevel, maxLevel );
   funcType                     aux( "aux", storage, minLevel, maxLevel );
   funcType                     error( "error", storage, minLevel, maxLevel + 1 );

   // operator for weak-form
   P2ElementwiseDivKGradOperator divKgrad( storage, minLevel, maxLevel );
   P2ConstantMassOperator        mass( storage, minLevel, maxLevel );

   // setup analytic solution for comparison
   for ( uint_t lvl = minLevel; lvl <= maxLevel + 1; lvl++ )
   {
      u_exact.interpolate( manufacturedSolution, lvl );
   }

   // run test cases
   std::vector< caseResult > results;
   for ( uint_t lvl = minLevel; lvl <= maxLevel; lvl++ )
   {
      // initialise functions
      u.interpolate( real_c( 0.0 ), lvl, Inner );
      u.interpolate( manufacturedSolution, lvl, DirichletBoundary );
      aux.interpolate( manufacturedRHS, lvl, All );

      // compute FEM right-hand side
      mass.apply( aux, rhs, lvl, All );

      // solve linear system
      PETScLUSolver< P2ElementwiseDivKGradOperator > solver( storage, lvl );
      solver.solve( divKgrad, u, rhs, lvl );

      // analyse error
      results.push_back( analyseCase( storage, lvl, u_exact, u, error ) );
      addToTable( results, lvl, minLevel, maxLevel );

      // output data for visualisation
      if ( outputVTK )
      {
         funcType k( "param func", storage, lvl, lvl );
         k.interpolate( parameterFunc, lvl, All );

         VTKOutput vtkOutput( "../../../output", "section7", storage );
         vtkOutput.add( u );
         vtkOutput.add( u_exact );
         vtkOutput.add( k );

         // have not computed error for this level, yet
         error.assign( {1.0, -1.0}, {u_exact, u}, lvl, All );
         vtkOutput.add( error );

         vtkOutput.write( lvl );
      }
   }
}
