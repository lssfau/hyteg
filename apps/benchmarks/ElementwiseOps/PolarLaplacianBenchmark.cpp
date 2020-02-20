/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Daniel Drzisga, Dominik Thoennes, Marcus Mohr.
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
#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/timing/Timer.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/ElementwiseOperatorPetsc.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p2_polar_laplacian.h"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

typedef enum
{
   POLARCOORDS_MAP,
   ANNULUS_MAP
} mapType;

template < typename opType >
void solve_using_geometry_map( mapType mType, uint_t minLevel, uint_t maxLevel, const char* fName );

template < typename opType >
void solve_using_pimped_form( uint_t minLevel, uint_t maxLevel, bool outputVTK );

typedef struct
{
   uint_t level;
   uint_t dofs;
   real_t maxNorm;
   real_t errNorm;
} caseResult;

void logSectionHeader( const char* opName, const char* approach )
{
   std::stringstream myStr;
   myStr << opName << " with " << approach;
   std::string hdr = myStr.str();
   size_t      len = hdr.length();
   std::string separator( len + 2, '-' );
   WALBERLA_LOG_INFO_ON_ROOT( separator << "\n " << hdr << ":\n" );
}

int main( int argc, char* argv[] )
{
   // Setup enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Going to use PETSc's LU solver
   PETScManager petscManager;

   // Set steering parameters
   uint_t minLevel = 2;
   uint_t maxLevel = 8;

   // Run tests with pimped form
   logSectionHeader( "P1ElementwiseOperator", "p1_polar_laplacian form" );
   solve_using_pimped_form< P1ElementwisePolarLaplaceOperator >( minLevel, maxLevel, false );

   logSectionHeader( "P2ElementwiseOperator", "p2_polar_laplacian form" );
   solve_using_pimped_form< P2ElementwisePolarLaplaceOperator >( minLevel, maxLevel, false );

   // Run tests with geometry map (polar coordinates)
   logSectionHeader( "P1ElementwiseOperator", "P1Form_laplace and PolarCoords map" );
   solve_using_geometry_map< P1ElementwiseBlendingLaplaceOperator >( POLARCOORDS_MAP, minLevel, maxLevel, "" );

   logSectionHeader( "P2ElementwiseOperator", "P2Form_laplace form and PolarCoords map" );
   solve_using_geometry_map< P2ElementwiseBlendingLaplaceOperator >( POLARCOORDS_MAP, minLevel, maxLevel, "" );

   // Run tests with geometry map (annulus)
   logSectionHeader( "P1ElementwiseOperator", "P1Form_laplace and AnnulusMap" );
   solve_using_geometry_map< P1ElementwiseBlendingLaplaceOperator >( ANNULUS_MAP, minLevel, maxLevel, "" );

   logSectionHeader( "P2ElementwiseOperator", "P2Form_laplace form and AnnulusMap" );
   solve_using_geometry_map< P2ElementwiseBlendingLaplaceOperator >( ANNULUS_MAP, minLevel, maxLevel, "" );

   return 0;
}

// ===============
//  generateTable
// ===============
void generateTable( const std::vector< caseResult >& results )
{
   std::stringstream tab;
   tab << "-------------------------------------------------------------------------------\n"
       << "| level |  #DoFs  |        disc. L2 norm        |         maximum norm        |\n"
       << "|       |         |     error    |   reduction  |     error    |   reduction  |\n"
       << "-------------------------------------------------------------------------------\n";
   real_t errOld = 0.0;
   real_t maxOld = 0.0;
   for ( auto entry : results )
   {
      tab << "| " << std::scientific << std::setw( 3 ) << entry.level << "   | " << std::setw( 7 ) << entry.dofs << " | "
          << entry.errNorm << " | " << errOld / entry.errNorm << " | " << entry.maxNorm << " | " << maxOld / entry.maxNorm
          << " |\n";
      errOld = entry.errNorm;
      maxOld = entry.maxNorm;
   }
   tab << "------------------------------------------------------------------------------\n";
   WALBERLA_LOG_INFO_ON_ROOT( "" << tab.str() );
}

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
   typedef typename FunctionTrait< typename funcType::template FunctionType< PetscInt > >::Tag funcTag;

   error.assign( {1.0, -1.0}, {u_exact, u}, level, All );
   real_t     npoints = static_cast< real_t >( numberOfGlobalDoFs< funcTag >( *storage, level ) );
   caseResult testCase;
   testCase.level   = level;
   testCase.dofs    = numberOfGlobalInnerDoFs< funcTag >( *storage, level );
   testCase.maxNorm = error.getMaxMagnitude( level );
   testCase.errNorm = std::sqrt( error.dotGlobal( error, level, All ) / npoints );

   return testCase;
}

// ==========================
//  solve_using_geometry_map
// ==========================
template < typename opType >
void solve_using_geometry_map( mapType mType, uint_t minLevel, uint_t maxLevel, const char* fName )
{
   typedef typename opType::srcType    funcType;
   std::shared_ptr< PrimitiveStorage > storage;
   std::string                         fileName( fName );

   // ----------------------------
   //  Perform the blending stuff
   // ----------------------------
   switch ( mType )
   {
   case POLARCOORDS_MAP:
   {
      // generate annulus mesh in polar coordinates
      real_t   rmin = 1.0;
      real_t   rmax = 2.0;
      Point2D  cornerLL( {rmin, 0.0} );
      Point2D  cornerUR( {rmax, 2.0 * pi} );
      MeshInfo meshInfo = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CROSS, 1, 6 );

      // set geometry mapping
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      for ( auto it : setupStorage.getFaces() )
      {
         setupStorage.setGeometryMap( it.second->getID(), std::make_shared< PolarCoordsMap >() );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         setupStorage.setGeometryMap( it.second->getID(), std::make_shared< PolarCoordsMap >() );
      }

      for ( auto it : setupStorage.getVertices() )
      {
         setupStorage.setGeometryMap( it.second->getID(), std::make_shared< PolarCoordsMap >() );
      }

      // prepare storage
      loadbalancing::roundRobin( setupStorage );
      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }
   break;

   case ANNULUS_MAP:
   {
      // generate coarse annulus mesh in cartesian coordinates
      MeshInfo meshInfo = MeshInfo::meshAnnulus( 1.0, 2.0, 0.0, 2.0 * pi, MeshInfo::CROSS, 6, 1 );

      // set geometry mapping
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(), std::make_shared< AnnulusMap >( face ) );
      }
      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(), std::make_shared< AnnulusMap >( edge, setupStorage ) );
      }

      // prepare storage
      loadbalancing::roundRobin( setupStorage );
      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }
   break;

   default:
      WALBERLA_ABORT( "Unknown value for 'method'! Please speak English!" );
   }

   // exact solution (in cartesian coordinates)
   std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
      real_t m   = 5.0;
      real_t rho = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t phi = std::atan2( x[1], x[0] ) + pi;
      return std::pow( 2, m ) / ( std::pow( 2, 2 * m ) + 1 ) * ( std::pow( rho, m ) + std::pow( rho, -m ) ) * std::sin( m * phi );
   };

   // generate functions
   funcType u_exact( "u_analytic", storage, minLevel, maxLevel );
   funcType u( "u_numeric", storage, minLevel, maxLevel );
   funcType rhs( "rhs", storage, minLevel, maxLevel );
   funcType error( "error", storage, minLevel, maxLevel );

   // Operator for weak-form
   opType lap( storage, minLevel, maxLevel );

   // run test cases
   std::vector< caseResult > results;
   for ( uint_t lvl = minLevel; lvl <= maxLevel; lvl++ )
   {
      // initialise functions
      u_exact.interpolate( solFunc, lvl );
      u.interpolate( real_c( 0.0 ), lvl, Inner );
      u.interpolate( solFunc, lvl, DirichletBoundary );
      rhs.interpolate( real_c( 0.0 ), lvl, All );

      // solve linear system
      PETScLUSolver< opType > solver( storage, lvl );
      solver.solve( lap, u, rhs, lvl );

      // analyse error
      results.push_back( analyseCase( storage, lvl, u_exact, u, error ) );

      // output data for visualisation
      if ( fileName.length() > 0 )
      {
         VTKOutput vtkOutput( "../../../output", fileName, storage );
         vtkOutput.add( u );
         vtkOutput.add( u_exact );
         vtkOutput.add( error );
         vtkOutput.write( lvl );
      }
   }

   generateTable( results );
}

// =========================
//  solve_using_pimped_form
// =========================
template < typename opType >
void solve_using_pimped_form( uint_t minLevel, uint_t maxLevel, bool outputVTK )
{
   // perform some template magic
   typedef typename opType::srcType funcType;
   // typedef typename FunctionTrait< typename opType::srcType::template FunctionType<PetscInt> >::Tag funcTag;

   // generate annulus mesh in polar coordinates
   real_t   rmin = 1.0;
   real_t   rmax = 2.0;
   Point2D  cornerLL( {rmin, 0.0} );
   Point2D  cornerUR( {rmax, 2.0 * pi} );
   MeshInfo meshInfo = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CROSS, 1, 6 );

   // prepare storage
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // exact solution (in polar coordinates)
   std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
      real_t m = 5.0;
      return pow( 2, m ) / ( pow( 2, 2 * m ) + 1 ) * ( pow( x[0], m ) + pow( x[0], -m ) ) * sin( m * x[1] );
   };

   // generate functions
   funcType u_exact( "u_analytic", storage, minLevel, maxLevel );
   funcType u( "u_numeric", storage, minLevel, maxLevel );
   funcType rhs( "rhs", storage, minLevel, maxLevel );
   funcType error( "error", storage, minLevel, maxLevel );

   // operator for weak-form
   opType lap( storage, minLevel, maxLevel );

   // run test cases
   std::vector< caseResult > results;
   for ( uint_t lvl = minLevel; lvl <= maxLevel; lvl++ )
   {
      // initialise functions
      u_exact.interpolate( solFunc, lvl );
      u.interpolate( real_c( 0.0 ), lvl, Inner );
      u.interpolate( solFunc, lvl, DirichletBoundary );
      rhs.interpolate( real_c( 0.0 ), lvl, All );

      // solve linear system
      PETScLUSolver< opType > solver( storage, lvl );
      solver.solve( lap, u, rhs, lvl );

      // analyse error
      results.push_back( analyseCase( storage, lvl, u_exact, u, error ) );

      // output data for visualisation
      if ( outputVTK )
      {
         VTKOutput vtkOutput( "../../../output", "elementwisePolar", storage );
         vtkOutput.add( u );
         vtkOutput.add( u_exact );
         vtkOutput.add( error );
         vtkOutput.write( lvl );
      }
   }

   generateTable( results );
}
