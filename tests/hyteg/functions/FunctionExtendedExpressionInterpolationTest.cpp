/*
 * Copyright (c) 2022 Marcus Mohr.
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
#include "core/Environment.h"
#include "core/debug/all.h"
#include "core/mpi/all.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using namespace hyteg;

namespace hyteg {

typedef enum
{
   CONST_FUNCS,
   POLY_RAT
} TestCase;

// some general settings for the tests
const uint_t level     = 2;
const bool   outputVTK = false;

template < typename myFuncType >
void run2DTest( TestCase testCase, std::shared_ptr< PrimitiveStorage > storage, const char* tag )
{
   // Setup HyTeG Functions
   myFuncType termOne( "term one in sum", storage, level, level );
   myFuncType termTwo( "term two in sum", storage, level, level );
   myFuncType termsCombined( "combined terms", storage, level, level );
   myFuncType difference( "error", storage, level, level );

   switch ( testCase )
   {
   case CONST_FUNCS: {
      WALBERLA_LOG_INFO_ON_ROOT( "[" << tag << "] Running CONST_FUNCS test in 2D:" );
      termOne.interpolate( real_c( 3.0 ), level, All );
      termTwo.interpolate( real_c( -2.0 ), level, All );

      // function object for extended expression interpolation
      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > extExpr =
          []( const Point3D&, const std::vector< real_t >& values ) { return values[0] + values[1]; };

      // perform extended expression interpolation
      termsCombined.interpolate( extExpr, { termOne, termTwo }, level, All );
   }
   break;

   case POLY_RAT: {
      WALBERLA_LOG_INFO_ON_ROOT( "[" << tag << "] Running POLY_RAT test in 2D:" );
      std::function< real_t( const Point3D& ) > polyFunc = []( const Point3D& x ) {
         return ( x[0] + real_c( 2.0 ) ) * ( x[1] - real_c( 3.0 ) );
      };
      std::function< real_t( const Point3D& ) > rat1Func = []( const Point3D& x ) {
         return real_c( 1.0 ) / ( x[0] + real_c( 2.0 ) );
      };
      termOne.interpolate( polyFunc, level, All );
      termTwo.interpolate( rat1Func, level, All );

      // function object for extended expression interpolation
      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > extExpr =
          []( const Point3D& x, const std::vector< real_t >& values ) {
             return values[0] * values[1] / ( x[1] - real_c( 3.0 ) );
          };

      // perform extended expression interpolation
      termsCombined.interpolate( extExpr, { termOne, termTwo }, level, All );
   }
   break;

   default:
      WALBERLA_ABORT( "Missing TestCase implementation!" );
   }

   // Exact result will always be constant to one everywhere,
   // compute difference
   difference.interpolate( real_c( 1.0 ), level, All );
   difference.assign( { real_c( 1.0 ), real_c( -1.0 ) }, { difference, termsCombined }, level, All );

   if ( outputVTK )
   {
      VTKOutput vtkOutput( "../../output", "FunctionExtendedExpressionInterpolationTest2D", storage );
      vtkOutput.add( termOne );
      vtkOutput.add( termTwo );
      vtkOutput.add( termsCombined );
      vtkOutput.add( difference );
      vtkOutput.write( level );
   }

   // ... and check result
   real_t error = difference.getMaxDoFMagnitude( level, All );
   WALBERLA_LOG_INFO_ON_ROOT( "--> Maximal magnitude of error = " << std::scientific << error );
   auto dp = std::is_same< real_t, double >();
   WALBERLA_CHECK_LESS( error, dp ? 1e-15 : 4e-7 );
}

template < typename myFuncType >
void run3DTest( TestCase testCase, std::shared_ptr< PrimitiveStorage > storage, const char* tag )
{
   // Setup HyTeG Functions
   myFuncType termOne( "term one in product", storage, level, level );
   myFuncType termTwo( "term two in product", storage, level, level );
   myFuncType termThr( "term three in product", storage, level, level );
   myFuncType termsCombined( "combined terms", storage, level, level );
   myFuncType difference( "error", storage, level, level );

   switch ( testCase )
   {
   case CONST_FUNCS: {
      WALBERLA_LOG_INFO_ON_ROOT( "[" << tag << "] Running CONST_FUNCS test in 3D:" );
      termOne.interpolate( real_c( 3.0 ), level, All );
      termTwo.interpolate( real_c( 1.0 / 6.0 ), level, All );
      termThr.interpolate( real_c( 2.0 ), level, All );

      // function object for extended expression interpolation
      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > extExpr =
          []( const Point3D&, const std::vector< real_t >& values ) { return values[0] * values[1] * values[2]; };

      // perform extended expression interpolation
      termsCombined.interpolate( extExpr, { termOne, termTwo, termThr }, level, All );
   }
   break;

   case POLY_RAT: {
      WALBERLA_LOG_INFO_ON_ROOT( "[" << tag << "] Running CONST_FUNCS test in 3D:" );
      std::function< real_t( const Point3D& ) > polyFunc = []( const Point3D& x ) {
         return ( x[0] + real_c( 2.0 ) ) * ( x[1] - real_c( 3.0 ) ) * ( x[2] + real_c( 5.0 ) );
      };
      std::function< real_t( const Point3D& ) > rat1Func = []( const Point3D& x ) {
         return real_c( 1.0 ) / ( x[0] + real_c( 2.0 ) );
      };
      std::function< real_t( const Point3D& ) > rat2Func = []( const Point3D& x ) {
         return real_c( 1.0 ) / ( x[1] - real_c( 3.0 ) );
      };
      termOne.interpolate( polyFunc, level, All );
      termTwo.interpolate( rat1Func, level, All );
      termThr.interpolate( rat2Func, level, All );

      // function object for extended expression interpolation
      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > extExpr =
          []( const Point3D& x, const std::vector< real_t >& values ) {
             return values[0] * values[1] * values[2] / ( x[2] + real_c( 5.0 ) );
          };

      // perform extended expression interpolation
      termsCombined.interpolate( extExpr, { termOne, termTwo, termThr }, level, All );
   }
   break;

   default:
      WALBERLA_ABORT( "Missing TestCase implementation!" );
   }

   // Exact result will always be constant to one everywhere,
   // compute difference
   difference.interpolate( real_c( 1.0 ), level, All );
   difference.assign( { real_c( 1.0 ), real_c( -1.0 ) }, { difference, termsCombined }, level, All );

   if ( outputVTK )
   {
      VTKOutput vtkOutput( "../../output", "FunctionExtendedExpressionInterpolationTest3D", storage );
      vtkOutput.add( termOne );
      vtkOutput.add( termTwo );
      vtkOutput.add( termThr );
      vtkOutput.add( termsCombined );
      vtkOutput.add( difference );
      vtkOutput.write( level );
   }

   // ... and check result
   real_t error = difference.getMaxDoFMagnitude( level, All );
   WALBERLA_LOG_INFO_ON_ROOT( "--> Maximal magnitude of error = " << std::scientific << error );
   auto dp = std::is_same< real_t, double >();
   WALBERLA_CHECK_LESS( error, dp ? 1e-15 : 4e-7 );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   // Setup enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::debug::enterTestMode();

   // ------------
   //  2D Testing
   // ------------

   // Generate mesh and primitive storage
   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( -1.0, -1.0 ), Point2D( 2.0, 1.1 ), MeshInfo::CRISSCROSS, 3, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // Run tests
   hyteg::run2DTest< P1Function< real_t > >( CONST_FUNCS, storage, "P1Function" );
   hyteg::run2DTest< P1Function< real_t > >( POLY_RAT, storage, "P1Function" );

   hyteg::run2DTest< EdgeDoFFunction< real_t > >( CONST_FUNCS, storage, "EdgeDoFFunction" );
   hyteg::run2DTest< EdgeDoFFunction< real_t > >( POLY_RAT, storage, "EdgeDoFFunction" );

   hyteg::run2DTest< P2Function< real_t > >( CONST_FUNCS, storage, "P2Function" );
   hyteg::run2DTest< P2Function< real_t > >( POLY_RAT, storage, "P2Function" );

   // NOT IMPLEMENTED FOR P0Function
   // hyteg::run2DTest< P0Function< real_t > >( CONST_FUNCS, storage, "P0Function" );
   // hyteg::run2DTest< P0Function< real_t > >( POLY_RAT, storage, "P0Function" );

   // NOT IMPLEMENTED FOR DGFunction (and those need more c'tor arguments)
   // hyteg::run2DTest< DGFunction< real_t > >( CONST_FUNCS, storage, "DGFunction" );
   // hyteg::run2DTest< DGFunction< real_t > >( POLY_RAT, storage, "DGFunction" );

   // ------------
   //  3D Testing
   // ------------
   // Generate mesh and primitive storage
   MeshInfo meshInfo3D = MeshInfo::meshSymmetricCuboid( Point3D( -1.0, -1.0, -1.0 ), Point3D( 2.0, 1.0, 3.0 ), 3, 1, 2 );
   SetupPrimitiveStorage setupStorage3D( meshInfo3D, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage3D );
   std::shared_ptr< PrimitiveStorage > storage3D = std::make_shared< PrimitiveStorage >( setupStorage3D );

   // Run tests
   hyteg::run3DTest< P1Function< real_t > >( CONST_FUNCS, storage3D, "P1Function" );
   hyteg::run3DTest< P1Function< real_t > >( POLY_RAT, storage3D, "P1Function" );

   hyteg::run3DTest< EdgeDoFFunction< real_t > >( CONST_FUNCS, storage3D, "EdgeDoFFunction" );
   hyteg::run3DTest< EdgeDoFFunction< real_t > >( POLY_RAT, storage3D, "EdgeDoFFunction" );

   hyteg::run3DTest< P2Function< real_t > >( CONST_FUNCS, storage3D, "P2Function" );
   hyteg::run3DTest< P2Function< real_t > >( POLY_RAT, storage3D, "P2Function" );

   return EXIT_SUCCESS;
}
