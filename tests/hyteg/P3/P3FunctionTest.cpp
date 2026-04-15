/*
 * Copyright (c) 2026 Marcus Mohr.
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

#include "hyteg/p3functionspace/P3Function.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p3functionspace/P3VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

void printSeparator()
{
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------------------" );
}

void testDoFCounting()
{
   printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "DoF Countingtion Test" );
   WALBERLA_LOG_INFO_ON_ROOT( "-> running 2D test" );
   {
      MeshInfo              mesh2D = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/penta_5el.msh" ) );
      SetupPrimitiveStorage setupStorage2D( mesh2D, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      std::shared_ptr< PrimitiveStorage > storage2D = std::make_shared< PrimitiveStorage >( setupStorage2D );

      const uint_t minLevel  = 0;
      const uint_t maxLevel  = 2;
      const uint_t numLevels = maxLevel - minLevel + 1;

      P3Function< real_t > p3Func( "P3Function in 2D", storage2D, minLevel, maxLevel );

      std::array< uint_t, numLevels > controlDoFs{ 31, 106, 391 };

      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         uint_t numDoFs = p3Func.getNumberOfGlobalDoFs( level );
         WALBERLA_LOG_INFO( "* level " << level << ": " << p3Func.getNumberOfGlobalDoFs( level ) );
         WALBERLA_LOG_INFO( "* level " << level << ": " << controlDoFs[level] << " (ctrl)" );
         WALBERLA_CHECK_EQUAL( numDoFs, controlDoFs[level] );
      }
   }
}

void testInstantiation()
{
   hyteg::printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "Basic Object Instantiation Test" );

   MeshInfo                            mesh2D = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
   SetupPrimitiveStorage               setupStorage2D( mesh2D, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage2D = std::make_shared< PrimitiveStorage >( setupStorage2D );

   const uint_t minLevel = 0;
   const uint_t maxLevel = 2;

   P3Function< float >   p3Float( "P3Function< float > in 2D", storage2D, minLevel, maxLevel );
   P3Function< double >  p3Double( "P3Function< double > in 2D", storage2D, minLevel, maxLevel );
   P3Function< int32_t > p3Int32( "P3Function< int32_t > in 2D", storage2D, minLevel, maxLevel );
   P3Function< int64_t > p3Int64( "P3Function< int64_t > in 2D", storage2D, minLevel, maxLevel );

   P3VectorFunction< float >   p3VectorFloat( "P3VectorFunction< float > in 2D", storage2D, minLevel, maxLevel );
   P3VectorFunction< double >  p3VectorDouble( "P3VectorFunction< double > in 2D", storage2D, minLevel, maxLevel );
   P3VectorFunction< int32_t > p3VectorInt32( "P3VectorFunction< int32_t > in 2D", storage2D, minLevel, maxLevel );
   P3VectorFunction< int64_t > p3VectorInt64( "P3VectorFunction< int64_t > in 2D", storage2D, minLevel, maxLevel );
}

void testInterpolation()
{
   hyteg::printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "Interpolation Test" );
   WALBERLA_LOG_INFO_ON_ROOT( "-> running 2D test" );
   {
      Point2D lowerLeft( real_c( 0.0 ), real_c( 0.0 ) );
      Point2D upperRight( real_c( 1.0 ), real_c( 1.0 ) );

      MeshInfo              mesh2D = MeshInfo::meshRectangle( lowerLeft, upperRight, MeshInfo::meshFlavour::CRISS, 1, 1 );
      SetupPrimitiveStorage setupStorage2D( mesh2D, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      std::shared_ptr< PrimitiveStorage > storage2D = std::make_shared< PrimitiveStorage >( setupStorage2D );

      const uint_t level = 2;

      P3Function< real_t > funcA( "P3Function 'A'", storage2D, level, level );
      P3Function< real_t > funcB( "P3Function 'B'", storage2D, level, level );

      WALBERLA_LOG_INFO_ON_ROOT( "* interpolation with constant expression" );

      // interpolation of constant function
      std::function< real_t( const hyteg::Point3D& ) > constExpr = []( const Point3D& ) -> real_t { return real_c( 2 ); };
      funcA.interpolate( constExpr, level, All );

      WALBERLA_CHECK_EQUAL( funcA.getMaxDoFValue( level ), real_c( 2.0 ) );
      WALBERLA_CHECK_EQUAL( funcA.getMinDoFValue( level ), real_c( 2.0 ) );
      WALBERLA_CHECK_EQUAL( funcA.getMaxDoFMagnitude( level ), real_c( 2.0 ) );

      WALBERLA_LOG_INFO_ON_ROOT( "* interpolation with non-constant expression" );

      // interpolation of 3rd order polynomial
      std::function< real_t( const hyteg::Point3D& ) > poly3 = []( const Point3D& x ) -> real_t {
         return real_c( 2 ) * x[0] * x[0] * x[0] + x[0] * x[1] - real_c( 0.5 );
      };
      funcB.interpolate( poly3, level, All );

      // testing P3Function as a whole
      real_t maxDoF  = funcB.getMaxDoFValue( level );
      real_t control = poly3( Point3D( real_c( 1.0 ), real_c( 1.0 ), real_c( 0.0 ) ) );
      WALBERLA_CHECK_EQUAL( maxDoF, control );

      WALBERLA_CHECK_EQUAL( funcB.getMinDoFValue( level ), real_c( -0.5 ) );

      // testing FaceDoF component
      control = poly3( Point3D( real_c( 11.0 / 12.0 ), real_c( 11.0 / 12.0 ), real_c( 0.0 ) ) );
      maxDoF  = funcB.getFaceDoFFunction().getMaxDoFValue( level );
      WALBERLA_CHECK_EQUAL( maxDoF, control );

      control       = poly3( Point3D( real_c( 1.0 / 12.0 ), real_c( 1.0 / 12.0 ), real_c( 0.0 ) ) );
      real_t minDoF = funcB.getFaceDoFFunction().getMinDoFValue( level );
      WALBERLA_CHECK_EQUAL( minDoF, control );

      // testing EdgeDoFBlue component
      control = poly3( Point3D( real_c( 1.0 ), real_c( 10.0 / 12.0 ), real_c( 0.0 ) ) );
      maxDoF  = funcB.getEdgeDoFFunctionBlue().getMaxDoFValue( level );
      WALBERLA_CHECK_EQUAL( maxDoF, control );

      // testing EdgeDoFRed component
      control = poly3( Point3D( real_c( 1.0 ), real_c( 11.0 / 12.0 ), real_c( 0.0 ) ) );
      maxDoF  = funcB.getEdgeDoFFunctionRed().getMaxDoFValue( level );
      WALBERLA_CHECK_EQUAL( maxDoF, control );
   }
}

void testAssignAndAdd()
{
   hyteg::printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "Assign + Add Test" );
   WALBERLA_LOG_INFO_ON_ROOT( "-> running 2D test" );
   {
      Point2D lowerLeft( real_c( 0.0 ), real_c( 0.0 ) );
      Point2D upperRight( real_c( 1.0 ), real_c( 1.0 ) );

      MeshInfo              mesh2D = MeshInfo::meshRectangle( lowerLeft, upperRight, MeshInfo::meshFlavour::CRISSCROSS, 1, 1 );
      SetupPrimitiveStorage setupStorage2D( mesh2D, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      std::shared_ptr< PrimitiveStorage > storage2D = std::make_shared< PrimitiveStorage >( setupStorage2D );

      const uint_t level = 2;

      P3Function< real_t > funcA( "P3Function 'A'", storage2D, level, level );
      P3Function< real_t > funcB( "P3Function 'B'", storage2D, level, level );
      P3Function< real_t > funcC( "P3Function 'C'", storage2D, level, level );

      bool onlyQuadratic = false;

      std::function< real_t( const hyteg::Point3D& ) > myExpr = [&onlyQuadratic]( const Point3D& x ) -> real_t {
         real_t quadPart = real_c( 3 ) * x[0] * ( x[0] + x[1] );
         return onlyQuadratic ? quadPart : quadPart + real_c( 2 ) * x[0] - real_c( 0.5 ) * x[1];
      };

      funcA.interpolate( myExpr, level );
      onlyQuadratic = true;
      funcB.interpolate( myExpr, level );

      WALBERLA_LOG_INFO_ON_ROOT( "* testing assign" );
      funcC.assign( { real_c( 1 ), real_c( -1 ) }, { funcA, funcB }, level, All );

      real_t maxDoF = funcC.getMaxDoFValue( level );
      real_t minDoF = funcC.getMinDoFValue( level );
      WALBERLA_CHECK_EQUAL( minDoF, real_c( -0.5 ) );
      WALBERLA_CHECK_EQUAL( maxDoF, real_c( +2.0 ) );

      WALBERLA_LOG_INFO_ON_ROOT( "* testing add" );
      funcC.add( { real_c( 0.5 ) }, { funcB }, level, All );
      maxDoF = funcC.getMaxDoFValue( level );
      minDoF = funcC.getMinDoFValue( level );
      WALBERLA_CHECK_EQUAL( minDoF, real_c( -0.5 ) );
      WALBERLA_CHECK_EQUAL( maxDoF, real_c( +4.5 ) );
   }
}

void testDotProduct()
{
   hyteg::printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "Dot Product Test" );
   WALBERLA_LOG_INFO_ON_ROOT( "-> running 2D test" );
   {
      Point2D lowerLeft( real_c( 0.0 ), real_c( 0.0 ) );
      Point2D upperRight( real_c( 1.0 ), real_c( 1.0 ) );

      MeshInfo              mesh2D = MeshInfo::meshRectangle( lowerLeft, upperRight, MeshInfo::meshFlavour::CRISSCROSS, 1, 1 );
      SetupPrimitiveStorage setupStorage2D( mesh2D, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      std::shared_ptr< PrimitiveStorage > storage2D = std::make_shared< PrimitiveStorage >( setupStorage2D );

      const uint_t level = 2;

      P3Function< real_t > funcA( "P3Function 'A'", storage2D, level, level );
      P3Function< real_t > funcB( "P3Function 'B'", storage2D, level, level );

      funcA.interpolate( real_c( 2 ), level );
      funcB.interpolate( real_c( 3 ), level );

      real_t dot = funcA.dotGlobal( funcB, level, All );

      uint_t numDoFs = funcA.getNumberOfGlobalDoFs( level );
      WALBERLA_LOG_INFO( "* dot product = " << dot );
      WALBERLA_CHECK_EQUAL( dot, real_c( 6 * numDoFs ) );
   }
}

void testEnumeration()
{
   hyteg::printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "Enumeration Test" );
   WALBERLA_LOG_INFO_ON_ROOT( "-> running 2D test" );
   {
      Point2D lowerLeft( real_c( 0.0 ), real_c( 0.0 ) );
      Point2D upperRight( real_c( 1.0 ), real_c( 1.0 ) );

      MeshInfo              mesh2D = MeshInfo::meshRectangle( lowerLeft, upperRight, MeshInfo::meshFlavour::CRISSCROSS, 1, 1 );
      SetupPrimitiveStorage setupStorage2D( mesh2D, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      std::shared_ptr< PrimitiveStorage > storage2D = std::make_shared< PrimitiveStorage >( setupStorage2D );

      const uint_t level = 2;

      P3Function< real_t > testFunc( "Test Function", storage2D, level, level );
      P3Function< idx_t >  numerator( "Enumerator", storage2D, level, level );

      numerator.enumerate( level );

      uint_t idxMin = numerator.getMinDoFValue( level );
      uint_t idxMax = numerator.getMaxDoFValue( level );
      WALBERLA_LOG_INFO( "* smallest DoF index = " << idxMin );
      WALBERLA_LOG_INFO( "* largest  DoF index = " << idxMax );
      WALBERLA_CHECK_EQUAL( idxMin, 0u );
      WALBERLA_CHECK_EQUAL( idxMax, 312u );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::testDoFCounting();
   hyteg::testInstantiation();
   hyteg::testInterpolation();
   hyteg::testAssignAndAdd();
   hyteg::testDotProduct();
   hyteg::testEnumeration();

   hyteg::printSeparator();
   return EXIT_SUCCESS;
}
