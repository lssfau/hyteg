/*
 * Copyright (c) 2021 Marcus Mohr.
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

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

// Using an annulus mesh we test that setting of BoundaryUIDs and
// that we can use them in the interpolation methods to set
// Dirichlet BC values.

template < typename func_t >
void runTest( bool useCentroids )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running BoundaryUIDTest for " << FunctionTrait< func_t >::getTypeName() );

   bool beVerbose = false;

   // -----------------------------------------
   //  Define markers for geometric boundaries
   // -----------------------------------------
   uint_t markerInnerBoundary = 11;
   uint_t markerOuterBoundary = 22;

   // --------------
   //  SetupStorage
   // --------------
   real_t innerRad    = 1.0;
   real_t outerRad    = 2.0;
   uint_t nLayers     = 2;
   real_t boundaryRad = 0.0;
   real_t tol         = real_c( 0.1 ) * ( outerRad - innerRad ) / real_c( nLayers );

   MeshInfo meshInfo = MeshInfo::meshAnnulus( innerRad, outerRad, MeshInfo::CROSS, 8, nLayers );

   // generate the setupStorage and associate blending map
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( setupStorage );

   // flag the inner and outer boundary by assigning different values
   auto onBoundary = [&boundaryRad, tol]( const Point3D& x ) {
      real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      return std::abs( boundaryRad - radius ) < tol;
   };

   if ( !useCentroids )
   {
      boundaryRad = outerRad;
      setupStorage.setMeshBoundaryFlagsByVertexLocation( markerOuterBoundary, onBoundary, true );

      boundaryRad = innerRad;
      setupStorage.setMeshBoundaryFlagsByVertexLocation( markerInnerBoundary, onBoundary, true );
   }
   else
   {
      boundaryRad = outerRad;
      setupStorage.setMeshBoundaryFlagsByCentroidLocation( markerOuterBoundary, onBoundary, true );

      boundaryRad = innerRad;
      setupStorage.setMeshBoundaryFlagsByCentroidLocation( markerInnerBoundary, onBoundary, true );
   }

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // report primitives and their flags
   if ( beVerbose )
   {
      std::stringstream sStr;
      setupStorage.toStream( sStr, true );
      WALBERLA_LOG_INFO_ON_ROOT( "" << sStr.str() );
   }

   // -----------------------
   //  Function Manipulation
   // -----------------------
   uint_t minLevel = 2;
   uint_t maxLevel = 3;
   func_t func1( "Test Func #1", storage, minLevel, maxLevel );
   func_t func2( "Test Func #2", storage, minLevel, maxLevel );

   // generate bc object and set different conditions on inside and outside
   BoundaryCondition bcs;
   BoundaryUID       outerBC = bcs.createDirichletBC( "Dirichlet on outer radius", markerOuterBoundary );
   BoundaryUID       innerBC = bcs.createDirichletBC( "Dirichlet on inner radius", markerInnerBoundary );

   real_t iValue = real_c( 30 ); // on inner boundary
   real_t mValue = real_c( 20 ); // in the interior
   real_t oValue = real_c( 10 ); // on outer boundary

   std::function< real_t( const Point3D& ) > DirichletInner = [iValue]( const Point3D& ) { return iValue; };
   std::function< real_t( const Point3D& ) > DirichletOuter = [oValue]( const Point3D& ) { return oValue; };

   func1.setBoundaryCondition( bcs );
   func2.setBoundaryCondition( bcs );

   // assign functions values
   func1.interpolate( mValue, maxLevel, All );
   func2.interpolate( mValue, maxLevel, All );
   if constexpr ( !std::is_base_of< CSFVectorFunction< func_t >, func_t >::value )
   {
      // use an expression
      func1.interpolate( DirichletInner, maxLevel, innerBC );
      func1.interpolate( DirichletOuter, maxLevel, outerBC );

      // use constant value
      func2.interpolate( iValue, maxLevel, innerBC );
      func2.interpolate( oValue, maxLevel, outerBC );
   }
   else
   {
      // use an expression
      func1.interpolate( { DirichletInner, DirichletInner }, maxLevel, innerBC );
      func1.interpolate( { DirichletOuter, DirichletOuter }, maxLevel, outerBC );

      // use constant value
      func2.interpolate( { iValue, iValue }, maxLevel, innerBC );
      func2.interpolate( { oValue, oValue }, maxLevel, outerBC );
   }

   // ------------------
   //  Control Function
   // ------------------
   std::function< real_t( const Point3D& ) > controlValues = [innerRad, outerRad, iValue, mValue, oValue]( const Point3D& x ) {
      real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t mytol  = real_c( std::is_same< real_t, double >() ? 1e-14 : 1e-6 );
      if ( std::abs( innerRad - radius ) < mytol )
      {
         return iValue;
      }
      else if ( std::abs( outerRad - radius ) < mytol )
      {
         return oValue;
      }
      return mValue;
   };

   func_t ctrl( "Test Func #3", storage, minLevel, maxLevel );
   func_t diff( "Test Func #4", storage, minLevel, maxLevel );
   for ( uint_t k = 1; k <= 2; ++k )
   {
      ctrl.interpolate( controlValues, maxLevel, All );
      if ( k == 1 )
      {
         diff.assign( { -1.0, 1.0 }, { func1, ctrl }, maxLevel, All );
      }
      else
      {
         diff.assign( { -1.0, 1.0 }, { func2, ctrl }, maxLevel, All );
      }
      real_t check = real_c( 0 );
      if constexpr ( std::is_base_of< CSFVectorFunction< func_t >, func_t >::value )
      {
         check = diff.getMaxComponentMagnitude( maxLevel, All );
      }
      else
      {
         check = diff.getMaxDoFMagnitude( maxLevel, All );
      }
      // WALBERLA_LOG_INFO_ON_ROOT( "k = " << k << ", check = " << check );
      WALBERLA_CHECK_FLOAT_EQUAL( check, real_c( 0 ) );
   }

   // ------------
   //  Output VTK
   // ------------
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "boundaryUIDTest";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( func1 );
      vtkOutput.add( func2 );
      vtkOutput.add( ctrl );
      vtkOutput.add( diff );
      vtkOutput.write( maxLevel );
   }
}

void captureTheFlags( bool useCentroids )
{
   WALBERLA_LOG_DETAIL_ON_ROOT( "Playing 'capture the flags'" );

   // -----------------------------------------
   //  Define markers for geometric boundaries
   // -----------------------------------------
   uint_t markerInnerBoundary = 11;
   uint_t markerOuterBoundary = 22;

   // --------------
   //  SetupStorage
   // --------------
   real_t innerRad    = 1.0;
   real_t outerRad    = 2.0;
   uint_t nLayers     = 2;
   real_t boundaryRad = 0.0;
   real_t tol         = real_c( 0.1 ) * ( outerRad - innerRad ) / real_c( nLayers );

   MeshInfo meshInfo = MeshInfo::meshAnnulus( innerRad, outerRad, MeshInfo::CROSS, 8, nLayers );

   // generate the setupStorage and associate blending map
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( setupStorage );

   // flag the inner and outer boundary by assigning different values
   auto onBoundary = [&boundaryRad, tol]( const Point3D& x ) {
      real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      return std::abs( boundaryRad - radius ) < tol;
   };

   if ( !useCentroids )
   {
      boundaryRad = outerRad;
      setupStorage.setMeshBoundaryFlagsByVertexLocation( markerOuterBoundary, onBoundary, true );

      boundaryRad = innerRad;
      setupStorage.setMeshBoundaryFlagsByVertexLocation( markerInnerBoundary, onBoundary, true );
   }
   else
   {
      boundaryRad = outerRad;
      setupStorage.setMeshBoundaryFlagsByCentroidLocation( markerOuterBoundary, onBoundary, true );

      boundaryRad = innerRad;
      setupStorage.setMeshBoundaryFlagsByCentroidLocation( markerInnerBoundary, onBoundary, true );
   }

   std::stringstream sStr;
   setupStorage.toStream( sStr, true );
   WALBERLA_LOG_DETAIL_ON_ROOT( "Here we go:\n" << sStr.str() );
}

// we mesh a rectangle and use a PolarCoordsMap to map it to a half annulus
template < typename func_t >
void centroidHardBlendingTest()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running BoundaryUIDTest::centroidHardBlendingTest for "
                              << FunctionTrait< func_t >::getTypeName() );

   bool beVerbose = false;

   // -----------------------------------------
   //  Define markers for geometric boundaries
   // -----------------------------------------
   uint_t markerInnerBoundary = 33;
   uint_t markerOuterBoundary = 66;
   uint_t markerSideBoundary  = 99;

   // --------------
   //  SetupStorage
   // --------------
   real_t innerRad    = real_c( 1 );
   real_t outerRad    = real_c( 2 );
   real_t phiMin      = real_c( 0 );
   real_t phiMax      = real_c( pi );
   real_t boundaryRad = 0.0;
   real_t tol         = real_c( 1e-5 );

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( innerRad, phiMin ), Point2D( outerRad, phiMax ), MeshInfo::CRISS, 3, 4 );
   // MeshInfo::meshRectangle( Point2D( innerRad, phiMin ), Point2D( outerRad, phiMax ), MeshInfo::DIAMOND, 3, 2 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   PolarCoordsMap::setMap( setupStorage );

   // ---------------
   //  BoundaryFlags
   // ---------------

   // oracle for curved boundaries
   auto onCurvedBoundary = [&boundaryRad, tol]( const Point3D& x ) {
      real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      return std::abs( boundaryRad - radius ) < tol;
   };

   // oracle for straight boundaries
   auto onStraightBoundary = [tol]( const Point3D& x ) { return std::abs( x[1] ) < tol; };

   setupStorage.setMeshBoundaryFlagsByCentroidLocation( markerSideBoundary, onStraightBoundary, true );

   boundaryRad = outerRad;
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( markerOuterBoundary, onCurvedBoundary, true );

   boundaryRad = innerRad;
   setupStorage.setMeshBoundaryFlagsByCentroidLocation( markerInnerBoundary, onCurvedBoundary, true );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // report primitives and their flags
   if ( beVerbose )
   {
      std::stringstream sStr;
      setupStorage.toStream( sStr, true );
      WALBERLA_LOG_INFO_ON_ROOT( "" << sStr.str() );
   }

   // -----------------------
   //  Function Manipulation
   // -----------------------
   uint_t minLevel = 2;
   uint_t maxLevel = 2;
   func_t test( "Test Func", storage, minLevel, maxLevel );
   func_t ctrl( "Ctrl Func", storage, minLevel, maxLevel );

   // generate bc object and set different conditions on inner, outer, and straight boundaries
   BoundaryCondition bcs;
   BoundaryUID       outerBC = bcs.createDirichletBC( "Dirichlet on outer radius", markerOuterBoundary );
   BoundaryUID       innerBC = bcs.createDirichletBC( "Dirichlet on inner radius", markerInnerBoundary );
   BoundaryUID       sideBC  = bcs.createDirichletBC( "Dirichlet on side boundary radius", markerSideBoundary );

   test.setBoundaryCondition( bcs );
   ctrl.setBoundaryCondition( bcs );

   // assign functions values
   real_t iValue = real_c( 30 ); // on inner boundary
   real_t mValue = real_c( 20 ); // in the interior
   real_t oValue = real_c( 10 ); // on outer boundary
   real_t sValue = real_c( -4 ); // on side boundary

   test.interpolate( mValue, maxLevel, All );
   test.interpolate( iValue, maxLevel, innerBC );
   test.interpolate( oValue, maxLevel, outerBC );
   test.interpolate( sValue, maxLevel, sideBC );

   // ------------------
   //  Control Function
   // ------------------
   std::function< real_t( const Point3D& ) > controlValues =
       [innerRad, outerRad, iValue, mValue, oValue, sValue]( const Point3D& x ) {
          real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] );
          real_t mytol  = real_c( std::is_same< real_t, double >() ? 1e-14 : 1e-6 );
          if ( std::abs( innerRad - radius ) < mytol )
          {
             return iValue;
          }
          else if ( std::abs( outerRad - radius ) < mytol )
          {
             return oValue;
          }
          else if ( std::abs( x[1] ) < mytol )
          {
             return sValue;
          }
          return mValue;
       };
   ctrl.interpolate( controlValues, maxLevel, All );

   // ------------
   //  Output VTK
   // ------------
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "centroidTest";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( test );
      vtkOutput.add( ctrl );
      vtkOutput.write( maxLevel );
   }

   // -----------------------
   //  check for differences
   // -----------------------
   func_t diff( "Diff Func", storage, minLevel, maxLevel );
   diff.assign( { -1.0, 1.0 }, { test, ctrl }, maxLevel, All );
   real_t check = diff.getMaxDoFMagnitude( maxLevel, All );
   WALBERLA_CHECK_FLOAT_EQUAL( check, real_c( 0 ) );
}

int main( int argc, char* argv[] )
{
   // -------------
   //  Basic Setup
   // -------------
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // -----------------------------------------------------------
   //  Run test with different function classes and flag setters
   // -----------------------------------------------------------

   bool useCentroids = false;
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "setupStorage::setMeshBoundaryFlagsByVertexLocation" );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   captureTheFlags( useCentroids );
   runTest< P1Function< real_t > >( useCentroids );
   runTest< P2Function< real_t > >( useCentroids );
   runTest< P1VectorFunction< real_t > >( useCentroids );
   runTest< P2VectorFunction< real_t > >( useCentroids );

   useCentroids = true;
   WALBERLA_LOG_INFO_ON_ROOT( "\n\n----------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "setupStorage::setMeshBoundaryFlagsByCentroidLocation" );
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------------------------------" );
   captureTheFlags( useCentroids );
   runTest< P1Function< real_t > >( useCentroids );
   runTest< P2Function< real_t > >( useCentroids );
   runTest< P1VectorFunction< real_t > >( useCentroids );
   runTest< P2VectorFunction< real_t > >( useCentroids );

   // use a strongly blended geometry
   centroidHardBlendingTest< P1Function< real_t > >();
   centroidHardBlendingTest< P2Function< real_t > >();

   return 0;
}
