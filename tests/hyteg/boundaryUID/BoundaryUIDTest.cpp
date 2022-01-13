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
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
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
void runTest()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running BoundaryUIDTest for " << FunctionTrait< func_t >::getTypeName() );

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
   real_t tol         = real_c( 0.5 ) * ( outerRad - innerRad ) / real_c( nLayers );

   MeshInfo meshInfo = MeshInfo::meshAnnulus( innerRad, outerRad, 0.0, 2.0 * pi, MeshInfo::CROSS, 8, nLayers );

   // flag the inner and outer boundary by assigning different values
   auto onBoundary = [&boundaryRad, tol]( const Point3D& x ) {
      real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      return std::abs( boundaryRad - radius ) < tol;
   };

   boundaryRad = outerRad;
   meshInfo.setMeshBoundaryFlagsByVertexLocation( markerOuterBoundary, onBoundary, true );

   boundaryRad = innerRad;
   meshInfo.setMeshBoundaryFlagsByVertexLocation( markerInnerBoundary, onBoundary, true );

   // generate the setupStorage and associate blending map
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // -----------------------
   //  Function Manipulation
   // -----------------------
   uint_t minLevel = 2;
   uint_t maxLevel = 3;
   func_t func( "Test Func #1", storage, minLevel, maxLevel );

   // generate bc object and set different conditions on inside and outside
   BoundaryCondition bcs;
   BoundaryUID       outerBC = bcs.createDirichletBC( "Dirichlet on outer radius", markerOuterBoundary );
   BoundaryUID       innerBC = bcs.createDirichletBC( "Dirichlet on inner radius", markerInnerBoundary );

   real_t iValue = real_c( 30 );  // on inner boundary
   real_t mValue = real_c( 20 );  // in the interior
   real_t oValue = real_c( 10 );  // on outer boundary

   std::function< real_t( const Point3D& ) > DirichletInner = [ iValue ]( const Point3D& ) { return iValue; };
   std::function< real_t( const Point3D& ) > DirichletOuter = [ oValue ]( const Point3D& ) { return oValue; };

   func.setBoundaryCondition( bcs );

   // assign functions values
   func.interpolate( mValue, maxLevel, All );
   if constexpr ( !std::is_base_of< CSFVectorFunction< func_t >, func_t >::value )
   {
      func.interpolate( DirichletInner, maxLevel, innerBC );
      func.interpolate( DirichletOuter, maxLevel, outerBC );

      // Would fail do to missing implementation for 1st argument being constant:
      // func.interpolate( iValue, maxLevel, innerBC );
      // func.interpolate( oValue, maxLevel, innerBC );
   }
   else
   {
      func.interpolate( {DirichletInner, DirichletInner}, maxLevel, innerBC );
      func.interpolate( {DirichletOuter, DirichletOuter}, maxLevel, outerBC );
   }

   // ------------------
   //  Control Function
   // ------------------
   std::function< real_t( const Point3D& ) > controlValues = [innerRad, outerRad, iValue, mValue, oValue]( const Point3D& x ) {
      real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t mytol  = 1e-14;
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

   func_t ctrl( "Test Func #2", storage, minLevel, maxLevel );
   func_t diff( "Test Func #3", storage, minLevel, maxLevel );
   ctrl.interpolate( controlValues, maxLevel, All );
   diff.assign( {-1.0, 1.0}, {func, ctrl}, maxLevel, All );
   real_t check = real_c( 0 );
   if constexpr ( std::is_base_of< CSFVectorFunction< func_t >, func_t >::value )
   {
      check = diff.getMaxComponentMagnitude( maxLevel, All );
   }
   else
   {
      check = diff.getMaxMagnitude( maxLevel, All );
   }
   WALBERLA_CHECK_FLOAT_EQUAL( check, real_c( 0 ) );

   // ------------
   //  Output VTK
   // ------------
   bool beVerbose = false;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "boundaryUIDTest";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( func );
      vtkOutput.add( ctrl );
      vtkOutput.add( diff );
      vtkOutput.write( maxLevel );
   }
}

int main( int argc, char* argv[] )
{
   // -------------
   //  Basic Setup
   // -------------
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // ------------------------------------------
   //  Run test with different function classes
   // ------------------------------------------
   runTest< P1Function< real_t > >();
   runTest< P2Function< real_t > >();
   runTest< P1VectorFunction< real_t > >();
   runTest< P2VectorFunction< real_t > >();

   return 0;
}
