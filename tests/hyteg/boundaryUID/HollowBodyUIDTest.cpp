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
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

// This test checks that the setting of the boundary flags for
// hollow bodies (full annulus and thick spherical shell) by
// the corresponding mesh generators works correctly and that
// the information is correctly passed through the primitive
// storage generation process.

template < typename func_t >
void runTest( std::shared_ptr< PrimitiveStorage >& storage, uint_t level, real_t innerRad, real_t outerRad, std::string fName )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Checking for " << FunctionTrait< func_t >::getTypeName() );

   bool beVerbose = false;

   // -----------------------
   //  Function Manipulation
   // -----------------------
   func_t func1( "Test Func #1", storage, level, level );
   func_t func2( "Test Func #2", storage, level, level );

   // generate bc object and set different conditions on inside and outside
   BoundaryCondition bcs;
   BoundaryUID       outerBC = bcs.createDirichletBC( "Dirichlet on outer radius", MeshInfo::flagOuterBoundary );
   BoundaryUID       innerBC = bcs.createDirichletBC( "Dirichlet on inner radius", MeshInfo::flagInnerBoundary );

   real_t iValue = real_c( 30 ); // on inner boundary
   real_t mValue = real_c( 20 ); // in the interior
   real_t oValue = real_c( 10 ); // on outer boundary

   func1.setBoundaryCondition( bcs );

   // assign functions values
   func1.interpolate( mValue, level, All );
   if constexpr ( !std::is_base_of< CSFVectorFunction< func_t >, func_t >::value )
   {
      func1.interpolate( iValue, level, innerBC );
      func1.interpolate( oValue, level, outerBC );
   }
   else
   {
      func1.interpolate( { iValue, iValue }, level, innerBC );
      func1.interpolate( { oValue, oValue }, level, outerBC );
   }

   func2.assign( { real_c( 1 ) }, { func1 }, level, All );

   // ------------------
   //  Control Function
   // ------------------
   std::function< real_t( const Point3D& ) > controlValues = [innerRad, outerRad, iValue, mValue, oValue]( const Point3D& x ) {
      real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
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

   func_t ctrl( "Test Func #3", storage, level, level );
   func_t diff( "Test Func #4", storage, level, level );
   for ( uint_t k = 1; k <= 2; ++k )
   {
      ctrl.interpolate( controlValues, level, All );
      if ( k == 1 )
      {
         diff.assign( { -1.0, 1.0 }, { func1, ctrl }, level, All );
      }
      else
      {
         diff.assign( { -1.0, 1.0 }, { func2, ctrl }, level, All );
      }
      real_t check = real_c( 0 );
      if constexpr ( std::is_base_of< CSFVectorFunction< func_t >, func_t >::value )
      {
         check = diff.getMaxComponentMagnitude( level, All );
      }
      else
      {
         check = diff.getMaxDoFMagnitude( level, All );
      }
      WALBERLA_LOG_INFO_ON_ROOT( "k = " << k << ", check = " << check );
      WALBERLA_CHECK_FLOAT_EQUAL( check, real_c( 0 ) );
   }

   // ------------
   //  Output VTK
   // ------------
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( func1 );
      vtkOutput.add( func2 );
      vtkOutput.add( ctrl );
      vtkOutput.add( diff );
      vtkOutput.write( level );
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

   real_t rmin{ real_c( 1.0 ) };
   real_t rmax{ real_c( 2.5 ) };

   // ----------------------
   //  Run test for ANNULUS
   // ----------------------

   std::string tag{ "HollowBodyUIDTest_" };
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Running HollowBodyUIDTest for an Annulus" );

      MeshInfo              meshInfo = MeshInfo::meshAnnulus( rmin, rmax, MeshInfo::CROSS, 8, 2 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      AnnulusMap::setMap( setupStorage );
      auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

      runTest< P1Function< real_t > >( storage, 2, rmin, rmax, tag + "AnnulusP1" );
      runTest< P2Function< real_t > >( storage, 2, rmin, rmax, tag + "AnnulusP2" );
   }

   // --------------------
   //  Run test for SHELL
   // --------------------

   {
      WALBERLA_LOG_INFO_ON_ROOT( "Running HollowBodyUIDTest for an Thick Spherical Shell" );

      uint_t                ntan     = 3;
      uint_t                nrad     = 3;
      MeshInfo              meshInfo = MeshInfo::meshSphericalShell( ntan, nrad, rmin, rmax );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      IcosahedralShellMap::setMap( setupStorage );
      auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

      runTest< P1Function< real_t > >( storage, 2, rmin, rmax, tag + "ShellP1" );
      runTest< P2Function< real_t > >( storage, 2, rmin, rmax, tag + "ShellP2" );
   }

   return 0;
}
