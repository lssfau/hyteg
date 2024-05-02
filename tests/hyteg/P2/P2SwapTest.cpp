/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
#include <functional>
#include <vector>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitives/PrimitiveID.hpp"
#include "hyteg/primitives/all.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;

namespace hyteg {

static void testP2Swap()
{
   const uint_t minLevel  = 2;
   const uint_t maxLevel  = 4;
   const uint_t testLevel = 3;

   MeshInfo mesh = MeshInfo::fromGmshFile( "../../meshes/quad_8el.msh" );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   auto                                timingTree = std::make_shared< walberla::WcTimingTree >();
   std::shared_ptr< PrimitiveStorage > storage    = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   P2Function< real_t > x( "x", storage, minLevel, maxLevel );
   P2Function< real_t > y( "y", storage, minLevel, maxLevel );
   P2Function< real_t > xCorrect( "xCorrect", storage, minLevel, maxLevel );
   P2Function< real_t > yCorrect( "yCorrect", storage, minLevel, maxLevel );
   P2Function< real_t > err( "err", storage, minLevel, maxLevel );

   // Interpolate

   std::function< real_t( const hyteg::Point3D& ) > funcX = []( const Point3D& xx ) -> real_t {
      return real_c( ( 1.0 + std::sin( xx[0] ) ) * ( 2.0 + xx[1] ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > funcY = []( const Point3D& xx ) -> real_t {
      return real_c( ( 1.0 + ( xx[0] / 5.0 ) ) * ( 42.0 + xx[1] ) );
   };

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      x.interpolate( funcX, level );
      y.interpolate( funcY, level );
      if ( level == testLevel )
      {
         xCorrect.interpolate( funcY, level );
         yCorrect.interpolate( funcX, level );
      }
      else
      {
         xCorrect.interpolate( funcX, level );
         yCorrect.interpolate( funcY, level );
      }
   }

   x.swap( y, testLevel );

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      err.assign( { 1.0, -1.0 }, { x, xCorrect }, level );
      const real_t errorX = err.dotGlobal( err, level );
      WALBERLA_LOG_DEVEL_ON_ROOT( errorX )
      WALBERLA_CHECK_LESS( errorX, 1e-16 );
      err.assign( { 1.0, -1.0 }, { y, yCorrect }, level );
      const real_t errorY = err.dotGlobal( err, level );
      WALBERLA_LOG_DEVEL_ON_ROOT( errorY )
      WALBERLA_CHECK_LESS( errorY, 1e-16 );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testP2Swap();

   return EXIT_SUCCESS;
}
