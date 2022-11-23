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
#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/timing/Timer.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

void usage()
{
   WALBERLA_LOG_INFO_ON_ROOT( "---------------\n 3DSurfaceDemo\n---------------\n"
                              << "This app demonstrates that in HyTeG we can\n"
                              << " - setup a surface mesh in 3D\n"
                              << " - use the resulting PrimitiveStorage for a Function object\n"
                              << " - interpolate an expression into that FE space\n"
                              << " - output the result for visualisation\n" );
}

int main( int argc, char* argv[] )
{
   // ---------------
   //  General Setup
   // ---------------

   // Setup enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Be talkative
   usage();

   // Generate a square mesh
   WALBERLA_LOG_INFO_ON_ROOT( " *** Generating square coarse mesh" );
   Point2D cornerLL( { real_c( -1 ), real_c( -1 ) } );
   Point2D cornerUR( { real_c( +1 ), real_c( +1 ) } );

   MeshInfo meshInfo = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CROSS, 2, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( " *** Modifying z-coordinates" );
   meshInfo.applyCoordinateMap( []( const Point3D& p ) {
      Point3D newPoint{ { p[0], p[1], p[0] * p[0] + real_c( 0.1 ) * p[1] } };
      return newPoint;
   } );

   WALBERLA_LOG_INFO_ON_ROOT( " *** Creating SetupPrimitiveStorage" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   WALBERLA_LOG_INFO_ON_ROOT( " *** Creating PrimitiveStorage" );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_LOG_INFO_ON_ROOT( " *** Performing FE interpolation" );
   uint_t                                    maxLevel = 3;
   std::function< real_t( const Point3D& ) > func     = []( const Point3D& x ) {
      // return real_c( 2 ) * x[0] * x[0] + real_c( -1 ) * x[1] + std::sin( x[2] );
      return std::sin( real_c( 10 ) * x[2] );
   };
   P2Function< real_t > u( "u", storage, maxLevel, maxLevel );
   u.interpolate( func, maxLevel );

   // output data for visualisation
   bool outputVTK = true;
   if ( outputVTK )
   {
      VTKOutput vtkOutput( "./output", "3DSurfaceDemo", storage );
      vtkOutput.add( u );
      vtkOutput.write( maxLevel );
   }

   return 0;
}
