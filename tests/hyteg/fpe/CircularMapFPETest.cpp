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

// check for floating-point exceptions in CircularMap()

#include <cfenv>
#include <core/Environment.h>

#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;


int main( int argc, char** argv )
{
#ifndef __APPLE__
   // should work with Intel, GCC, Clang and even MSVC compiler /nope not MSVC
   #ifndef _MSC_VER
      feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
   #endif
#endif
   // environment stuff
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Set mesh and primitives
   MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   Point3D circleCenter{ { -0.5, 0.5, 0 } };
   real_t  circleRadius = 0.25;
   Face& face = *setupStorage.getFace( 14 );
   CircularMap myMap( face, setupStorage, circleCenter, circleRadius );

   std::array< Point3D, 3 >coords = face.getCoordinates();
   for( uint_t k = 0; k < 3; k++ ) {
     WALBERLA_LOG_INFO_ON_ROOT( "vertex #" << k << " has coords: " << coords[k] );
   }

   // works
   WALBERLA_LOG_INFO_ON_ROOT( "Testing with point on edge" );
   Matrix2r jacMat;
   Point3D myVertex1( { 0.0, 0.95, 0.0 } );
   myMap.evalDF( myVertex1, jacMat );

   // fails
   WALBERLA_LOG_INFO_ON_ROOT( "Testing with location of upper vertex" );
   Point3D myVertex2( { 0.0, 1.0, 0.0 } );
   myMap.evalDF( myVertex2, jacMat );

   WALBERLA_LOG_INFO_ON_ROOT( "* Check! No floating point exceptions detected!" );

   return EXIT_SUCCESS;
}
