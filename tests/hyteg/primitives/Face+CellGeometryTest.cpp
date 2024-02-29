/*
* Copyright (c) 2024 Marcus Mohr.
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

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;
using walberla::uint_c;

using namespace hyteg;

void runFaceTest()
{
   // define a given triangle by its vertices and specify its properties
   Point2D a( { real_c( -1.5 ), real_c( +1.0 ) } );
   Point2D b( { real_c( +2.3 ), real_c( +3.0 ) } );
   Point2D c( { real_c( +4.0 ), real_c( -0.9 ) } );

   const real_t checkValArea   = real_c( 9.109999999999999 );
   const real_t checkValRadius = real_c( 1.268137586698838e+00 );

   // turn triangle into a PrimitiveStorage object
   MeshInfo meshInfo = MeshInfo::singleTriangle( a, b, c );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   PrimitiveStorage      storage( setupStorage );

   // extract single face and check methods that return geometric properties
   for ( const auto& faceIter : storage.getFaces() )
   {
      auto face = faceIter.second;

      real_t area = face->getArea();
      WALBERLA_LOG_INFO_ON_ROOT( "-> area = " << area );
      WALBERLA_CHECK_FLOAT_EQUAL( area, checkValArea );

      real_t radius = face->getIncircleRadius();
      WALBERLA_LOG_INFO_ON_ROOT( "-> incircle radius = " << radius );
      WALBERLA_CHECK_FLOAT_EQUAL( radius, checkValRadius );
   };
}

void runCellTest()
{
   WALBERLA_LOG_INFO_ON_ROOT( "using a regular tetrahedron" );
   {
      // define vertices of a regular tetrahedron and specify its properties
      real_t aux1 = real_c( std::sqrt( 3.0 ) / 2.0 );
      real_t aux2 = real_c( std::sqrt( 2.0 / 3.0 ) );

      Point3D a( { real_c( -0.5 ), real_c( 0 ), real_c( 0 ) } );
      Point3D b( { real_c( +0.5 ), real_c( 0 ), real_c( 0 ) } );
      Point3D c( { real_c( 0 ), aux1, real_c( 0 ) } );
      Point3D d( { real_c( 0 ), aux1 / real_c( 3 ), aux2 } );

      const real_t checkValVolume = std::sqrt( real_c( 2 ) ) / real_c( 12 );
      const real_t checkValRadius = std::sqrt( real_c( 6 ) ) / real_c( 12 );

      // turn tetrahedron into a PrimitiveStorage object
      MeshInfo meshInfo = MeshInfo::singleTetrahedron( { a, b, c, d } );

      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      PrimitiveStorage      storage( setupStorage );

      // extract single cell and check methods that return geometric properties
      for ( const auto& cellIter : storage.getCells() )
      {
         auto cell = cellIter.second;

         real_t volume = cell->getVolume();
         WALBERLA_LOG_INFO_ON_ROOT( "-> volume = " << volume );
         WALBERLA_CHECK_FLOAT_EQUAL( volume, checkValVolume );

         real_t radius = cell->getInsphereRadius();
         WALBERLA_LOG_INFO_ON_ROOT( "-> insphere radius = " << radius );
         WALBERLA_CHECK_FLOAT_EQUAL( radius, checkValRadius );
      };
   }

   WALBERLA_LOG_INFO_ON_ROOT( "using an arbitrary tetrahedron" );
   {
      // define vertices of a regular tetrahedron and specify its properties
      Point3D a( { real_c( 0.00 ), real_c( 0.0 ), real_c( +0.0 ) } );
      Point3D b( { real_c( 2.00 ), real_c( 0.0 ), real_c( +0.0 ) } );
      Point3D c( { real_c( 1.75 ), real_c( 1.5 ), real_c( -1.0 ) } );
      Point3D d( { real_c( 0.20 ), real_c( 0.3 ), real_c( +1.0 ) } );

      const real_t checkValVolume = real_c( 6.000000000000000e-01 );
      const real_t checkValRadius = real_c( 2.964715904083913e-01 );

      // turn tetrahedron into a PrimitiveStorage object
      MeshInfo meshInfo = MeshInfo::singleTetrahedron( { a, b, c, d } );

      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      PrimitiveStorage      storage( setupStorage );

      // extract single cell and check methods that return geometric properties
      for ( const auto& cellIter : storage.getCells() )
      {
         auto cell = cellIter.second;

         real_t volume = cell->getVolume();
         WALBERLA_LOG_INFO_ON_ROOT( "-> volume = " << volume );
         WALBERLA_CHECK_FLOAT_EQUAL( volume, checkValVolume );

         real_t radius = cell->getInsphereRadius();
         WALBERLA_LOG_INFO_ON_ROOT( "-> insphere radius = " << radius );
         WALBERLA_CHECK_FLOAT_EQUAL( radius, checkValRadius );
      };
   }
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "Testing geometric property methods of Face:" );
   runFaceTest();

   WALBERLA_LOG_INFO_ON_ROOT( "Testing geometric property methods of Cell:" );
   runCellTest();

   return EXIT_SUCCESS;
}
