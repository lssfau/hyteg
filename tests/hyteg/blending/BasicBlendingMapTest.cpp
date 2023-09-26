/*
 * Copyright (c) 2017-2023 Marcus Mohr.
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

// This test checks that we can generate instances for each blending map,
// that querying for affineness and identity works correctly and that
// we can (de)serialize a map and obtain a new map of the same type

#include <core/Environment.h>

#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AffineMap3D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/IdentityMap.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/geometry/ThinShellMap.hpp"
#include "hyteg/geometry/TokamakMap.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

bool blendingMapIsIdentity( const GeometryMap& map )
{
   return dynamic_cast< const IdentityMap* >( &map ) != nullptr;
}

bool blendingMapIsAffine( const GeometryMap& map )
{
   return ( dynamic_cast< const IdentityMap* >( &map ) != nullptr ) ||
          ( dynamic_cast< const AffineMap2D* >( &map ) != nullptr ) || ( dynamic_cast< const AffineMap3D* >( &map ) != nullptr );
}

struct MapProperties
{
   bool              isIdentity;
   bool              isAffine;
   GeometryMap::Type type;
};

auto genMap( const std::string& variant )
{
   std::shared_ptr< GeometryMap > map;

   if ( variant == "IdentityMap" )
   {
      map = std::make_shared< IdentityMap >();
   }

   else if ( variant == "PolarCoordsMap" )
   {
      map = std::make_shared< PolarCoordsMap >();
   }

   else if ( variant == "AffineMap2D" )
   {
      map = std::make_shared< AffineMap2D >( Matrix2r(), Point2D() );
   }

   else if ( variant == "AffineMap3D" )
   {
      map = std::make_shared< AffineMap3D >( Matrix3r(), Point3D() );
   }

   else if ( variant == "AnnulusMap" )
   {
      MeshInfo              meshInfo = MeshInfo::meshAnnulus( real_c( 1 ), real_c( 2 ), MeshInfo::CRISS, 6, 1 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      Face                  face = *( setupStorage.getFaces().begin()->second );
      map                        = std::make_shared< AnnulusMap >( face );
   }

   else if ( variant == "CircularMap" )
   {
      MeshInfo meshInfo =
          MeshInfo::singleTriangle( { real_c( 0 ), real_c( 0 ) }, { real_c( -1 ), real_c( 0 ) }, { real_c( 1 ), real_c( 0 ) } );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      Face                  face = *( setupStorage.getFaces().begin()->second );
      Point3D               center{ real_c( 0 ), real_c( 0 ), real_c( 0 ) };
      map = std::make_shared< CircularMap >( face, setupStorage, center, real_c( 0.5 ) );
   }

   else if ( variant == "IcosahedralShellMap" )
   {
      MeshInfo              meshInfo = MeshInfo::meshSphericalShell( 2, 2, real_c( 1.0 ), real_c( 2.0 ) );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      Cell                  cell = *( setupStorage.getCells().begin()->second );
      map                        = std::make_shared< IcosahedralShellMap >( cell, setupStorage );
   }

   else if ( variant == "ThinShellMap" )
   {
      real_t                radius{ real_c( 2 ) };
      MeshInfo              meshInfo = MeshInfo::meshThinSphericalShell( 2, radius );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      Face                  face = *( setupStorage.getFaces().begin()->second );
      map                        = std::make_shared< ThinShellMap >( face, radius );
   }

   else if ( variant == "TokamakMap" )
   {
      const uint_t                toroidalResolution         = 3;
      const uint_t                poloidalResolution         = 2;
      const real_t                radiusOriginToCenterOfTube = real_c( 6.2 );
      const std::vector< real_t > tubeLayerRadii             = { real_c( 1.2 ), real_c( 2.2 ), real_c( 3 ) };
      const real_t                torodialStartAngle         = real_c( 0 );
      const real_t                polodialStartAngle         = real_c( 0 );

      real_t delta = std::sin( real_c( 0.33 ) );
      real_t r1    = real_c( 2.0 );
      real_t r2    = real_c( 3.7 );

      MeshInfo meshInfo = MeshInfo::meshTorus( toroidalResolution,
                                               poloidalResolution,
                                               radiusOriginToCenterOfTube,
                                               tubeLayerRadii,
                                               torodialStartAngle,
                                               polodialStartAngle );

      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      Cell                  cell = *( setupStorage.getCells().begin()->second );

      map = std::make_shared< TokamakMap >( cell,
                                            setupStorage,
                                            toroidalResolution,
                                            poloidalResolution,
                                            radiusOriginToCenterOfTube,
                                            tubeLayerRadii,
                                            torodialStartAngle,
                                            polodialStartAngle,
                                            delta,
                                            r1,
                                            r2 );
   }

   else
   {
      WALBERLA_ABORT( "Map variant '" << variant << "' not supported!" );
   }

   return map;
}

int main( int argc, char* argv[] )
{
   // basic setup
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // generate a list of maps to check
   std::vector< std::string > mapVariants{ "IdentityMap",
                                           "PolarCoordsMap",
                                           "AffineMap2D",
                                           "AffineMap3D",
                                           "AnnulusMap",
                                           "CircularMap",
                                           "IcosahedralShellMap",
                                           "ThinShellMap",
                                           "TokamakMap" };

   // run tests
   const auto badGuy = genMap( "TokamakMap" );

   for ( const auto& variant : mapVariants )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "*** Testing " << variant << " ***" );
      auto candidate = genMap( variant );

      // test property query methods
      WALBERLA_LOG_INFO_ON_ROOT( "isIdentity = " << std::boolalpha << candidate->isIdentity() );
      WALBERLA_CHECK_EQUAL( candidate->isIdentity(), blendingMapIsIdentity( *candidate ) );

      WALBERLA_LOG_INFO_ON_ROOT( "isAffine   = " << std::boolalpha << candidate->isAffine() );
      WALBERLA_CHECK_EQUAL( candidate->isAffine(), blendingMapIsAffine( *candidate ) );

      // test (de)serialisation
      if ( variant != "TokamakMap" )
      {
         walberla::mpi::SendBuffer sendBuffer;

         WALBERLA_LOG_INFO_ON_ROOT( "serializing map" );
         GeometryMap::serialize( candidate, sendBuffer );

         WALBERLA_LOG_INFO_ON_ROOT( "deserializing map" );
         walberla::mpi::RecvBuffer      recvBuffer( sendBuffer );
         std::shared_ptr< GeometryMap > clone = GeometryMap::deserialize( recvBuffer );

         // check that clone has identical type
         WALBERLA_LOG_INFO_ON_ROOT( "checking type of deserialized map" );
         const GeometryMap& candidateRef = *candidate;
         const GeometryMap& cloneRef     = *clone;
         const GeometryMap& badGuyRef    = *badGuy;
         WALBERLA_CHECK( typeid( candidateRef ) == typeid( cloneRef ) );
         WALBERLA_CHECK( !( typeid( cloneRef ) == typeid( badGuyRef ) ) ); // avoid problems with false positives
      }

      WALBERLA_LOG_INFO_ON_ROOT( "" );
   }

   return 0;
}
