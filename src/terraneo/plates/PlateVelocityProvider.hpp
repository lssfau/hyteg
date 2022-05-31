/*
 * Copyright (c) 2022 Berta Vilacis, Marcus Mohr.
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

#pragma once

#include "terraneo/dataimport/io.hpp"
#include "terraneo/plates/PlateRotationProvider.hpp"
#include "terraneo/plates/PlateStorage.hpp"
#include "terraneo/plates/SmoothingStrategies.hpp"

// preserve ordering of includes
#include "terraneo/plates/functionsForPlates.hpp"

namespace terraneo {
namespace plates {

/// API class for computation of velocities from plate reconstructions
class PlateVelocityProvider
{
 public:
   PlateVelocityProvider( std::string nameOfTopologiesFile, std::string nameOfRotationsFile )
   : plateTopologies_( nameOfTopologiesFile,
                       real_c( 1 ),
                       []( const std::string& filename ) { return terraneo::io::readJsonFile( filename ); } )
   , plateRotations_( nameOfRotationsFile,
                      []( const std::string& filename ) { return terraneo::io::readRotationsFile( filename ); } )
   {}

   vec3D getPointVelocity( const vec3D& point, const real_t age )
   {
      return getPointVelocity( point, age, LinearDistanceSmoother{0.015} );
   }

   uint_t findPlateID( const vec3D& point, const real_t age )
   {
      uint_t plateID{0};
      bool   plateFound{false};
      real_t distance{real_c( -1 )};

      std::tie( plateFound, plateID, distance ) = findPlateAndDistance( age, plateTopologies_, point );
      if ( plateFound )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Point found on plate with ID = " << plateID << ", distance to boundary = " << distance );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Point not found on any plate!" );
      }
      return plateID;
   }

   template < typename SmoothingStrategy >
   vec3D getPointVelocity( const vec3D& point, const real_t age, SmoothingStrategy computeSmoothing )
   {
      uint_t plateID{0};
      bool   plateFound{false};
      real_t distance{real_c( -1 )};

      std::tie( plateFound, plateID, distance ) = findPlateAndDistance( age, plateTopologies_, point );
      if ( plateFound )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Point found on plate with ID = " << plateID << ", distance to boundary = " << distance );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Point not found on any plate!" );
      }

      real_t smoothingFactor = computeSmoothing( distance );

      WALBERLA_LOG_INFO_ON_ROOT( "Smoothing Factor: " << smoothingFactor << "\n" );
      WALBERLA_LOG_INFO_ON_ROOT( "Plate ID: " << plateID << "\n" );
      if ( plateID == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "No plate ID assigned (Plate ID =0). Velocity set to (0.0,0.0,0.0)" );
         return {real_c( 0 ), real_c( 0 ), real_c( 0 )};
      }
      else
      {
         return computeCartesianVelocityVector( plateRotations_, plateID, age, point, smoothingFactor );
      }
   };

 private:
   PlateStorage          plateTopologies_;
   PlateRotationProvider plateRotations_;
};

} // namespace plates
} // namespace terraneo
