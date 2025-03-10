/*
 * Copyright (c) 2022-2025 Berta Vilacis, Marcus Mohr, Nils Kohl.
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

#include "terraneo/dataimport/FileIO.hpp"
#include "terraneo/helpers/conversions.hpp"
#include "terraneo/plates/PlateNotFoundHandlers.hpp"
#include "terraneo/plates/PlateRotationProvider.hpp"
#include "terraneo/plates/PlateStorage.hpp"
#include "terraneo/plates/SmoothingStrategies.hpp"

// preserve ordering of includes
#include "LocalAveragingPointWeightProvider.hpp"
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

   /// Returns the plate ID for a point at a given age stage
   ///
   /// This is basically an auxilliary function for testing plate detection and allows to
   /// generate data to visualise plate movement.
   uint_t findPlateID( const vec3D& point, const real_t age ) const
   {
      uint_t plateID{ idWhenNoPlateFound };
      bool   plateFound{ false };
      real_t distance{ real_c( -1 ) };

      // Transform the point to Lon, Lat, Radius - to preform all caculations
      // We use the Lon, Lat coordinates
      const vec3D pointLonLat = terraneo::conversions::cart2sph( point );

      std::tie( plateFound, plateID, distance ) = findPlateAndDistance( age, plateTopologies_, pointLonLat, idWhenNoPlateFound );
      return plateID;
   }

   /// Returns velocity vector for a point determined from the velocity of the associated plate at given age stage
   ///
   /// This is the convenience (non-expert) version of the method which uses a
   /// - LinearDistanceSmoother{ 0.015 }
   /// - DefaultPlateNotFoundHandler{}
   vec3D getPointVelocity( const vec3D& point, const real_t age )
   {
      return getPointVelocity( point, age, LinearDistanceSmoother{ 0.015 }, DefaultPlateNotFoundHandler{} );
   }

   /// Returns velocity vector for a point determined from the velocity of the associated plate at given age stage
   ///
   /// This is the expert version of the method which allows to explicitly set a SmoothingStrategy and
   /// a PlateNotFoundStrategy.
   template < typename SmoothingStrategy, typename PlateNotFoundStrategy >
   vec3D getPointVelocity( const vec3D&            point,
                           const real_t            age,
                           SmoothingStrategy       computeSmoothing,
                           PlateNotFoundStrategy&& errorHandler )
   {
      uint_t plateID{ 0 };
      bool   plateFound{ false };
      real_t distance{ real_c( -1 ) };

      // Transform the point to Lon, Lat, Radius - to preform all calculations
      // We use the Lon, Lat coordinates
      vec3D pointLonLat = terraneo::conversions::cart2sph( point );

      std::tie( plateFound, plateID, distance ) = findPlateAndDistance( age, plateTopologies_, pointLonLat, idWhenNoPlateFound );

      if ( !plateFound )
      {
         return errorHandler( point, age );
      }
      // else
      // {
      //    WALBERLA_LOG_DETAIL_ON_ROOT( "Point found on plate with ID = " << plateID << ", distance to boundary = " << distance );
      // }

      const real_t smoothingFactor = computeSmoothing( distance );

      WALBERLA_LOG_DETAIL_ON_ROOT( "Smoothing Factor: " << smoothingFactor << "\n" );
      WALBERLA_LOG_DETAIL_ON_ROOT( "Plate ID: " << plateID << "\n" );
      return computeCartesianVelocityVector( plateRotations_, plateID, age, pointLonLat, smoothingFactor );
   }

   /// Computes a weighted average of the velocity around the given point using the provided point and weights.
   ///
   /// Averaging is only applied if at least one of the provided points is located on a different plate.
   template < typename PlateNotFoundStrategy >
   vec3D getLocallyAveragedPointVelocity( const vec3D&                             point,
                                          const real_t                             age,
                                          const LocalAveragingPointWeightProvider& pointWeightProvider,
                                          PlateNotFoundStrategy&&                  errorHandler )
   {
      uint_t plateID{ 0 };
      bool   plateFound{ false };
      real_t distance{ real_c( -1 ) };

      // Transform the point to Lon, Lat, Radius - to perform all calculations
      // We use the Lon, Lat coordinates
      const vec3D pointLonLat = conversions::cart2sph( point );

      std::tie( plateFound, plateID, distance ) = findPlateAndDistance( age, plateTopologies_, pointLonLat, idWhenNoPlateFound );

      if ( !plateFound )
      {
         return errorHandler( point, age );
      }

      vec3D  avgVelCart( 0, 0, 0 );
      real_t weightSum = 0;

      if ( pointWeightProvider.maxDistance( pointLonLat ) < distance )
      {
         // We do not apply averaging since all points that would be used for averaging are on the same plate.
         WALBERLA_LOG_DETAIL_ON_ROOT( "No averaging." );
         WALBERLA_LOG_DETAIL_ON_ROOT( "Plate ID: " << plateID << "\n" );
         return computeCartesianVelocityVector( plateRotations_, plateID, age, pointLonLat, 1.0 );
      }

      const auto pointsAndWeights = pointWeightProvider.samplePointsAndWeightsLonLat( pointLonLat );

      // We average since at least some of the samples are possibly on at least one other plate.
      for ( const auto& [samplePointSphLonLat, weight] : pointsAndWeights )
      {
         uint_t avgPointPlateID{ 0 };
         bool   avgPointPlateFound{ false };
         real_t avgPointDistance{ real_c( -1 ) };

         std::tie( avgPointPlateFound, avgPointPlateID, avgPointDistance ) =
             findPlateAndDistance( age, plateTopologies_, samplePointSphLonLat, idWhenNoPlateFound );

         if ( avgPointPlateFound )
         {
            // This is possibly slightly inaccurate since we are averaging over the cartesian velocity vectors and then projecting
            // out the normal component. It would be better to average in the "lonlat-space" and then convert and return the
            // cartesian vector. On the other hand, averaging the plate velocities is already a somewhat arbitrary and physically
            // meaningless approximation in the first place, so this might just work.
            avgVelCart +=
                weight * computeCartesianVelocityVector( plateRotations_, avgPointPlateID, age, samplePointSphLonLat, 1.0 );
            weightSum += weight;
         }
      }

      avgVelCart /= weightSum;

      const auto n                = point.normalized();
      const auto dot              = n.dot( avgVelCart );
      const auto normalComponent  = dot * n;
      const auto tangentComponent = avgVelCart - normalComponent;

      return tangentComponent;
   }

   /// Query function to obtain a vector of plate stages available in the datafiles
   const std::vector< real_t >& getListOfPlateStages() const { return plateTopologies_.getListOfPlateStages(); }

   /// Plate ID to be used when no associated plate was found for a point
   const uint_t idWhenNoPlateFound{ 0 };

 private:
   PlateStorage          plateTopologies_;
   PlateRotationProvider plateRotations_;
};

} // namespace plates
} // namespace terraneo
