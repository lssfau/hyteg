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

#include <limits>

// CGAL: compilation with -lCGAL -lmpfr -lgmp -DCGAL_DISABLE_ROUNDING_MATH_CHECK
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>

#include "terraneo/helpers/conversions.hpp"
#include "terraneo/plates/functionsForGeometry.hpp"
#include "terraneo/plates/functionsForRotations.hpp"
#include "terraneo/plates/types.hpp"

namespace terraneo {
namespace plates {

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;

/// Determine to which plate a point belongs
///
/// The function returns a bool to indicate whether any plate matched, the plate's ID and
/// the distance from this plate's boundary
std::tuple< bool, uint_t, real_t > findPlateAndDistance( const real_t age, const PlateStorage& plateStore, const vec3D& point )
{
   // query all plates for given age stage
   auto& plates = plateStore.getPlatesForStage( age );

   // be pessimistic
   bool   plateFound{ false };
   uint_t plateID{ 0 };
   real_t distance{ std::numeric_limits< real_t >::max() };

   for ( auto& currentPlate : plates )
   {
      // rotate point using same rotation to xy-plane as was applied to this plate
      mat3D          rotMtx     = terraneo::plates::getRotationMatrixPolygon( currentPlate.center );
      vec3D          rotPoint   = rotMtx * point;
      const Polygon& bdrPolygon = currentPlate.boundary;
      Point_2        polyPoints[bdrPolygon.size()];
      for ( int index = 0; index < bdrPolygon.size(); ++index )
      {
         polyPoints[index] = Point_2( bdrPolygon[index]( 0 ), bdrPolygon[index]( 1 ) );
      }

      if ( rotPoint[2] > real_c( 0 ) )
      {
         real_t a;

         // check if the point belongs to the polygon
         switch (
             CGAL::bounded_side_2( polyPoints, polyPoints + bdrPolygon.size(), Point_2( rotPoint[0], rotPoint[1] ), Kernel() ) )
         {
         case CGAL::ON_BOUNDED_SIDE:
            for ( int k = 0; k < bdrPolygon.size() - 1; ++k )
            {
               a = terraneo::plates::getDistanceLinePoint(
                   rotPoint,
                   { bdrPolygon[k]( 0 ), bdrPolygon[k]( 1 ), bdrPolygon[k]( 2 ) },
                   { bdrPolygon[k + 1]( 0 ), bdrPolygon[k + 1]( 1 ), bdrPolygon[k + 1]( 2 ) } );
               if ( bdrPolygon[k] != bdrPolygon[k + 1] and a <= distance )
               {
                  distance = a;
               }
            }
            WALBERLA_LOG_INFO_ON_ROOT( "Distance (km): " << distance );
            plateFound = true;
            break;

         case CGAL::ON_BOUNDARY:
            WALBERLA_LOG_INFO_ON_ROOT( " is on the polygon boundary." );
            distance   = real_c( 0 );
            plateFound = true;
            break;

         case CGAL::ON_UNBOUNDED_SIDE:
            break;
         }
      }
      else
      {
         std::string exception{ "How to handle z value <= 0?" };
         WALBERLA_LOG_INFO_ON_ROOT( exception );
         // throw( exception );
      }

      // plate found then leave loop
      if ( plateFound )
      {
         plateID = currentPlate.id;
         break;
      }
   }

   return std::make_tuple( plateFound, plateID, distance );
}

/// From the Euler vector compute the velocity in xyz
vec3D eulerVectorToVelocity( const vec3D& point, vec3D& wXYZ, const real_t smoothing )
{
   real_t earthRadius = plates::constants::earthRadiusInKm * real_c( 1e3 );
   real_t toms        = real_c( 3600 * 24 * 365 ); // conversions factor cm/yr -> m/s

   vec3D eVector;
   vec3D pxyz;

   eVector = terraneo::conversions::degToRad( wXYZ ) * real_c( 1e-6 );

   // Transform to the point to the xyz in a sphere of earthRadius;
   pxyz = terraneo::conversions::cart2sph( point );
   WALBERLA_LOG_INFO_ON_ROOT( "Point (lon, lat): [" << pxyz[0] << ", " << pxyz[1] << "]" );
   pxyz = terraneo::conversions::sph2cart( { pxyz[0], pxyz[1] }, earthRadius );
   WALBERLA_LOG_INFO_ON_ROOT( "Point (x,y,z): [" << point[0] << ", " << point[1] << ", " << point[2] << "]\n" );
   vec3D v = eVector.cross( pxyz );

   WALBERLA_LOG_INFO_ON_ROOT( "Velocity (x,y,z) in cm/yr: [" << v[0] << ", " << v[1] << ", " << v[2] << "] " );

   v *= smoothing / toms;

   WALBERLA_LOG_INFO_ON_ROOT( "Velocity (x,y,z) in m/s:   [" << v[0] << ", " << v[1] << ", " << v[2] << "] " );

   WALBERLA_LOG_INFO_ON_ROOT( "Velocity magnitude: " << sqrt( v[0] * v[0] + v[1] * v[1] + v[2] * v[2] ) * toms * 100 );
   return v;
}

/// Get the velocity in given the plate id, create the reconstruction path, get
/// the rotations and calculate the velocity
vec3D computeCartesianVelocityVector( const PlateRotationProvider& rotData,
                                      const int                    plateID,
                                      const real_t                 age,
                                      const vec3D&                 point,
                                      const real_t                 smoothing )
{
   // age of the euler pole is defined by ((age1 +age2)/2)
   std::array< real_t, 2 >     time{ age, age + 1 };
   std::vector< RotationInfo > recTree;

   std::vector< FiniteRotation >      FinRot;
   const std::vector< RotationInfo >& rotations = rotData.getRotations();

   int pID = plateID;

   while ( pID != 0 )
   {
      rotIter_t rangeBegin;
      rotIter_t rangeEnd;

      for ( rotIter_t it = rotations.begin(); it != rotations.end(); ++it )
      {
         if ( it->plateID == pID )
         {
            rangeBegin = it;
            rangeEnd   = rangeBegin + 1;
            while ( rangeEnd->plateID == pID )
            {
               ++rangeEnd;
            }
            // append to list of finite rotations
            pID = terraneo::plates::determineSeriesOfFiniteRotations( rangeBegin, rangeEnd, time, FinRot );
            break;
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Looping ... (pID = " << pID << ")" );
   }

   std::array< FiniteRotation, 2 > finNahs = terraneo::plates::combineSeriesOfFiniteRotations( FinRot );

   // compute Euler Vector
   vec3D lonlatang = terraneo::plates::stagePoleF( finNahs[0].lonLatAng, finNahs[1].lonLatAng );
   lonlatang[2]    = lonlatang[2] / ( finNahs[1].time - finNahs[0].time );
   vec3D wXYZ      = terraneo::conversions::sph2cart( { lonlatang[0], lonlatang[1] }, lonlatang[2] );

   return eulerVectorToVelocity( point, wXYZ, smoothing );
}

} // namespace plates
} // namespace terraneo
