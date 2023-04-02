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

#include "terraneo/helpers/conversions.hpp"
#include "terraneo/plates/types.hpp"

namespace terraneo::plates {

/// Given two points on a sphere with the radius of the Earth, it returs the
/// distance between them in km using the Haversine formula
inline real_t distancePointPoint( const vec3D& lonlat1, const vec3D& lonlat2 )
{
   real_t phi1 = conversions::degToRad( lonlat1[1] );
   real_t phi2 = conversions::degToRad( lonlat2[1] );
   real_t dphi = conversions::degToRad( lonlat2[1] - lonlat1[1] );
   real_t dlam = conversions::degToRad( lonlat2[0] - lonlat1[0] );

   real_t a = std::sin( dphi * real_c( 0.5 ) ) * std::sin( dphi * real_c( 0.5 ) ) +
              std::cos( phi1 ) * std::cos( phi2 ) * sin( dlam * real_c( 0.5 ) ) * sin( dlam * real_c( 0.5 ) );
   real_t c = real_c( 2 ) * std::atan2( std::sqrt( a ), std::sqrt( real_c( 1 ) - a ) );
   return c * plates::constants::earthRadiusInKm;
}

/// The following function that calculates the intersection between a point
/// and a line is inspired by the Python code from geopy, see:
/// https://github.com/geopy/geopy/blob/master/geopy/distance.py
vec3D intersectPointWithLine( const vec3D& point, const vec3D& lineStart, const vec3D& lineEnd )
{
   WALBERLA_ASSERT( lineStart != lineEnd, "Line in intersectPointWithLine() should not be degenerate." );

   // find the intersect point (the projected point) on the line
   vec3D  lineVec    = ( lineEnd - lineStart );
   real_t lineLength = lineVec.norm();

   // calculate the magnitude of the interesection point
   real_t mag = 0.0;
   mag        = ( point - lineStart ).dot( lineVec );
   mag        = mag / ( lineLength * lineLength );
   vec3D intersection;

   // if closest point does not fall within the line segment, take the shorter distance to an endpoint
   real_t eps = real_c( 1e-5 );
   if ( mag < eps || mag > real_c( 1 ) )
   {
      real_t ix = ( point - lineStart ).norm();
      real_t iy = ( point - lineEnd ).norm();
      if ( ix > iy )
      {
         intersection = lineEnd;
      }
      else
      {
         intersection = lineStart;
      }
   }
   else
   {
      intersection = lineStart + mag * lineVec;
   }

   return intersection;
}

/// Calls all the functions to return the distance between the point and the line segment that
/// constructs the polygon, the returned value is in km
real_t getDistanceLinePoint( const vec3D& point, const vec3D& pini, const vec3D& pend )
{
   vec3D intersectPoint = intersectPointWithLine( point, pini, pend );
   vec3D pintersect     = conversions::cart2sph( intersectPoint );
   vec3D ppoint         = conversions::cart2sph( point );
   return distancePointPoint( ppoint, pintersect );
}

} // namespace terraneo::plates
