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

namespace terraneo {
namespace plates {

using terraneo::plates::FiniteRotation;
using terraneo::plates::rotIter_t;

/// Computes a matrix for rotation around a given axis vector by given angle
mat3D xyzRotationMatrix( const vec3D& vector, real_t angle )
{
   mat3D  rotMat;
   real_t x       = vector( 0 );
   real_t y       = vector( 1 );
   real_t z       = vector( 2 );
   rotMat( 0, 0 ) = cos( angle ) + x * x * ( real_c( real_c( 1 ) ) - cos( angle ) );
   rotMat( 0, 1 ) = x * y * ( real_c( 1 ) - cos( angle ) ) - z * sin( angle );
   rotMat( 0, 2 ) = x * z * ( real_c( 1 ) - cos( angle ) ) + y * sin( angle );
   rotMat( 1, 0 ) = y * x * ( real_c( 1 ) - cos( angle ) ) + z * sin( angle );
   rotMat( 1, 1 ) = cos( angle ) + y * y * ( real_c( 1 ) - cos( angle ) );
   rotMat( 1, 2 ) = y * z * ( real_c( 1 ) - cos( angle ) ) - x * sin( angle );
   rotMat( 2, 0 ) = z * x * ( real_c( 1 ) - cos( angle ) ) - y * sin( angle );
   rotMat( 2, 1 ) = z * y * ( real_c( 1 ) - cos( angle ) ) + x * sin( angle );
   rotMat( 2, 2 ) = cos( angle ) + z * z * ( real_c( 1 ) - cos( angle ) );
   return rotMat;
}

/// Obtains the rotation matrix to move a plate polygon with a center to (0,0,0)
mat3D getRotationMatrixPolygon( const vec3D& axis )
{
   real_t xm    = axis( 0 );
   real_t ym    = axis( 1 );
   real_t zm    = axis( 2 );
   real_t norm  = std::sqrt( xm * xm + ym * ym );
   real_t angle = atan2( norm, zm );
   return xyzRotationMatrix( { ym / norm, -xm / norm, real_c( 0 ) }, angle );
}

/// Determine rotation matrix for a finite rotation
mat3D rotationMatrix( const real_t lon, const real_t lat, const real_t angleInDegree )
{
   real_t angleInRadians = conversions::degToRad( angleInDegree );
   vec3D  eXYZ           = conversions::sph2cart( { lon, lat } );
   return xyzRotationMatrix( eXYZ, angleInRadians );
}

/// Calculate longitude, latitude and angle of rotation for the reconstruction path from the data of the rotational file
vec3D rotMatrix2LonLatW( const mat3D& Rot )
{
   vec3D  lonlatang;
   real_t sqrtRot = sqrt( ( Rot( 2, 1 ) - Rot( 1, 2 ) ) * ( Rot( 2, 1 ) - Rot( 1, 2 ) ) +
                          ( Rot( 0, 2 ) - Rot( 2, 0 ) ) * ( Rot( 0, 2 ) - Rot( 2, 0 ) ) +
                          ( Rot( 1, 0 ) - Rot( 0, 1 ) ) * ( Rot( 1, 0 ) - Rot( 0, 1 ) ) );
   // latitude
   if ( sqrtRot == real_c( 0 ) ) // <- direct comparision of FP value for zero?
                                 // should probably add a tolerance
   {
      lonlatang[1] = real_c( 0 );
   }
   else
   {
      lonlatang[1] = asin( ( Rot( 1, 0 ) - Rot( 0, 1 ) ) / sqrtRot );
   }

   // longitude
   lonlatang[0] = atan2( Rot( 0, 2 ) - Rot( 2, 0 ), Rot( 2, 1 ) - Rot( 1, 2 ) );

   // angle
   lonlatang[2] = atan2( sqrtRot, Rot( 0, 0 ) + Rot( 1, 1 ) + Rot( 2, 2 ) - real_c( 1 ) );
   lonlatang    = conversions::radToDeg( lonlatang );

   return lonlatang;
}

/// Calculate stage Pole for a pair of rotations
vec3D stagePoleF( const vec3D& rot1, const vec3D& rot2 )
{
   mat3D R1     = rotationMatrix( rot1[0], rot1[1], rot1[2] );
   mat3D R2     = rotationMatrix( rot2[0], rot2[1], rot2[2] * real_c( -1 ) );
   mat3D stageR = R1 * R2;
   return rotMatrix2LonLatW( stageR );
}

/// Get the intermediate rotations of the reconstruction tree
vec3D intermediateRotations( const terraneo::plates::RotationInfo& rot1, const terraneo::plates::RotationInfo& rot2, real_t time )
{
   real_t dT        = ( rot2.time - time ) / ( rot2.time - rot1.time );
   mat3D  R2        = rotationMatrix( rot2.longitude, rot2.latitude, rot2.angle );
   vec3D  lonlatang = stagePoleF( { rot1.longitude, rot1.latitude, rot1.angle }, { rot2.longitude, rot2.latitude, rot2.angle } );
   mat3D  stgR      = rotationMatrix( lonlatang[0], lonlatang[1], lonlatang[2] * dT );
   mat3D  intR      = stgR * R2;

   return rotMatrix2LonLatW( intR );
}

/// Given the input data from the rotational find the lines of the specific
/// time asked and get the rotations
int determineSeriesOfFiniteRotations( const rotIter_t&               rangeBegin,
                                      const rotIter_t&               rangeEnd,
                                      const std::array< real_t, 2 >& time,
                                      std::vector< FiniteRotation >& finRot )
{
   std::vector< terraneo::plates::RotationInfo > fin0;
   int                                           pID = -1;

   terraneo::plates::RotationInfo aux;
   aux.time        = real_c( 0 );
   aux.longitude   = real_c( 0 );
   aux.latitude    = real_c( 90 );
   aux.angle       = real_c( 0 );
   aux.plateID     = 0;
   aux.conjugateID = 0;
   fin0.push_back( aux );

   for ( rotIter_t rot = rangeBegin; rot != rangeEnd; ++rot )
   {
      fin0.push_back( *rot );
   }

   for ( int ii = 0; ii < time.size(); ++ii )
   {
      real_t t = time[ii];

      // determine in which position of the asked time to grab the corresponding lines
      //
      // NOTE: if loop does not return a hit, the 0 remains unchanged!
      // Starting from 1, since fin0 has that leading extra entry
      int jj = 0;
      for ( int idx = 1; idx < fin0.size(); ++idx )
      {
         if ( fin0[idx].time - t >= real_c( 0 ) )
         {
            jj = idx;
            break;
         }
      }

      if ( fin0.size() - 2 > jj )
      {
         if ( int( time[0] ) == int( fin0[jj].time ) && fin0[jj].time == fin0[jj + 1].time )
            jj = jj + 2;
      }

      vec3D lonlatAng = intermediateRotations( fin0[jj - 1], fin0[jj], t );
      finRot.push_back( { t, lonlatAng } );
      pID = fin0[jj].conjugateID;
   }

   return pID;
}

/// Combine a sequence of finite rotations
std::array< FiniteRotation, 2 > combineSeriesOfFiniteRotations( const std::vector< FiniteRotation >& FinRot )
{
   std::array< FiniteRotation, 2 > absFin;

   for ( int jj = 0; jj < 2; ++jj )
   {
      absFin[jj] = { real_c( 0 ), { real_c( 0 ), real_c( 90 ), real_c( 0 ) } };

      for ( int ii = jj; ii < FinRot.size(); ii = ii + 2 )
      {
         mat3D Rot1           = rotationMatrix( absFin[jj].lonLatAng[0], absFin[jj].lonLatAng[1], absFin[jj].lonLatAng[2] );
         mat3D Rot2           = rotationMatrix( FinRot[ii].lonLatAng[0], FinRot[ii].lonLatAng[1], FinRot[ii].lonLatAng[2] );
         Rot1                 = Rot2 * Rot1;
         absFin[jj].lonLatAng = rotMatrix2LonLatW( Rot1 );
         absFin[jj].time      = FinRot[ii].time;
      }
   }

   return absFin;
}

} // namespace plates
} // namespace terraneo
