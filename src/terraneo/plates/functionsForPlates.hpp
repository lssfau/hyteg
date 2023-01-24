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

#include <boost/geometry.hpp>

#include "terraneo/helpers/conversions.hpp"
#include "terraneo/plates/functionsForGeometry.hpp"
#include "terraneo/plates/functionsForRotations.hpp"
#include "terraneo/plates/types.hpp"

namespace terraneo {
namespace plates {

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::spherical_equatorial<boost::geometry::degree>> spherical_point;
typedef boost::geometry::model::polygon<spherical_point> polygon_on_sphere;

/// Determine to which plate a point belongs
///
/// The function returns a bool to indicate whether any plate matched, the plate's ID and
/// the distance from this plate's boundary

inline std::tuple< bool, uint_t, real_t >
    findPlateAndDistance( const real_t age, const PlateStorage& plateStore, const vec3D& point, uint_t idWhenNoPlateFound )
{
   // query all plates for given age stage
   auto& plates = plateStore.getPlatesForStage( std::ceil( age ) );

   // be pessimistic
   bool   plateFound{ false };
   uint_t plateID{ idWhenNoPlateFound };
   real_t distance{ std::numeric_limits< real_t >::max() };

   // Create the point in the surface of a sphere from the library boost::geometry 
   spherical_point pntSph(point[0], point[1]);

   //loop over the plates available
   for ( auto& currentPlate : plates )
   {
      // create the polygon on the surface of a sphere and populate with the plate boundary coordinates
      const Polygon& bdrPolygon = currentPlate.boundary;
      polygon_on_sphere polygonOnSphere;
      for ( int index = 0; index < bdrPolygon.size(); ++index )
      {
         boost::geometry::append(polygonOnSphere.outer(), spherical_point(bdrPolygon[index]( 0 ), bdrPolygon[index]( 1 )));
      }
      // Correct the geometry
      boost::geometry::correct(polygonOnSphere);

      // check if the point belongs to the polygon
      if (boost::geometry::within(pntSph, polygonOnSphere))
      {
         // calculate the distace from the polygon to the point
         boost::geometry::for_each_segment(polygonOnSphere, [&distance, &pntSph](const auto& segment){
            distance = std::min<real_t>(distance, boost::geometry::distance(segment, pntSph));
         });
         // distance at the surface of the Earth
         distance = distance * 6371.0; 
         // plate is found 
         plateFound = true;
      }
      else{
         // ;
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

/// From the Euler vector compute the surface velocity in xyz
inline vec3D eulerVectorToVelocity( const vec3D& point, vec3D& wXYZ, const real_t smoothing )
{
   real_t earthRadius = plates::constants::earthRadiusInKm * real_c( 1e3 );
   real_t toms        = real_c( 3600 * 24 * 365 ); // conversions factor cm/yr -> m/s

   vec3D eVector;
   vec3D pxyz;

   eVector = terraneo::conversions::degToRad( wXYZ ) * real_c( 1e-6 );

   // Transform to the point to the xyz in a sphere of earthRadius;
   pxyz    = terraneo::conversions::sph2cart( { point[0], point[1] }, earthRadius );
   vec3D v = eVector.cross( pxyz );

   v *= smoothing / toms;

   return v;
}

/// Get the velocity in given the plate id, create the reconstruction path, get
/// the rotations and calculate the velocity
inline vec3D computeCartesianVelocityVector( const PlateRotationProvider& rotData,
                                             const int                    plateID,
                                             const real_t                 age,
                                             const vec3D&                 point,
                                             const real_t                 smoothing )
{
   // age of the euler pole is defined by ((age1 + age2)/2)
   // This is valid when the velocities are calculated every 1 Myrs. 
   // Otherwise this needs to be changed to the desired time step
   // taking into conserdation the time resolution of the plate boundaries available.
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
      // WALBERLA_LOG_DETAIL_ON_ROOT( "Looping ... (pID = " << pID << ")" );
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
