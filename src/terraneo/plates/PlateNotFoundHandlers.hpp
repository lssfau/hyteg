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

#pragma once

#include <algorithm>
#include <vector>

namespace terraneo::plates {

/// \defgroup PlateNotFoundHandlers PlateNotFoundHandlers
///
/// Functors for handling the case the we fail to find a plate associated with a point.
/// These will be called by PlateVelocityProvider::getPointVelocity() in this case
/// and get passed the point as vec3D object and the age as realt.
/// @{

/// Default strategy
///
/// The default strategy when we find no plate for a point is to set the velocity to zero.
class DefaultPlateNotFoundHandler
{
 public:
   vec3D operator()( const vec3D& point, const real_t age )
   {
      WALBERLA_UNUSED( point );
      WALBERLA_UNUSED( age );
      return {real_c( 0 ), real_c( 0 ), real_c( 0 )};
   }
};

/// Statistics strategy
///
/// The statistics strategy will record the points and ages for which we fail to find plates
/// and can be queried for a report.
class StatisticsPlateNotFoundHandler
{
 public:
   vec3D operator()( const vec3D& point, const real_t age )
   {
      pointsNotFound_[age].push_back( point );
      return {real_c( 0 ), real_c( 0 ), real_c( 0 )};
   }

   void generateReport()
   {
      // remove duplicates from list of points:
      // we can only compare vec3D for equality, but for sorting we'd need "<"
      // std::sort( pointsNotFound_.begin(), pointsNotFound_.end() );
      // auto last = std::unique( pointsNotFound_.begin(), pointsNotFound_.end() );
      // pointsNotFound_.erase( last, pointsNotFound_.end() );

      for ( const auto& elem : pointsNotFound_ )
      {
         real_t                      age          = elem.first;
         const std::vector< vec3D >& pointsForAge = elem.second;
         WALBERLA_LOG_INFO( "List of points w/o plateID for age = " << age << ":" );
         for ( uint k = 0; k < pointsForAge.size(); ++k )
         {
            WALBERLA_LOG_INFO( "" << std::setw( 3 ) << k << ": (" << std::scientific << std::showpos << pointsForAge[k]( 0 )
                                  << ", " << pointsForAge[k]( 1 ) << ", " << pointsForAge[k]( 2 ) << ")" );
         }
      }
   }

 private:
   std::map< real_t, std::vector< vec3D > > pointsNotFound_;
};

/// @}

} // namespace terraneo::plates
