/*
 * Copyright (c) 2024 Eugenio D'Ascoli, Ponsuganth Ilangovan, Nils Kohl.
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

#include "hyteg/types/PointND.hpp"

#include "terraneo/dataimport/ParameterIO.hpp"

namespace terraneo {
real_t interpolateDataValues( const hyteg::Point3D&        x,
                              const std::vector< real_t >& radius,
                              const std::vector< real_t >& values,
                              real_t                       rMin,
                              real_t                       rMax )
{
   real_t pointRadius = x.norm();
   real_t retVal      = 1.0;

   // Check if radius and viscosity std::vector are filled and at least 2 entries for radius are present
   if ( radius.size() != values.size() || radius.size() < 2 )
   {
      WALBERLA_ABORT( "Value- and radius vector must be of the same size and contain at least two elements." );
   }
   // Check if the min and max radius of the input profile are reasonable

   if ( radius[radius.size() - 1] < rMin || radius[0] > rMax )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Inconsistent radial profile loaded!" );
      WALBERLA_LOG_INFO_ON_ROOT( "Min radius of radial profile: " << radius[radius.size() - 1] );
      WALBERLA_LOG_INFO_ON_ROOT( "Max radius of radial profile: " << radius[0] );
      WALBERLA_ABORT( "Cancel simulation run" );
   }

   // Loop over radius vector and find a datapoint that lies in between two given values
   // If true: perform interpolation to estimate the value
   // Check boundaries and set values accordingly

   if ( pointRadius >= radius[0] )
   {
      retVal = values[0];
   }
   else if ( pointRadius <= radius[radius.size() - 1] )
   {
      retVal = values[radius.size() - 1];
   }
   else
   {
      uint_t count = 0;
      while ( pointRadius < radius[count] )
      {
         ++count;
      }
      real_t interpolFactor = ( pointRadius - radius[count] ) / ( radius[count - 1] - radius[count] );
      retVal                = ( interpolFactor * ( values[count - 1] - values[count] ) ) + values[count];
   }
   return retVal;
}
} // namespace terraneo