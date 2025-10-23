/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include "../Data/DataTools.hpp"
#include "TemperatureDependentDensityModel.hpp"

using walberla::real_c;
using walberla::real_t;

namespace MantleConvection {

// expects range to be sorted and to be of the same length as values
class LinearDensityInterpolation : public TemperatureDependentDensityModel< real_t >
{
 public:
   using TemperatureDependentDensityModel< real_t >::prefix_;

   LinearDensityInterpolation( NondimensionalisationParameters& nondimensionalisation,
                               walberla::Config::BlockHandle&   parameters,
                               std::vector< real_t >&           range,
                               std::vector< real_t >&           values,
                               std::string                      prefix = "" )
   : TemperatureDependentDensityModel< real_t >( prefix )
   , range_( range )
   , values_( values )
   {
      WALBERLA_UNUSED(nondimensionalisation);
      WALBERLA_UNUSED(parameters);
      if ( values.size() != range.size() )
      {
         WALBERLA_LOG_WARNING( "The values and range vectors of LinearDensityInterpolation should be of the same size!" )
      }
   }

   real_t evaluate( const hyteg::Point3D& x, real_t temp ) override
   {
      WALBERLA_UNUSED( temp );
      return MantleConvection::linearInterpolation( x.norm(), range_, values_ );
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "####################################################"                                    << "\n";
      os << std::string( offset, ' ') << "########### Linear Density Interpolation ###########"                                    << "\n";
      os << std::string( offset, ' ') << "####################################################"                                    << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 11 ) << std::left << "prefix_: "  << prefix_                              ;
      // clang-format on

      return os;
   }

 public:
   std::vector< real_t > range_;  // vector containing the range of the interpolation values
   std::vector< real_t > values_; // vector containing the interpolation values
};

inline std::ostream& operator<<( std::ostream& os, const LinearDensityInterpolation& ldi )
{
   return ldi.print( os );
}
} // namespace MantleConvection