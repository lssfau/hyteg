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
#include "../Pressure/PressureModel.hpp"
#include "TemperatureDependentDensityModel.hpp"

using walberla::real_c;
using walberla::real_t;

namespace MantleConvection {

// expects rangeP, rangeT to be sorted and to be of the same length as the respective number of row and columns of values
class PressureProfileDensityBilinearInterpolation : public TemperatureDependentDensityModel< real_t >
{
 public:
   using TemperatureDependentDensityModel< real_t >::prefix_;

   PressureProfileDensityBilinearInterpolation( NondimensionalisationParameters&                  nondimensionalisation,
                                                walberla::Config::BlockHandle&                    parameters,
                                                const std::shared_ptr< PressureModel< real_t > >& PressureProfile,
                                                std::vector< real_t >&                            rangeP,
                                                std::vector< real_t >&                            rangeT,
                                                std::vector< std::vector< real_t > >&             values,
                                                std::string                                       prefix = "" )
   : TemperatureDependentDensityModel< real_t >( prefix )
   , rangeP_( rangeP )
   , rangeT_( rangeT )
   , values_( values )
   , PressureProfile_( PressureProfile )
   {
      WALBERLA_UNUSED(nondimensionalisation);
      WALBERLA_UNUSED(parameters);
      if ( ( values.size() != rangeP.size() ) || ( values.front().size() != rangeT.size() ) )
      {
         WALBERLA_LOG_WARNING(
             "The values and range vector P dimension mismatch (rows) for PressureProfileDensityBilinearInterpolation!" );
      }
      for ( uint_t i = 0; i < values.size(); i++ )
      {
         if ( ( values.at( i ).size() != rangeT.size() ) )
         {
            WALBERLA_LOG_WARNING(
                "The values and range vector T dimension mismatch (colums) for PressureProfileDensityBilinearInterpolation in row "
                << i << "!" );
         }
      }
   }

   real_t evaluate( const hyteg::Point3D& x, real_t temp ) override
   {
      return MantleConvection::bilinearInterpolation( PressureProfile_->evaluate( x ), temp, rangeP_, rangeT_, values_ );
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#######################################################################"                 << "\n";
      os << std::string( offset, ' ') << "########### Pressure Profile Density Bilinear Interpolation ###########"                 << "\n";
      os << std::string( offset, ' ') << "#######################################################################"                 << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 11 ) << std::left << "prefix_: "  << prefix_                              ;
      // clang-format on

      return os;
   }

 public:
   std::vector< real_t >                rangeP_; // vector containing the pressure range of the interpolation values
   std::vector< real_t >                rangeT_; // vector containing the temperature range of the interpolation values
   std::vector< std::vector< real_t > > values_; // vector containing the interpolation values

   std::shared_ptr< PressureModel< real_t > >
       PressureProfile_; // space dependent pressure profile used to evaluate the given table at point (p(x),T)
};

inline std::ostream& operator<<( std::ostream& os, const PressureProfileDensityBilinearInterpolation& ppdbi )
{
   return ppdbi.print( os );
}
} // namespace MantleConvection