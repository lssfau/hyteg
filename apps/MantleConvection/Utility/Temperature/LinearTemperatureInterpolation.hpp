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
#include "TemperatureModel.hpp"

using walberla::real_t;

namespace MantleConvection {

// expects range to be sorted and to be of the same length as values
// if range and values are between 0 and 1 set unitRangeAndValues to true
class LinearTemperatureInterpolation : public TemperatureModel< real_t >
{
 public:
   using TemperatureModel< real_t >::prefix_;

   LinearTemperatureInterpolation( NondimensionalisationParameters&                 nondimensionalisation,
                                   walberla::Config::BlockHandle&                   parameters,
                                   std::function< bool( const hyteg::Point3D& x ) > surfaceFct,
                                   std::function< bool( const hyteg::Point3D& x ) > CMBFct,
                                   std::vector< real_t >&                           range,
                                   std::vector< real_t >&                           values,
                                   bool                                             unitRangeAndValues = false,
                                   std::string                                      prefix             = "" )
   : TemperatureModel< real_t >( prefix )
   , surfaceFct_( surfaceFct )
   , CMBFct_( CMBFct )
   , range_( range )
   , values_( values )
   , unitRangeAndValues_( unitRangeAndValues )
   {
      WALBERLA_UNUSED(parameters);
      temperatureSurface_ = nondimensionalisation.temperatureSurface_;
      temperatureCMB_     = nondimensionalisation.temperatureCMB_;
      radiusCMB_          = nondimensionalisation.radiusCMB_;

      if ( values.size() != range.size() )
      {
         WALBERLA_LOG_WARNING( "The values and range vectors of LinearTemperatureInterpolation should be of the same size!" )
      }
   }

   real_t evaluateNoCheck( const hyteg::Point3D& x ) override
   {
      real_t radius = x.norm();
      if ( unitRangeAndValues_ )
      {
         radius -= radiusCMB_;
      }

      real_t res = MantleConvection::linearInterpolation( radius, range_, values_ );

      if ( unitRangeAndValues_ )
      {
         return temperatureSurface_ + res;
      }
      return res;
   }

   real_t evaluate( const hyteg::Point3D& x ) override
   {
      if ( checkSurface( x ) )
      {
         return temperatureSurface_;
      }
      else if ( checkCMB( x ) )
      {
         return temperatureCMB_;
      }

      return evaluateNoCheck( x );
   }

   bool checkSurface( const hyteg::Point3D& x ) override { return surfaceFct_( x ); }
   bool checkCMB( const hyteg::Point3D& x ) override { return CMBFct_( x ); }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "######################################################"                                    << "\n";
      os << std::string( offset, ' ') << "########## Linear Temperature Interpolation ##########"                                    << "\n";
      os << std::string( offset, ' ') << "######################################################"                                    << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                           << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "prefix_: "             << prefix_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "temperatureSurface_: " << temperatureSurface_  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "temperatureCMB_: "     << temperatureCMB_      << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "radiusCMB_: "          << radiusCMB_           << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "unitRangeAndValues_: " << unitRangeAndValues_         ;
      // clang-format on

      return os;
   }

 public:
   std::function< bool( const hyteg::Point3D& x ) > surfaceFct_;         // returns true if x is on the surface
   std::function< bool( const hyteg::Point3D& x ) > CMBFct_;             // returns true if x is on the cmb
   real_t                                           temperatureSurface_; // nondimensional temperature surface
   real_t                                           temperatureCMB_;     // nondimensional temperature CMB
   real_t                                           radiusCMB_;          // nondimensional radius cmb, rMin
   std::vector< real_t >                            range_;  // vector containing the range of the interpolation values
   std::vector< real_t >                            values_; // vector containing the interpolation values
   bool                                             unitRangeAndValues_; // true if range and values are between 0 and 1
};

inline std::ostream& operator<<( std::ostream& os, const LinearTemperatureInterpolation& lti )
{
   return lti.print( os );
}
} // namespace MantleConvection