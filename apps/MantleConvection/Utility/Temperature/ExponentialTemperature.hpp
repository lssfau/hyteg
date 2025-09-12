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

#include "TemperatureModel.hpp"

using walberla::real_t;

namespace MantleConvection {

class ExponentialTemperature : public TemperatureModel< real_t >
{
 public:
   using TemperatureModel< real_t >::prefix_;

   ExponentialTemperature( NondimensionalisationParameters&                 nondimensionalisation,
                           walberla::Config::BlockHandle&                   parameters,
                           std::function< bool( const hyteg::Point3D& x ) > surfaceFct,
                           std::function< bool( const hyteg::Point3D& x ) > CMBFct,
                           std::string                                      prefix = "" )
   : TemperatureModel< real_t >( prefix )
   , surfaceFct_( surfaceFct )
   , CMBFct_( CMBFct )
   {
      temperatureSurface_ = nondimensionalisation.temperatureSurface_;
      temperatureCMB_     = nondimensionalisation.temperatureCMB_;
      radiusSurface_      = nondimensionalisation.radiusSurface_;
      Di_                 = nondimensionalisation.Di_;
      deltaT_             = nondimensionalisation.deltaT_;

      temperatureAdiabaticSurfaceDimensional_ =
          parameters.getParameter< real_t >( prefix + std::string( "temperatureAdiabaticSurfaceDimensional" ) );

      temperatureAdiabaticSurface_ = temperatureAdiabaticSurfaceDimensional_ / nondimensionalisation.deltaT_;
   }

   real_t evaluateNoCheck( const hyteg::Point3D& x ) override
   {
      return temperatureAdiabaticSurface_ * std::exp( Di_ * ( radiusSurface_ - x.norm() ) );
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
      os << std::string( offset, ' ') << "#####################################################"                                                                                  << "\n";
      os << std::string( offset, ' ') << "############## Exponential Temperature ##############"                                                                                  << "\n";
      os << std::string( offset, ' ') << "#####################################################"                                                                                  << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                                                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 42 ) << std::left << "prefix_: "                                   << prefix_                                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 42 ) << std::left << "temperatureSurface_: "                       << temperatureSurface_                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 42 ) << std::left << "temperatureCMB_: "                           << temperatureCMB_                             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 42 ) << std::left << "radiusSurface_: "                            << radiusSurface_                              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 42 ) << std::left << "Di_: "                                       << Di_                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 42 ) << std::left << "deltaT_: "                                   << deltaT_                                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 42 ) << std::left << "temperatureAdiabaticSurfaceDimensional_: "   << temperatureAdiabaticSurfaceDimensional_     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 42 ) << std::left << "temperatureAdiabaticSurface_: "              << temperatureAdiabaticSurface_                       ;
      // clang-format on

      return os;
   }

 public:
   std::function< bool( const hyteg::Point3D& x ) > surfaceFct_;         // returns true if x is on the surface
   std::function< bool( const hyteg::Point3D& x ) > CMBFct_;             // returns true if x is on the cmb
   real_t                                           temperatureSurface_; // nondimensional temperature surface
   real_t                                           temperatureCMB_;     // nondimensional temperature CMB
   real_t                                           radiusSurface_;      // nondimensional radius surface, rMax
   real_t                                           Di_;                 // Dissipation number
   real_t                                           deltaT_;             // temperature difference between surface and CMB [K]

   real_t temperatureAdiabaticSurfaceDimensional_; // dimensional adiabatic surface temperature
   real_t temperatureAdiabaticSurface_;            // nondimensional adiabatic surface temperature
};

inline std::ostream& operator<<( std::ostream& os, const ExponentialTemperature& et )
{
   return et.print( os );
}

} // namespace MantleConvection