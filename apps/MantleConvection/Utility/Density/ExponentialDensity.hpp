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

#include "TemperatureDependentDensityModel.hpp"

using walberla::real_c;
using walberla::real_t;

namespace MantleConvection {

class ExponentialDensity : public TemperatureDependentDensityModel< real_t >
{
 public:
   using TemperatureDependentDensityModel< real_t >::prefix_;

   ExponentialDensity( NondimensionalisationParameters& nondimensionalisation,
                       walberla::Config::BlockHandle&   parameters,
                       std::string                      prefix = "" )
   : TemperatureDependentDensityModel< real_t >( prefix )
   {
      radiusSurface_ = nondimensionalisation.radiusSurface_;
      GammaRef_      = nondimensionalisation.GammaRef_;
      Di_            = nondimensionalisation.Di_;

      rhoSurfDimensional_ = parameters.getParameter< real_t >( "rhoSurf" );
      rhoSurf_            = rhoSurfDimensional_ / nondimensionalisation.rhoRef_;
   }

   real_t evaluate( const hyteg::Point3D& x, real_t temp ) override
   {
      WALBERLA_UNUSED( temp );
      return rhoSurf_ * std::exp( Di_ / GammaRef_ * ( radiusSurface_ - x.norm() ) );
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#####################################################"                                           << "\n";
      os << std::string( offset, ' ') << "################ Exponential Density ################"                                           << "\n";
      os << std::string( offset, ' ') << "#####################################################"                                           << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "prefix_: "                << prefix_                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "rhoSurfDimensional_: "    << rhoSurfDimensional_     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "radiusSurface_: "         << radiusSurface_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "GammaRef_: "              << GammaRef_               << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "Di_: "                    << Di_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "rhoSurf_: "               << rhoSurf_                       ;
      // clang-format on

      return os;
   }

 public:
   real_t rhoSurfDimensional_; // dimensional density at the surface
   real_t radiusSurface_;      // nondimensional radius surface, rMax
   real_t GammaRef_;           // grueneisen parameter nondim. reference value
   real_t Di_;                 // Dissipation number
   real_t rhoSurf_;            // nondimensional surface density
};

inline std::ostream& operator<<( std::ostream& os, const ExponentialDensity& ed )
{
   return ed.print( os );
}

} // namespace MantleConvection