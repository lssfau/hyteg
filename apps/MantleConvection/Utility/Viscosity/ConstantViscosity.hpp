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

#include "TemperatureDependentViscosityModel.hpp"

using walberla::real_c;
using walberla::real_t;

namespace MantleConvection {

class ConstantViscosity : public TemperatureDependentViscosityModel< real_t >
{
 public:
   using TemperatureDependentViscosityModel< real_t >::prefix_;

   ConstantViscosity( NondimensionalisationParameters& nondimensionalisation,
                      walberla::Config::BlockHandle&   parameters,
                      real_t                           scaling = real_c( 1.0 ),
                      std::string                      prefix  = "" )
   : TemperatureDependentViscosityModel< real_t >( prefix )
   , scaling_( scaling )
   {
      WALBERLA_UNUSED(nondimensionalisation);
      WALBERLA_UNUSED(parameters);
   }

   real_t evaluate( const hyteg::Point3D& x, real_t temp ) override
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( temp );
      return scaling_;
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "######################################################"             << "\n";
      os << std::string( offset, ' ') << "################# Constant Viscosity #################"             << "\n";
      os << std::string( offset, ' ') << "######################################################"             << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 11 ) << std::left << "prefix_: "  << prefix_  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 11 ) << std::left << "scaling_: " << scaling_        ;

      // clang-format on

      return os;
   }

 public:
   real_t scaling_; // constant scaling
};

inline std::ostream& operator<<( std::ostream& os, const ConstantViscosity& cv )
{
   return cv.print( os );
}

} // namespace MantleConvection