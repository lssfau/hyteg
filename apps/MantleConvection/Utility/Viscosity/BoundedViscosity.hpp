/*
 * Copyright (c) 2025 Andreas Burkhart.
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

// Bounded viscosity
// Modifies an existing by making sure that the viscosity lies within a given range.
// This is enforced via std::min, std::max.

class BoundedViscosity : public TemperatureDependentViscosityModel< real_t >
{
 public:
   using TemperatureDependentViscosityModel< real_t >::prefix_;

   BoundedViscosity( NondimensionalisationParameters&                                       nondimensionalisation,
                     walberla::Config::BlockHandle&                                         parameters,
                     const std::shared_ptr< TemperatureDependentViscosityModel< real_t > >& ModifiedModel,
                     std::string                                                            prefix = "" )
   : TemperatureDependentViscosityModel< real_t >( prefix )
   , ModifiedModel_( ModifiedModel )
   {
      minVisc_ = parameters.getParameter< real_t >( prefix + std::string( "minVisc" ) );
      maxVisc_ = parameters.getParameter< real_t >( prefix + std::string( "maxVisc" ) );
   }

   real_t evaluate( const hyteg::Point3D& x, real_t temp ) override
   {
      return std::max( std::min( ModifiedModel_->evaluate( x, temp ), maxVisc_ ), minVisc_ );
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "###################################################"                                                                 << "\n";
      os << std::string( offset, ' ') << "################ Bounded Viscosity ################"                                                                 << "\n";
      os << std::string( offset, ' ') << "###################################################"                                                                 << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 11 ) << std::left << "prefix_: "                           << prefix_                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 11 ) << std::left << "minVisc_: "                          << minVisc_                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 11 ) << std::left << "maxVisc_: "                          << maxVisc_                         << "\n";
      // clang-format on

      return ModifiedModel_->print( os, offset + 3 );
   }

 public:
   std::shared_ptr< TemperatureDependentViscosityModel< real_t > > ModifiedModel_;

   real_t minVisc_; // nondimensional minimum viscosity
   real_t maxVisc_; // nondimensional maximum viscosity
};

inline std::ostream& operator<<( std::ostream& os, const BoundedViscosity& bv )
{
   return bv.print( os );
}

} // namespace MantleConvection