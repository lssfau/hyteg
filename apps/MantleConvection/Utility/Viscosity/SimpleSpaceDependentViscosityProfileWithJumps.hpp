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

// This is equivalent to the viscosityTerraComp.csv profile in the data folder
class SimpleSpaceDependentViscosityProfileWithJumps : public TemperatureDependentViscosityModel< real_t >
{
 public:
   using TemperatureDependentViscosityModel< real_t >::prefix_;

   SimpleSpaceDependentViscosityProfileWithJumps( NondimensionalisationParameters& nondimensionalisation,
                                                  walberla::Config::BlockHandle&   parameters,
                                                  real_t                           scaling = real_c( 1.0 ),
                                                  std::string                      prefix  = "" )
   : TemperatureDependentViscosityModel< real_t >( prefix )
   , scaling_( scaling )
   {
      WALBERLA_UNUSED(parameters);
      etaRef_    = nondimensionalisation.etaRef_;
      radiusCMB_ = nondimensionalisation.radiusCMB_;
   }

   real_t evaluate( const hyteg::Point3D& x, real_t temp ) override
   {
      WALBERLA_UNUSED( temp );

      real_t pos = ( x.norm() - radiusCMB_ );

      return ( ( pos <= 0.75 ) ?
                   ( 1.25e+23 ) :
                   ( ( pos <= 0.75840517241379302 ) ?
                         ( 3.54619628745821e+26 * std::pow( pos, 2 ) - 5.43063135866058e+26 * pos + 2.0794881073001899e+26 ) :
                         ( ( pos <= 0.76400862068965503 ) ?
                               ( 1.8656863432129701e+26 * std::pow( pos, 2 ) - 2.8816164906437698e+26 * pos +
                                 1.11289507706839e+26 ) :
                               ( ( pos <= 0.77241379310344804 ) ?
                                     ( 1.06929133385375e+26 * std::pow( pos, 2 ) - 1.66471118539444e+26 * pos +
                                       6.4803200518165701e+25 ) :
                                     ( ( pos <= 0.85172413999999996 ) ?
                                           ( 1.5e+22 ) :
                                           ( ( pos <= 0.85948275862069001 ) ?
                                                 ( 3.8579740318342099e+25 * std::pow( pos, 2 ) - 6.6965235842778496e+25 * pos +
                                                   2.90638521562008e+25 ) :
                                                 ( ( pos <= 0.86465517241379297 ) ?
                                                       ( 2.2385551110116e+25 * std::pow( pos, 2 ) - 3.9127983014155496e+25 * pos +
                                                         1.7101032729417599e+25 ) :
                                                       ( ( pos <= 0.87241379310344802 ) ?
                                                             ( 1.38463307937214e+25 * std::pow( pos, 2 ) -
                                                               2.4361020984252302e+25 * pos + 1.0716867679420699e+25 ) :
                                                             ( ( pos <= 0.93793103 ) ?
                                                                   ( 2.5e+21 ) :
                                                                   ( ( pos <= 0.95344827586206904 ) ?
                                                                         ( 1.9551462059537398e+25 * std::pow( pos, 2 ) -
                                                                           3.6535322574369698e+25 * pos + 1.70704057747143e+25 ) :
                                                                         ( ( pos <= 0.96379310344827596 ) ?
                                                                               ( 5.10752322498446e+25 * std::pow( pos, 2 ) -
                                                                                 9.6647891247610802e+25 * pos +
                                                                                 4.57275182542852e+25 ) :
                                                                               ( ( pos <= 0.97413793103448298 ) ?
                                                                                     ( 1.28735188889197e+26 * std::pow( pos, 2 ) -
                                                                                       2.4634415249381102e+26 * pos +
                                                                                       1.17865630354825e+26 ) :
                                                                                     ( ( pos <= 0.98448275862069001 ) ?
                                                                                           ( 2.7542030108487401e+26 *
                                                                                                 std::pow( pos, 2 ) -
                                                                                             5.3212721590952799e+26 * pos +
                                                                                             2.5706169141506601e+26 ) :
                                                                                           ( 8.7854567127168094e+26 *
                                                                                                 std::pow( pos, 2 ) -
                                                                                             1.7196602723807901e+27 * pos +
                                                                                             8.4161460110911101e+26 ) ) ) ) ) ) ) ) ) ) ) ) ) ) /
             etaRef_ * scaling_;
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "###########################################################"                              << "\n";
      os << std::string( offset, ' ') << "### Simple Space Dependent Viscosity Profile With Jumps ###"                              << "\n";
      os << std::string( offset, ' ') << "###########################################################"                              << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 13 ) << std::left            << "prefix_: "    << prefix_           << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 13 ) << std::left            << "etaRef_: "    << etaRef_           << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 13 ) << std::left            << "radiusCMB_: " << radiusCMB_        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 13 ) << std::left            << "scaling_: "   << scaling_                 ;
      // clang-format on

      return os;
   }

 public:
   real_t scaling_;   // scaling for the profile
   real_t etaRef_;    // viscosity nondim. reference value                  [Pa s]
   real_t radiusCMB_; // nondimensional radius cmb, rMin
};

inline std::ostream& operator<<( std::ostream& os, const SimpleSpaceDependentViscosityProfileWithJumps& spvp )
{
   return spvp.print( os );
}

} // namespace MantleConvection