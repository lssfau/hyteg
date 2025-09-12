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

// Basic Frank–Kamenetskii flow, see e.g.
// D.A. May and J. Brown and L. Le Pourhiet,
// A scalable, matrix-free multigrid preconditioner for finite element discretizations of heterogeneous Stokes flow,
// Computer Methods in Applied Mechanics and Engineering, Volume 290, p.496-523, 2015,
// https://doi.org/10.1016/j.cma.2015.03.014
//
// or
//
// Lin, Yi-An and Colli, Lorenzo and Wu, Jonny
// NW Pacific-Panthalassa Intra-Oceanic Subduction During Mesozoic Times From Mantle Convection and Geoid Models
// Geochemistry, Geophysics, Geosystems, Volume 23, Number 11, Pages e2022GC010514,
// https://doi.org/10.1029/2022GC010514
//
// Note: The viscosity function should be defined for the full temperature T = T_d + T_s
class ExponentialViscosity : public TemperatureDependentViscosityModel< real_t >
{
 public:
   using TemperatureDependentViscosityModel< real_t >::prefix_;

   ExponentialViscosity( NondimensionalisationParameters&                                       nondimensionalisation,
                         walberla::Config::BlockHandle&                                         parameters,
                         const std::shared_ptr< TemperatureDependentViscosityModel< real_t > >& ModifiedModel,
                         std::string                                                            prefix = "" )
   : TemperatureDependentViscosityModel< real_t >( prefix )
   , ModifiedModel_( ModifiedModel )
   {
      temperatureSurface_ = nondimensionalisation.temperatureSurface_;
      radiusSurface_      = nondimensionalisation.radiusSurface_;

      rockChemicalCompositionParameter_ =
          parameters.getParameter< real_t >( prefix + std::string( "rockChemicalCompositionParameter" ) );
      depthDependency_ = parameters.getParameter< real_t >( prefix + std::string( "depthDependency" ) );
      additiveOffSet_  = parameters.getParameter< real_t >( prefix + std::string( "additiveOffSet" ) );
   }

   real_t evaluate( const hyteg::Point3D& x, real_t temp ) override
   {
      // get the temperature to be between 0 and 1 for the exp function
      temp -= temperatureSurface_;

      // if ( temp < real_c( 0 ) )
      // {
      //    WALBERLA_LOG_WARNING_ON_ROOT( "Temperature " << temp
      //                                                 << " smaller than surface temperature detected in viscosity evaluation!" );
      // }
      // else if ( temp > real_c( 1 ) )
      // {
      //    WALBERLA_LOG_WARNING_ON_ROOT( "Temperature " << temp
      //                                                 << " greater than cmb temperature detected in viscosity evaluation!" );
      // }

      real_t pos = ( radiusSurface_ - x.norm() );

      return ModifiedModel_->evaluate( x, temp ) * std::exp( -rockChemicalCompositionParameter_ * temp + depthDependency_ * pos + additiveOffSet_ );
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#######################################################"                                                                 << "\n";
      os << std::string( offset, ' ') << "################ Exponential Viscosity ################"                                                                 << "\n";
      os << std::string( offset, ' ') << "#######################################################"                                                                 << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 36 ) << std::left << "prefix_: "                           << prefix_                              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 36 ) << std::left << "temperatureSurface_: "               << temperatureSurface_                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 36 ) << std::left << "radiusSurface_: "                    << radiusSurface_                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 36 ) << std::left << "rockChemicalCompositionParameter_: " << rockChemicalCompositionParameter_    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 36 ) << std::left << "depthDependency_: "                  << depthDependency_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 36 ) << std::left << "additiveOffSet_: "                   << additiveOffSet_                      << "\n";
      // clang-format on

      return ModifiedModel_->print( os, offset + 3 );
   }

 public:
   std::shared_ptr< TemperatureDependentViscosityModel< real_t > > ModifiedModel_;

   real_t temperatureSurface_; // nondimensional temperature surface
   real_t radiusSurface_;      // nondimensional radius surface, rMax

   real_t rockChemicalCompositionParameter_; // chemical composition parameter of the rock, activation energy
   real_t depthDependency_;                  // depth dependency factor for the viscosity, activation volume
   real_t
       additiveOffSet_; // additive offset inside the exp function (set to 0.5*rockChemicalCompositionParameter_ for Frank–Kamenetskii type 2 flow)
};

inline std::ostream& operator<<( std::ostream& os, const ExponentialViscosity& ev )
{
   return ev.print( os );
}

} // namespace MantleConvection