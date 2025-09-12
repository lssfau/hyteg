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

#include "core/DataTypes.h"
#include "core/config/Config.h"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace MantleConvection {

struct NondimensionalisationParameters
{
 public:
   NondimensionalisationParameters()
   {}

   NondimensionalisationParameters( walberla::Config::BlockHandle& parameters, std::string prefix = "" )
   : prefix_( prefix )
   {
      // independent parameters
      temperatureSurfaceDimensional_ = parameters.getParameter< real_t >( prefix_ + std::string( "temperatureSurface" ) );
      temperatureCMBDimensional_     = parameters.getParameter< real_t >( prefix_ + std::string( "temperatureCMB" ) );
      radiusSurfaceDimensional_      = parameters.getParameter< real_t >( prefix_ + std::string( "radiusSurface" ) );
      radiusCMBDimensional_          = parameters.getParameter< real_t >( prefix_ + std::string( "radiusCMB" ) );
      etaRef_                        = parameters.getParameter< real_t >( prefix_ + std::string( "etaRef" ) );
      rhoRef_                        = parameters.getParameter< real_t >( prefix_ + std::string( "rhoRef" ) );
      C_pRef_                        = parameters.getParameter< real_t >( prefix_ + std::string( "C_pRef" ) );
      alphaRef_                      = parameters.getParameter< real_t >( prefix_ + std::string( "alphaRef" ) );
      gRef_                          = parameters.getParameter< real_t >( prefix_ + std::string( "gRef" ) );
      GammaRef_                      = parameters.getParameter< real_t >( prefix_ + std::string( "GammaRef" ) );
      uRef_                          = parameters.getParameter< real_t >( prefix_ + std::string( "uRef" ) );
      kRef_                          = parameters.getParameter< real_t >( prefix_ + std::string( "kRef" ) );

      // dependent parameters
      d_             = radiusSurfaceDimensional_ - radiusCMBDimensional_;
      radiusSurface_ = radiusSurfaceDimensional_ / d_;
      radiusCMB_     = radiusCMBDimensional_ / d_;

      deltaT_             = temperatureCMBDimensional_ - temperatureSurfaceDimensional_;
      temperatureSurface_ = temperatureSurfaceDimensional_ / deltaT_;
      temperatureCMB_     = temperatureCMBDimensional_ / deltaT_;

      tRef_       = d_ / uRef_;
      pRef_       = etaRef_ * uRef_ / d_;
      HRef_       = uRef_ * deltaT_ * C_pRef_ / d_;
      Q_LRef_     = rhoRef_ * uRef_ * deltaT_ * C_pRef_ / d_;
      K_TRef_     = GammaRef_ * rhoRef_ * C_pRef_ / alphaRef_;
      kappa_TRef_ = alphaRef_ / GammaRef_ / rhoRef_ / C_pRef_;
      kappaRef_   = kRef_ / rhoRef_ / C_pRef_;

      Ra_    = rhoRef_ * alphaRef_ * gRef_ * deltaT_ * std::pow( d_, 3 ) / kappaRef_ / etaRef_;
      Di_    = alphaRef_ * gRef_ * d_ / C_pRef_;
      gamma_ = Di_ / GammaRef_;
      Pe_    = uRef_ * d_ / kappaRef_;
      xi_    = alphaRef_ * deltaT_;
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const
   {
      // clang-format off
      os << std::string( offset, ' ') << "######################################################"                                                            << "\n";
      os << std::string( offset, ' ') << "########## Nondimensionalisation Parameters ##########"                                                            << "\n";
      os << std::string( offset, ' ') << "######################################################"                                                            << "\n";
      os << std::string( offset, ' ') << "   " << "------Independent Parameters------"                                                                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "prefix_: "                           << prefix_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "temperatureSurfaceDimensional_: "    << temperatureSurfaceDimensional_ << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "temperatureCMBDimensional_: "        << temperatureCMBDimensional_     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "radiusSurfaceDimensional_: "         << radiusSurfaceDimensional_      << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "radiusCMBDimensional_: "             << radiusCMBDimensional_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "etaRef_: "                           << etaRef_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "rhoRef_: "                           << rhoRef_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "C_pRef_: "                           << C_pRef_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "alphaRef_: "                         << alphaRef_                      << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "gRef_: "                             << gRef_                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "GammaRef_: "                         << GammaRef_                      << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "uRef_: "                             << uRef_                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "kRef_: "                             << kRef_                          << "\n";
      os << std::string( offset, ' ') << "   " << "-------Dependent parameters-------"                                                                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "d_: "                                << d_                             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "radiusSurface_: "                    << radiusSurface_                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "radiusCMB_: "                        << radiusCMB_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "deltaT_: "                           << deltaT_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "temperatureSurface_: "               << temperatureSurface_            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "temperatureCMB_: "                   << temperatureCMB_                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "tRef_: "                             << tRef_                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "pRef_: "                             << pRef_                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "HRef_: "                             << HRef_                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "Q_LRef_: "                           << Q_LRef_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "K_TRef_: "                           << K_TRef_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "kappa_TRef_: "                       << kappa_TRef_                    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "kappaRef_: "                         << kappaRef_                      << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "Ra_: "                               << Ra_                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "Di_: "                               << Di_                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "Pe_: "                               << Pe_                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "gamma_: "                            << gamma_                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 33 ) << std::left << "xi_: "                               << xi_                                   ;
      // clang-format on

      return os;
   }

 public:
   // clang-format off
   std::string prefix_;

   // independent parameters
   real_t temperatureSurfaceDimensional_; // dimensional temperature surface                    [K]
   real_t temperatureCMBDimensional_;     // dimensional temperature CMB                        [K]
   real_t radiusSurfaceDimensional_;      // dimensional radius surface                         [m]
   real_t radiusCMBDimensional_;          // dimensional radius CMB                             [m]
   real_t etaRef_;                        // viscosity nondim. reference value                  [Pa s]
   real_t rhoRef_;                        // density nondim. reference value                    [kg / m^3]
   real_t C_pRef_;                        // specific heat capacity nondim. reference value     [J / K / kg]
   real_t alphaRef_;                      // thermal expansivity nondim. reference value        [1 / K]
   real_t gRef_;                          // gravity nondim. reference value                    [m / s^2]
   real_t GammaRef_;                      // grueneisen parameter nondim. reference value       []
   real_t uRef_;                          // velocity nondim. reference value                   [m / s]
   real_t kRef_;                          // thermal conductivity nondim. reference value       [W / m / K]

   // dependent parameters
   real_t d_;                           // dimensional mantle thickness                         [m]
   real_t radiusSurface_;               // nondimensional radius surface, rMax
   real_t radiusCMB_;                   // nondimensional radius cmb, rMin

   real_t deltaT_;                      // dim. temperature difference between surface and CMB  [K]
   real_t temperatureSurface_;          // nondimensional temperature surface
   real_t temperatureCMB_;              // nondimensional temperature CMB
   
   real_t tRef_;                        // time nondim. reference value                         [s]
   real_t pRef_;                        // pressure nondim. reference value                     [Pa]
   real_t HRef_;                        // internal heating nondim. reference value             [W / kg]
   real_t Q_LRef_;                      // latent heat,phase transition nondim. reference value [W / m^3]
   real_t K_TRef_;                      // isothermal bulk modulus nondim. reference value      [Pa]
   real_t kappa_TRef_;                  // isothermal compressibility nondim. reference value   [1 / Pa]
   real_t kappaRef_;                    // thermal diffusivity nondim. reference value          [m^2 / s]

   real_t Ra_;                          // Rayleigh number
   real_t Di_;                          // Dissipation number
   real_t Pe_;                          // Peclet number
   real_t gamma_;                       // mantle compressibility
   real_t xi_;                          // driving term

   // clang-format on
};

inline std::ostream& operator<<( std::ostream& os, const NondimensionalisationParameters& nd )
{
   return nd.print( os );
}

} // namespace MantleConvection