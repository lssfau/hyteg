/*
 * Copyright (c) 2024 Eugenio D'Ascoli, Ponsuganth Ilangovan, Nils Kohl.
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

#include "hyteg/types/PointND.hpp"

#include "terraneo/helpers/TerraNeoParameters.hpp"

namespace terraneo {
real_t viscosityFunction( const hyteg::Point3D& x, real_t Temperature, const TerraNeoParameters& TN )
{
   real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
   real_t retVal = 1.0;

   // If a viscosity profile is provided, use it, otherwise use a constant background viscosity
   if ( TN.simulationParameters.haveViscosityProfile )
   {
      retVal = terraneo::interpolateDataValues( x,
                                                TN.physicalParameters.radius,
                                                TN.physicalParameters.viscosityProfile,
                                                TN.domainParameters.rMin,
                                                TN.domainParameters.rMax );
   }
   else
   {
      retVal = TN.physicalParameters.viscosity;
   }
   //scale background viscosity by temperature- and depth-dependent factors
   //depth-dependent factor counteracts the decrease in viscosity due to increasing temperature with depth
   if ( TN.simulationParameters.tempDependentViscosity )
   {
      // Account for non-dim temperature to be between 0-1
      Temperature -= TN.physicalParameters.surfaceTemp / ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );

      switch ( TN.simulationParameters.tempDependentViscosityType )
      {
      //Frank–Kamenetskii type 1
      case 0: {
         retVal *= std::exp( -TN.physicalParameters.activationEnergy * ( Temperature ) +
                             TN.physicalParameters.depthViscosityFactor * ( TN.domainParameters.rMax - radius ) /
                                 ( TN.domainParameters.rMax - TN.domainParameters.rMin ) );
         break;
      }
      //Frank–Kamenetskii type 2
      case 1: {
         retVal *= std::exp( TN.physicalParameters.activationEnergy * ( real_c( 0.5 ) - Temperature ) +
                             TN.physicalParameters.depthViscosityFactor * ( TN.domainParameters.rMax - radius ) /
                                 ( TN.domainParameters.rMax - TN.domainParameters.rMin ) );
         break;
      }

      //with respect to mean
      case 2: {
         uint_t shell = static_cast< uint_t >(
             std::round( real_c( TN.simulationParameters.numLayers ) *
                         ( ( radius - TN.domainParameters.rMin ) / ( TN.domainParameters.rMax - TN.domainParameters.rMin ) ) ) );

         retVal *= std::exp( -TN.physicalParameters.activationEnergy *
                             ( Temperature - TN.physicalParameters.temperatureProfile.at( shell ) ) );

         break;
      }
      //Arrhenius type
      case 3: {
         retVal *= std::exp( TN.physicalParameters.activationEnergy *
                                 ( ( real_c( 1 ) / ( Temperature + real_c( 0.25 ) ) ) - real_c( 1.45 ) ) +
                             TN.physicalParameters.depthViscosityFactor * ( TN.domainParameters.rMax - radius ) /
                                 ( TN.domainParameters.rMax - TN.domainParameters.rMin ) );

         break;
      }
      //Frank–Kamenetskii type 1
      default: {
         retVal *= std::exp( -TN.physicalParameters.activationEnergy * ( Temperature ) +
                             TN.physicalParameters.depthViscosityFactor * ( TN.domainParameters.rMax - radius ) /
                                 ( TN.domainParameters.rMax - TN.domainParameters.rMin ) );
         break;
      }
      }

      //impose min viscosity
      if ( retVal < TN.physicalParameters.viscosityLowerBound )
      {
         retVal = TN.physicalParameters.viscosityLowerBound;
      }

      //impose max viscosity
      if ( retVal > TN.physicalParameters.viscosityUpperBound )
      {
         retVal = TN.physicalParameters.viscosityUpperBound;
      }
   }

   retVal /= TN.physicalParameters.referenceViscosity;

   return retVal;
}
} // namespace terraneo