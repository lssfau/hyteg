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

#include "Convection.hpp"
#include "terraneo/dataimport/ParameterIO.hpp"

namespace terraneo {

////////////////////////
// Public functions  //
///////////////////////

const SimulationParameters& ConvectionSimulation::getSimulationParams()
{
   return TN.simulationParameters;
}

real_t ConvectionSimulation::viscosityFunction( const Point3D& x, real_t Temperature )
{
   real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
   real_t retVal = 1.0;

   // If a viscosity profile is provided, use it, otherwise use a constant background viscosity

   if ( TN.simulationParameters.haveViscosityProfile )
   {
      // Check if radius and viscosity std::vector are filled and at least 2 entries for radius are present
      if ( TN.physicalParameters.radius.size() != TN.physicalParameters.viscosityProfile.size() ||
           TN.physicalParameters.radius.size() < 2 )
      {
         WALBERLA_ABORT( "Viscosity and radius vector must be of the same size and contain at least two elements." );
      }

      // Loop over radius vector and find a datapoint that lies in between two given values
      // If true: perform interpolation to estimate the value
      // Check boundaries and set viscosity values accordingly

      if ( radius >= TN.physicalParameters.radius[0] )
      {
         retVal = TN.physicalParameters.viscosityProfile[0];
      }
      else if ( radius <= TN.physicalParameters.radius[TN.physicalParameters.radius.size()] )
      {
         retVal = TN.physicalParameters.viscosityProfile[TN.physicalParameters.radius.size()];
      }
      else
      {
         uint_t count = 0;
         while ( radius < TN.physicalParameters.radius[count] )
         {
            ++count;
         }
         real_t interpolFactor = ( radius - TN.physicalParameters.radius[count] ) /
                                 ( TN.physicalParameters.radius[count - 1] - TN.physicalParameters.radius[count] );
         retVal = ( interpolFactor *
                    ( TN.physicalParameters.viscosityProfile[count - 1] - TN.physicalParameters.viscosityProfile[count] ) ) +
                  TN.physicalParameters.viscosityProfile[count];
      }
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

real_t ConvectionSimulation::densityFunction( const Point3D& x )
{
   auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
   real_t retVal;

   //implement adiabatic compression, determined by dissipation number and gruneisen parameter
   real_t rho = TN.physicalParameters.surfaceDensity *
                std::exp( TN.physicalParameters.dissipationNumber * ( TN.domainParameters.rMax - radius ) /
                          TN.physicalParameters.grueneisenParameter );

   retVal = rho / TN.physicalParameters.referenceDensity;

   return retVal;
}

real_t ConvectionSimulation::diffPreFactorFunction( const Point3D& x )
{
   return ( real_c( 1.0 ) ) / ( densityFunction( x ) * TN.physicalParameters.pecletNumber );
}

void ConvectionSimulation::updatePlateVelocities( StokesFunction& U )
{
   uint_t coordIdx = 0;

   //function to return plate velocities, copied and adapted from PlateVelocityDemo.cpp.
   std::function< real_t( const Point3D& ) > Velocity = [&coordIdx]( const Point3D& x ) {
      terraneo::vec3D coords{ x[0], x[1], x[2] };
      //get velocity at current plate age (intervals of 1Ma)
      terraneo::vec3D velocity = oracle->getPointVelocity(
          coords,
          TN.simulationParameters.plateAge,
          terraneo::plates::LinearDistanceSmoother{ real_c( 1 ) / TN.simulationParameters.plateSmoothingDistance },
          terraneo::plates::DefaultPlateNotFoundHandler{} );

      return velocity[int_c( coordIdx )] /
             ( TN.physicalParameters.characteristicVelocity *
               TN.simulationParameters.plateVelocityScaling ); //non-dimensionalise by dividing by characteristic velocity
   };

   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; ++l )
   {
      for ( coordIdx = 0; coordIdx < 3; ++coordIdx )
      {
         //interpolate current plate velocities at the surface
         U.uvw()[coordIdx].interpolate( Velocity, l, idSurface );

         //just used for setting up p2p1StokesFunctionContainer["StokesLHSPrev"] in the initialisation
         if ( TN.simulationParameters.timeStep == 0 )
         {
            p2p1StokesFunctionContainer["StokesLHSPrev"]->uvw()[coordIdx].interpolate( Velocity, l, idSurface );
         }
      }
   }
}

void ConvectionSimulation::updateViscosity()
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > viscosityInit =
       [&]( const Point3D& x, const std::vector< real_t >& Temperature ) { return viscosityFunction( x, Temperature[0] ); };

   auto viscosityInitInv = [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
      return ( 1.0 / ( viscosityFunction( x, Temperature[0] ) ) );
   };

   // Before interpolation: Ensure the new reference viscosity is set to the min current viscosity if a viscosity profile
   // or a temperature dependent viscosity is utilized.
   // This is will ensure that that the minimum non-dimensionalised value for the viscosity = 1.
   p2ScalarFunctionContainer["ViscosityFE"]->interpolate( viscosityInit, { *(p2ScalarFunctionContainer[std::string("TemperatureFE")]) }, TN.domainParameters.maxLevel, All );

   real_t maxViscosity = p2ScalarFunctionContainer["ViscosityFE"]->getMaxValue( TN.domainParameters.maxLevel ) * TN.physicalParameters.referenceViscosity;
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Max viscosity [Pa s]: " << maxViscosity );
   real_t minRefViscosity = p2ScalarFunctionContainer["ViscosityFE"]->getMinValue( TN.domainParameters.maxLevel ) * TN.physicalParameters.referenceViscosity;
   WALBERLA_LOG_INFO_ON_ROOT( "New update reference viscosity [Pa s]: " << minRefViscosity );

   if ( TN.simulationParameters.tempDependentViscosity || TN.simulationParameters.haveViscosityProfile )
   {
      // Update reference viscosity with new min Viscosity value
      TN.physicalParameters.referenceViscosity = minRefViscosity;

      TN.physicalParameters.rayleighNumber =
          ( TN.physicalParameters.referenceDensity * TN.physicalParameters.gravity * TN.physicalParameters.thermalExpansivity *
            std::pow( TN.physicalParameters.mantleThickness, 3 ) *
            ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) ) /
          ( TN.physicalParameters.referenceViscosity * TN.physicalParameters.thermalDiffusivity );
   }

   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      p2ScalarFunctionContainer["ViscosityFE"]->interpolate( viscosityInit, { *(p2ScalarFunctionContainer[std::string("TemperatureFE")]) }, l, All );
      p2ScalarFunctionContainer["ViscosityFEInv"]->interpolate( viscosityInitInv, { *(p2ScalarFunctionContainer[std::string("TemperatureFE")]) }, l, All );
   }
}

void ConvectionSimulation::normalFunc( const Point3D& p, Point3D& n )
{
   real_t radius = p.norm();
   if ( std::abs( radius - TN.domainParameters.rMax ) < std::abs( radius - TN.domainParameters.rMin ) )
   {
      n = Point3D( { p[0] / radius, p[1] / radius, p[2] / radius } );
   }
   else
   {
      n = Point3D( { -p[0] / radius, -p[1] / radius, -p[2] / radius } );
   }
}

//returns a reference adiabat relevant for the Earth, commonly implemented in TALA
real_t ConvectionSimulation::referenceTemperatureFunction( const Point3D& x )
{
   auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

   if ( ( radius - TN.domainParameters.rMin ) < real_c( 1e-10 ) )
   {
      return ( TN.physicalParameters.cmbTemp ) / ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );
   }
   else if ( ( TN.domainParameters.rMax - radius ) < real_c( 1e-10 ) )
   {
      return ( TN.physicalParameters.surfaceTemp ) / ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );
   }

   real_t temp = TN.physicalParameters.adiabatSurfaceTemp *
                 std::exp( ( TN.physicalParameters.dissipationNumber * ( TN.domainParameters.rMax - radius ) ) );

   real_t retVal = ( temp ) / ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );

   return retVal;
}

void ConvectionSimulation::outputTimingTree()
{
   auto timer = storage->getTimingTree();
   writeTimingTreeJSON( *timer, TN.outputParameters.outputDirectory + "/" + "TimingTree.json" );
}

} // namespace terraneo