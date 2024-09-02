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

real_t ConvectionSimulation::densityFunction( const Point3D& x )
{
   real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
   real_t retVal;

   if ( TN.simulationParameters.radialProfile )
   {
      updateNonDimParameters( x );
   }

   //implement adiabatic compression, determined by dissipation number and gruneisen parameter
   real_t rho = TN.physicalParameters.surfaceDensity *
                std::exp( TN.physicalParameters.dissipationNumber * ( TN.domainParameters.rMax - radius ) /
                          TN.physicalParameters.grueneisenParameter );

   retVal = rho / TN.physicalParameters.referenceDensity;

   return retVal;
}

real_t ConvectionSimulation::diffPreFactorFunction( const Point3D& x )
{
   if ( TN.simulationParameters.radialProfile )
   {
      updateNonDimParameters( x );
   }
   return ( real_c( 1.0 ) ) / ( densityFunction( x ) * TN.physicalParameters.pecletNumber );
}

real_t ConvectionSimulation::adiabaticCoefficientFunction( const Point3D& x )
{
   if ( TN.simulationParameters.radialProfile )
   {
      updateNonDimParameters( x );
   }
   return TN.physicalParameters.dissipationNumber;
}

real_t ConvectionSimulation::constantEnergyCoefficientFunction( const Point3D& x )
{
   real_t intHeatingFactor = 1.0;
   if ( TN.simulationParameters.radialProfile )
   {
      updateNonDimParameters( x );
   }
   return TN.physicalParameters.hNumber * intHeatingFactor;
}

real_t ConvectionSimulation::surfaceTempCoefficientFunction( const Point3D& )
{
   return 0.0;
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

         //just used for setting up p2p1StokesFunctionContainer["VelocityFEPrev"] in the initialisation
         if ( TN.simulationParameters.timeStep == 0 )
         {
            p2p1StokesFunctionContainer["VelocityFEPrev"]->uvw()[coordIdx].interpolate( Velocity, l, idSurface );
         }
      }
   }
}

void ConvectionSimulation::updateViscosity()
{
   if ( p2InjectionOperator == nullptr )
   {
      p2InjectionOperator =
          std::make_shared< P2toP2QuadraticInjection >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
   }

   for(uint_t level = TN.domainParameters.maxLevel; level > TN.domainParameters.minLevel; level--)
   {
      p2InjectionOperator->restrict(*( p2ScalarFunctionContainer[std::string( "TemperatureFE" )] ), level, All);
   }

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > viscosityInit =
       [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
          updateNonDimParameters( x );
          return terraneo::viscosityFunction( x, Temperature[0], TN );
       };

   auto viscosityInitInv = [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
      updateNonDimParameters( x );
      return ( 1.0 / ( terraneo::viscosityFunction( x, Temperature[0], TN ) ) );
   };

   // Before interpolation: Ensure the new reference viscosity is set to the min current viscosity if a viscosity profile
   // or a temperature dependent viscosity is utilized.
   // This is will ensure that that the minimum non-dimensionalised value for the viscosity = 1.
   p2ScalarFunctionContainer["ViscosityFE"]->interpolate(
       viscosityInit, { *( p2ScalarFunctionContainer[std::string( "TemperatureFE" )] ) }, TN.domainParameters.maxLevel, All );

   real_t maxViscosity = p2ScalarFunctionContainer["ViscosityFE"]->getMaxValue( TN.domainParameters.maxLevel ) *
                         TN.physicalParameters.referenceViscosity;
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Max viscosity [Pa s]: " << maxViscosity );
   real_t minRefViscosity = p2ScalarFunctionContainer["ViscosityFE"]->getMinValue( TN.domainParameters.maxLevel ) *
                            TN.physicalParameters.referenceViscosity;
   WALBERLA_LOG_INFO_ON_ROOT( "New update reference viscosity [Pa s]: " << minRefViscosity );

   if ( TN.simulationParameters.tempDependentViscosity || TN.simulationParameters.haveViscosityProfile )
   {
      // Update reference viscosity with new min Viscosity value
      TN.physicalParameters.referenceViscosity = minRefViscosity;
   }

   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      p2ScalarFunctionContainer["ViscosityFE"]->interpolate(
          viscosityInit, { *( p2ScalarFunctionContainer[std::string( "TemperatureFE" )] ) }, l, All );
      p2ScalarFunctionContainer["ViscosityFEInv"]->interpolate(
          viscosityInitInv, { *( p2ScalarFunctionContainer[std::string( "TemperatureFE" )] ) }, l, All );
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

void ConvectionSimulation::oppositeGravity( const Point3D& p, Point3D& n )
{
   real_t radius = p.norm();
   n             = Point3D( { p[0] / radius, p[1] / radius, p[2] / radius } );
}

//returns a reference adiabat relevant for the Earth, commonly implemented in TALA
real_t ConvectionSimulation::referenceTemperatureFunction( const Point3D& x )
{
   if ( TN.simulationParameters.radialProfile )
   {
      updateNonDimParameters( x );
   }

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

void ConvectionSimulation::updateNonDimParameters( const Point3D& x )
{
   // This function updates the Rayleigh Number, Dissipation number, Peclet number or
   // h Number if any corresponding radial profile for thermal expansivity, specific heat capacity,
   // grueneisen Parameter or viscosity is provided.

   if ( TN.simulationParameters.haveThermalExpProfile )
   {
      // update Thermal Expansivity parameter
      TN.physicalParameters.thermalExpansivity = terraneo::interpolateDataValues( x,
                                                                                  TN.physicalParameters.radiusAlpha,
                                                                                  TN.physicalParameters.thermalExpansivityProfile,
                                                                                  TN.domainParameters.rMin,
                                                                                  TN.domainParameters.rMax );
   }

   if ( TN.simulationParameters.haveSpecificHeatCapProfile )
   {
      // update Specific Heat Capacity (at const. pressure)
      TN.physicalParameters.specificHeatCapacity =
          terraneo::interpolateDataValues( x,
                                           TN.physicalParameters.radiusCp,
                                           TN.physicalParameters.specificHeatCapacityProfile,
                                           TN.domainParameters.rMin,
                                           TN.domainParameters.rMax );
   }

   if ( TN.simulationParameters.haveGrueneisenProfile )
   {
      // update Grueneisen parameter
      TN.physicalParameters.grueneisenParameter = terraneo::interpolateDataValues( x,
                                                                                   TN.physicalParameters.radiusGamma,
                                                                                   TN.physicalParameters.grueneisenProfile,
                                                                                   TN.domainParameters.rMin,
                                                                                   TN.domainParameters.rMax );
   }

   TN.physicalParameters.thermalDiffusivity =
       TN.physicalParameters.thermalConductivity /
       ( TN.physicalParameters.referenceDensity * TN.physicalParameters.specificHeatCapacity );

   TN.physicalParameters.rayleighNumber =
       ( TN.physicalParameters.referenceDensity * TN.physicalParameters.gravity * TN.physicalParameters.thermalExpansivity *
         TN.physicalParameters.mantleThickness * TN.physicalParameters.mantleThickness * TN.physicalParameters.mantleThickness *
         ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) ) /
       ( TN.physicalParameters.referenceViscosity * TN.physicalParameters.thermalDiffusivity );

   TN.physicalParameters.pecletNumber = ( TN.physicalParameters.characteristicVelocity * TN.physicalParameters.mantleThickness ) /
                                        TN.physicalParameters.thermalDiffusivity;
   TN.physicalParameters.dissipationNumber =
       ( TN.physicalParameters.thermalExpansivity * TN.physicalParameters.gravity * TN.physicalParameters.mantleThickness ) /
       TN.physicalParameters.specificHeatCapacity;

   TN.physicalParameters.hNumber = ( TN.physicalParameters.internalHeatingRate * TN.physicalParameters.mantleThickness ) /
                                   ( TN.physicalParameters.specificHeatCapacity * TN.physicalParameters.characteristicVelocity *
                                     ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) );
}

void ConvectionSimulation::outputTimingTree()
{
   auto timer = storage->getTimingTree();
   writeTimingTreeJSON( *timer, TN.outputParameters.outputDirectory + "/" + "TimingTree.json" );
}

} // namespace terraneo