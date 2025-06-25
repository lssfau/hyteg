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
   real_t rho;

   if ( TN.simulationParameters.haveDensityProfile && TN.simulationParameters.compressible )
   {
      // Interpolate density values from given profile
      rho = interpolateDataValues( x,
                                   TN.physicalParameters.radiusDensity,
                                   TN.physicalParameters.densityProfile,
                                   TN.domainParameters.rMin,
                                   TN.domainParameters.rMax );
   }
   else
   {
      if ( !TN.simulationParameters.compressible )
      {
         // If incompressible density is rho / rho_ref = 1
         rho = TN.physicalParameters.referenceDensity;
      }
      else
      {
         // Implement adiabatic compression, determined by dissipation number and gruneisen parameter
         rho = TN.physicalParameters.surfaceDensity *
               std::exp( TN.physicalParameters.dissipationNumber * ( TN.domainParameters.rMax - radius ) /
                         TN.physicalParameters.grueneisenParameter );
      }
   }
   retVal = rho / TN.physicalParameters.referenceDensity;

   return retVal;
}

real_t ConvectionSimulation::diffPreFactorFunction( const Point3D& x )
{
   if ( TN.simulationParameters.haveSpecificHeatCapProfile )
   {
      TN.physicalParameters.specificHeatCapacityRadial =
          terraneo::interpolateDataValues( x,
                                           TN.physicalParameters.radiusCp,
                                           TN.physicalParameters.specificHeatCapacityProfile,
                                           TN.domainParameters.rMin,
                                           TN.domainParameters.rMax );
      return ( ( real_c( 1.0 ) / ( densityFunction( x ) * TN.physicalParameters.pecletNumber ) *
                 ( TN.physicalParameters.specificHeatCapacity / TN.physicalParameters.specificHeatCapacityRadial ) ) );
   }
   else
   {
      return ( real_c( 1.0 ) ) / ( densityFunction( x ) * TN.physicalParameters.pecletNumber );
   }
}

real_t ConvectionSimulation::adiabaticCoefficientFunction( const Point3D& x )
{
   if ( TN.simulationParameters.haveSpecificHeatCapProfile && TN.simulationParameters.haveThermalExpProfile )
   {
      TN.physicalParameters.thermalExpansivityRadial =
          terraneo::interpolateDataValues( x,
                                           TN.physicalParameters.radiusAlpha,
                                           TN.physicalParameters.thermalExpansivityProfile,
                                           TN.domainParameters.rMin,
                                           TN.domainParameters.rMax );
      TN.physicalParameters.specificHeatCapacityRadial =
          terraneo::interpolateDataValues( x,
                                           TN.physicalParameters.radiusCp,
                                           TN.physicalParameters.specificHeatCapacityProfile,
                                           TN.domainParameters.rMin,
                                           TN.domainParameters.rMax );

      return ( TN.physicalParameters.dissipationNumber *
               ( TN.physicalParameters.thermalExpansivityRadial / TN.physicalParameters.thermalExpansivity ) *
               ( TN.physicalParameters.specificHeatCapacity / TN.physicalParameters.specificHeatCapacityRadial ) );
   }
   else
   {
      return TN.physicalParameters.dissipationNumber;
   }
}

real_t ConvectionSimulation::constantEnergyCoefficientFunction( const Point3D& x )
{
   real_t intHeatingFactor = 1.0;
   if ( TN.simulationParameters.haveSpecificHeatCapProfile )
   {
      TN.physicalParameters.specificHeatCapacityRadial =
          terraneo::interpolateDataValues( x,
                                           TN.physicalParameters.radiusCp,
                                           TN.physicalParameters.specificHeatCapacityProfile,
                                           TN.domainParameters.rMin,
                                           TN.domainParameters.rMax );
      return ( TN.physicalParameters.hNumber * intHeatingFactor *
               ( TN.physicalParameters.specificHeatCapacity / TN.physicalParameters.specificHeatCapacityRadial ) );
   }
   else
   {
      return TN.physicalParameters.hNumber * intHeatingFactor;
   }
}

real_t ConvectionSimulation::surfaceTempCoefficientFunction( const Point3D& )
{
   return 0.0;
}

void ConvectionSimulation::updatePlateVelocities( StokesFunction& U )
{
   uint_t                                           coordIdx = 0;
   terraneo::plates::StatisticsPlateNotFoundHandler handlerWithStatistics;

   // callback function for computing the velocity components
   std::function< real_t( const Point3D& ) > Velocity = [&coordIdx, &handlerWithStatistics]( const Point3D& point ) {
      vec3D coords{ point[0], point[1], point[2] };
      vec3D velocity = oracle->getLocallyAveragedPointVelocity(
          coords, TN.simulationParameters.plateAge, *avgPointProvider, handlerWithStatistics );
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

   for ( uint_t level = TN.domainParameters.maxLevel; level > TN.domainParameters.minLevel; level-- )
   {
      p2InjectionOperator->restrict( *( p2ScalarFunctionContainer[std::string( "TemperatureFE" )] ), level, All );
   }

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > viscosityInit =
       [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
          return terraneo::viscosityFunction( x, Temperature[0], TN );
       };

   auto viscosityInitInv = [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
      return ( 1.0 / ( terraneo::viscosityFunction( x, Temperature[0], TN ) ) );
   };

   // Before interpolation: Ensure the new reference viscosity is set to the min current viscosity if a viscosity profile
   // or a temperature dependent viscosity is utilized.
   // This is will ensure that that the minimum non-dimensionalised value for the viscosity = 1.
   p2ScalarFunctionContainer["ViscosityFE"]->interpolate(
       viscosityInit, { *( p2ScalarFunctionContainer[std::string( "TemperatureFE" )] ) }, TN.domainParameters.maxLevel, All );

   real_t maxViscosity = p2ScalarFunctionContainer["ViscosityFE"]->getMaxDoFValue( TN.domainParameters.maxLevel ) *
                         TN.physicalParameters.referenceViscosity;

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Max viscosity [Pa s]: " << maxViscosity );

   real_t minRefViscosity = p2ScalarFunctionContainer["ViscosityFE"]->getMinDoFValue( TN.domainParameters.maxLevel ) *
                            TN.physicalParameters.referenceViscosity;
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
   auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

   if ( ( radius - TN.domainParameters.rMin ) < real_c( 1e-10 ) )
   {
      return ( TN.physicalParameters.cmbTemp ) / ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );
   }
   else if ( ( TN.domainParameters.rMax - radius ) < real_c( 1e-10 ) )
   {
      return ( TN.physicalParameters.surfaceTemp ) / ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );
   }

   if ( TN.simulationParameters.adaptiveRefTemp && TN.simulationParameters.timeStep > 0 )
   {
      uint_t shell = static_cast< uint_t >(
          std::round( real_c( TN.simulationParameters.numLayers ) *
                      ( ( radius - TN.domainParameters.rMin ) / ( TN.domainParameters.rMax - TN.domainParameters.rMin ) ) ) );
      WALBERLA_ASSERT( shell < TN.physicalParameters.temperatureProfile.size() );
      real_t retVal = TN.physicalParameters.temperatureProfile.at( shell );
      return retVal;
   }
   if ( !TN.simulationParameters.compressible && TN.simulationParameters.volAvrgTemperatureDev &&
        TN.simulationParameters.timeStep > 0 )
   {
      return TN.simulationParameters.avrgTemperatureVol;
   }
   else
   {
      if ( TN.simulationParameters.haveTemperatureProfile )
      {
         real_t temp = interpolateDataValues( x,
                                              TN.physicalParameters.radiusT,
                                              TN.physicalParameters.temperatureInputProfile,
                                              TN.domainParameters.rMin,
                                              TN.domainParameters.rMax );

         real_t retVal = ( temp ) / ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );
         return retVal;
      }
      else
      {
         real_t temp = TN.physicalParameters.adiabatSurfaceTemp *
                       std::exp( ( TN.physicalParameters.dissipationNumber * ( TN.domainParameters.rMax - radius ) ) );

         real_t retVal = ( temp ) / ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );

         return retVal;
      }
   }
}

void ConvectionSimulation::outputTimingTree()
{
   auto timer = storage->getTimingTree();
   writeTimingTreeJSON( *timer, TN.outputParameters.outputDirectory + "/" + "TimingTree.json" );
}

// SQL database for runtime analysis
void ConvectionSimulation::initTimingDB()
{
   auto gitHash = hyteg::buildinfo::gitSHA1();
   db           = std::make_shared< hyteg::FixedSizeSQLDB >( TN.outputParameters.fileNameSQLdb, true );

   db->setConstantEntry( "gitHash", gitHash );
   db->setConstantEntry( "cfl_max", TN.simulationParameters.cflMax );
   db->setConstantEntry( "unknownsTemperature", TN.simulationParameters.unknownsTemperature );
   db->setConstantEntry( "unknownsStokes", TN.simulationParameters.unknownsStokes );
   db->setConstantEntry( "h_min", TN.simulationParameters.hMin );
   db->setConstantEntry( "h_max", TN.simulationParameters.hMax );
   db->setConstantEntry( "numProcessors", TN.domainParameters.numProcessors );
   db->setConstantEntry( "nRad", TN.domainParameters.nRad );
   db->setConstantEntry( "nTAn", TN.domainParameters.nTan );
   db->setConstantEntry( "maxLevel", TN.domainParameters.maxLevel );
   db->setConstantEntry( "minLevel", TN.domainParameters.minLevel );
}
// This function calculates the average heatflow through the layer above the CMB and out of the Earth's surface
// It takes the mean temperature of the corresponding layers as input and calculates the heat flow
// as q = -k *  grad T * A [TW].
void ConvectionSimulation::calculateHeatflow( const std::shared_ptr< RadialProfile >& temperatureProfile )
{
   const real_t pi                  = walberla::math::pi;
   const real_t redimTemp           = TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp;
   const real_t thermalConductivity = TN.physicalParameters.thermalConductivity;
   const real_t mantleThickness     = TN.physicalParameters.mantleThickness;

   const std::vector< real_t >& radius    = temperatureProfile->shellRadii;
   const std::vector< real_t >& meanTemp  = temperatureProfile->mean;
   const uint_t                 numLayers = radius.size();

   const real_t areaSurface = 4 * pi * std::pow( radius[numLayers - 1], 2 );
   const real_t areaCMB     = 4 * pi * std::pow( radius[0], 2 );

   // Calculate heat flows
   const real_t heatFlowFactor = -thermalConductivity * mantleThickness * 1e-12;

   real_t heatFlowCMB = heatFlowFactor * ( meanTemp[1] - meanTemp[0] ) * redimTemp / ( radius[1] - radius[0] ) * areaCMB;

   real_t heatFlowSurface = heatFlowFactor * ( meanTemp[numLayers - 1] - meanTemp[numLayers - 2] ) * redimTemp /
                            ( radius[numLayers - 1] - radius[numLayers - 2] ) * areaSurface;

   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "Average heatflow CMB: " << heatFlowCMB << " TW" );
   WALBERLA_LOG_INFO_ON_ROOT( "Average heatflow Surface: " << heatFlowSurface << " TW" );
   if ( TN.outputParameters.createTimingDB )
   {
      db->setVariableEntry( "avrg_Heatflow_CMB_TW", heatFlowCMB );
      db->setVariableEntry( "avrg_Heatflow_Surface_TW", heatFlowSurface );
   }
}

// This function calculates the average heatflow through the layer above the CMB and out of the Earth's surface
// It takes the mean temperature of the corresponding layers as input and calculates the heat flow
// as q = -k *  grad T * A [TW].
void ConvectionSimulation::calculateHeatflowIntegral( const std::shared_ptr< RadialProfile >& temperatureProfile )
{
   const real_t pi                  = walberla::math::pi;
   const real_t redimTemp           = TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp;
   const real_t thermalConductivity = TN.physicalParameters.thermalConductivity;
   const real_t mantleThickness     = TN.physicalParameters.mantleThickness;

   const std::vector< real_t >& radius    = temperatureProfile->shellRadii;
   const std::vector< real_t >& meanTemp  = temperatureProfile->mean;
   const uint_t                 numLayers = radius.size();

   // Calculate heat flows
   const real_t heatFlowFactor = thermalConductivity * mantleThickness * 1e-12;

   real_t dTdr = ( meanTemp[1] - meanTemp[0] ) / ( radius[1] - radius[0] );

   uint_t nuSamples   = 101u;
   real_t hGradient   = 1e-2;
   real_t epsBoundary = 1e-7;

   real_t dTdrIntegralOuter = nusseltcalc::calculateNusseltNumberSphere3D( *( p2ScalarFunctionContainer["TemperatureFE"] ),
                                                                           TN.domainParameters.maxLevel,
                                                                           hGradient,
                                                                           TN.domainParameters.rMax,
                                                                           epsBoundary,
                                                                           nuSamples );

   real_t dTdrIntegralInner =
       nusseltcalc::calculateNusseltNumberSphere3D( *( p2ScalarFunctionContainer["TemperatureFE"] ),
                                                    TN.domainParameters.maxLevel,
                                                    hGradient,
                                                    TN.domainParameters.rMin + hGradient + 2.0 * epsBoundary,
                                                    epsBoundary,
                                                    nuSamples );

   real_t heatFlowSurface = heatFlowFactor * dTdrIntegralOuter * redimTemp;
   real_t heatFlowCMB     = heatFlowFactor * dTdrIntegralInner * redimTemp;

   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "Average heatflow CMB: " << heatFlowCMB << " TW" );
   WALBERLA_LOG_INFO_ON_ROOT( "Average heatflow Surface: " << heatFlowSurface << " TW" );
   if ( TN.outputParameters.createTimingDB )
   {
      db->setVariableEntry( "avrg_Heatflow_CMB_TW", heatFlowCMB );
      db->setVariableEntry( "avrg_Heatflow_Surface_TW", heatFlowSurface );
   }
}

} // namespace terraneo
