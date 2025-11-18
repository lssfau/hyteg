/*
 * Copyright (c) 2024 Eugenio D'Ascoli, Ponsuganth Ilangovan, Nils Kohl, Isabel Papanagnou.
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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
const SimulationParameters& ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::getSimulationParams()
{
   return TN.simulationParameters;
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::densityFunction( const Point3D& x )
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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::hydrostaticPressureFunction( const Point3D& x )
{
   auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
   real_t retVal;
   real_t rho_PDA;
   real_t hydrostaticPressure;
   
   // Adiabatic reference density
   if ( TN.simulationParameters.havePressureProfile )
   {
      hydrostaticPressure = interpolateDataValues( x,
                                      TN.physicalParameters.radiusPressure,
                                      TN.physicalParameters.pressureProfile,
                                      TN.domainParameters.rMin,
                                      TN.domainParameters.rMax );
   }
   else
   {
      if ( TN.simulationParameters.haveDensityProfile )
      {
         rho_PDA = interpolateDataValues( x,
                                       TN.physicalParameters.radiusDensity,
                                       TN.physicalParameters.densityProfile,
                                       TN.domainParameters.rMin,
                                       TN.domainParameters.rMax );
      }
      else
      {
         rho_PDA = TN.physicalParameters.surfaceDensity *
                  std::exp( TN.physicalParameters.dissipationNumber * ( TN.domainParameters.rMax - radius ) / TN.physicalParameters.grueneisenParameter );
       
      } 

      // Hydrostatic pressure
      hydrostaticPressure = rho_PDA * TN.physicalParameters.gravity * ( ( TN.domainParameters.rMax - radius ) * TN.physicalParameters.mantleThickness ) ;
   }

   //non-dimensionalise
   retVal = hydrostaticPressure * ( TN.physicalParameters.mantleThickness / (TN.physicalParameters.referenceViscosity * TN.physicalParameters.characteristicVelocity) );

   return retVal;
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::diffPreFactorFunction( const Point3D& x )
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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::adiabaticCoefficientFunction( const Point3D& x )
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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::constantEnergyCoefficientFunction( const Point3D& x )
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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::surfaceTempCoefficientFunction( const Point3D& )
{
   return 0.0;
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::updatePlateVelocities( P2P1StokesFunction_T& U )
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
         if ( !TN.simulationParameters.simulationInitialised )
         {
            p2p1StokesFunctionContainer["VelocityFEPrev"]->uvw()[coordIdx].interpolate( Velocity, l, idSurface );
         }
      }
   }
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::updateViscosity()
{
   std::shared_ptr< P2ScalarFunction_T >& temperatureP2  = p2ScalarFunctionContainer.at( "TemperatureFE" );
   std::shared_ptr< P2ScalarFunction_T >& viscosityP2    = p2ScalarFunctionContainer.at( "ViscosityFE" );
   std::shared_ptr< P2ScalarFunction_T >& viscosityInvP2 = p2ScalarFunctionContainer.at( "ViscosityFEInv" );

   std::shared_ptr< P0ScalarFunction_T >& viscosityP0 = p0ScalarFunctionContainer.at( "ViscosityFEP0" );

   if ( p2InjectionOperator_ == nullptr )
   {
      p2InjectionOperator_ =
          std::make_shared< P2toP2QuadraticInjection >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
   }

   for ( uint_t level = TN.domainParameters.maxLevel; level > TN.domainParameters.minLevel; level-- )
   {
      p2InjectionOperator_->restrict( *temperatureP2, level, All );
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
   viscosityP2->interpolate( viscosityInit, { *temperatureP2 }, TN.domainParameters.maxLevel, All );

   communication::syncFunctionBetweenPrimitives( viscosityP2->getVertexDoFFunction(), TN.domainParameters.maxLevel );

   P1toP0Conversion(
       viscosityP2->getVertexDoFFunction(), *viscosityP0, TN.domainParameters.maxLevel, AveragingType::ARITHMETIC_QP );

   P0toP0AveragedInjection p0toP0AveragedInjection( AveragingType::ARITHMETIC, false );

   p0toP0AveragedInjection.restrictToAllLowerLevels( *viscosityP0, TN.domainParameters.maxLevel );

   real_t maxViscosity = viscosityP2->getMaxDoFValue( TN.domainParameters.maxLevel ) * TN.physicalParameters.referenceViscosity;

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Max viscosity [Pa s]: " << maxViscosity );

   real_t minRefViscosity =
       viscosityP2->getMinDoFValue( TN.domainParameters.maxLevel ) * TN.physicalParameters.referenceViscosity;
   WALBERLA_LOG_INFO_ON_ROOT( "New update reference viscosity [Pa s]: " << minRefViscosity );

   if ( TN.simulationParameters.tempDependentViscosity || TN.simulationParameters.haveViscosityProfile )
   {
      // Update reference viscosity with new min Viscosity value
      TN.physicalParameters.referenceViscosity = minRefViscosity;
      // TN.physicalParameters.rayleighNumber = ( TN.physicalParameters.referenceDensity * TN.physicalParameters.gravity * TN.physicalParameters.thermalExpansivity *
      // std::pow( TN.physicalParameters.mantleThickness, 3 ) *
      // ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) ) /
      // ( TN.physicalParameters.referenceViscosity * TN.physicalParameters.thermalDiffusivity );

      if ( TN.simulationParameters.pdaForm == true )
      {
         TN.physicalParameters.rayleighNumberPDA =
          ( TN.physicalParameters.referenceHydDensity * TN.physicalParameters.gravity * 
            std::pow( TN.physicalParameters.mantleThickness, 3 ) ) /
          ( TN.physicalParameters.referenceViscosity * TN.physicalParameters.thermalDiffusivity );

      }
      else
      {
         TN.physicalParameters.rayleighNumber =
             ( TN.physicalParameters.referenceDensity * TN.physicalParameters.gravity * TN.physicalParameters.thermalExpansivity *
               std::pow( TN.physicalParameters.mantleThickness, 3 ) *
               ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) ) /
             ( TN.physicalParameters.referenceViscosity * TN.physicalParameters.thermalDiffusivity );
      }
   }

   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      viscosityP2->interpolate( viscosityInit, { *temperatureP2 }, l, All );
      viscosityInvP2->interpolate( viscosityInitInv, { *temperatureP2 }, l, All );
   }
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::normalFunc( const Point3D& p, Point3D& n )
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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::oppositeGravity( const Point3D& p, Point3D& n )
{
   real_t radius = p.norm();
   n             = Point3D( { p[0] / radius, p[1] / radius, p[2] / radius } );
}

//returns a reference adiabat relevant for the Earth, commonly implemented in TALA
template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::referenceTemperatureFunction( const Point3D& x )
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

   if ( TN.simulationParameters.adaptiveRefTemp && TN.simulationParameters.simulationInitialised )
   {
      uint_t shell = uint_c(
          std::round( real_c( TN.simulationParameters.numLayers ) *
                      ( ( radius - TN.domainParameters.rMin ) / ( TN.domainParameters.rMax - TN.domainParameters.rMin ) ) ) );
      WALBERLA_ASSERT( shell < TN.physicalParameters.temperatureProfile.size() );
      real_t retVal = TN.physicalParameters.temperatureProfile.at( shell );
      return retVal;
   }
   if ( !TN.simulationParameters.compressible && TN.simulationParameters.volAvrgTemperatureDev &&
        TN.simulationParameters.simulationInitialised )
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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::outputTimingTree()
{
   auto timer = storage->getTimingTree();
   writeTimingTreeJSON( *timer, TN.outputParameters.outputDirectory + "/" + "TimingTree.json" );
}

// SQL database for runtime analysis
template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::initTimingDB()
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
template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::calculateHeatflow(
    const std::shared_ptr< RadialProfile >& temperatureProfile )
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
template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::calculateHeatflowIntegral(
    const std::shared_ptr< RadialProfile >& temperatureProfile )
{
   const real_t pi                  = walberla::math::pi;
   const real_t redimTemp           = TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp;
   const real_t thermalConductivity = TN.physicalParameters.thermalConductivity;
   const real_t mantleThickness     = TN.physicalParameters.mantleThickness;

   const std::vector< real_t >& radius    = temperatureProfile->shellRadii;
   const std::vector< real_t >& meanTemp  = temperatureProfile->mean;
   const uint_t                 numLayers = radius.size();

   std::shared_ptr< P2ScalarFunction_T >& temperatureP2 = p2ScalarFunctionContainer.at( "TemperatureFE" );

   // Calculate heat flows
   const real_t heatFlowFactor = thermalConductivity * mantleThickness * 1e-12;

   real_t dTdr = ( meanTemp[1] - meanTemp[0] ) / ( radius[1] - radius[0] );

   uint_t nuSamples   = 101u;
   real_t hGradient   = 1e-2;
   real_t epsBoundary = 1e-7;

   real_t dTdrIntegralOuter = nusseltcalc::calculateNusseltNumberSphere3D(
       *temperatureP2, TN.domainParameters.maxLevel, hGradient, TN.domainParameters.rMax, epsBoundary, nuSamples );

   real_t dTdrIntegralInner =
       nusseltcalc::calculateNusseltNumberSphere3D( *temperatureP2,
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

std::function< real_t( const Point3D&, const std::vector< real_t >& ) > DoFCounter( real_t threshold )
{
   return [threshold]( const Point3D&, const std::vector< real_t >& values ) { return ( values[0] < threshold ) ? 1.0 : 0.0; };
}
   
//template < typename TemperatureFunction_T, typename ViscosityFunction_T >
//real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::getTableValue( std::vector< std::vector< real_t > >& paramPTS , real_t& fullPressure , real_t& fullTemp )
real_t getTableValue( std::vector< std::vector< real_t > >& paramPTS , real_t& fullPressure , real_t& fullTemp )
{
   // Search for values in the pressure and temperature vectors that are closest to the given P,T pair, 
   // return the corresponding indices iP and iTS, respectively, and retrieve the associated parameter value

   uint_t searchParam = 50;
   uint_t nTsteps = TN.simulationParameters.numTsteps / searchParam;
   uint_t nPsteps = TN.simulationParameters.numPsteps / searchParam;
   real_t retVal = 0;

   // Ensuring that values remain within the bounds of the thermodynamic table
   if ( uint_t( fullTemp ) > TN.physicalParameters.temperatureRangePT[0][TN.simulationParameters.numTsteps - 1] )
   {
      fullTemp = TN.physicalParameters.temperatureRangePT[0][TN.simulationParameters.numTsteps - 1];
   }
   if ( uint_t( fullTemp ) < TN.physicalParameters.temperatureRangePT[0][0] )
   {
      fullTemp = TN.physicalParameters.temperatureRangePT[0][0];
   }
   
   for ( uint_t i = 0; i < nPsteps; i++ )
   {      
      uint_t n_ = i * searchParam;
      uint_t n = ( i + 1 ) * searchParam;

      if ( TN.physicalParameters.pressureRangePT[n][0] >= fullPressure )
      {
         for ( uint_t k = n_; k <= n; k++ )
         {
            if ( TN.physicalParameters.pressureRangePT[k][0] >= fullPressure )
            {
               if ( k != 0 )
               {
                  TN.simulationParameters.iP_min = k - 1;
                  TN.simulationParameters.iP_max = k;
               }
               else
               {
                  TN.simulationParameters.iP_min = k;
                  TN.simulationParameters.iP_max = k + 1;
               }               
               real_t diffP_max = std::abs( TN.physicalParameters.pressureRangePT[TN.simulationParameters.iP_max][0] - fullPressure );
               real_t diffP_min = std::abs( TN.physicalParameters.pressureRangePT[TN.simulationParameters.iP_min][0] - fullPressure );
                           
               if ( diffP_max <= diffP_min )
               {
                  TN.simulationParameters.iP = TN.simulationParameters.iP_max;
               }
               else
               {
                  TN.simulationParameters.iP = TN.simulationParameters.iP_min;
               }

               for ( uint_t j = 0; j < nTsteps; j++ )
               {
                  uint_t m_ = j * searchParam;
                  uint_t m = ( j + 1 ) * searchParam;

                  if ( TN.physicalParameters.temperatureRangePT[0][m] >= fullTemp )
                  { 
                     for ( uint_t o = m_; o <= m; o++ )
                     {                        
                        if ( TN.physicalParameters.temperatureRangePT[0][o] >= fullTemp )
                        {      
                           if ( o != 0 )
                           {
                              TN.simulationParameters.iTS_min = o - 1;
                              TN.simulationParameters.iTS_max = o;
                           }
                           else
                           {
                              TN.simulationParameters.iTS_min = o;
                              TN.simulationParameters.iTS_max = o + 1;
                           } 
                           real_t diffTS_max = std::abs( TN.physicalParameters.temperatureRangePT[0][TN.simulationParameters.iTS_max] - fullTemp );
                           real_t diffTS_min = std::abs( TN.physicalParameters.temperatureRangePT[0][TN.simulationParameters.iTS_min] - fullTemp );
                           
                           if ( diffTS_max <= diffTS_min )
                           {
                              TN.simulationParameters.iTS = TN.simulationParameters.iTS_max;
                           }
                           else
                           {
                              TN.simulationParameters.iTS = TN.simulationParameters.iTS_min;
                           }
                           retVal = paramPTS[TN.simulationParameters.iP][TN.simulationParameters.iTS];
                           break; 
                        }
                     }
                     break;
                  }
               }
               break; 
            }
         }
         break;
      }
   }
   return retVal;
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::updatePDARHSOperators()
{
   std::shared_ptr< P2ScalarFunction_T >& hydDensityP2                  = p2ScalarFunctionContainer.at( "nondimHydrostaticDensityPTFE" );
   std::shared_ptr< P1ScalarFunction_T >& hydDensityP1Prev              = p1ScalarFunctionContainer.at( "P1nondimHydrostaticDensityPTFEPrev" );
   std::shared_ptr< P1ScalarFunction_T >& hydDensityP1Prev2             = p1ScalarFunctionContainer.at( "P1nondimHydrostaticDensityPTFEPrev2" );
   std::shared_ptr< P1ScalarFunction_T >& hydDensityCoeffTimDerivative  = p1ScalarFunctionContainer.at( "P1timeDerivativeHydDensityCoeff" );

   updateThermodynamicParams( All );

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > getTimeDerivativeHydDensityBDF2 = 
      [=]( const Point3D& , const std::vector< real_t >& HydrostaticDensityFE ) {
         real_t timeDerivativeHydDensity = ( ( ( 
            ( real_c( 2 ) * TN.simulationParameters.dt + TN.simulationParameters.dtPrev ) * HydrostaticDensityFE[0] / ( TN.simulationParameters.dt + TN.simulationParameters.dtPrev ) ) 
            - ( ( TN.simulationParameters.dt + TN.simulationParameters.dtPrev ) * HydrostaticDensityFE[1] / ( TN.simulationParameters.dtPrev ) ) 
            + ( ( TN.simulationParameters.dt * TN.simulationParameters.dt ) * HydrostaticDensityFE[2] / ( TN.simulationParameters.dtPrev * ( TN.simulationParameters.dt + TN.simulationParameters.dtPrev ) ) ) ) 
            / ( TN.simulationParameters.dt * HydrostaticDensityFE[0] ) );
         return timeDerivativeHydDensity;
      };

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > getTimeDerivativeHydDensity = 
      [=]( const Point3D& , const std::vector< real_t >& HydrostaticDensityFE ) {
         real_t timeDerivativeHydDensity = ( ( HydrostaticDensityFE[0] - HydrostaticDensityFE[1] ) / ( TN.simulationParameters.dt * HydrostaticDensityFE[0] ) );
         return timeDerivativeHydDensity;
      };

   if ( TN.simulationParameters.timeStep > 1 )
   {
      hydDensityCoeffTimDerivative->interpolate( getTimeDerivativeHydDensityBDF2, { hydDensityP2->getVertexDoFFunction(),
                                                   *( hydDensityP1Prev ),
                                                   *( hydDensityP1Prev2 ) },
                                                   TN.domainParameters.maxLevel, All );
   }
   else
   {
      hydDensityCoeffTimDerivative->interpolate( getTimeDerivativeHydDensity, { hydDensityP2->getVertexDoFFunction(), *( hydDensityP1Prev ) }, TN.domainParameters.maxLevel, All );
   }

   P1frozenVelocityPDA = std::make_shared< P1FrozenVelocityFullOperator_T >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, hydDensityP2->getVertexDoFFunction() );
   P1timeDerivativePDA = std::make_shared< P1TimeDerivativeHydDensityOperator_T >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *( hydDensityCoeffTimDerivative ) );
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::ConvectionSimulation::updateThermodynamicParams( hyteg::DoFType flag )
{

   std::shared_ptr< P2ScalarFunction_T >& temperatureP2             = p2ScalarFunctionContainer.at( "TemperatureFE" );
   std::shared_ptr< P2ScalarFunction_T >& dimHydDensityP2           = p2ScalarFunctionContainer.at( "dimHydrostaticDensityPTFE" );
   std::shared_ptr< P2ScalarFunction_T >& hydDensityP2              = p2ScalarFunctionContainer.at( "nondimHydrostaticDensityPTFE" );

   // Hydrostatic pressure (non-dimensionalised value)
   std::function<real_t( const Point3D& ) > getHydrostaticPressure = [&](const Point3D& x )
   {
      real_t retVal =  hydrostaticPressureFunction( x );
      return retVal;
   };

   // Hydrostatic density 
   std::function<real_t( const Point3D&, const std::vector< real_t >& ) > getHydrostaticDensity = [&](const Point3D& x, const std::vector< real_t >& Temperature )
   {
      // Redimensionalising temperature and (full) pressure for the thermodynamic table search
      real_t fullHydrostaticPressure = TN.physicalParameters.referenceViscosity * TN.physicalParameters.characteristicVelocity * hydrostaticPressureFunction( x ) / TN.physicalParameters.mantleThickness;
      real_t retVal;

      // Getting parameter from thermodynamic table in dimensional form
      real_t fullTemp = ( ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) * Temperature[0] );
      retVal = getTableValue( TN.physicalParameters.materialDensityPT, fullHydrostaticPressure, fullTemp );
      
      return retVal;
   };

   // Obtaining finite element fields for relevant thermodynamic parameters for the current (updated) temperature and pressure fields
   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      // Dimensional parameter FE fields
      // Dimensional thermodynamic parameters (P,T)
      dimHydDensityP2->interpolate( getHydrostaticDensity, { *( temperatureP2 ) }, l, flag );
   }
   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      // Non-dimensional parameter FE fields
      // Non-dimensionalsed thermodynamic parameters (P,T) - P2
      hydDensityP2->assign( { real_c( 1 ) / TN.physicalParameters.referenceHydDensity }, {*( dimHydDensityP2 )}, l, flag );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Updated thermodynamic parameters" );
   
   if ( TN.simulationParameters.updateThermodynamicReferenceValues == true )
   {
      // Updating reference values for thermodynamic parameters
      TN.physicalParameters.referenceHydDensity = dimHydDensityP2->getMinDoFValue( TN.domainParameters.maxLevel );
      
      // Updating the relevant non-dimensional numbers and 'helper quantities' with the updated reference values 
      TN.physicalParameters.rayleighNumberPDA =
         ( TN.physicalParameters.referenceHydDensity * TN.physicalParameters.gravity *  std::pow( TN.physicalParameters.mantleThickness, 3 ) ) /
         ( TN.physicalParameters.referenceViscosity * TN.physicalParameters.thermalDiffusivity );

      //WALBERLA_LOG_INFO_ON_ROOT( "Thermodynamic reference values and non-dimensional numbers updated" );
      WALBERLA_LOG_INFO_ON_ROOT( "PDA Rayleigh number           : " << TN.physicalParameters.rayleighNumberPDA );
      WALBERLA_LOG_INFO_ON_ROOT( "Hydrostatic density reference : " << TN.physicalParameters.referenceHydDensity );
   };
}

} // namespace terraneo
