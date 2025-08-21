/*
 * Copyright (c) 2024 Eugenio D'Ascoli.
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

#include <core/DataTypes.h>
#include <core/config/Config.h>
#include <core/logging/Logging.h>
#include <core/mpi/MPIManager.h>

#include "hyteg/Levelinfo.hpp"

#include "terraneo/dataimport/FileIO.hpp"
#include "terraneo/helpers/TerraNeoParameters.hpp"
#include "terraneo/plates/PlateVelocityProvider.hpp"
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

namespace terraneo {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

struct ParameterFileVersion
{
   uint_t major_ = 0;
   uint_t minor_ = 1;

   ParameterFileVersion( uint_t major, uint_t minor )
   : major_( major )
   , minor_( minor )
   {}

   ParameterFileVersion( const std::string& version )
   {
      size_t dotPos = version.find( '.' );

      if ( dotPos == std::string::npos )
      {
         WALBERLA_ABORT( "Invalid Parameter File Version string: " << version );
      }

      try
      {
         major_ = std::stoi( version.substr( 0, dotPos ) );
         minor_ = std::stoi( version.substr( dotPos + 1 ) );
      } catch ( const std::invalid_argument& e )
      {
         WALBERLA_ABORT( "Invalid Parameter File Version string: " << version << ". Error: " << e.what() );
      } catch ( const std::out_of_range& e )
      {
         WALBERLA_ABORT( "Parameter File Version number out of range: " << version << ". Error: " << e.what() );
      }
   }

   bool operator==( const ParameterFileVersion& other ) const { return major_ == other.major_ && minor_ == other.minor_; }

   bool operator!=( const ParameterFileVersion& other ) const { return !( *this == other ); }

   bool operator<( const ParameterFileVersion& other ) const
   {
      return ( major_ < other.major_ ) || ( major_ == other.major_ && minor_ < other.minor_ );
   }

   bool operator<=( const ParameterFileVersion& other ) const { return *this < other || *this == other; }

   bool operator>( const ParameterFileVersion& other ) const { return !( *this <= other ); }

   bool operator>=( const ParameterFileVersion& other ) const { return !( *this < other ); }
};

// Helper function to convert vector to string for logging
template < typename T >
std::string vectorToString( const std::vector< T >& vec )
{
   std::ostringstream oss;
   for ( size_t i = 0; i < vec.size(); ++i )
   {
      oss << vec[i];
      if ( i != vec.size() - 1 )
         oss << ", ";
   }
   return oss.str();
}

/**
 * @brief Loads the radial profiles if provided via the main configuration block.
 *
 * This function reads and extracts radial profile data from the main configuration block.
 * It parses the radial profiles and stores the values in the correpsonding std::vectors. 
 * These vectors are later populated to the application where the corresponding value can be 
 * extracted for each radial layer and be further processed. 
 */

template < typename T >
inline void loadRadialProfile( const walberla::Config::BlockHandle& mainConf,
                               const std::string&                   profileKey,
                               T&                                   radiusData,
                               T&                                   profileData,
                               bool&                                haveProfileFlag,
                               const real_t                         mantleThickness )
{
   SimulationParameters simulationParam;

   if ( mainConf.isDefined( profileKey ) )
   {
      const std::string filePath = mainConf.getParameter< std::string >( profileKey );
      auto              jsonData = io::readJsonFile( filePath );

      const auto radiusKey = "Radius (m)";
      const auto valueKey  = simulationParam.getKeyValue( jsonData );

      WALBERLA_CHECK_GREATER( jsonData.count( radiusKey ), 0, "No key '" << radiusKey << "' in " << profileKey << " file." );
      WALBERLA_CHECK_GREATER( jsonData.count( valueKey ), 0, "No key '" << valueKey << "' in " << profileKey << " file." );

      radiusData  = jsonData[radiusKey].get< std::vector< real_t > >();
      profileData = jsonData[valueKey].get< std::vector< real_t > >();

      WALBERLA_CHECK_EQUAL(
          radiusData.size(), profileData.size(), "Mismatch in size for " << profileKey << " radius and values." );

      for ( auto& radius : radiusData )
      {
         radius /= mantleThickness; // Non-dimensionalize radius
      }

      haveProfileFlag = true;
   }
}

/**
 * @brief Parses the configuration parameters from the main configuration block.
 *
 * This function reads and extracts various domain, model, simulation, and initialization parameters from the main configuration block.
 * It populates the corresponding variables and performs necessary calculations for non-dimensional numbers.
 */

// define plate velocity oracle
std::shared_ptr< terraneo::plates::PlateVelocityProvider > oracle;
// define plate velocity averaging point provider
std::shared_ptr< terraneo::plates::UniformCirclesPointWeightProvider > avgPointProvider;

inline TerraNeoParameters parseConfig( const walberla::Config::BlockHandle& mainConf )
{
   /* ParameterFileVersion */
   const ParameterFileVersion parameterFileVerion =
       mainConf.isDefined( "version" ) ? ParameterFileVersion( mainConf.getParameter< std::string >( "version" ) ) :
                                         ParameterFileVersion( 0, 1 );

   DomainParameters                             domainParam;
   SolverParameters                             solverParam;
   OutputParameters                             outputParam;
   SimulationParameters                         simulationParam;
   PhysicalParameters                           physicalParam;
   TemperatureDeivationInitialisationParameters initialisationParam;

   /*############ DOMAIN PARAMETERS ############*/

   domainParam.rCMB     = mainConf.getParameter< real_t >( "rCMB" );
   domainParam.rSurface = mainConf.getParameter< real_t >( "rSurface" );
   domainParam.nTan     = mainConf.getParameter< uint_t >( "nTan" );
   domainParam.nRad     = mainConf.getParameter< uint_t >( "nRad" );
   domainParam.minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   domainParam.maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   //calculate non-dimensional radii such that mantle thickness = 1
   physicalParam.mantleThickness = domainParam.rSurface - domainParam.rCMB;
   domainParam.rMin              = domainParam.rCMB / ( domainParam.rSurface - domainParam.rCMB );
   domainParam.rMax              = domainParam.rSurface / ( domainParam.rSurface - domainParam.rCMB );

   // Read macro layer radii from parameter file
   if ( mainConf.isDefined( "inputMacroLayers" ) )
   {
      std::string macroLayersStr;
      std::string layerRadius;

      macroLayersStr = mainConf.getParameter< std::string >( "inputMacroLayers" );
      std::stringstream ss( macroLayersStr );

      // Parse input string to vector
      while ( std::getline( ss, layerRadius, ',' ) )
      {
         // Trim leading/trailing spaces
         layerRadius.erase( 0, layerRadius.find_first_not_of( " \t" ) );
         layerRadius.erase( layerRadius.find_last_not_of( " \t" ) + 1 );
         domainParam.macroLayers.push_back( real_c( std::stof( layerRadius ) ) );
      }
      // non-dimensionalize the boundary radii
      for ( real_t& entry : domainParam.macroLayers )
      {
         entry = entry / ( domainParam.rSurface - domainParam.rCMB );
      }

      // Sort macroLayers to ensure ordering from bottom to top boundary
      std::sort( domainParam.macroLayers.begin(), domainParam.macroLayers.end() );

      // Add outer boundaries
      domainParam.macroLayers.insert( domainParam.macroLayers.begin(), domainParam.rMin );
      domainParam.macroLayers.push_back( domainParam.rMax );

      // Override nRad
      domainParam.nRad = static_cast< uint_t >( domainParam.macroLayers.size() );
   }
   // else determine equally spaced layer radii from nRad
   else
   {
      real_t layerRadius;

      for ( uint_t layer = 0; layer < domainParam.nRad; layer++ )
      {
         layerRadius =
             domainParam.rMin + real_c( layer ) * ( domainParam.rMax - domainParam.rMin ) / real_c( domainParam.nRad - 1 );
         domainParam.macroLayers.push_back( layerRadius );
      }
   }

   /*############ MODEL PARAMETERS ############*/
   if ( mainConf.isDefined( "temperatureInputProfile" ) )
   {
      simulationParam.fileTemperatureInputProfile = mainConf.getParameter< std::string >( "temperatureInputProfile" );
      const std::string profileKey                = "temperatureInputProfile";

      loadRadialProfile< std::vector< real_t > >( mainConf,
                                                  profileKey,
                                                  physicalParam.radiusT,
                                                  physicalParam.temperatureInputProfile,
                                                  simulationParam.haveTemperatureProfile,
                                                  physicalParam.mantleThickness );
   }

   if ( mainConf.isDefined( "viscosityProfile" ) )
   {
      simulationParam.fileViscosityProfile = mainConf.getParameter< std::string >( "viscosityProfile" );
      const std::string profileKey         = "viscosityProfile";

      loadRadialProfile< std::vector< real_t > >( mainConf,
                                                  profileKey,
                                                  physicalParam.radius,
                                                  physicalParam.viscosityProfile,
                                                  simulationParam.haveViscosityProfile,
                                                  physicalParam.mantleThickness );
   }

   // Profile for thermal expansivity
   if ( mainConf.isDefined( "thermalExpansivityProfile" ) )
   {
      simulationParam.fileThermalExpProfile = mainConf.getParameter< std::string >( "thermalExpansivityProfile" );
      const std::string profileKey          = "thermalExpansivityProfile";

      loadRadialProfile< std::vector< real_t > >( mainConf,
                                                  profileKey,
                                                  physicalParam.radiusAlpha,
                                                  physicalParam.thermalExpansivityProfile,
                                                  simulationParam.haveThermalExpProfile,
                                                  physicalParam.mantleThickness );
   }

   // Radial Profile for specific heat capacity at constant pressur

   if ( mainConf.isDefined( "specificHeatCapacityProfile" ) )
   {
      simulationParam.fileSpecificHeatCap = mainConf.getParameter< std::string >( "specificHeatCapacityProfile" );
      const std::string profileKey        = "specificHeatCapacityProfile";

      loadRadialProfile< std::vector< real_t > >( mainConf,
                                                  profileKey,
                                                  physicalParam.radiusCp,
                                                  physicalParam.specificHeatCapacityProfile,
                                                  simulationParam.haveSpecificHeatCapProfile,
                                                  physicalParam.mantleThickness );
   }

   if ( mainConf.isDefined( "densityProfile" ) )
   {
      simulationParam.fileDensityProfile = mainConf.getParameter< std::string >( "densityProfile" );
      const std::string profileKey       = "densityProfile";

      loadRadialProfile< std::vector< real_t > >( mainConf,
                                                  profileKey,
                                                  physicalParam.radiusDensity,
                                                  physicalParam.densityProfile,
                                                  simulationParam.haveDensityProfile,
                                                  physicalParam.mantleThickness );
   }

   physicalParam.cmbTemp                = mainConf.getParameter< real_t >( "cmbTemp" );
   physicalParam.surfaceTemp            = mainConf.getParameter< real_t >( "surfaceTemp" );
   physicalParam.thermalExpansivity     = mainConf.getParameter< real_t >( "thermalExpansivity" );
   physicalParam.thermalConductivity    = mainConf.getParameter< real_t >( "thermalConductivity" );
   physicalParam.grueneisenParameter    = mainConf.getParameter< real_t >( "grueneisenParameter" );
   physicalParam.specificHeatCapacity   = mainConf.getParameter< real_t >( "specificHeatCapacity" );
   physicalParam.internalHeatingRate    = mainConf.getParameter< real_t >( "internalHeatingRate" );
   physicalParam.characteristicVelocity = mainConf.getParameter< real_t >( "characteristicVelocity" );
   physicalParam.surfaceDensity         = mainConf.getParameter< real_t >( "surfaceDensity" );
   physicalParam.referenceDensity       = mainConf.getParameter< real_t >( "referenceDensity" );
   physicalParam.viscosity              = mainConf.getParameter< real_t >( "viscosity" );
   physicalParam.referenceViscosity     = physicalParam.viscosity;

   // Set all radial varying parameters to input reference values to avoid inconsistent calculations on non-dim numbers
   physicalParam.specificHeatCapacityRadial = physicalParam.specificHeatCapacity;
   physicalParam.thermalExpansivityRadial   = physicalParam.thermalExpansivity;
   physicalParam.characteristicVelocity =
       physicalParam.thermalConductivity /
       ( physicalParam.referenceDensity * physicalParam.specificHeatCapacity * physicalParam.mantleThickness );
   //used to calculate non-D numbers
   physicalParam.thermalDiffusivity =
       physicalParam.thermalConductivity / ( physicalParam.referenceDensity * physicalParam.specificHeatCapacity );

   //calculating non-D numbers

   physicalParam.rayleighNumber =
       ( physicalParam.referenceDensity * physicalParam.gravity * physicalParam.thermalExpansivity *
         std::pow( physicalParam.mantleThickness, 3 ) * ( physicalParam.cmbTemp - physicalParam.surfaceTemp ) ) /
       ( physicalParam.referenceViscosity * physicalParam.thermalDiffusivity );

   physicalParam.pecletNumber =
       ( physicalParam.characteristicVelocity * physicalParam.mantleThickness ) / physicalParam.thermalDiffusivity;

   physicalParam.dissipationNumber =
       ( physicalParam.thermalExpansivity * physicalParam.gravity * physicalParam.mantleThickness ) /
       physicalParam.specificHeatCapacity;

   physicalParam.hNumber = ( physicalParam.internalHeatingRate * physicalParam.mantleThickness ) /
                           ( physicalParam.specificHeatCapacity * physicalParam.characteristicVelocity *
                             ( physicalParam.cmbTemp - physicalParam.surfaceTemp ) );

   /*############ SIMULATION PARAMETERS ############*/

   simulationParam.cflMax                 = mainConf.getParameter< real_t >( "cflMax" );
   simulationParam.fixedTimestep          = mainConf.getParameter< bool >( "fixedTimestep" );
   simulationParam.maxTimestepSize        = mainConf.getParameter< real_t >( "maxTimestepSize" );
   simulationParam.dtConstant             = mainConf.getParameter< real_t >( "dtConstant" );
   simulationParam.predictorCorrector     = mainConf.getParameter< bool >( "predictorCorrector" );
   simulationParam.maxNumTimesteps        = mainConf.getParameter< uint_t >( "maxNumTimesteps" );
   simulationParam.adaptiveRefTemp        = mainConf.getParameter< bool >( "adaptiveRefTemp" );
   simulationParam.volAvrgTemperatureDev  = mainConf.getParameter< bool >( "volAvrgTemperatureDev" );
   simulationParam.tempDependentViscosity = mainConf.getParameter< bool >( "tempDependentViscosity" );
   simulationParam.simulationType         = mainConf.getParameter< std::string >( "simulationType" );
   simulationParam.compressible           = mainConf.getParameter< bool >( "compressible" );
   simulationParam.frozenVelocity         = mainConf.getParameter< bool >( "frozenVelocity" );
   simulationParam.shearHeating           = mainConf.getParameter< bool >( "shearHeating" );
   if ( simulationParam.shearHeating )
   {
      if ( parameterFileVerion < ParameterFileVersion( 0, 2 ) )
      {
         simulationParam.lithosphereShearHeatingScaling = mainConf.getParameter< real_t >( "shearHeatingScaling" );
      }
      else
      {
         simulationParam.lithosphereShearHeatingScaling = mainConf.getParameter< real_t >( "lithosphereShearHeatingScaling" );
      }
      simulationParam.lithosphereThickness = mainConf.getParameter< real_t >( "lithosphereThickness" );
   }
   simulationParam.adiabaticHeating = mainConf.getParameter< bool >( "adiabaticHeating" );
   simulationParam.internalHeating  = mainConf.getParameter< bool >( "internalHeating" );
   simulationParam.boundaryCond     = mainConf.getParameter< uint_t >( "boundaryCond" );
   simulationParam.timingAnalysis   = mainConf.getParameter< bool >( "timingAnalysis" );

   simulationParam.checkTemperatureConsistency     = mainConf.getParameter< bool >( "checkTemperatureConsistency" );
   simulationParam.temperatureConsistencyThreshold = mainConf.getParameter< real_t >( "temperatureConsistencyThreshold" );

   simulationParam.particleLocationRadius    = mainConf.getParameter< real_t >( "particleLocationRadius" );
   simulationParam.projectPointsBackToDomain = mainConf.getParameter< bool >( "projectPointsBackToDomain" );
   simulationParam.cautionedEvaluate         = mainConf.getParameter< bool >( "cautionedEvaluate" );

   if ( simulationParam.tempDependentViscosity || simulationParam.haveViscosityProfile )
   {
      simulationParam.tempDependentViscosityType = mainConf.getParameter< uint_t >( "tempDependentViscosityType" );
      physicalParam.activationEnergy             = mainConf.getParameter< real_t >( "activationEnergy" );
      physicalParam.depthViscosityFactor         = mainConf.getParameter< real_t >( "depthViscosityFactor" );
      physicalParam.viscosityLowerBound          = mainConf.getParameter< real_t >( "viscosityLowerBound" );
      physicalParam.viscosityUpperBound          = mainConf.getParameter< real_t >( "viscosityUpperBound" );
   }

   simulationParam.finalAge = mainConf.getParameter< real_t >( "finalAge" );
   //simulation parameters for circulation models only:
   if ( simulationParam.simulationType == "CirculationModel" )
   {
      if ( mainConf.isDefined( "fnameTopologies" ) && mainConf.isDefined( "fnameReconstructions" ) )
      {
         simulationParam.fnameTopologies        = mainConf.getParameter< std::string >( "fnameTopologies" );
         simulationParam.fnameReconstructions   = mainConf.getParameter< std::string >( "fnameReconstructions" );
         simulationParam.initialAge             = mainConf.getParameter< real_t >( "initialAge" );
         simulationParam.ageMa                  = simulationParam.initialAge;
         simulationParam.agePrev                = simulationParam.initialAge;
         simulationParam.plateAge               = simulationParam.initialAge;
         simulationParam.plateVelocityScaling   = mainConf.getParameter< real_t >( "plateVelocityScaling" );
         simulationParam.plateSmoothingDistance = mainConf.getParameter< real_t >( "plateSmoothingDistance" );
      }

      else
      {
         WALBERLA_ABORT(
             "To run a circulation model, both a plate topology file and rotation file need to be specified in the parameter file. Aborting." );
      }
   }

   /*############ INITIALISATION PARAMETERS ############*/

   // Check parameter file ParameterFileVersion
   uint_t numHarmonics;
   if ( parameterFileVerion < ParameterFileVersion( 0, 2 ) )
   {
      bool temperatureNoise = mainConf.getParameter< bool >( "temperatureNoise" );
      bool superposition    = mainConf.getParameter< bool >( "superposition" );

      if ( temperatureNoise )
      {
         initialisationParam.initialTemperatureDeviationMethod = INITIAL_TEMPERATURE_DEVIATION_METHOD::WHITE_NOISE;
      }
      else if ( superposition )
      {
         initialisationParam.initialTemperatureDeviationMethod = INITIAL_TEMPERATURE_DEVIATION_METHOD::RANDOM_SUPERPOSITION_SPH;
      }
      else
      {
         initialisationParam.initialTemperatureDeviationMethod = INITIAL_TEMPERATURE_DEVIATION_METHOD::SINGLE_SPH;
      }

      switch ( initialisationParam.initialTemperatureDeviationMethod )
      {
      case INITIAL_TEMPERATURE_DEVIATION_METHOD::RANDOM_SUPERPOSITION_SPH:
         initialisationParam.tempInit                    = mainConf.getParameter< uint_t >( "tempInit" );
         initialisationParam.initialTemperatureSteepness = mainConf.getParameter< real_t >( "initialTemperatureSteepness" );
         initialisationParam.lmax                        = mainConf.getParameter< uint_t >( "lmax" );
         initialisationParam.lmin                        = mainConf.getParameter< uint_t >( "lmin" );
         initialisationParam.buoyancyFactor              = mainConf.getParameter< real_t >( "buoyancyFactor" );

         initialisationParam.sphTool = std::make_shared< SphericalHarmonicsTool >( initialisationParam.lmax );

         numHarmonics = ( ( initialisationParam.lmax + 1 ) * ( initialisationParam.lmax + 1 ) ) -
                        ( initialisationParam.lmin ) * ( initialisationParam.lmin );
         initialisationParam.superpositionRand.reserve( numHarmonics );
         walberla::math::seedRandomGenerator( initialisationParam.superpositionSPHRandomSeed );

         for ( uint_t i = 0; i < numHarmonics; i++ )
         {
            initialisationParam.superpositionRand.push_back( walberla::math::realRandom( real_c( -1 ), real_c( 1 ) ) );
         }
         break;
      case INITIAL_TEMPERATURE_DEVIATION_METHOD::SINGLE_SPH:
         initialisationParam.tempInit                    = mainConf.getParameter< uint_t >( "tempInit" );
         initialisationParam.initialTemperatureSteepness = mainConf.getParameter< real_t >( "initialTemperatureSteepness" );
         initialisationParam.deg                         = mainConf.getParameter< uint_t >( "degree" );
         initialisationParam.ord                         = mainConf.getParameter< int >( "order" );
         initialisationParam.buoyancyFactor              = mainConf.getParameter< real_t >( "buoyancyFactor" );
         initialisationParam.sphTool                     = std::make_shared< SphericalHarmonicsTool >( initialisationParam.deg );
         break;
      case INITIAL_TEMPERATURE_DEVIATION_METHOD::WHITE_NOISE:
         initialisationParam.buoyancyFactor = mainConf.getParameter< real_t >( "noiseFactor" );
         break;
      }
   }
   else
   {
      uint_t initialTemperatureDeviationMethodUint_t = mainConf.getParameter< uint_t >( "initialTemperatureDeviationMethod" );
      initialisationParam.initialTemperatureDeviationMethod =
          static_cast< INITIAL_TEMPERATURE_DEVIATION_METHOD >( initialTemperatureDeviationMethodUint_t );

      initialisationParam.buoyancyFactor = mainConf.getParameter< real_t >( "buoyancyFactor" );

      switch ( initialisationParam.initialTemperatureDeviationMethod )
      {
      case INITIAL_TEMPERATURE_DEVIATION_METHOD::RANDOM_SUPERPOSITION_SPH:
         initialisationParam.tempInit                    = mainConf.getParameter< uint_t >( "tempInit" );
         initialisationParam.initialTemperatureSteepness = mainConf.getParameter< real_t >( "initialTemperatureSteepness" );
         initialisationParam.lmax                        = mainConf.getParameter< uint_t >( "lmax" );
         initialisationParam.lmin                        = mainConf.getParameter< uint_t >( "lmin" );
         initialisationParam.superpositionSPHRandomSeed  = mainConf.getParameter< uint_t >( "superpositionSPHRandomSeed" );

         initialisationParam.sphTool = std::make_shared< SphericalHarmonicsTool >( initialisationParam.lmax );

         numHarmonics = ( ( initialisationParam.lmax + 1 ) * ( initialisationParam.lmax + 1 ) ) -
                        ( initialisationParam.lmin ) * ( initialisationParam.lmin );
         initialisationParam.superpositionRand.reserve( numHarmonics );
         walberla::math::seedRandomGenerator( initialisationParam.superpositionSPHRandomSeed );

         for ( uint_t i = 0; i < numHarmonics; i++ )
         {
            initialisationParam.superpositionRand.push_back( walberla::math::realRandom( real_c( -1 ), real_c( 1 ) ) );
         }
         break;
      case INITIAL_TEMPERATURE_DEVIATION_METHOD::SINGLE_SPH:
         initialisationParam.tempInit                    = mainConf.getParameter< uint_t >( "tempInit" );
         initialisationParam.initialTemperatureSteepness = mainConf.getParameter< real_t >( "initialTemperatureSteepness" );
         initialisationParam.deg                         = mainConf.getParameter< uint_t >( "degree" );
         initialisationParam.ord                         = mainConf.getParameter< int >( "order" );
         initialisationParam.sphTool                     = std::make_shared< SphericalHarmonicsTool >( initialisationParam.deg );
         break;
      default:
         WALBERLA_LOG_WARNING( "Handling for this INITIAL_TEMPERATURE_DEVIATION_METHOD is not implemented" )
      }
   }

   /*############ SOLVER PARAMETERS ############*/
   solverParam.rotFactor = mainConf.getParameter< real_t >( "rotFactor" );

   solverParam.solverFlag  = mainConf.getParameter< uint_t >( "solverFlag" );
   solverParam.solverPETSc = mainConf.getParameter< uint_t >( "PETScFlag" );

   solverParam.numPowerIterations                = mainConf.getParameter< uint_t >( "numPowerIterations" );
   solverParam.ChebyshevOrder                    = mainConf.getParameter< uint_t >( "ChebyshevOrder" );
   solverParam.ChebyshevSpectralRadiusUpperLimit = mainConf.getParameter< real_t >( "ChebyshevSpectralRadiusUpperLimit" );
   solverParam.ChebyshevSpectralRadiusLowerLimit = mainConf.getParameter< real_t >( "ChebyshevSpectralRadiusLowerLimit" );

   solverParam.FGMRESOuterIterations = mainConf.getParameter< uint_t >( "FGMRESOuterIterations" );
   solverParam.FGMRESRestartLength   = mainConf.getParameter< uint_t >( "FGMRESRestartLength" );
   solverParam.FGMRESTolerance       = mainConf.getParameter< real_t >( "FGMRESTolerance" );

   solverParam.uzawaIterations    = mainConf.getParameter< uint_t >( "uzawaIterations" );
   solverParam.uzawaOmega         = mainConf.getParameter< real_t >( "uzawaOmega" );
   solverParam.estimateUzawaOmega = mainConf.getParameter< bool >( "estimateUzawaOmega" );

   solverParam.ABlockMGIterations         = mainConf.getParameter< uint_t >( "ABlockMGIterations" );
   solverParam.ABlockMGTolerance          = mainConf.getParameter< real_t >( "ABlockMGTolerance" );
   solverParam.ABlockMGPreSmooth          = mainConf.getParameter< uint_t >( "ABlockMGPreSmooth" );
   solverParam.ABlockMGPostSmooth         = mainConf.getParameter< uint_t >( "ABlockMGPostSmooth" );
   solverParam.ABlockCoarseGridIterations = mainConf.getParameter< uint_t >( "ABlockCoarseGridIterations" );
   solverParam.ABlockCoarseGridTolerance  = mainConf.getParameter< real_t >( "ABlockCoarseGridTolerance" );

   solverParam.SchurMGIterations         = mainConf.getParameter< uint_t >( "SchurMGIterations" );
   solverParam.SchurMGTolerance          = mainConf.getParameter< real_t >( "SchurMGTolerance" );
   solverParam.SchurMGPreSmooth          = mainConf.getParameter< uint_t >( "SchurMGPreSmooth" );
   solverParam.SchurMGPostSmooth         = mainConf.getParameter< uint_t >( "SchurMGPostSmooth" );
   solverParam.SchurCoarseGridIterations = mainConf.getParameter< uint_t >( "SchurCoarseGridIterations" );
   solverParam.SchurCoarseGridTolerance  = mainConf.getParameter< real_t >( "SchurCoarseGridTolerance" );

   solverParam.stokesKillTolerance = mainConf.getParameter< real_t >( "stokesKillTolerance" );
   solverParam.useRotationWrapper  = mainConf.getParameter< bool >( "useRotationWrapper" );

   solverParam.stokesMaxNumIterations           = mainConf.getParameter< uint_t >( "stokesMaxNumIterations" );
   solverParam.stokesRelativeResidualUTolerance = mainConf.getParameter< real_t >( "stokesRelativeResidualUTolerance" );
   solverParam.stokesAbsoluteResidualUTolerance = mainConf.getParameter< real_t >( "stokesAbsoluteResidualUTolerance" );

   solverParam.diffusionMaxNumIterations           = mainConf.getParameter< uint_t >( "diffusionMaxNumIterations" );
   solverParam.diffusionAbsoluteResidualUTolerance = mainConf.getParameter< real_t >( "diffusionAbsoluteResidualUTolerance" );

   solverParam.stokesUzawaCoarseGridIter       = mainConf.getParameter< uint_t >( "stokesUzawaCoarseGridIter" );
   solverParam.stokesUzawaCoarseGridTol        = mainConf.getParameter< real_t >( "stokesUzawaCoarseGridTol" );
   solverParam.stokesSmoothIncrementCoarseGrid = mainConf.getParameter< uint_t >( "stokesSmoothIncrementCoarseGrid" );

   outputParam.outputDirectory     = mainConf.getParameter< std::string >( "outputDirectory" );
   outputParam.modelBaseName       = mainConf.getParameter< std::string >( "modelBaseName" );
   outputParam.overrideModelFolder = mainConf.getParameter< bool >( "overrideModelFolder" );
   outputParam.dataOutput          = mainConf.getParameter< bool >( "dataOutput" );
   outputParam.vtk                 = mainConf.getParameter< bool >( "vtk" );

   outputParam.OutputVelocity         = mainConf.getParameter< bool >( "OutputVelocity" );
   outputParam.OutputTemperature      = mainConf.getParameter< bool >( "OutputTemperature" );
   outputParam.OutputInterval         = mainConf.getParameter< uint_t >( "OutputInterval" );
   outputParam.OutputProfilesInterval = mainConf.getParameter< uint_t >( "OutputProfilesInterval" );
   outputParam.outputVertexDoFs       = mainConf.getParameter< bool >( "OutputVertexDoFs" );

   outputParam.ADIOS2ParamKey     = mainConf.getParameter< std::string >( "ADIOS2ParamKey" );
   outputParam.ADIOS2Value        = mainConf.getParameter< std::string >( "ADIOS2Value" );
   outputParam.ADIOS2OutputConfig = mainConf.getParameter< std::string >( "ADIOS2OutputConfig" );

   outputParam.ADIOS2StoreCheckpointPath     = mainConf.getParameter< std::string >( "ADIOS2StoreCheckpointPath" );
   outputParam.ADIOS2StoreCheckpointFilename = mainConf.getParameter< std::string >( "ADIOS2StoreCheckpointFilename" );

   outputParam.ADIOS2StartCheckpointPath     = mainConf.getParameter< std::string >( "ADIOS2StartCheckpointPath" );
   outputParam.ADIOS2StartCheckpointFilename = mainConf.getParameter< std::string >( "ADIOS2StartCheckpointFilename" );

   outputParam.ADIOS2StartFromCheckpoint = mainConf.getParameter< bool >( "ADIOS2StartFromCheckpoint" );
   outputParam.ADIOS2StoreCheckpoint     = mainConf.getParameter< bool >( "ADIOS2StoreCheckpoint" );

   outputParam.ADIOS2StoreCheckpointFrequency = mainConf.getParameter< uint_t >( "ADIOS2StoreCheckpointFrequency" );

   outputParam.outputProfiles = mainConf.getParameter< bool >( "outputProfiles" );

   outputParam.outputMyr = mainConf.getParameter< bool >( "outputMyr" );

   if ( outputParam.outputMyr )
   {
      outputParam.outputIntervalMyr = mainConf.getParameter< uint_t >( "outputIntervalMyr" );
      outputParam.OutputInterval    = 1;
   }

   outputParam.fileNameSQLdb  = mainConf.getParameter< std::string >( "SQLdatabaseFileName" );
   outputParam.createTimingDB = mainConf.getParameter< bool >( "createTimingDB" );

   //Start from timesteps and model age != 0 when starting from checkpoint and set output counters accordingly
   if ( outputParam.ADIOS2StartFromCheckpoint )
   {
      simulationParam.timeStep0      = mainConf.getParameter< uint_t >( "checkpointStartTimestep" );
      simulationParam.modelRunTimeMa = mainConf.getParameter< real_t >( "checkpointStartTimeMa" );
      simulationParam.modelTime      = ( simulationParam.modelRunTimeMa * physicalParam.characteristicVelocity *
                                    simulationParam.plateVelocityScaling * simulationParam.secondsPerMyr ) /
                                  physicalParam.mantleThickness;
      if ( outputParam.outputMyr )
      {
         outputParam.checkpointCount =
             static_cast< uint_t >( simulationParam.modelRunTimeMa / outputParam.ADIOS2StoreCheckpointFrequency ) + 1;
         outputParam.dataOutputCount =
             static_cast< uint_t >( simulationParam.modelRunTimeMa / outputParam.outputIntervalMyr ) + 1;
      }
      simulationParam.timeStep = simulationParam.timeStep0;
      WALBERLA_LOG_INFO_ON_ROOT( "Starting from checkpoint at timestep " << simulationParam.timeStep0 << " and model runtime "
                                                                         << simulationParam.modelRunTimeMa << " Ma." );
      if ( outputParam.outputMyr )
      {
         WALBERLA_LOG_INFO_ON_ROOT(
             "Next output at " << outputParam.dataOutputCount * outputParam.outputIntervalMyr << " Ma and next checkpoint at "
                               << outputParam.checkpointCount * outputParam.ADIOS2StoreCheckpointFrequency << " ." );
      }
   }

   simulationParam.verbose = mainConf.getParameter< bool >( "verbose" );
   //number of radial layers at max level (x2 for P2 elements)
   simulationParam.numLayers =
       2 * ( domainParam.nRad - 1 ) * ( hyteg::levelinfo::num_microvertices_per_edge( domainParam.maxLevel ) - 1 );

   TerraNeoParameters terraNeoParameters;
   terraNeoParameters.domainParameters         = domainParam;
   terraNeoParameters.solverParameters         = solverParam;
   terraNeoParameters.outputParameters         = outputParam;
   terraNeoParameters.simulationParameters     = simulationParam;
   terraNeoParameters.physicalParameters       = physicalParam;
   terraNeoParameters.initialisationParameters = initialisationParam;

   return terraNeoParameters;
}

/**
 * @brief Prints the configuration parameters to the log file.
 *
 * This function prints the domain, model, physical, non-dimensional, initialisation, and simulation parameters
 * to the root process and to the log file.
 */

inline void printConfig( const TerraNeoParameters& terraNeoParameters, std::string paramsPath = "" )
{
   const auto outputParam         = terraNeoParameters.outputParameters;
   const auto domainParam         = terraNeoParameters.domainParameters;
   const auto physicalParam       = terraNeoParameters.physicalParameters;
   const auto simulationParam     = terraNeoParameters.simulationParameters;
   const auto initialisationParam = terraNeoParameters.initialisationParameters;
   const auto solverParam         = terraNeoParameters.solverParameters;

   if ( paramsPath != "" )
   {
      //logging for model info
      WALBERLA_ROOT_SECTION()
      {
         walberla::logging::Logging::instance()->includeLoggingToFile( paramsPath );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Domain Parameters    ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "Min radius   : " << domainParam.rMin );
   WALBERLA_LOG_INFO_ON_ROOT( "Max radius   : " << domainParam.rMax );
   WALBERLA_LOG_INFO_ON_ROOT( "nTan         : " << domainParam.nTan );
   WALBERLA_LOG_INFO_ON_ROOT( "nRad         : " << domainParam.nRad );
   WALBERLA_LOG_INFO_ON_ROOT( "minLevel     : " << domainParam.minLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "maxLevel     : " << domainParam.maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "Domain Volume: " << domainParam.domainVolume() );
   WALBERLA_LOG_INFO_ON_ROOT( "Radii of macro layers: " << vectorToString( domainParam.macroLayers ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Physical Parameters    ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "Surface temperature          : " << physicalParam.surfaceTemp );
   WALBERLA_LOG_INFO_ON_ROOT( "Temperature CMB              : " << physicalParam.cmbTemp );
   WALBERLA_LOG_INFO_ON_ROOT( "Thermal conductivity         : " << physicalParam.thermalConductivity );
   WALBERLA_LOG_INFO_ON_ROOT( "Grueneisen parameter         : " << physicalParam.grueneisenParameter );
   WALBERLA_LOG_INFO_ON_ROOT( "Internal heating rate        : " << physicalParam.internalHeatingRate );
   WALBERLA_LOG_INFO_ON_ROOT( "Thermal diffusivity          : " << physicalParam.thermalDiffusivity );
   WALBERLA_LOG_INFO_ON_ROOT( "Characteristic velocity      : " << physicalParam.characteristicVelocity );

   if ( simulationParam.haveTemperatureProfile )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Temperature profile name     : " << simulationParam.fileTemperatureInputProfile );
   }
   if ( simulationParam.haveViscosityProfile )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Viscosity profile name       : " << simulationParam.fileViscosityProfile );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Viscosity                    : " << physicalParam.viscosity );
   }
   if ( simulationParam.haveThermalExpProfile )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Thermal exp. profile name    : " << simulationParam.fileThermalExpProfile );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Thermal expansivity          : " << physicalParam.thermalExpansivity );
   }

   if ( simulationParam.haveSpecificHeatCapProfile )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Specific heat cap. profile   : " << simulationParam.fileSpecificHeatCap );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Specific heat capacity       : " << physicalParam.specificHeatCapacity );
   }

   if ( simulationParam.haveDensityProfile )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Density profile              : " << simulationParam.fileDensityProfile );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Reference density            : " << physicalParam.referenceDensity );
   }

   if ( simulationParam.haveViscosityProfile || simulationParam.tempDependentViscosity )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Viscosity lower bound        : " << physicalParam.viscosityLowerBound );
      WALBERLA_LOG_INFO_ON_ROOT( "Viscosity upper bound        : " << physicalParam.viscosityUpperBound );
      if ( simulationParam.tempDependentViscosity )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "T dependent viscosity type: " << simulationParam.tempDependentViscosityType );
         WALBERLA_LOG_INFO_ON_ROOT( "Activation energy         : " << physicalParam.activationEnergy );
         WALBERLA_LOG_INFO_ON_ROOT( "Depth viscosity factor    : " << physicalParam.depthViscosityFactor );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Non-dimensial Numbers    ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "Rayleigh Number   : " << physicalParam.rayleighNumber );
   WALBERLA_LOG_INFO_ON_ROOT( "Peclet Number     : " << physicalParam.pecletNumber );
   WALBERLA_LOG_INFO_ON_ROOT( "Dissipation Number: " << physicalParam.dissipationNumber );
   WALBERLA_LOG_INFO_ON_ROOT( "hNumber           : " << physicalParam.hNumber );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Init Parameters   ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   switch ( initialisationParam.initialTemperatureDeviationMethod )
   {
   case INITIAL_TEMPERATURE_DEVIATION_METHOD::SINGLE_SPH:
      WALBERLA_LOG_INFO_ON_ROOT( "Single SPH" );
      WALBERLA_LOG_INFO_ON_ROOT( "tempInit                     : " << initialisationParam.tempInit );
      WALBERLA_LOG_INFO_ON_ROOT( "degree                       : " << initialisationParam.deg );
      WALBERLA_LOG_INFO_ON_ROOT( "order                        : " << initialisationParam.ord );
      WALBERLA_LOG_INFO_ON_ROOT( "Initial temperature steepness: " << initialisationParam.initialTemperatureSteepness );
      break;
   case INITIAL_TEMPERATURE_DEVIATION_METHOD::RANDOM_SUPERPOSITION_SPH:
      WALBERLA_LOG_INFO_ON_ROOT( "Random superposition SPH" );
      WALBERLA_LOG_INFO_ON_ROOT( "tempInit                     : " << initialisationParam.tempInit );
      WALBERLA_LOG_INFO_ON_ROOT( "lmin                         : " << initialisationParam.lmin );
      WALBERLA_LOG_INFO_ON_ROOT( "lmax                         : " << initialisationParam.lmax );
      WALBERLA_LOG_INFO_ON_ROOT( "Random Seed                  : " << initialisationParam.superpositionSPHRandomSeed );
      WALBERLA_LOG_INFO_ON_ROOT( "Initial temperature steepness: " << initialisationParam.initialTemperatureSteepness );
      break;
   case INITIAL_TEMPERATURE_DEVIATION_METHOD::WHITE_NOISE:
   default:
      WALBERLA_LOG_INFO_ON_ROOT( "White Noise" );
      break;
   }
   WALBERLA_LOG_INFO_ON_ROOT( "Buoyancy factor              : " << initialisationParam.buoyancyFactor );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Simulation Parameters    ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "Simulation type         : " << simulationParam.simulationType );
   WALBERLA_LOG_INFO_ON_ROOT( "Boundary condition      : " << simulationParam.boundaryCond );
   WALBERLA_LOG_INFO_ON_ROOT( "Unknowns temperature    : " << simulationParam.unknownsTemperature );
   WALBERLA_LOG_INFO_ON_ROOT( "Unknowns Stokes         : " << simulationParam.unknownsStokes );
   WALBERLA_LOG_INFO_ON_ROOT( "hMin                    : " << simulationParam.hMin );
   WALBERLA_LOG_INFO_ON_ROOT( "hMax                    : " << simulationParam.hMax );
   WALBERLA_LOG_INFO_ON_ROOT( "Fixed timestep          : " << ( simulationParam.fixedTimestep ? "true" : "false" ) );
   if ( simulationParam.fixedTimestep )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "dtConstant              : " << simulationParam.dtConstant );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "cflMax                  : " << simulationParam.cflMax );
   WALBERLA_LOG_INFO_ON_ROOT( "MaxNumTimesteps         : " << simulationParam.maxNumTimesteps );
   WALBERLA_LOG_INFO_ON_ROOT( "Compressible            : " << ( simulationParam.compressible ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Shear heating           : " << ( simulationParam.shearHeating ? "true" : "false" ) );
   if ( simulationParam.shearHeating )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Shear heating scaling factor : " << simulationParam.lithosphereShearHeatingScaling );
      WALBERLA_LOG_INFO_ON_ROOT( "Lithosphere thickness [km]   : " << simulationParam.lithosphereThickness );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Adiabatic heating       : " << ( simulationParam.adiabaticHeating ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Internal heating        : " << ( simulationParam.internalHeating ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Frozen velocity         : " << ( simulationParam.frozenVelocity ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "T-dependent viscosity   : " << ( simulationParam.tempDependentViscosity ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "adaptive Ref Temp.      : " << ( simulationParam.adaptiveRefTemp ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT(
       "volumetric avrg Ref Temp.      : " << ( simulationParam.volAvrgTemperatureDev ? "true" : "false" ) );

   if ( simulationParam.simulationType == "CirculationModel" )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Filename topologies     : " << simulationParam.fnameTopologies );
      WALBERLA_LOG_INFO_ON_ROOT( "Filename reconstructions: " << simulationParam.fnameReconstructions );
      WALBERLA_LOG_INFO_ON_ROOT( "Initital age            : " << simulationParam.initialAge );
      WALBERLA_LOG_INFO_ON_ROOT( "Final age               : " << simulationParam.finalAge );
      WALBERLA_LOG_INFO_ON_ROOT( "Plate velocity scaling  : " << simulationParam.plateVelocityScaling );
      WALBERLA_LOG_INFO_ON_ROOT( "Plate smoothing distance: " << simulationParam.plateSmoothingDistance );
   }
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Output Parameters    ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "Output directory  : " << outputParam.outputDirectory );
   WALBERLA_LOG_INFO_ON_ROOT( "Output base name  : " << outputParam.outputBaseName );
   WALBERLA_LOG_INFO_ON_ROOT( "data Output       : " << ( outputParam.dataOutput ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "vtk               : " << ( outputParam.vtk ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Output velocity   : " << ( outputParam.OutputVelocity ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Output interval   : " << outputParam.OutputInterval );
   WALBERLA_LOG_INFO_ON_ROOT( "Output Vertex DoFs: " << ( outputParam.outputVertexDoFs ? "true" : "false" ) );
   if ( outputParam.outputProfiles && simulationParam.tempDependentViscosity )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Output temperature & viscosity profiles: "
                                 << "true" );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Output temperature profiles: " << ( outputParam.outputProfiles ? "true" : "false" ) );
   }

   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Solver Parameters    ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   if ( solverParam.solverPETSc == 1u )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Use PETSc solver for coarse grid       : " << "true" );
   }
   if ( solverParam.solverFlag == 0u )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "FGMRES solver outer iterations         : " << solverParam.FGMRESOuterIterations );
      WALBERLA_LOG_INFO_ON_ROOT( "FGMRES solver tolerance                : " << solverParam.FGMRESTolerance );
      WALBERLA_LOG_INFO_ON_ROOT( "Uzawa smoother iterations              : " << solverParam.uzawaIterations );

      WALBERLA_LOG_INFO_ON_ROOT( "A-Block multigrid iterations           : " << solverParam.ABlockMGIterations );
      WALBERLA_LOG_INFO_ON_ROOT( "A-Block multigrid solver tolerance     : " << solverParam.ABlockMGTolerance );
      WALBERLA_LOG_INFO_ON_ROOT( "A-Block multigrid pre-smoothing steps  : " << solverParam.ABlockMGPreSmooth );
      WALBERLA_LOG_INFO_ON_ROOT( "A-Block multigrid post-smoothing steps : " << solverParam.ABlockMGPostSmooth );
      WALBERLA_LOG_INFO_ON_ROOT( "Schur multigrid iterations             : " << solverParam.SchurMGIterations );
      WALBERLA_LOG_INFO_ON_ROOT( "Schur multigrid solver tolerance       : " << solverParam.SchurMGTolerance );
      WALBERLA_LOG_INFO_ON_ROOT( "Schur multigrid pre-smoothing steps    : " << solverParam.SchurMGPreSmooth );
      WALBERLA_LOG_INFO_ON_ROOT( "Schur multigrid post-smoothing steps   : " << solverParam.SchurMGPostSmooth );
      WALBERLA_LOG_INFO_ON_ROOT( "Diffusion max num iterations           : " << solverParam.diffusionMaxNumIterations );
      WALBERLA_LOG_INFO_ON_ROOT( "Diffusion absolute residual U tolerance: " << solverParam.diffusionAbsoluteResidualUTolerance );
      WALBERLA_LOG_INFO_ON_ROOT( "Stokes kill-tolerance                  : " << solverParam.stokesKillTolerance );
      WALBERLA_LOG_INFO_ON_ROOT( " " );
      WALBERLA_LOG_INFO_ON_ROOT( "A-Block coarse grid iterations         : " << solverParam.ABlockCoarseGridIterations );
      WALBERLA_LOG_INFO_ON_ROOT( "A-Block coarse grid tolerance          : " << solverParam.ABlockCoarseGridTolerance );
      WALBERLA_LOG_INFO_ON_ROOT( "Schur coarse grid iterations           : " << solverParam.SchurCoarseGridIterations );
      WALBERLA_LOG_INFO_ON_ROOT( "Schur coarse grid tolerance            : " << solverParam.SchurCoarseGridTolerance );
   }
   else if ( solverParam.solverFlag == 1u )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Max Stokes coarse grid iterations      : " << solverParam.stokesMaxNumIterations );
      WALBERLA_LOG_INFO_ON_ROOT( "Stokes relative residual U tolerance   : " << solverParam.stokesRelativeResidualUTolerance );
      WALBERLA_LOG_INFO_ON_ROOT( "Stokes relative absolute U tolerance   : " << solverParam.stokesAbsoluteResidualUTolerance );
   }

   WALBERLA_ROOT_SECTION()
   {
      walberla::logging::Logging::instance()->stopLoggingToFile();
   }
}

} // namespace terraneo
