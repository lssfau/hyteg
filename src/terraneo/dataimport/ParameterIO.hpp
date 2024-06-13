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

namespace terraneo {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

/**
 * @brief Parses the configuration parameters from the main configuration block.
 *
 * This function reads and extracts various domain, model, simulation, and initialization parameters from the main configuration block.
 * It populates the corresponding variables and performs necessary calculations for non-dimensional numbers.
 */

// define plate velocity oracle
std::shared_ptr< terraneo::plates::PlateVelocityProvider > oracle;

inline TerraNeoParameters parseConfig( const walberla::Config::BlockHandle& mainConf )
{
   DomainParameters         domainParam;
   SolverParameters         solverParam;
   OutputParameters         outputParam;
   SimulationParameters     simulationParam;
   PhysicalParameters       physicalParam;
   InitialisationParameters initialisationParam;

   /*############ DOMAIN PARAMETERS ############*/

   domainParam.rCMB     = mainConf.getParameter< real_t >( "rCMB" );
   domainParam.rSurface = mainConf.getParameter< real_t >( "rSurface" );
   domainParam.nTan     = mainConf.getParameter< uint_t >( "nTan" );
   domainParam.nRad     = mainConf.getParameter< uint_t >( "nRad" );
   domainParam.minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   domainParam.maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   //calculate non-dimensional radii such that mantle thickness = 1
   domainParam.rMin = domainParam.rCMB / ( domainParam.rSurface - domainParam.rCMB );
   domainParam.rMax = domainParam.rSurface / ( domainParam.rSurface - domainParam.rCMB );

   /*############ MODEL PARAMETERS ############*/

   if ( mainConf.isDefined( "viscosityProfile" ) )
   {
      simulationParam.fileViscosityProfile = mainConf.getParameter< std::string >( "viscosityProfile" );
      auto viscosityJson                   = io::readJsonFile( simulationParam.fileViscosityProfile );

      const auto radiusKey    = "Radius (m)";
      const auto viscosityKey = "Viscosity (Pa s)";

      WALBERLA_CHECK_GREATER( viscosityJson.count( radiusKey ), 0, "No key '" << radiusKey << "' in viscosity profile file." )
      WALBERLA_CHECK_GREATER(
          viscosityJson.count( viscosityKey ), 0, "No key '" << viscosityKey << "' in viscosity profile file." )

      physicalParam.radius           = viscosityJson[radiusKey].get< std::vector< real_t > >();
      physicalParam.viscosityProfile = viscosityJson[viscosityKey].get< std::vector< real_t > >();

      WALBERLA_CHECK_EQUAL( physicalParam.radius.size(), physicalParam.viscosityProfile.size() )

      simulationParam.haveViscosityProfile = true;

      for ( uint_t i = 0; i < physicalParam.radius.size(); i++ )
      {
         // non-dimensionalise radius
         physicalParam.radius[i] /= physicalParam.mantleThickness;
      }

      // ¿Why is there a commented std::reverse?
      // std::reverse( physicalParam.viscosityProfile.begin(), physicalParam.viscosityProfile.end() );
   }

   else
   {
      physicalParam.viscosity = mainConf.getParameter< real_t >( "viscosity" );
   }

   physicalParam.cmbTemp                     = mainConf.getParameter< real_t >( "cmbTemp" );
   physicalParam.surfaceTemp                 = mainConf.getParameter< real_t >( "surfaceTemp" );
   physicalParam.initialTemperatureSteepness = mainConf.getParameter< real_t >( "initialTemperatureSteepness" );
   physicalParam.thermalExpansivity          = mainConf.getParameter< real_t >( "thermalExpansivity" );
   physicalParam.thermalConductivity         = mainConf.getParameter< real_t >( "thermalConductivity" );
   physicalParam.grueneisenParameter         = mainConf.getParameter< real_t >( "grueneisenParameter" );
   physicalParam.specificHeatCapacity        = mainConf.getParameter< real_t >( "specificHeatCapacity" );
   physicalParam.internalHeatingRate         = mainConf.getParameter< real_t >( "internalHeatingRate" );
   physicalParam.characteristicVelocity      = mainConf.getParameter< real_t >( "characteristicVelocity" );
   physicalParam.surfaceDensity              = mainConf.getParameter< real_t >( "surfaceDensity" );
   physicalParam.referenceDensity            = mainConf.getParameter< real_t >( "referenceDensity" );

   //used to calculate non-D numbers
   physicalParam.mantleThickness = domainParam.rSurface - domainParam.rCMB;
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
   simulationParam.dtConstant             = mainConf.getParameter< real_t >( "dtConstant" );
   simulationParam.maxNumTimesteps        = mainConf.getParameter< uint_t >( "maxNumTimesteps" );
   simulationParam.adaptiveRefTemp        = mainConf.getParameter< bool >( "adaptiveRefTemp" );
   simulationParam.tempDependentViscosity = mainConf.getParameter< bool >( "tempDependentViscosity" );
   simulationParam.simulationType         = mainConf.getParameter< std::string >( "simulationType" );
   simulationParam.compressible           = mainConf.getParameter< bool >( "compressible" );
   simulationParam.shearHeating           = mainConf.getParameter< bool >( "shearHeating" );
   simulationParam.adiabaticHeating       = mainConf.getParameter< bool >( "adiabaticHeating" );
   simulationParam.internalHeating        = mainConf.getParameter< bool >( "internalHeating" );
   simulationParam.boundaryCond           = mainConf.getParameter< uint_t >( "boundaryCond" );
   simulationParam.timingAnalysis         = mainConf.getParameter< bool >( "timingAnalysis" );

   if ( simulationParam.tempDependentViscosity )
   {
      simulationParam.tempDependentViscosityType = mainConf.getParameter< uint_t >( "tempDependentViscosityType" );
      simulationParam.resetSolver                = mainConf.getParameter< bool >( "resetSolver" );
      simulationParam.resetSolverFrequency       = mainConf.getParameter< uint_t >( "resetSolverFrequency" );
      physicalParam.activationEnergy             = mainConf.getParameter< real_t >( "activationEnergy" );
      physicalParam.depthViscosityFactor         = mainConf.getParameter< real_t >( "depthViscosityFactor" );
      physicalParam.viscosityLowerBound          = mainConf.getParameter< real_t >( "viscosityLowerBound" );
      physicalParam.viscosityUpperBound          = mainConf.getParameter< real_t >( "viscosityUpperBound" );
   }

   //simulation parameters for circulation models only:
   if ( simulationParam.simulationType == "CirculationModel" )
   {
      if ( mainConf.isDefined( "fnameTopologies" ) && mainConf.isDefined( "fnameReconstructions" ) )
      {
         simulationParam.fnameTopologies        = mainConf.getParameter< std::string >( "fnameTopologies" );
         simulationParam.fnameReconstructions   = mainConf.getParameter< std::string >( "fnameReconstructions" );
         simulationParam.initialAge             = mainConf.getParameter< real_t >( "initialAge" );
         simulationParam.finalAge               = mainConf.getParameter< real_t >( "finalAge" );
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

   initialisationParam.tempInit       = mainConf.getParameter< uint_t >( "tempInit" );
   initialisationParam.deg            = mainConf.getParameter< uint_t >( "degree" );
   initialisationParam.ord            = mainConf.getParameter< int >( "order" );
   initialisationParam.lmax           = mainConf.getParameter< uint_t >( "lmax" );
   initialisationParam.lmin           = mainConf.getParameter< uint_t >( "lmin" );
   initialisationParam.superposition  = mainConf.getParameter< bool >( "superposition" );
   initialisationParam.buoyancyFactor = mainConf.getParameter< real_t >( "buoyancyFactor" );

   initialisationParam.temperatureNoise             = mainConf.getParameter< bool >( "temperatureNoise" );
   initialisationParam.temperatureSphericalHarmonic = mainConf.getParameter< bool >( "temperatureSphericalHarmonic" );
   initialisationParam.noiseFactor                  = mainConf.getParameter< real_t >( "noiseFactor" );

   /*############ SOLVER PARAMETERS ############*/

   solverParam.numPowerIterations    = mainConf.getParameter< uint_t >( "numPowerIterations" );
   solverParam.FGMRESOuterIterations = mainConf.getParameter< uint_t >( "FGMRESOuterIterations" );
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
   solverParam.SchurMGTolerance          = mainConf.getParameter< uint_t >( "SchurMGTolerance" );
   solverParam.SchurMGPreSmooth          = mainConf.getParameter< uint_t >( "SchurMGPreSmooth" );
   solverParam.SchurMGPostSmooth         = mainConf.getParameter< uint_t >( "SchurMGPostSmooth" );
   solverParam.SchurCoarseGridIterations = mainConf.getParameter< uint_t >( "SchurCoarseGridIterations" );
   solverParam.SchurCoarseGridTolerance  = mainConf.getParameter< real_t >( "SchurCoarseGridTolerance" );

   solverParam.stokesKillTolerance = mainConf.getParameter< uint_t >( "stokesKillTolerance" );

   solverParam.diffusionMaxNumIterations           = mainConf.getParameter< uint_t >( "diffusionMaxNumIterations" );
   solverParam.diffusionAbsoluteResidualUTolerance = mainConf.getParameter< real_t >( "diffusionAbsoluteResidualUTolerance" );

   outputParam.outputDirectory = mainConf.getParameter< std::string >( "outputDirectory" );
   outputParam.outputBaseName  = mainConf.getParameter< std::string >( "outputBaseName" );
   outputParam.dataOutput      = mainConf.getParameter< bool >( "dataOutput" );
   outputParam.vtk             = mainConf.getParameter< bool >( "vtk" );

   outputParam.OutputVelocity      = mainConf.getParameter< bool >( "OutputVelocity" );
   outputParam.OutputTemperature   = mainConf.getParameter< bool >( "OutputTemperature" );
   outputParam.OutputInterval      = mainConf.getParameter< uint_t >( "OutputInterval" );
   outputParam.vtkOutputVertexDoFs = mainConf.getParameter< bool >( "OutputVertexDoFs" );

   outputParam.ADIOS2ParamKey           = mainConf.getParameter< std::string >( "ADIOS2ParamKey" );
   outputParam.ADIOS2Value              = mainConf.getParameter< std::string >( "ADIOS2Value" );
   outputParam.ADIOS2OutputConfig       = mainConf.getParameter< std::string >( "ADIOS2OutputConfig" );
   outputParam.ADIOS2CheckpointPath     = mainConf.getParameter< std::string >( "ADIOS2CheckpointPath" );
   outputParam.ADIOS2CheckpointFilename = mainConf.getParameter< std::string >( "ADIOS2CheckpointFilename" );

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

inline void printConfig( const TerraNeoParameters& terraNeoParameters )
{
   const auto outputParam         = terraNeoParameters.outputParameters;
   const auto domainParam         = terraNeoParameters.domainParameters;
   const auto physicalParam       = terraNeoParameters.physicalParameters;
   const auto simulationParam     = terraNeoParameters.simulationParameters;
   const auto initialisationParam = terraNeoParameters.initialisationParameters;
   const auto solverParam         = terraNeoParameters.solverParameters;

   //logging for model info
   WALBERLA_ROOT_SECTION()
   {
      walberla::logging::Logging::instance()->includeLoggingToFile( outputParam.outputDirectory + "/" +
                                                                    outputParam.outputBaseName + "_params.out" );
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
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Physical Parameters    ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "Surface Temperature          : " << physicalParam.surfaceTemp );
   WALBERLA_LOG_INFO_ON_ROOT( "Temperature CMB              : " << physicalParam.cmbTemp );
   WALBERLA_LOG_INFO_ON_ROOT( "Initial Temperature steepness: " << physicalParam.initialTemperatureSteepness );
   WALBERLA_LOG_INFO_ON_ROOT( "Thermal Expansivity          : " << physicalParam.thermalExpansivity );
   WALBERLA_LOG_INFO_ON_ROOT( "Thermal Conductivity         : " << physicalParam.thermalConductivity );
   WALBERLA_LOG_INFO_ON_ROOT( "Specific Heat Capacity       : " << physicalParam.specificHeatCapacity );
   WALBERLA_LOG_INFO_ON_ROOT( "Internal Heating Rate        : " << physicalParam.internalHeatingRate );
   WALBERLA_LOG_INFO_ON_ROOT( "Thermal Diffusivity          : " << physicalParam.thermalDiffusivity );
   WALBERLA_LOG_INFO_ON_ROOT( "Characteristic Velocity      : " << physicalParam.characteristicVelocity );

   WALBERLA_LOG_INFO_ON_ROOT( "Reference density            : " << physicalParam.referenceDensity );

   if ( simulationParam.haveViscosityProfile )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Viscosity profile name       : " << simulationParam.fileViscosityProfile );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "viscosity                    : " << physicalParam.viscosity );
   }

   if ( simulationParam.tempDependentViscosity )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "T dependent Viscosity type: " << simulationParam.tempDependentViscosityType );
      WALBERLA_LOG_INFO_ON_ROOT( "Activation Energy         : " << physicalParam.activationEnergy );
      WALBERLA_LOG_INFO_ON_ROOT( "Depth Viscosity Factor    : " << physicalParam.depthViscosityFactor );
      WALBERLA_LOG_INFO_ON_ROOT( "Viscosity lower bound     : " << physicalParam.viscosityLowerBound );
      WALBERLA_LOG_INFO_ON_ROOT( "Viscosity upper bound     : " << physicalParam.viscosityUpperBound );
      WALBERLA_LOG_INFO_ON_ROOT( "reset Solver              : " << ( simulationParam.resetSolver ? "true" : "false" ) );
      WALBERLA_LOG_INFO_ON_ROOT( "reset Solver Frequency    : " << ( simulationParam.resetSolverFrequency ) );
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
   WALBERLA_LOG_INFO_ON_ROOT( "tempInit       : " << initialisationParam.tempInit );
   WALBERLA_LOG_INFO_ON_ROOT( "degree         : " << initialisationParam.deg );
   WALBERLA_LOG_INFO_ON_ROOT( "order          : " << initialisationParam.ord );
   WALBERLA_LOG_INFO_ON_ROOT( "lmin           : " << initialisationParam.lmin );
   WALBERLA_LOG_INFO_ON_ROOT( "lmax           : " << initialisationParam.lmax );
   WALBERLA_LOG_INFO_ON_ROOT( "Superposition  : " << initialisationParam.superposition );
   WALBERLA_LOG_INFO_ON_ROOT( "Buoyancy factor: " << initialisationParam.buoyancyFactor );
   WALBERLA_LOG_INFO_ON_ROOT( "Noise Factor   : " << initialisationParam.noiseFactor );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Simulation Parameters    ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "Simulation type         : " << simulationParam.simulationType );
   WALBERLA_LOG_INFO_ON_ROOT( "Unknowns Temperature    : " << simulationParam.unknownsTemperature );
   WALBERLA_LOG_INFO_ON_ROOT( "Unknowns Stokes         : " << simulationParam.unknownsStokes );
   WALBERLA_LOG_INFO_ON_ROOT( "hMin                    : " << simulationParam.hMin );
   WALBERLA_LOG_INFO_ON_ROOT( "hMax                    : " << simulationParam.hMax );
   WALBERLA_LOG_INFO_ON_ROOT( "Fixed timestep          : " << ( simulationParam.fixedTimestep ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "dtConstant              : " << simulationParam.dtConstant );
   WALBERLA_LOG_INFO_ON_ROOT( "Timestep                : " << simulationParam.timestep );
   WALBERLA_LOG_INFO_ON_ROOT( "MaxNumTimesteps         : " << simulationParam.maxNumTimesteps );
   WALBERLA_LOG_INFO_ON_ROOT( "Compressible            : " << ( simulationParam.compressible ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Shear heating           : " << ( simulationParam.shearHeating ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Shear heating           : " << ( simulationParam.adiabaticHeating ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Shear heating           : " << ( simulationParam.internalHeating ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "T-dependent Viscosity   : " << ( simulationParam.tempDependentViscosity ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "adaptive Ref Temp.      : " << ( simulationParam.adaptiveRefTemp ? "true" : "false" ) );

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
   WALBERLA_LOG_INFO_ON_ROOT( "Output Base Name  : " << outputParam.outputBaseName );
   WALBERLA_LOG_INFO_ON_ROOT( "data Output       : " << ( outputParam.dataOutput ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "vtk               : " << ( outputParam.vtk ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Output Velocity   : " << ( outputParam.OutputVelocity ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Output Interval   : " << outputParam.OutputInterval );
   WALBERLA_LOG_INFO_ON_ROOT( "Output Vertex DoFs: " << ( outputParam.vtkOutputVertexDoFs ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Output Vertex DoFs: " << ( outputParam.outputProfiles ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----    Solver Parameters    ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "FGMRES solver outer iterations         : " << solverParam.FGMRESOuterIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "FGMRES solver tolerance                : " << solverParam.FGMRESTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "Uzawa smoother iterations              : " << solverParam.uzawaIterations );

   WALBERLA_LOG_INFO_ON_ROOT( "A-Block multigrid solver tolerance     : " << solverParam.ABlockMGTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "A-Block multigrid iterations           : " << solverParam.ABlockMGIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "A-Block multigrid pre-smoothing steps  : " << solverParam.ABlockMGPreSmooth );
   WALBERLA_LOG_INFO_ON_ROOT( "A-Block multigrid post-smoothing steps : " << solverParam.ABlockMGPostSmooth );

   WALBERLA_LOG_INFO_ON_ROOT( "Schur multigrid solver tolerance       : " << solverParam.SchurMGTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "Schur multigrid iterations             : " << solverParam.SchurMGIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "Schur multigrid pre-smoothing steps    : " << solverParam.SchurMGPreSmooth );
   WALBERLA_LOG_INFO_ON_ROOT( "Schur multigrid post-smoothing steps   : " << solverParam.SchurMGPostSmooth );

   WALBERLA_LOG_INFO_ON_ROOT( "Diffusion max num iterations           : " << solverParam.diffusionMaxNumIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "Diffusion absolute residual U-tolerance: " << solverParam.diffusionAbsoluteResidualUTolerance );

   WALBERLA_LOG_INFO_ON_ROOT( "Stokes kill-tolerance                  : " << solverParam.stokesKillTolerance );

   WALBERLA_ROOT_SECTION() { walberla::logging::Logging::instance()->stopLoggingToFile(); }
}

} // namespace terraneo