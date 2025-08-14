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

#include <chrono>
#include <cmath>
#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/math/Random.h>
#include <filesystem>
#include <sstream>
#include <thread>
#include <vector>

#include "core/DataTypes.h"

#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

namespace terraneo {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

// Domain information

struct DomainParameters
{
   //geometric information

   real_t rCMB     = real_c( 3480000 );
   real_t rSurface = real_c( 6371000 );

   //calculate non-dimensional radii such that mantle thickness = 1
   real_t rMin = rCMB / ( rSurface - rCMB );
   real_t rMax = rSurface / ( rSurface - rCMB );

   uint_t nTan          = 2;
   uint_t nRad          = 2;
   uint_t minLevel      = 0;
   uint_t maxLevel      = 1;
   uint_t numProcessors = 1;

   real_t domainVolume() const
   {
      return ( real_c( 4.0 ) / real_c( 3.0 ) ) * walberla::math::pi * rMax * rMax * rMax -
             ( real_c( 4.0 ) / real_c( 3.0 ) ) * walberla::math::pi * rMin * rMin * rMin;
   }
};

// Solver parameters

struct SolverParameters
{
   uint_t solverFlag  = 0u;
   uint_t solverPETSc = 0u;

   // Stokes solver parameters

   uint_t numPowerIterations                = 25;
   uint_t ChebyshevOrder                    = 3u;
   real_t ChebyshevSpectralRadiusUpperLimit = 1.2;
   real_t ChebyshevSpectralRadiusLowerLimit = 0.1;

   uint_t FGMRESOuterIterations      = 5;
   uint_t FGMRESRestartLength        = 5;
   real_t FGMRESTolerance            = 1e-6;
   uint_t uzawaIterations            = 5;
   real_t uzawaOmega                 = real_c( 0.3 );
   bool   estimateUzawaOmega         = false;
   uint_t ABlockMGIterations         = 5;
   real_t ABlockMGTolerance          = 1e-6;
   uint_t ABlockMGPreSmooth          = 3;
   uint_t ABlockMGPostSmooth         = 3;
   uint_t ABlockCoarseGridIterations = 5;
   real_t ABlockCoarseGridTolerance  = 1e-6;
   uint_t SchurMGIterations          = 5;
   real_t SchurMGTolerance           = 1e-6;
   uint_t SchurMGPreSmooth           = 3;
   uint_t SchurMGPostSmooth          = 3;
   uint_t SchurCoarseGridIterations  = 5;
   real_t SchurCoarseGridTolerance   = 1e-6;
   real_t stokesKillTolerance        = real_c( 1000 );

   bool useRotationWrapper = true;

   // Uzawa type multigrid solver
   real_t initialResidualU                 = real_c( 0 );
   real_t vCycleResidualUPrev              = real_c( 0 );
   uint_t numVCycles                       = 0;
   real_t averageResidualReductionU        = real_c( 0 );
   uint_t stokesMaxNumIterations           = 5;
   real_t stokesRelativeResidualUTolerance = 1e-6;
   real_t stokesAbsoluteResidualUTolerance = 1e-6;
   uint_t stokesUzawaCoarseGridIter        = 10;
   real_t stokesUzawaCoarseGridTol         = 1e-6;
   uint_t stokesSmoothIncrementCoarseGrid  = 2;

   // Diffusion solver parameters

   uint_t diffusionMaxNumIterations           = 10000;
   real_t diffusionAbsoluteResidualUTolerance = real_c( 10000 );

   real_t gmresApproximationToleranceTransport = real_c( 1e-5 );

   real_t rotFactor = 0.0;
};
// Output parameters
struct OutputParameters
{
   std::string outputDirectory = std::string( "output" );
   std::string modelBaseName   = std::string( "conv_sim" );

   std::string ADIOS2OutputConfig = std::string( "ADIOS2config.xml" );
   std::string ADIOS2ParamKey     = std::string( "NumAggregators" );
   std::string ADIOS2Value        = std::string( "16" );

   std::string ADIOS2StoreCheckpointPath     = std::string( "output" );
   std::string ADIOS2StoreCheckpointFilename = std::string( "conv_sim_store" );

   std::string ADIOS2StartCheckpointPath     = std::string( "output" );
   std::string ADIOS2StartCheckpointFilename = std::string( "conv_sim_store" );

   std::string ADIOS2CheckpointPath     = std::string( "output" );
   std::string ADIOS2CheckpointFilename = std::string( "conv_sim" );

   std::string fileNameSQLdb = std::string( "SQLdatabaseFileName" );

   bool ADIOS2StartFromCheckpoint = false;
   bool ADIOS2StoreCheckpoint     = false;

   uint_t ADIOS2StoreCheckpointFrequency = 100U;

   bool dataOutput = true;
   bool vtk        = true;

   uint_t OutputInterval         = 1;
   uint_t OutputProfilesInterval = 1;
   uint_t checkpointCount        = 1;

   bool OutputVelocity    = true;
   bool OutputViscosity   = true;
   bool OutputDensity     = true;
   bool OutputTemperature = true;
   bool OutputVerbose     = false;
   bool OutputDimensional = false;

   bool   outputMyr         = false;
   uint_t outputIntervalMyr = 1;
   real_t prevOutputTime    = real_c( 0 );

   bool outputVertexDoFs = true;
   bool outputProfiles   = false;

   bool createTimingDB = false;
};
// Simulation parameters
struct SimulationParameters
{
   std::string getKeyValue( nlohmann::json& jsonData )
   {
      const auto  keyRadius = "Radius (m)";
      std::string keyValue  = "Value";

      // Iterate over the keys in the JSON object
      for ( auto& [key, value] : jsonData.items() )
      {
         if ( key != keyRadius )
         {
            keyValue = key;
            break;
         }
      }

      // Check that both radius and value data are available within the provided profile
      WALBERLA_CHECK_GREATER( jsonData.count( keyRadius ), 0, "No key '" << keyRadius << "' in profile file." );
      WALBERLA_CHECK_GREATER( jsonData.count( keyValue ), 0, "No key '" << keyValue << "' in profile file." );

      return keyValue;
   }

   //Parameters derived from other parameters
   uint_t unknownsTemperature = 0;
   uint_t unknownsStokes      = 0;
   real_t hMin                = real_c( 0 );
   real_t hMax                = real_c( 0 );
   uint_t numLayers           = 0;

   //Parameters given via config file
   bool        fixedTimestep              = false;
   uint_t      timeStep                   = 0;
   real_t      maxTimestepSize            = 0;
   real_t      avrgTemperatureVol         = real_c( 0 );
   real_t      modelTime                  = real_c( 0 );
   real_t      dtPrev                     = real_c( 0 );
   real_t      dt                         = real_c( 0 );
   real_t      dtConstant                 = real_c( 0 );
   real_t      cflMax                     = real_c( 1 );
   uint_t      timestep                   = 0;
   std::string simulationType             = std::string( "ConvectionModel" );
   uint_t      maxNumTimesteps            = 100;
   bool        adaptiveRefTemp            = false;
   bool        tempDependentViscosity     = false;
   bool        volAvrgTemperatureDev      = false;
   uint_t      tempDependentViscosityType = 0;

   //circulation model parameters
   real_t      initialAge     = real_c( 100 ); //initial age for circulation models
   real_t      finalAge       = real_c( 0 );   //final age for circulation models
   real_t      ageMa          = initialAge;    //currrent age during circulation model
   real_t      agePrev        = initialAge;    //age of previous timestep during circulation model
   real_t      plateAge       = initialAge;    //current age of plates being implemented (intervals of 1Myr)
   real_t      modelRunTimeMa = real_c( 0 );
   std::string fnameTopologies;
   std::string fnameReconstructions;
   real_t      plateVelocityScaling        = real_c( 1 );
   real_t      plateSmoothingDistance      = real_c( 110 );
   bool        compressible                = true; // default: Compressible foŕmulation
   bool        frozenVelocity              = false; // default: Non-frozen
   bool        shearHeating                = true; //default: include shear heating
   bool        adiabaticHeating            = true; //default: include adiabatic heating
   bool        internalHeating             = true; //default: include internal heating
   uint_t      boundaryCond                = 1;    // default: No-Slip/No-Slip
   bool        boundaryCondFreeSlip        = false;
   bool        verbose                     = false;
   bool        haveTemperatureProfile      = false;
   bool        haveViscosityProfile        = false;
   bool        haveThermalExpProfile       = false;
   bool        haveSpecificHeatCapProfile  = false;
   bool        haveDensityProfile          = false;
   bool        predictorCorrector          = false;
   std::string fileTemperatureInputProfile = std::string( "TemperatureInputProfile.json" );
   std::string fileViscosityProfile        = std::string( "ViscosityProfile.json" );
   std::string fileThermalExpProfile       = std::string( "ThermalExpProfile.json" );
   std::string fileSpecificHeatCap         = std::string( "SpecificHeatCapProfile.json" );
   std::string fileDensityProfile          = std::string( "DensityProfile.json" );

   //needed for conversions in the simulation
   real_t secondsPerMyr = real_c( 3.154e7 * 1e6 );

   // Shear heating scaling for mantle ciruclation model with
   // predifned Lithosphere thickness in km
   real_t lithosphereShearHeatingScaling = 1e-5;
   real_t lithosphereThickness           = real_c( 100 );

   // Needed for timing analysis of the simulation run
   bool timingAnalysis = true;

   // Failsafe parameters
   bool checkTemperatureConsistency = false;     // Check if radially averaged temperatures goes
                                                 // out of the expected range [surfaceTemp, cmbTemp]
   real_t temperatureConsistencyThreshold = 1.0; // Kelvin

   real_t initialTimestepSize                    = 1.0; // Ma
   uint_t initialNStepsForTimestepLinearIncrease = 20u;

   // MMOC Parameters
   real_t particleLocationRadius    = real_c( 1e-2 );
   bool   projectPointsBackToDomain = false;
   bool   cautionedEvaluate         = false;
};

enum class INITIAL_TEMPERATURE_DEVIATION_METHOD : uint_t
{
   WHITE_NOISE              = 0,
   SINGLE_SPH               = 1,
   RANDOM_SUPERPOSITION_SPH = 2,
};
struct TemperatureDeivationInitialisationParameters
{
   uint_t                               tempInit                          = 3;
   INITIAL_TEMPERATURE_DEVIATION_METHOD initialTemperatureDeviationMethod = INITIAL_TEMPERATURE_DEVIATION_METHOD::WHITE_NOISE;
   real_t                               buoyancyFactor                    = real_c( 0.01 );
   // Single SPH
   uint_t deg = 4;
   int    ord = 2;
   // Superposition SPH
   uint_t                                    lmax                        = 4;
   uint_t                                    lmin                        = 0;
   uint_t                                    superpositionSPHRandomSeed  = 42;
   real_t                                    initialTemperatureSteepness = real_c( 10 );
   std::vector< real_t >                     superpositionRand;
   std::shared_ptr< SphericalHarmonicsTool > sphTool;
};

struct PhysicalParameters
{
   // Either use profiles or constant values (decision for each individually, e.g. temp profile + constant density + viscosity profile is possible)
   // profiles
   std::vector< real_t > radius;
   std::vector< real_t > radiusT;
   std::vector< real_t > radiusCp;
   std::vector< real_t > radiusAlpha;
   std::vector< real_t > radiusDensity;
   std::vector< real_t > temperatureProfile;
   std::vector< real_t > velocityProfile;
   std::vector< real_t > viscosityProfile;
   std::vector< real_t > temperatureInputProfile;
   std::vector< real_t > thermalExpansivityProfile;
   std::vector< real_t > specificHeatCapacityProfile;
   std::vector< real_t > densityProfile;

   //temperature
   //physical versions used to calculate non-D numbers, others used in simulation
   //non-dimensionalisation is set up so that cmb_temp=1 and surface_temp=1 for any inputted physical temperatures
   real_t surfaceTemp = real_c( 300 );
   real_t cmbTemp     = real_c( 4200 );

   //material parameters
   real_t thermalExpansivity         = real_c( 2.238 * 1e-5 );
   real_t thermalExpansivityRadial   = real_c( 2.238 * 1e-5 );
   real_t thermalConductivity        = real_c( 3 );
   real_t specificHeatCapacity       = real_c( 1260 );
   real_t specificHeatCapacityRadial = real_c( 1260 );
   real_t internalHeatingRate        = real_c( 1e-12 );
   real_t density                    = real_c( 0 );
   real_t referenceDensity           = real_c( 4500 );
   real_t surfaceDensity             = real_c( 3300 );
   real_t referenceViscosity         = real_c( 1e22 );
   real_t viscosity                  = real_c( 1e22 );
   real_t grueneisenParameter        = real_c( 1.1 );
   real_t adiabatSurfaceTemp         = real_c( 1600 );
   real_t activationEnergy           = real_c( 5 );
   real_t depthViscosityFactor       = real_c( 3 );
   real_t viscosityLowerBound        = real_c( 1e19 );
   real_t viscosityUpperBound        = real_c( 1e24 );

   //gravity

   real_t gravity = real_c( 9.81 );

   //numbers required to get non-D numbers

   real_t mantleThickness        = real_c( 2900000 );
   real_t thermalDiffusivity     = thermalConductivity / ( referenceDensity * specificHeatCapacity );
   real_t characteristicVelocity = thermalConductivity / ( referenceDensity * specificHeatCapacity * mantleThickness );

   //non-D numbers derived from other parameters

   real_t rayleighNumber = ( referenceDensity * gravity * thermalExpansivity * mantleThickness * mantleThickness *
                             mantleThickness * ( cmbTemp - surfaceTemp ) ) /
                           ( referenceViscosity * thermalDiffusivity );
   real_t pecletNumber      = ( characteristicVelocity * mantleThickness ) / thermalDiffusivity;
   real_t dissipationNumber = ( thermalExpansivity * gravity * mantleThickness ) / specificHeatCapacity;
   real_t hNumber =
       ( internalHeatingRate * mantleThickness ) / ( specificHeatCapacity * characteristicVelocity * ( cmbTemp - surfaceTemp ) );
};

struct TerraNeoParameters
{
   DomainParameters                             domainParameters;
   SolverParameters                             solverParameters;
   OutputParameters                             outputParameters;
   SimulationParameters                         simulationParameters;
   TemperatureDeivationInitialisationParameters initialisationParameters;
   PhysicalParameters                           physicalParameters;
};

} // namespace terraneo