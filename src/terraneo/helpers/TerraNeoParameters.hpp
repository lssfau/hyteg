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

namespace terraneo {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

// Domain information

struct DomainParameters
{
   //geometric information

   real_t rCMB     = real_c( 3471000 );
   real_t rSurface = real_c( 6371000 );

   //calculate non-dimensional radii such that mantle thickness = 1
   real_t rMin = rCMB / ( rSurface - rCMB );
   real_t rMax = rSurface / ( rSurface - rCMB );

   uint_t nTan     = 2;
   uint_t nRad     = 2;
   uint_t minLevel = 0;
   uint_t maxLevel = 1;

   real_t domainVolume() const
   {
      return ( real_c( 4.0 ) / real_c( 3.0 ) ) * walberla::math::pi * rMax * rMax * rMax -
             ( real_c( 4.0 ) / real_c( 3.0 ) ) * walberla::math::pi * rMin * rMin * rMin;
   }
};

// Solver parameters

struct SolverParameters
{
   uint_t solverFlag = 0u;
   // Stokes solver parameters
   uint_t numPowerIterations         = 25;
   uint_t FGMRESOuterIterations      = 5;
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

   // Uzawa type multigrid solver
   real_t initialResidualU                 = real_c( 0 );
   real_t vCycleResidualUPrev              = real_c( 0 );
   uint_t numVCycles                       = 0;
   real_t averageResidualReductionU        = real_c( 0 );
   uint_t stokesMaxNumIterations           = 5;
   real_t stokesRelativeResidualUTolerance = 1e-6;
   real_t stokesAbsoluteResidualUTolerance = 1e-6;

   // Diffusion solver parameters

   uint_t diffusionMaxNumIterations           = 10000;
   real_t diffusionAbsoluteResidualUTolerance = real_c( 10000 );

   real_t gmresApproximationToleranceTransport = real_c( 1e-5 );
};
// Output parameters
struct OutputParameters
{
   std::string outputDirectory = std::string( "output" );
   std::string outputBaseName  = std::string( "conv_sim" );

   std::string ADIOS2OutputConfig       = std::string( "ADIOS2config.xml" );
   std::string ADIOS2ParamKey           = std::string( "NumAggregators" );
   std::string ADIOS2Value              = std::string( "16" );
   std::string ADIOS2CheckpointPath     = std::string( "output" );
   std::string ADIOS2CheckpointFilename = std::string( "conv_sim" );

   bool ADIOS2StartFromCheckpoint = false;
   bool ADIOS2StoreCheckpoint     = false;

   uint_t ADIOS2StoreCheckpointFrequency = 100U;

   bool   dataOutput        = true;
   bool   vtk               = true;
   bool   OutputVelocity    = true;
   bool   OutputViscosity   = true;
   bool   OutputTemperature = true;
   uint_t OutputInterval    = 1;

   bool   outputMyr         = false;
   uint_t outputIntervalMyr = 1;
   real_t prevOutputTime    = real_c( 0 );

   bool vtkOutputVertexDoFs = true;
   bool outputProfiles      = false;
};
// Simulation parameters
struct SimulationParameters
{
   //Parameters derived from other parameters
   uint_t unknownsTemperature = 0;
   uint_t unknownsStokes      = 0;
   real_t hMin                = real_c( 0 );
   real_t hMax                = real_c( 0 );
   uint_t numLayers           = 0;

   //Parameters given via config file
   bool        fixedTimestep              = false;
   uint_t      timeStep                   = 0;
   real_t      modelTime                  = real_c( 0 );
   real_t      dtPrev                     = real_c( 0 );
   real_t      dt                         = real_c( 0 );
   real_t      dtConstant                 = real_c( 0 );
   real_t      cflMax                     = real_c( 1 );
   uint_t      timestep                   = 0;
   std::string simulationType             = std::string( "ConvectionModel" );
   uint_t      maxNumTimesteps            = 100;
   bool        resetSolver                = false;
   uint_t      resetSolverFrequency       = 100;
   bool        adaptiveRefTemp            = false;
   bool        tempDependentViscosity     = false;
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
   real_t      plateVelocityScaling   = real_c( 1 );
   real_t      plateSmoothingDistance = real_c( 110 );
   bool        compressible           = true; // default: Compressible foŕmulation
   bool        shearHeating           = true; //default: include shear heating
   bool        adiabaticHeating       = true; //default: include adiabatic heating
   bool        internalHeating        = true; //default: include internal heating
   uint_t      boundaryCond           = 1;    // default: No-Slip/No-Slip
   bool        boundaryCondFreeSlip   = false;
   bool        verbose                = false;
   bool        haveViscosityProfile   = false;
   std::string fileViscosityProfile   = std::string( "ViscosityProfile.txt" );

   //needed for conversions in the simulation
   real_t secondsPerMyr = real_c( 3.154e7 * 1e6 );

   // Shear heating scaling for mantle ciruclation model with
   // predifned Lithosphere thickness in km
   real_t shearHeatingScaling  = 1e-5;
   real_t lithosphereThickness = real_c( 100 );

   // Needed for timing analysis of the simulation run
   bool timingAnalysis = true;
};

struct InitialisationParameters
{
   uint_t                tempInit                     = 3;
   uint_t                deg                          = 4;
   int                   ord                          = 2;
   uint_t                lmax                         = 4;
   uint_t                lmin                         = 0;
   bool                  superposition                = false;
   bool                  temperatureNoise             = true;
   bool                  temperatureSphericalHarmonic = false;
   real_t                noiseFactor                  = real_c( 0.1 );
   std::vector< real_t > superpositionRand;
   real_t                buoyancyFactor = real_c( 0.01 );
};

struct PhysicalParameters
{
   // Either use profiles or constant values (decision for each individually, e.g. temp profile + constant density + viscosity profile is possible)
   // profiles
   std::vector< real_t > radius;
   std::vector< real_t > temperatureProfile;
   std::vector< real_t > viscosityProfile;
   real_t                initialTemperatureSteepness = real_c( 10 );

   //temperature
   //physical versions used to calculate non-D numbers, others used in simulation
   //non-dimensionalisation is set up so that cmb_temp=1 and surface_temp=1 for any inputted physical temperatures
   real_t surfaceTemp = real_c( 300 );
   real_t cmbTemp     = real_c( 4200 );

   //material parameters
   real_t thermalExpansivity   = real_c( 2.238 * 1e-5 );
   real_t thermalConductivity  = real_c( 3 );
   real_t specificHeatCapacity = real_c( 1260 );
   real_t internalHeatingRate  = real_c( 1e-12 );
   real_t referenceDensity     = real_c( 4500 );
   real_t surfaceDensity       = real_c( 3300 );
   real_t referenceViscosity   = real_c( 1e22 );
   real_t viscosity            = real_c( 1e22 );
   real_t grueneisenParameter  = real_c( 1.1 );
   real_t adiabatSurfaceTemp   = real_c( 1600 );
   real_t activationEnergy     = real_c( 5 );
   real_t depthViscosityFactor = real_c( 3 );
   real_t viscosityLowerBound  = real_c( 1e19 );
   real_t viscosityUpperBound  = real_c( 1e24 );

   //gravity

   real_t gravity = real_c( 9.81 );

   //numbers required to get non-D numbers

   real_t characteristicVelocity = real_c( 5e-9 );

   real_t mantleThickness    = real_c( 2900000 );
   real_t thermalDiffusivity = thermalConductivity / ( referenceDensity * specificHeatCapacity );

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
   DomainParameters         domainParameters;
   SolverParameters         solverParameters;
   OutputParameters         outputParameters;
   SimulationParameters     simulationParameters;
   InitialisationParameters initialisationParameters;
   PhysicalParameters       physicalParameters;
};

} // namespace terraneo