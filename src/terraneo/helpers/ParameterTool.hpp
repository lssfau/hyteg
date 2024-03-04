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

#include "core/extern/json.hpp"

#include "hyteg/Levelinfo.hpp"

#include "terraneo/helpers/TerraNeoDataStructures.hpp"
#include "terraneo/plates/PlateVelocityProvider.hpp"

namespace terraneo {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using json = nlohmann::json;

/**
 * @brief Reads data from a file and populates a 2D vector with the values.
 *
 * This function reads data from a file (JSON or txt/csv) specified by the filename parameter (i.e. a radial viscosity or temperature profile). 
 * The function reads the data and populates the data_vector parameter, which is a 2D vector of real_t values. 
 * The num_columns parameter specifies the number of columns expected in each row.
 *
 * @param filename The name of the file to read the data from.
 * @param data_vector The 2D vector to populate with the data.
 * @param num_columns The number of columns expected in each row.
 * @return True if the data was successfully read and populated in the data_vector, false otherwise.
 */

inline bool readDataFile( const std::string& filename, std::vector< std::vector< real_t > >& data_vector, uint_t num_columns )
{
   std::string           line;
   std::ifstream         myfile( filename );
   uint_t                counter = 0;
   real_t                tmp_var;
   std::vector< real_t > tmp_vec;
   tmp_vec.reserve( num_columns );

   if ( myfile.is_open() )
   {
      bool gotData = false;

      while ( getline( myfile, line ) )
      {
         tmp_vec.clear();
         std::stringstream sstream( line );

         // Check if the file is in JSON format
         if ( filename.substr( filename.find_last_of( "." ) + 1 ) == "json" )
         {
            json j;
            try
            {
               j = json::parse( line );
            } catch ( const std::exception& e )
            {
               WALBERLA_LOG_WARNING_ON_ROOT( "Failed to parse JSON in file " << filename << ": " << e.what() );
               return gotData;
            }
            // Check if the number of arrays in the JSON file is 2
            if ( j.size() != 2 && j.size() == num_columns )
            {
               WALBERLA_LOG_WARNING_ON_ROOT( "Invalid number of arrays in JSON file " << filename
                                                                                      << ": Expected: " << num_columns );
               WALBERLA_LOG_WARNING_ON_ROOT( "Number of arrays found: " << j.size() );
               return gotData;
            }
            // Read values from JSON array
            for ( const auto& item : j.items() )
            {
               const std::string& key   = item.key();
               const json&        value = item.value();

               // Check if the key exists in the JSON object
               if ( !value.is_array() )
               {
                  WALBERLA_LOG_WARNING_ON_ROOT( "Invalid JSON format in file " << filename << ": Expected an array for key "
                                                                               << key );
                  return gotData;
               }
               for ( const auto& val : value )
               {
                  if ( !val.is_number() )
                  {
                     WALBERLA_LOG_WARNING_ON_ROOT( "Invalid value type for key " << key << " in JSON file " << filename
                                                                                 << ": Expected a number" );
                     return gotData;
                  }
                  tmp_vec.push_back( val.get< real_t >() );
               }
            }

            // Add the temporary vector to the data vector
            data_vector.push_back( tmp_vec );
            // Increment the counter
            counter++;
            gotData = true;
         }
         else
         {
            // File is in .txt format
            // Cycle through columns and read values into tmp_var
            for ( uint_t i = 0; i < num_columns; ++i )
            {
               if ( sstream.eof() )
               {
                  if ( i < num_columns - 1 )
                  {
                     WALBERLA_LOG_WARNING_ON_ROOT( "Not enough elements in row " << counter << " of file " << filename );
                     return gotData;
                  }
                  else
                  {
                     break; // Skip the warning if it's the last column in the row
                  }
               }

               sstream >> tmp_var;
               tmp_vec.push_back( tmp_var );
            }
         }

         // Append row to vector
         data_vector.push_back( tmp_vec );
         gotData = true;
         // Move to next row
         counter++;
      }
      myfile.close();
   }
   else
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Unable to open file " << filename );
      return false;
   }

   return true;
}

/**
 * @brief Performs linear interpolation between data points in a 2D vector.
 *
 * This function takes a 2D vector of real_t values and performs linear interpolation to estimate the value at a given independent variable.
 * The independent variable is compared with the values in the first column of the data_vector to determine the appropriate data points for interpolation.
 * Linear interpolation is then performed between these data points to estimate the value at the given independent variable.
 *
 * @param data_vector The 2D vector containing the data points.
 * @param independent_var The independent variable for which the value is to be estimated.
 * @return The estimated value at the given independent variable.
 */

real_t linearInterpolateBetween( std::vector< std::vector< real_t > >& data_vector, const real_t& independent_var )
{
   if ( independent_var <= data_vector[0][0] )
      return data_vector[0][1];

   if ( independent_var > data_vector[data_vector.size() - 1][0] )
      return data_vector[data_vector.size() - 1][1];

   uint_t cur_row = 0;
   while ( independent_var > data_vector[cur_row][0] )
      ++cur_row;

   real_t normalized_factor =
       ( independent_var - data_vector[cur_row - 1][0] ) / ( data_vector[cur_row][0] - data_vector[cur_row - 1][0] );
   real_t retVal =
       ( normalized_factor * ( data_vector[cur_row][1] - data_vector[cur_row - 1][1] ) ) + data_vector[cur_row - 1][1];
   return retVal;
}

/**
 * @brief Parses the configuration parameters from the main configuration block.
 *
 * This function reads and extracts various domain, model, simulation, and initialization parameters from the main configuration block.
 * It populates the corresponding variables and performs necessary calculations for non-dimensional numbers.
 *
 * @param mainConf The main configuration block containing the parameters.
 */

DomainParameters         domainParam;
SolverParameters         solverParam;
OutputParameters         outputParam;
SimulationParameters     simulationParam;
PhysicalParameters       physicalParam;
InitialisationParameters initialisationParam;

inline void parseConfig( const walberla::Config::BlockHandle& mainConf )
{
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
      if ( readDataFile( simulationParam.fileViscosityProfile, physicalParam.viscosityProfile, 2 ) )
      {
         simulationParam.haveViscosityProfile = true;

         for ( uint_t i = 0; i < physicalParam.viscosityProfile.size(); i++ )
         {
            //non-dimensionalise radius
            physicalParam.viscosityProfile[i][0] /= physicalParam.mantleThickness;
         }

         std::reverse( physicalParam.viscosityProfile.begin(), physicalParam.viscosityProfile.end() );
      }
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
   simulationParam.boundaryCond           = mainConf.getParameter< uint_t >( "boundaryCond" );
   solverParam.solverType                 = mainConf.getParameter< uint_t >( "SolverType" );
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

         auto oracle = std::make_shared< terraneo::plates::PlateVelocityProvider >( simulationParam.fnameTopologies,
                                                                                    simulationParam.fnameReconstructions );
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

   solverParam.stokesMaxNumIterations           = mainConf.getParameter< uint_t >( "stokesMaxNumIterations" );
   solverParam.stokesAbsoluteResidualUTolerance = mainConf.getParameter< real_t >( "stokesAbsoluteResidualUTolerance" );
   solverParam.stokesRelativeResidualUTolerance = mainConf.getParameter< real_t >( "stokesRelativeResidualUTolerance" );

   solverParam.coarseGridAbsoluteResidualTolerance = mainConf.getParameter< real_t >( "coarseGridAbsoluteResidualTolerance" );
   solverParam.coarseGridRelativeResidualTolerance = mainConf.getParameter< real_t >( "coarseGridRelativeResidualTolerance" );

   solverParam.uzawaOmega                      = mainConf.getParameter< real_t >( "uzawaOmega" );
   solverParam.uzawaInnerIterations            = mainConf.getParameter< uint_t >( "uzawaInnerIterations" );
   solverParam.uzawaPreSmooth                  = mainConf.getParameter< uint_t >( "uzawaPreSmooth" );
   solverParam.uzawaPostSmooth                 = mainConf.getParameter< uint_t >( "uzawaPostSmooth" );
   solverParam.numVCyclesPerLevel              = mainConf.getParameter< uint_t >( "numVCyclesPerLevel" );
   solverParam.estimateUzawaOmega              = mainConf.getParameter< bool >( "estimateUzawaOmega" );
   solverParam.fullMultigrid                   = mainConf.getParameter< bool >( "fullMultigrid" );
   solverParam.preComputeStokesElementMatrices = mainConf.getParameter< bool >( "preComputeStokesElementMatrices" );
   solverParam.chebyshevIterations             = mainConf.getParameter< uint_t >( "chebyshevIterations" );
   solverParam.stokesKillTolerance             = mainConf.getParameter< real_t >( "stokesKillTolerance" );

   solverParam.diffusionMaxNumIterations           = mainConf.getParameter< uint_t >( "diffusionMaxNumIterations" );
   solverParam.diffusionAbsoluteResidualUTolerance = mainConf.getParameter< real_t >( "diffusionAbsoluteResidualUTolerance" );

   solverParam.gmresApproximationToleranceTransport = mainConf.getParameter< real_t >( "gmresApproximationToleranceTransport" );

   solverParam.solverType = mainConf.getParameter< uint_t >( "SolverType" );

   outputParam.outputDirectory     = mainConf.getParameter< std::string >( "outputDirectory" );
   outputParam.outputBaseName      = mainConf.getParameter< std::string >( "outputBaseName" );
   outputParam.dataOutput          = mainConf.getParameter< bool >( "dataOutput" );
   outputParam.ADIOS2ParamKey      = mainConf.getParameter< std::string >( "ADIOS2ParamKey" );
   outputParam.ADIOS2Value         = mainConf.getParameter< std::string >( "ADIOS2Value" );
   outputParam.vtk                 = mainConf.getParameter< bool >( "vtk" );
   outputParam.outputConfig        = mainConf.getParameter< std::string >( "outputConfig" );
   outputParam.OutputVelocity      = mainConf.getParameter< bool >( "OutputVelocity" );
   outputParam.OutputTemperature   = mainConf.getParameter< bool >( "OutputTemperature" );
   outputParam.OutputInterval      = mainConf.getParameter< uint_t >( "OutputInterval" );
   outputParam.vtkOutputVertexDoFs = mainConf.getParameter< bool >( "OutputVertexDoFs" );

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
}

/**
 * @brief Prints the configuration parameters to the log file.
 *
 * This function prints the domain, model, physical, non-dimensional, initialisation, and simulation parameters
 * to the root process and to the log file.
 */

inline void printConfig()
{
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
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes Max Num Iterations              : " << solverParam.stokesMaxNumIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes absolute residual U-tolerance   : " << solverParam.stokesAbsoluteResidualUTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes relative residual U-tolerance   : " << solverParam.stokesRelativeResidualUTolerance );

   WALBERLA_LOG_INFO_ON_ROOT( "Coarse grid absolute residual tolerance: " << solverParam.coarseGridAbsoluteResidualTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "Coarse grid relative residual tolerance: " << solverParam.coarseGridRelativeResidualTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "Uzawa inner iterations                 : " << solverParam.uzawaInnerIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "Uzawa pre-smooth                       : " << solverParam.uzawaPreSmooth );
   WALBERLA_LOG_INFO_ON_ROOT( "Uzawa post-smooth                      : " << solverParam.uzawaPostSmooth );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of V-Cycles per level           : " << solverParam.numVCyclesPerLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "Uzawa Omega                            : " << solverParam.uzawaOmega );
   WALBERLA_LOG_INFO_ON_ROOT( "Diffusion max num iterations           : " << solverParam.diffusionMaxNumIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "Diffusion absolute residual U-tolerance: " << solverParam.diffusionAbsoluteResidualUTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes kill-tolerance                  : " << solverParam.stokesKillTolerance );

   WALBERLA_LOG_INFO_ON_ROOT( "Full-multigrid                         : " << ( solverParam.fullMultigrid ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT(
       "pre-compute Stokes Element Matrices    : " << ( solverParam.preComputeStokesElementMatrices ? "true" : "false" ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Chebyshev iterations                   : " << solverParam.chebyshevIterations );

   WALBERLA_ROOT_SECTION()
   {
      walberla::logging::Logging::instance()->stopLoggingToFile();
   }
}

} // namespace terraneo