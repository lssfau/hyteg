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

namespace terraneo {
void ConvectionSimulation::setupOutput()
{
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup Output -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   const std::string& outputDirectory = TN.outputParameters.outputDirectory;
   const std::string& modelBaseName   = TN.outputParameters.modelBaseName;

   if ( outputDirectoriesCreated == false )
   {
      modelPath = walberla::format( "%s/%s", outputDirectory.c_str(), modelBaseName.c_str() );

      bool modelPathExists = std::filesystem::is_directory( modelPath );

      WALBERLA_ROOT_SECTION()
      {
         if ( modelPathExists )
         {
            WALBERLA_ABORT( "Path with model name already exists, use override switch to overwrite data" );
         }
         else
         {
            std::filesystem::create_directory( modelPath );
         }
      }

      modelParaviewPath       = walberla::format( "%s/output", modelPath.c_str() );
      modelCheckpointPath     = walberla::format( "%s/checkpoint", modelPath.c_str() );
      modelRadialProfilesPath = walberla::format( "%s/radialProfiles", modelPath.c_str() );

      WALBERLA_ROOT_SECTION()
      {
         std::filesystem::create_directory( modelParaviewPath );
         std::filesystem::create_directory( modelCheckpointPath );
         std::filesystem::create_directory( modelRadialProfilesPath );
      }

      outputDirectoriesCreated = true;
   }

   if ( TN.outputParameters.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Output format: VTK" );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      vtkOutput =
          std::make_shared< hyteg::VTKOutput >( modelParaviewPath, modelBaseName, storage, TN.outputParameters.OutputInterval );
      vtkOutput->setVTKDataFormat( hyteg::vtk::DataFormat::BINARY );
   }
   else
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      WALBERLA_LOG_INFO_ON_ROOT( "Output format: ADIOS2 - BP5" );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      _output =
          std::make_shared< AdiosWriter >( modelParaviewPath, modelBaseName, TN.outputParameters.ADIOS2OutputConfig, storage );
      _output->setParameter( TN.outputParameters.ADIOS2ParamKey, TN.outputParameters.ADIOS2Value );

      if ( TN.outputParameters.ADIOS2StartFromCheckpoint )
      {
         checkpointImporter = std::make_shared< AdiosCheckpointImporter >( TN.outputParameters.ADIOS2StartCheckpointPath,
                                                                           TN.outputParameters.ADIOS2StartCheckpointFilename,
                                                                           TN.outputParameters.ADIOS2OutputConfig );
      }

#else
      WALBERLA_LOG_INFO_ON_ROOT( "No submodule ADIOS2 enabled! " );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "No data output! " );

      WALBERLA_ABORT( "TerraNeo not compiled with ADIOS2 but output is requested with ADIOS2" );
#endif
   }

   std::function< real_t( const Point3D& ) > scaleTRef = [this]( const Point3D& x ) {
      real_t refTemp = referenceTemperatureFct( x );
      return ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) * refTemp;
   };

   // Setup output of Temperature and viscosity in SI units [K] and [Pa s], respectively
   p2ScalarFunctionContainer["ReferenceTemperature[K]"]->interpolate( scaleTRef, TN.domainParameters.maxLevel, All );
   p2ScalarFunctionContainer["Viscosity[Pas]"]->assign( { ( TN.physicalParameters.referenceViscosity ) },
                                                        { *( p2ScalarFunctionContainer["ViscosityFE"] ) },
                                                        TN.domainParameters.maxLevel,
                                                        All );

   if ( TN.outputParameters.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " VTK is NOT recommended for speed and memory management " );
      if ( TN.outputParameters.OutputTemperature )
      {
         if ( TN.outputParameters.outputVertexDoFs )
         {
            vtkOutput->add( p2ScalarFunctionContainer["Temperature[K]"]->getVertexDoFFunction() );
            vtkOutput->add( p2ScalarFunctionContainer["TemperatureDev"]->getVertexDoFFunction() );
         }
         else
         {
            vtkOutput->add( *( p2ScalarFunctionContainer["Temperature[K]"] ) );
            vtkOutput->add( *( p2ScalarFunctionContainer["TemperatureDev"] ) );
         }
      }

      if ( TN.outputParameters.OutputVelocity )
      {
         if ( TN.outputParameters.outputVertexDoFs )
         {
            vtkOutput->add( p2p1StokesFunctionContainer["VelocityFE"]->uvw()[0].getVertexDoFFunction() );
            vtkOutput->add( p2p1StokesFunctionContainer["VelocityFE"]->uvw()[1].getVertexDoFFunction() );
            vtkOutput->add( p2p1StokesFunctionContainer["VelocityFE"]->uvw()[2].getVertexDoFFunction() );
         }
         else
         {
            vtkOutput->add( p2p1StokesFunctionContainer["VelocityFE"]->uvw() );
         }
      }
      else
      {
         if ( TN.outputParameters.outputVertexDoFs )
         {
            vtkOutput->add( p2ScalarFunctionContainer["VelocityMagnitudeSquared"]->getVertexDoFFunction() );
         }
         else
         {
            vtkOutput->add( *( p2ScalarFunctionContainer["VelocityMagnitudeSquared"] ) );
         }
      }
   }
   else
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      if ( TN.outputParameters.OutputTemperature )
      {
         _output->add( *( p2ScalarFunctionContainer["Temperature[K]"] ) );
         _output->add( *( p2ScalarFunctionContainer["TemperatureDev"] ) );
      }

      if ( TN.outputParameters.OutputVelocity )
      {
         if ( TN.outputParameters.outputVertexDoFs )
         {
            _output->add( p2p1StokesFunctionContainer["VelocityFE"]->uvw()[0].getVertexDoFFunction() );
            _output->add( p2p1StokesFunctionContainer["VelocityFE"]->uvw()[1].getVertexDoFFunction() );
            _output->add( p2p1StokesFunctionContainer["VelocityFE"]->uvw()[2].getVertexDoFFunction() );
         }
         else
         {
            _output->add( p2p1StokesFunctionContainer["VelocityFE"]->uvw() );
         }
         _output->add( ( p2p1StokesFunctionContainer["StokesRHS"] )->p() );
      }
      else
      {
         _output->add( *( p2ScalarFunctionContainer["VelocityMagnitudeSquared"] ) );
      }

      _output->add( p2ScalarFunctionContainer["DensityFE"]->getVertexDoFFunction() );
      if ( TN.outputParameters.outputVertexDoFs )
      {
         _output->add( p2ScalarFunctionContainer["Viscosity[Pas]"]->getVertexDoFFunction() );
      }
      else

      {
         _output->add( *( p2ScalarFunctionContainer["Viscosity[Pas]"] ) );
      }
      // Add attributes to adios2 output
      // There must be a nicer way to collect these attributes

      _output->setAllAttributes( attributeList_ );
#else
      WALBERLA_LOG_INFO_ON_ROOT( "No valid output format specified!" );
#endif
   }
}

////////////////////////
// Private functions //
///////////////////////

void ConvectionSimulation::dataOutput()
{
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "**************************" );
   WALBERLA_LOG_INFO_ON_ROOT( "****   Write Output   ****" );
   WALBERLA_LOG_INFO_ON_ROOT( "**************************" );

   const std::string& modelBaseName = TN.outputParameters.modelBaseName;

   if ( !TN.simulationParameters.adaptiveRefTemp || TN.simulationParameters.timeStep == 0 )
   {
      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > temperatureDevFunction =
          [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
             auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
             real_t retVal;
             uint_t shell = static_cast< uint_t >( std::round(
                 real_c( TN.simulationParameters.numLayers ) *
                 ( ( radius - TN.domainParameters.rMin ) / ( TN.domainParameters.rMax - TN.domainParameters.rMin ) ) ) );
             retVal       = ( Temperature[0] - temperatureProfiles->mean.at( shell ) );
             retVal *= ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ); // Redimensionalise Temperature [K]
             return retVal;
          };

      for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; ++l )
      {
         p2ScalarFunctionContainer["TemperatureDev"]->interpolate(
             temperatureDevFunction, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, l, All );
      }
   }
   else
   {
      p2ScalarFunctionContainer["TemperatureDev"]->assign(
          { ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) },
          { *( p2ScalarFunctionContainer["TemperatureDev"] ) },
          TN.domainParameters.maxLevel,
          All );
   }
   p2ScalarFunctionContainer["Temperature[K]"]->assign( { ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) },
                                                        { *( p2ScalarFunctionContainer["TemperatureFE"] ) },
                                                        TN.domainParameters.maxLevel,
                                                        All );

   if ( TN.simulationParameters.tempDependentViscosity )
   {
      p2ScalarFunctionContainer["Viscosity[Pas]"]->assign( { ( TN.physicalParameters.referenceViscosity ) },
                                                           { *( p2ScalarFunctionContainer["ViscosityFE"] ) },
                                                           TN.domainParameters.maxLevel,
                                                           All );
   }

   uint_t outputTime;

   //For circulation model, output with plate age in filename
   if ( TN.simulationParameters.simulationType == "CirculationModel" )
   {
      outputTime = uint_c( std::round( TN.simulationParameters.ageMa ) );
   }

   //For convection model, output with number of timesteps
   else
   {
      outputTime = uint_c( std::round( TN.simulationParameters.timeStep ) );
   }

   if ( TN.outputParameters.dataOutput && TN.outputParameters.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "****    Write Output VTK   ****" );
      storage->getTimingTree()->start( "VTK data output" );
      vtkOutput->write( TN.domainParameters.maxLevel, outputTime );
      storage->getTimingTree()->stop( "VTK data output" );
   }
   else if ( TN.outputParameters.dataOutput )
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      if ( TN.simulationParameters.timeStep % TN.outputParameters.OutputInterval == 0u )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "****   Write Output ADIOS2 ****" );
         storage->getTimingTree()->start( "Adios2 data output" );
         _output->write( TN.domainParameters.maxLevel, TN.simulationParameters.timeStep );
         storage->getTimingTree()->stop( "Adios2 data output" );
      }

      if ( !TN.outputParameters.outputMyr && TN.simulationParameters.timeStep > 0U && TN.outputParameters.ADIOS2StoreCheckpoint &&
           TN.simulationParameters.timeStep % TN.outputParameters.ADIOS2StoreCheckpointFrequency == 0U )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );

         checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
         checkpointExporter->registerFunction(
             *( p2ScalarFunctionContainer["TemperatureFE"] ), TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
         checkpointExporter->storeCheckpoint( modelCheckpointPath, modelBaseName, attributeList_ );
      }
      else if ( TN.outputParameters.outputMyr && TN.outputParameters.ADIOS2StoreCheckpoint &&
                ( ( TN.simulationParameters.modelRunTimeMa ) >=
                  real_c( TN.outputParameters.checkpointCount * TN.outputParameters.ADIOS2StoreCheckpointFrequency ) ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );
         checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
         checkpointExporter->registerFunction(
             *( p2ScalarFunctionContainer["TemperatureFE"] ), TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
         checkpointExporter->storeCheckpoint(
             TN.outputParameters.ADIOS2StoreCheckpointPath, TN.outputParameters.ADIOS2StoreCheckpointFilename, attributeList_ );
      }
      else if ( TN.outputParameters.outputMyr && TN.outputParameters.ADIOS2StoreCheckpoint &&
                ( ( TN.simulationParameters.modelRunTimeMa ) >=
                  real_c( TN.outputParameters.checkpointCount * TN.outputParameters.ADIOS2StoreCheckpointFrequency ) ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );
         checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
         checkpointExporter->registerFunction(
             *( p2ScalarFunctionContainer["TemperatureFE"] ), TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
         checkpointExporter->storeCheckpoint( TN.outputParameters.ADIOS2StoreCheckpointPath,
                                              TN.outputParameters.ADIOS2StoreCheckpointFilename );
         TN.outputParameters.checkpointCount += 1;
      }
#else
      WALBERLA_LOG_INFO_ON_ROOT( "No valid output format specified! " );
#endif
   }

   if ( TN.outputParameters.outputProfiles )
   {
      modelRadialProfilesPath = walberla::format( "%s/radialProfiles", modelPath.c_str() );
      // Redimensionalise tempeature to SI unit [K]
      for ( uint_t i = 0; i < temperatureProfiles->mean.size(); i++ )
      {
         temperatureProfiles->mean[i] *= ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );
         temperatureProfiles->max[i] *= ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );
         temperatureProfiles->min[i] *= ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp );
      }

      temperatureProfiles->logToFile( walberla::format( "%s/%s_TempProfile_%s.dat",
                                                        modelRadialProfilesPath.c_str(),
                                                        modelBaseName.c_str(),
                                                        std::to_string( outputTime ).c_str() ),
                                      "temperature" );

      real_t cmPerYear = 3.15e9;
      for ( uint_t i = 0; i < velocityProfiles->rms.size(); i++ )
      {
         velocityProfiles->mean[i] *= ( TN.physicalParameters.characteristicVelocity * cmPerYear );
         velocityProfiles->max[i] *= ( TN.physicalParameters.characteristicVelocity * cmPerYear );
         velocityProfiles->min[i] *= ( TN.physicalParameters.characteristicVelocity * cmPerYear );
         velocityProfiles->rms[i] *= ( TN.physicalParameters.characteristicVelocity * cmPerYear );
      }

      velocityProfiles->logToFile( walberla::format( "%s/%s_VelocityProfile_%s.dat",
                                                     modelRadialProfilesPath.c_str(),
                                                     modelBaseName.c_str(),
                                                     std::to_string( outputTime ).c_str() ),
                                   "viscosity" );

      if ( TN.simulationParameters.tempDependentViscosity )
      {
         // Redimensionalise viscosity to SI unit [Pa s]
         for ( uint_t i = 0; i < viscosityProfiles->mean.size(); i++ )
         {
            viscosityProfiles->mean[i] *= TN.physicalParameters.referenceViscosity;
            viscosityProfiles->max[i] *= TN.physicalParameters.referenceViscosity;
            viscosityProfiles->min[i] *= TN.physicalParameters.referenceViscosity;
         }

         viscosityProfiles->logToFile( walberla::format( "%s/%s_ViscProfile_%s.dat",
                                                         modelRadialProfilesPath.c_str(),
                                                         modelBaseName.c_str(),
                                                         std::to_string( outputTime ).c_str() ),
                                       "viscosity" );
      }
   }

   TN.outputParameters.prevOutputTime = std::round( TN.simulationParameters.modelRunTimeMa );
}

void ConvectionSimulation::outputCheckpoint()
{
#ifdef HYTEG_BUILD_WITH_ADIOS2
   if ( !TN.outputParameters.outputMyr && TN.simulationParameters.timeStep > 0U &&
        TN.simulationParameters.timeStep % TN.outputParameters.ADIOS2StoreCheckpointFrequency == 0U )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );
      checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
      checkpointExporter->registerFunction(
          *( p2ScalarFunctionContainer["TemperatureFE"] ), TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
      checkpointExporter->storeCheckpoint( TN.outputParameters.ADIOS2StoreCheckpointPath,
                                           TN.outputParameters.ADIOS2StoreCheckpointFilename );
   }
   else if ( TN.outputParameters.outputMyr &&
             ( ( TN.simulationParameters.modelRunTimeMa ) >=
               real_c( TN.outputParameters.checkpointCount * TN.outputParameters.ADIOS2StoreCheckpointFrequency ) ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );
      checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
      checkpointExporter->registerFunction(
          *( p2ScalarFunctionContainer["TemperatureFE"] ), TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
      checkpointExporter->storeCheckpoint( TN.outputParameters.ADIOS2StoreCheckpointPath,
                                           TN.outputParameters.ADIOS2StoreCheckpointFilename );
      TN.outputParameters.checkpointCount += 1;
   }

#else
   WALBERLA_LOG_INFO_ON_ROOT( "No valid output format for checkpoint data specified! " );
#endif
}

} // namespace terraneo
