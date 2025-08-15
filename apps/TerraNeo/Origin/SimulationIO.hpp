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
template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::setupOutput()
{
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup Output -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureFE = p2p1StokesFunctionContainer.at( "VelocityFE" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesResidual     = p2p1StokesFunctionContainer.at( "StokesResidual" );

   std::shared_ptr< P2ScalarFunction_T >& referenceTemperature     = p2ScalarFunctionContainer.at( "ReferenceTemperature" );
   std::shared_ptr< P2ScalarFunction_T >& viscosityP2              = p2ScalarFunctionContainer.at( "ViscosityFE" );
   std::shared_ptr< P2ScalarFunction_T >& temperatureP2            = p2ScalarFunctionContainer.at( "TemperatureFE" );
   std::shared_ptr< P2ScalarFunction_T >& temperatureK             = p2ScalarFunctionContainer.at( "Temperature[K]" );
   std::shared_ptr< P2ScalarFunction_T >& temperatureDev           = p2ScalarFunctionContainer.at( "TemperatureDev" );
   std::shared_ptr< P2ScalarFunction_T >& velocityMagnitudeSquared = p2ScalarFunctionContainer.at( "VelocityMagnitudeSquared" );
   std::shared_ptr< P2ScalarFunction_T >& densityP2                = p2ScalarFunctionContainer.at( "DensityFE" );

   std::shared_ptr< P2ScalarFunction_T >& shearHeatingCoeff = p2ScalarFunctionContainer.at( "ShearHeatingTermCoeff" );

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
      vtkOutput_ =
          std::make_shared< hyteg::VTKOutput >( modelParaviewPath, modelBaseName, storage, TN.outputParameters.OutputInterval );
      vtkOutput_->setVTKDataFormat( hyteg::vtk::DataFormat::BINARY );
   }
   else
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      WALBERLA_LOG_INFO_ON_ROOT( "Output format: ADIOS2 - BP5" );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      adiosOutput_ =
          std::make_shared< AdiosWriter >( modelParaviewPath, modelBaseName, TN.outputParameters.ADIOS2OutputConfig, storage );
      adiosOutput_->setParameter( TN.outputParameters.ADIOS2ParamKey, TN.outputParameters.ADIOS2Value );

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

   // std::function< real_t( const Point3D& ) > scaleTRef = [this]( const Point3D& x ) {
   //    real_t refTemp = referenceTemperatureFct( x );
   //    return ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) * refTemp;
   // };

   // Setup output of Temperature and viscosity in SI units [K] and [Pa s], respectively
   referenceTemperature->interpolate( referenceTemperatureFct, TN.domainParameters.maxLevel, All );
   // viscosityPas->assign( { ( TN.physicalParameters.referenceViscosity ) }, { *viscosityP2 }, TN.domainParameters.maxLevel, All );

   auto addFunctionsForOutput = [&, this]< typename Output_T >( const std::shared_ptr< Output_T >& output ) {
      if ( TN.outputParameters.OutputTemperature )
      {
         if ( TN.outputParameters.outputVertexDoFs )
         {
            output->add( temperatureP2->getVertexDoFFunction() );
         }
         else
         {
            output->add( *temperatureP2 );
         }
      }

      if ( TN.outputParameters.OutputVelocity )
      {
         if ( TN.outputParameters.outputVertexDoFs )
         {
            output->add( velocityPressureFE->uvw()[0].getVertexDoFFunction() );
            output->add( velocityPressureFE->uvw()[1].getVertexDoFFunction() );
            output->add( velocityPressureFE->uvw()[2].getVertexDoFFunction() );
         }
         else
         {
            output->add( velocityPressureFE->uvw() );
         }

         output->add( velocityPressureFE->p() );
      }

      if ( TN.outputParameters.OutputViscosity )
      {
         if ( TN.outputParameters.outputVertexDoFs )
         {
            output->add( viscosityP2->getVertexDoFFunction() );
         }
         else
         {
            output->add( *viscosityP2 );
         }
      }

      if ( TN.outputParameters.OutputDensity )
      {
         if ( TN.outputParameters.outputVertexDoFs )
         {
            output->add( densityP2->getVertexDoFFunction() );
         }
         else
         {
            output->add( *densityP2 );
         }
      }

      if ( TN.outputParameters.OutputVerbose )
      {
         if ( TN.outputParameters.outputVertexDoFs )
         {
            output->add( shearHeatingCoeff->getVertexDoFFunction() );
            output->add( temperatureDev->getVertexDoFFunction() );
            output->add( referenceTemperature->getVertexDoFFunction() );
         }
         else
         {
            output->add( *shearHeatingCoeff );
            output->add( *temperatureDev );
            output->add( *referenceTemperature );
         }
      }
   };

   if ( TN.outputParameters.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " VTK is NOT recommended for speed and memory management " );
      addFunctionsForOutput( vtkOutput_ );
   }
   else
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      addFunctionsForOutput( adiosOutput_ );

      // Add attributes to adios2 output
      // There must be a nicer way to collect these attributes
      adiosOutput_->setAllAttributes( attributeList_ );
#else
      WALBERLA_LOG_INFO_ON_ROOT( "No valid output format specified!" );
#endif
   }
}

////////////////////////
// Private functions //
///////////////////////

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::dataOutput()
{
   WALBERLA_LOG_INFO_ON_ROOT( " " );
   WALBERLA_LOG_INFO_ON_ROOT( "**************************" );
   WALBERLA_LOG_INFO_ON_ROOT( "****   Write Output   ****" );
   WALBERLA_LOG_INFO_ON_ROOT( "**************************" );

   std::shared_ptr< P2ScalarFunction_T >& temperatureDev = p2ScalarFunctionContainer.at( "TemperatureDev" );
   std::shared_ptr< P2ScalarFunction_T >& temperatureP2  = p2ScalarFunctionContainer.at( "TemperatureFE" );
   std::shared_ptr< P2ScalarFunction_T >& temperatureK   = p2ScalarFunctionContainer.at( "Temperature[K]" );
   std::shared_ptr< P2ScalarFunction_T >& viscosityPas   = p2ScalarFunctionContainer.at( "Viscosity[Pas]" );
   std::shared_ptr< P2ScalarFunction_T >& viscosityP2    = p2ScalarFunctionContainer.at( "ViscosityFE" );

   const std::string& modelBaseName = TN.outputParameters.modelBaseName;

   if ( TN.outputParameters.OutputTemperature && TN.outputParameters.OutputDimensional )
   {
      temperatureP2->assign( { ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) },
                             { *temperatureP2 },
                             TN.domainParameters.maxLevel,
                             All );
   }

   if ( TN.simulationParameters.tempDependentViscosity && TN.outputParameters.OutputViscosity &&
        TN.outputParameters.OutputDimensional )
   {
      viscosityP2->assign(
          { ( TN.physicalParameters.referenceViscosity ) }, { *viscosityP2 }, TN.domainParameters.maxLevel, All );
   }

   uint_t outputTime = 0u;

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
      vtkOutput_->write( TN.domainParameters.maxLevel, outputTime );
      storage->getTimingTree()->stop( "VTK data output" );
   }
   else if ( TN.outputParameters.dataOutput )
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      if ( TN.simulationParameters.timeStep % TN.outputParameters.OutputInterval == 0u )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "****   Write Output ADIOS2 ****" );
         storage->getTimingTree()->start( "Adios2 data output" );
         adiosOutput_->write( TN.domainParameters.maxLevel, TN.simulationParameters.timeStep );
         storage->getTimingTree()->stop( "Adios2 data output" );
      }

      if ( !TN.outputParameters.outputMyr && TN.simulationParameters.timeStep > 0U && TN.outputParameters.ADIOS2StoreCheckpoint &&
           TN.simulationParameters.timeStep % TN.outputParameters.ADIOS2StoreCheckpointFrequency == 0U )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );

         checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
         checkpointExporter->registerFunction( *temperatureP2, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
         checkpointExporter->storeCheckpoint( modelCheckpointPath, modelBaseName, attributeList_ );
      }
      else if ( TN.outputParameters.outputMyr && TN.outputParameters.ADIOS2StoreCheckpoint &&
                ( ( TN.simulationParameters.modelRunTimeMa ) >=
                  real_c( TN.outputParameters.checkpointCount * TN.outputParameters.ADIOS2StoreCheckpointFrequency ) ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );
         checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
         checkpointExporter->registerFunction( *temperatureP2, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
         checkpointExporter->storeCheckpoint(
             TN.outputParameters.ADIOS2StoreCheckpointPath, TN.outputParameters.ADIOS2StoreCheckpointFilename, attributeList_ );
      }
      else if ( TN.outputParameters.outputMyr && TN.outputParameters.ADIOS2StoreCheckpoint &&
                ( ( TN.simulationParameters.modelRunTimeMa ) >=
                  real_c( TN.outputParameters.checkpointCount * TN.outputParameters.ADIOS2StoreCheckpointFrequency ) ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );
         checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
         checkpointExporter->registerFunction( *temperatureP2, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
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

   if ( TN.outputParameters.OutputTemperature && TN.outputParameters.OutputDimensional )
   {
      temperatureP2->assign( { real_c( 1.0 ) / ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) },
                             { *temperatureP2 },
                             TN.domainParameters.maxLevel,
                             All );
   }

   if ( TN.simulationParameters.tempDependentViscosity && TN.outputParameters.OutputViscosity &&
        TN.outputParameters.OutputDimensional )
   {
      viscosityP2->assign(
          { ( real_c( 1.0 ) / TN.physicalParameters.referenceViscosity ) }, { *viscosityP2 }, TN.domainParameters.maxLevel, All );
   }

   TN.outputParameters.prevOutputTime = std::round( TN.simulationParameters.modelRunTimeMa );
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::outputCheckpoint()
{
   std::shared_ptr< P2ScalarFunction_T >& temperatureP2 = p2ScalarFunctionContainer.at( "TemperatureFE" );

#ifdef HYTEG_BUILD_WITH_ADIOS2
   if ( !TN.outputParameters.outputMyr && TN.simulationParameters.timeStep > 0U &&
        TN.simulationParameters.timeStep % TN.outputParameters.ADIOS2StoreCheckpointFrequency == 0U )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );
      checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
      checkpointExporter->registerFunction( *temperatureP2, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
      checkpointExporter->storeCheckpoint( TN.outputParameters.ADIOS2StoreCheckpointPath,
                                           TN.outputParameters.ADIOS2StoreCheckpointFilename );
   }
   else if ( TN.outputParameters.outputMyr &&
             ( ( TN.simulationParameters.modelRunTimeMa ) >=
               real_c( TN.outputParameters.checkpointCount * TN.outputParameters.ADIOS2StoreCheckpointFrequency ) ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "****   Write Checkpoint ADIOS2 ****" );
      checkpointExporter = std::make_shared< AdiosCheckpointExporter >( TN.outputParameters.ADIOS2OutputConfig );
      checkpointExporter->registerFunction( *temperatureP2, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
      checkpointExporter->storeCheckpoint( TN.outputParameters.ADIOS2StoreCheckpointPath,
                                           TN.outputParameters.ADIOS2StoreCheckpointFilename );
      TN.outputParameters.checkpointCount += 1;
   }
#else
   WALBERLA_LOG_INFO_ON_ROOT( "No valid output format for checkpoint data specified! " );
#endif
}

} // namespace terraneo
