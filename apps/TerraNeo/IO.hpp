#pragma once

#include "Convection.hpp"

namespace terraneo {
void ConvectionSimulation::setupOutput()
{
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup Output -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   if ( TN.outputParameters.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Output format: VTK" );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      vtkOutput = std::make_shared< hyteg::VTKOutput >(
          TN.outputParameters.outputDirectory, TN.outputParameters.outputBaseName, storage, TN.outputParameters.OutputInterval );
      vtkOutput->setVTKDataFormat( hyteg::vtk::DataFormat::BINARY );
   }
   else
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      WALBERLA_LOG_INFO_ON_ROOT( "Output format: ADIOS2 - BP4" );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      _output = std::make_shared< AdiosWriter >(
          TN.outputParameters.outputDirectory, TN.outputParameters.outputBaseName, TN.outputParameters.outputConfig, storage );
      _output->setParameter( TN.outputParameters.ADIOS2ParamKey, TN.outputParameters.ADIOS2Value );
#else
      WALBERLA_LOG_INFO_ON_ROOT( " No submodule ADIOS2 enabled. No data output " );
#endif
   }

   if ( TN.outputParameters.vtk )
   {
      if ( TN.outputParameters.OutputTemperature )
      {
         if ( TN.outputParameters.vtkOutputVertexDoFs )
         {
            vtkOutput->add( temperature->getVertexDoFFunction() );
            vtkOutput->add( temperatureDev->getVertexDoFFunction() );
         }
         else
         {
            vtkOutput->add( *temperature );
            vtkOutput->add( *temperatureDev );
         }
      }

      if ( TN.outputParameters.OutputVelocity )
      {
         if ( TN.outputParameters.vtkOutputVertexDoFs )
         {
            vtkOutput->add( stokesLHS->uvw()[0].getVertexDoFFunction() );
            vtkOutput->add( stokesLHS->uvw()[1].getVertexDoFFunction() );
            vtkOutput->add( stokesLHS->uvw()[2].getVertexDoFFunction() );
         }
         else
         {
            vtkOutput->add( stokesLHS->uvw()[0] );
            vtkOutput->add( stokesLHS->uvw()[1] );
            vtkOutput->add( stokesLHS->uvw()[2] );
         }
      }
      else
      {
         if ( TN.outputParameters.vtkOutputVertexDoFs )
         {
            vtkOutput->add( velocityMagnitudeSquared->getVertexDoFFunction() );
         }
         else
         {
            vtkOutput->add( *velocityMagnitudeSquared );
         }
      }
   }
   else
   {
#ifdef HYTEG_BUILD_WITH_ADIOS2
      if ( TN.outputParameters.OutputTemperature )
      {
         _output->add( *temperature );
         _output->add( *temperatureDev );
      }

      if ( TN.outputParameters.OutputVelocity )
      {
         _output->add( stokesLHS->uvw()[0] );
         _output->add( stokesLHS->uvw()[1] );
         _output->add( stokesLHS->uvw()[2] );

         // stokes RHS velocity

         _output->add( stokesRHS->uvw()[0] );
         _output->add( stokesRHS->uvw()[1] );
         _output->add( stokesRHS->uvw()[2] );

         // stokes RHS pressure field
         _output->add( stokesRHS->p() );
      }
      else
      {
         _output->add( *velocityMagnitudeSquared );
      }

      _output->add( *temperatureReference );

      _output->add( *diffusionFE );
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

   if ( !TN.simulationParameters.adaptiveRefTemp )
   {
      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > temperatureDevFunction =
          [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
             auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
             real_t retVal;

             uint_t shell = static_cast< uint_t >( std::round(
                 real_c( TN.simulationParameters.numLayers ) *
                 ( ( radius - TN.domainParameters.rMin ) / ( TN.domainParameters.rMax - TN.domainParameters.rMin ) ) ) );
             retVal       = ( Temperature[0] - temperatureProfiles->mean.at( shell ) );
             return retVal;
          };

      for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; ++l )
      {
         temperatureDev->interpolate( temperatureDevFunction, { *temperature }, l, All );
      }
   }

   uint_t outputTime;

   //with plates, output with plate age in filename
   if ( TN.simulationParameters.simulationType == "CirculationModel" )
   {
      outputTime = uint_c( std::round( TN.simulationParameters.ageMa ) );
   }

   //otherwise, output with time since beginning of simulation
   else
   {
      outputTime = uint_c( std::round( TN.simulationParameters.modelRunTimeMa ) );
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
      WALBERLA_LOG_INFO_ON_ROOT( "****   Write Output ADIOS2 ****" );
      storage->getTimingTree()->start( "Adios2 data output" );
      _output->write( TN.domainParameters.maxLevel, TN.simulationParameters.timeStep );
      storage->getTimingTree()->stop( "Adios2 data output" );
#else
      WALBERLA_LOG_INFO_ON_ROOT( "No valid output format specified! " );
#endif
   }

   if ( TN.outputParameters.outputProfiles )
   {
      temperatureProfiles->logToFile( TN.outputParameters.outputDirectory + "/" + "Profiles" + "/" +
                                          TN.outputParameters.outputBaseName + "_TempProfile_" + std::to_string( outputTime ) +
                                          ".dat",
                                      "temperature" );

      if ( TN.simulationParameters.tempDependentViscosity )
      {
         viscosityProfiles->logToFile( TN.outputParameters.outputDirectory + "/" + "Profiles" + "/" +
                                           TN.outputParameters.outputBaseName + "_ViscProfile_" + std::to_string( outputTime ) +
                                           ".dat",
                                       "viscosity" );
      }
   }

   TN.outputParameters.prevOutputTime = std::round( TN.simulationParameters.modelRunTimeMa );
}

} // namespace terraneo