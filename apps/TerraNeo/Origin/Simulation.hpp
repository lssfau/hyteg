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
#pragma once

#include "Convection.hpp"
#include "terraneo/dataimport/ParameterIO.hpp"

////////////////////////
//   Initialisation   //
////////////////////////

#include "ModelInit.hpp"
#include "SimulationHelpers.hpp"
#include "SimulationIO.hpp"

namespace terraneo {

/////////////////////////
//   Simulation Step   //
/////////////////////////

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::step()
{
   walberla::WcTimer localTimerStep;
   if ( TN.simulationParameters.timeStep == 0 )
   {
      if ( TN.outputParameters.ADIOS2StartFromCheckpoint )
      {
#ifdef HYTEG_BUILD_WITH_ADIOS2
         WALBERLA_LOG_INFO_ON_ROOT( "ADIOS2 load Temperature field from Checkpoint" );
         // Get Information about checkpoint data
         std::vector< AdiosCheckpointImporter::FunctionDescription > checkpointDataInfo =
             checkpointImporter->getFunctionDetails();
         // Import temperature field from checkpoint
         // restore function will only work out of the box at the same level -> if new level is higher -> prolongate, if lower -> restrict
         // Check if current simulation max Level is compatible with checkpoint data maxLevel
         if ( TN.domainParameters.maxLevel > checkpointDataInfo[0].maxLevel )
         {
            WALBERLA_LOG_INFO_ON_ROOT(
                "Restore Temperature field from checkpoint data maxLevel: " << checkpointDataInfo[0].maxLevel );
            checkpointImporter->restoreFunction( *( p2ScalarFunctionContainer["TemperatureFE"] ),
                                                 TN.domainParameters.minLevel,
                                                 checkpointDataInfo[0].maxLevel,
                                                 0,
                                                 true );

            WALBERLA_LOG_INFO_ON_ROOT( "Prolongate Temperature field checkpoint data from level: "
                                       << checkpointDataInfo[0].maxLevel
                                       << " to new maxLevel: " << TN.domainParameters.maxLevel );
            for ( uint_t l = checkpointDataInfo[0].maxLevel; l < TN.domainParameters.maxLevel; l++ )
            {
               p2ProlongationOperator_->prolongate( *( p2ScalarFunctionContainer["TemperatureFE"] ), l, All );
               WALBERLA_LOG_INFO_ON_ROOT( "Temperatuer field prolongated from level:  " << l << " to level: " << l + 1 );
            }
         }
         else
         {
            WALBERLA_LOG_INFO_ON_ROOT(
                "Restore Temperature field from checkpoint data at level: " << checkpointDataInfo[0].maxLevel );
            checkpointImporter->restoreFunction( *( p2ScalarFunctionContainer["TemperatureFE"] ),
                                                 TN.domainParameters.minLevel,
                                                 TN.domainParameters.maxLevel,
                                                 0,
                                                 true );
         }
#else
         WALBERLA_ABORT( " No submodule ADIOS2 enabled! Loading Checkpoint not possible - Abort simulation " );
#endif
         solveStokes();
         dataOutput();
      }
      else
      {
         dataOutput();
         solveStokes();
      }

      ++TN.simulationParameters.timeStep;
      p2p1StokesFunctionContainer["VelocityFEPrev"]->assign(
          { real_c( 1 ) }, { *( p2p1StokesFunctionContainer["VelocityFE"] ) }, TN.domainParameters.maxLevel, All );
   } //end timestep0 stokes

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "-------- Time step: " << TN.simulationParameters.timeStep << " --------" );

   real_t vMax = velocityMaxMagnitude( p2p1StokesFunctionContainer["VelocityFE"]->uvw(),
                                       p2p1StokesFunctionContainer["StokesTmp1"]->uvw().component( 0u ),
                                       *( p2ScalarFunctionContainer["VelocityMagnitudeSquared"] ),
                                       TN.domainParameters.maxLevel,
                                       All );

   TN.simulationParameters.dtPrev = TN.simulationParameters.dt;

   if ( TN.simulationParameters.fixedTimestep )
   {
      TN.simulationParameters.dt = TN.simulationParameters.dtConstant *
                                   ( TN.physicalParameters.characteristicVelocity * TN.simulationParameters.plateVelocityScaling *
                                     TN.simulationParameters.secondsPerMyr ) /
                                   TN.physicalParameters.mantleThickness;
   }

   else
   {
      // Within the first 20 timesteps gradually increase the timestep size from 1000 a to 100000 a
      if ( TN.simulationParameters.timeStep < 20 )
      {
         // Initial Step size is 1000 a
         real_t stepSize = 0.001;
         TN.simulationParameters.dt =
             ( TN.simulationParameters.timeStep * 10 * stepSize * TN.physicalParameters.characteristicVelocity *
               TN.simulationParameters.plateVelocityScaling * TN.simulationParameters.secondsPerMyr ) /
             TN.physicalParameters.mantleThickness;
      }
      else
      {
         TN.simulationParameters.dt = ( TN.simulationParameters.cflMax / vMax ) * TN.simulationParameters.hMin;
      }
      // TN.simulationParameters.dt = ( TN.simulationParameters.cflMax / vMax ) * TN.simulationParameters.hMin;
   }

   real_t stepSize = ( TN.simulationParameters.dt * TN.physicalParameters.mantleThickness ) /
                     ( TN.physicalParameters.characteristicVelocity * TN.simulationParameters.plateVelocityScaling *
                       TN.simulationParameters.secondsPerMyr );
   // Limit the timestep size to a maximum value (in Ma).
   if ( stepSize > ( TN.simulationParameters.maxTimestepSize ) )
   {
      TN.simulationParameters.dt = ( ( TN.simulationParameters.maxTimestepSize ) *
                                     ( TN.physicalParameters.characteristicVelocity *
                                       TN.simulationParameters.plateVelocityScaling * TN.simulationParameters.secondsPerMyr ) ) /
                                   ( TN.physicalParameters.mantleThickness );
      stepSize = TN.simulationParameters.maxTimestepSize;
   }

   TN.simulationParameters.modelTime += TN.simulationParameters.dt;

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Step size: " << stepSize << " Ma" );
   WALBERLA_LOG_INFO_ON_ROOT( "Velocity max magnitude: " << vMax << " " );

   // convert vMax to SI-unit velocity in [cm/a]
   real_t MaxVelocityMagSI = ( ( 100.0 * 365.0 * 24.0 * 3600.0 ) * ( vMax * TN.physicalParameters.characteristicVelocity ) );

   if ( TN.simulationParameters.simulationType == "CirculationModel" )
   {
      //update current and previous age
      TN.simulationParameters.agePrev = TN.simulationParameters.ageMa;
      TN.simulationParameters.ageMa -= ( TN.simulationParameters.dt * TN.physicalParameters.mantleThickness ) /
                                       ( TN.physicalParameters.characteristicVelocity *
                                         TN.simulationParameters.plateVelocityScaling * TN.simulationParameters.secondsPerMyr );

      //if the time goes below finalAge Ma we return from current step and finish the circulation model in main
      if ( TN.simulationParameters.ageMa <= TN.simulationParameters.finalAge )
      {
         return;
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Model age: " << TN.simulationParameters.ageMa << " Ma" )
   }

   TN.simulationParameters.modelRunTimeMa =
       ( TN.simulationParameters.modelTime * TN.physicalParameters.mantleThickness ) /
       ( TN.physicalParameters.characteristicVelocity * TN.simulationParameters.plateVelocityScaling *
         TN.simulationParameters.secondsPerMyr );

   WALBERLA_LOG_INFO_ON_ROOT( "Model runtime: " << TN.simulationParameters.modelRunTimeMa << " Ma" )

   p2ScalarFunctionContainer["TemperatureFE"]->assign(
       { real_c( 1 ) }, { *( p2ScalarFunctionContainer["TemperaturePrev"] ) }, TN.domainParameters.maxLevel, All );

   //######################################################//
   //                  ENERGY EQUATION                     //
   //######################################################//

   /*############ ADVECTION STEP ############*/

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Advection Step ------" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   const std::shared_ptr< TemperatureFunction_T >& temperatureFE = [&] {
      if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
      {
         return p2ScalarFunctionContainer.at( "TemperatureFE" );
      }
      else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
      {
         return p1ScalarFunctionContainer.at( "TemperatureFEP1" );
      }
      else
      {
         WALBERLA_ABORT( "Unknown type" );
      }
   }();

   if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
   {
      temperatureMMOCOperator_->step( *temperatureFE,
                                      p2p1StokesFunctionContainer["VelocityFE"]->uvw(),
                                      p2p1StokesFunctionContainer["VelocityFEPrev"]->uvw(),
                                      TN.domainParameters.maxLevel,
                                      hyteg::Inner,
                                      TN.simulationParameters.dt,
                                      1u,
                                      true,
                                      true,
                                      false );
   }
   else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
   {
      P2toP1Conversion( p2p1StokesFunctionContainer["VelocityFEPrev"]->uvw(),
                        *( p1VectorFunctionContainer["VelocityFEPrevP1"] ),
                        TN.domainParameters.maxLevel + 1,
                        All );

      P2toP1Conversion( p2p1StokesFunctionContainer["VelocityFE"]->uvw(),
                        *( p1VectorFunctionContainer["VelocityFEP1"] ),
                        TN.domainParameters.maxLevel + 1,
                        All );

      P2toP1Conversion( *( p2ScalarFunctionContainer["TemperatureFE"] ), *temperatureFE, TN.domainParameters.maxLevel + 1, All );

      temperatureMMOCOperator_->step( *temperatureFE,
                                      *( p1VectorFunctionContainer["VelocityFEP1"] ),
                                      *( p1VectorFunctionContainer["VelocityFEPrevP1"] ),
                                      TN.domainParameters.maxLevel + 1,
                                      hyteg::Inner,
                                      TN.simulationParameters.dt,
                                      1u,
                                      true,
                                      false,
                                      false );
   }
   else
   {
      WALBERLA_ABORT( "Unsupported function type" );
   }

   // Reset temperature on boundary to initial values

   // p2ScalarFunctionContainer["TemperaturePrev"]->assign(
   //     { real_c( 1 ) }, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, TN.domainParameters.maxLevel, All );

   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      std::function< real_t( const Point3D& ) > initTemperature;
      switch ( TN.initialisationParameters.initialTemperatureDeviationMethod )
      {
      case INITIAL_TEMPERATURE_DEVIATION_METHOD::SINGLE_SPH:
         initTemperature = temperatureSingleSPH( *temperatureInitParams, *temperatureReferenceFct );
         break;
      case INITIAL_TEMPERATURE_DEVIATION_METHOD::RANDOM_SUPERPOSITION_SPH:
         initTemperature = temperatureRandomSuperpositioneSPH( *temperatureInitParams, *temperatureReferenceFct );
         break;
      case INITIAL_TEMPERATURE_DEVIATION_METHOD::WHITE_NOISE:
         initTemperature = temperatureWhiteNoise( *temperatureInitParams, *temperatureReferenceFct );
         break;
      default:
         WALBERLA_ABORT( "Unknown initial temperature deviation method" );
      }

      p2ScalarFunctionContainer["TemperatureFE"]->interpolate( initTemperature, l, DirichletBoundary );
      p1ScalarFunctionContainer["TemperatureFEP1"]->interpolate( initTemperature, l, DirichletBoundary );
   }

   /*############ ENERGY STEP ############*/

   // setupEnergyRHS();
   localTimerStep.start();
   storage->getTimingTree()->start( "TerraNeo solve Energy equation" );
   solveEnergy();
   storage->getTimingTree()->stop( "TerraNeo solve Energy equation" );
   localTimerStep.end();
   real_t timerSolveEnergy = localTimerStep.last();

   //######################################################//
   //                  STOKES EQUATIONS                    //
   //######################################################//

   // update velocity field storing the velocity field of the Prev timestep
   p2p1StokesFunctionContainer["VelocityFEPrev"]->assign(
       { real_c( 1 ) }, { *( p2p1StokesFunctionContainer["VelocityFE"] ) }, TN.domainParameters.maxLevel, All );

   if ( TN.simulationParameters.simulationType == "CirculationModel" )
   {
      //update plate velocity boundary condition prior to stokes solve
      //update if >=1Myr has passed since Prev update
      if ( ( ( TN.simulationParameters.plateAge - TN.simulationParameters.ageMa ) >= real_c( 1 ) ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "" );
         WALBERLA_LOG_INFO_ON_ROOT( "Update plates" );

         //save the age of the current update (rounded to 1Myr intervals)
         TN.simulationParameters.plateAge = std::round( TN.simulationParameters.ageMa );
         updatePlateVelocities( *( p2p1StokesFunctionContainer["VelocityFE"] ) );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "Plate age: " << TN.simulationParameters.plateAge << " Ma" )
   }

   //update ref temp vector based on new temperature field

   // P1toP2Conversion( *( p1ScalarFunctionContainer["TemperatureFEP1"] ),
   //                   *( p2ScalarFunctionContainer["TemperatureFE"] ),
   //                   TN.domainParameters.maxLevel - 1,
   //                   All );

   // P2toP2QuadraticProlongation prolongationP2;

   // prolongationP2.prolongate(*( p2ScalarFunctionContainer["TemperatureFE"] ), TN.domainParameters.maxLevel - 1, All );

   temperatureProfiles = std::make_shared< RadialProfile >( computeRadialProfile( *( p2ScalarFunctionContainer["TemperatureFE"] ),
                                                                                  TN.domainParameters.rMin,
                                                                                  TN.domainParameters.rMax,
                                                                                  TN.domainParameters.nRad,
                                                                                  TN.domainParameters.maxLevel ) );

   TN.physicalParameters.temperatureProfile = temperatureProfiles->mean;
   // calculateHeatflow( temperatureProfiles );
   // calculateHeatflowIntegral( temperatureProfiles );

   if ( TN.simulationParameters.checkTemperatureConsistency )
   {
      // Consistency check for out of range Temperatures
      for ( uint_t i = 0; i < temperatureProfiles->min.size(); i++ )
      {
         // Redimensionalise temperature
         real_t dimFactor = TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp;

         if ( ( temperatureProfiles->min[i] ) <
                  TN.physicalParameters.surfaceTemp - TN.simulationParameters.temperatureConsistencyThreshold ||
              ( temperatureProfiles->max[i] ) >
                  TN.physicalParameters.cmbTemp + TN.simulationParameters.temperatureConsistencyThreshold )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Out_of_Range Temperature: " << temperatureProfiles->min[i] * dimFactor
                                                                    << " detected at shell radii: "
                                                                    << temperatureProfiles->shellRadii[i] );
            WALBERLA_LOG_INFO_ON_ROOT( "Dump data" );
            dataOutput();
            WALBERLA_ABORT( "Aborting simulation run" );
         }
      }
   }

   solveStokes();
   // Compute rms velocity radially
   velocityProfiles = std::make_shared< RadialProfile >( computeRadialProfile( p2p1StokesFunctionContainer["VelocityFE"]->uvw(),
                                                                               TN.domainParameters.rMin,
                                                                               TN.domainParameters.rMax,
                                                                               TN.domainParameters.nRad,
                                                                               TN.domainParameters.maxLevel ) );
   TN.physicalParameters.velocityProfile = velocityProfiles->rms;

   // update viscosity Profiles for logging
   if ( TN.simulationParameters.tempDependentViscosity )
   {
      viscosityProfiles = std::make_shared< RadialProfile >( computeRadialProfile( *( p2ScalarFunctionContainer["ViscosityFE"] ),
                                                                                   TN.domainParameters.rMin,
                                                                                   TN.domainParameters.rMax,
                                                                                   TN.domainParameters.nRad,
                                                                                   TN.domainParameters.maxLevel ) );
   }

   //########################################################################################//
   // THIS IS AN ATTEMPT TO IMPROVE THE COUPLING
   // THE ENERGY AND STOKES ARE SOLVED AGAIN!
   //
   // IN THE FUTURE THIS COULD BE DONE IN A LOOP
   // TO ENSURE EVEN TIGHTER COUPLING
   //########################################################################################//

   bool predictorCorrector = TN.simulationParameters.predictorCorrector;

   if ( predictorCorrector )
   {
      //######################################################//
      //                ENERGY EQUATION AGAIN                 //
      //######################################################//

      /*############ ADVECTION STEP ############*/

      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( "------- Advection Step ------" );
      WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
      {
         temperatureMMOCOperator_->step( *temperatureFE,
                                         p2p1StokesFunctionContainer["VelocityFE"]->uvw(),
                                         p2p1StokesFunctionContainer["VelocityFEPrev"]->uvw(),
                                         TN.domainParameters.maxLevel,
                                         hyteg::Inner,
                                         TN.simulationParameters.dt,
                                         1u,
                                         true,
                                         true,
                                         false );
      }
      else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
      {
         P2toP1Conversion( p2p1StokesFunctionContainer["VelocityFEPrev"]->uvw(),
                           *( p1VectorFunctionContainer["VelocityFEPrevP1"] ),
                           TN.domainParameters.maxLevel + 1,
                           All );

         P2toP1Conversion( p2p1StokesFunctionContainer["VelocityFE"]->uvw(),
                           *( p1VectorFunctionContainer["VelocityFEP1"] ),
                           TN.domainParameters.maxLevel + 1,
                           All );

         P2toP1Conversion(
             *( p2ScalarFunctionContainer["TemperatureFE"] ), *temperatureFE, TN.domainParameters.maxLevel + 1, All );

         temperatureMMOCOperator_->step( *temperatureFE,
                                         *( p1VectorFunctionContainer["VelocityFEP1"] ),
                                         *( p1VectorFunctionContainer["VelocityFEPrevP1"] ),
                                         TN.domainParameters.maxLevel + 1,
                                         hyteg::Inner,
                                         TN.simulationParameters.dt,
                                         1u,
                                         true,
                                         false,
                                         false );

         P1toP2Conversion( *( p1ScalarFunctionContainer["TemperatureFEP1"] ),
                           *( p2ScalarFunctionContainer["TemperatureFE"] ),
                           TN.domainParameters.maxLevel,
                           All );
      }
      else
      {
         WALBERLA_ABORT( "Unsupported function type" );
      }

      // p2ScalarFunctionContainer["TemperaturePrev"]->assign(
      //     { real_c( 1 ) }, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, TN.domainParameters.maxLevel, All );

      for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
      {
         std::function< real_t( const Point3D& ) > initTemperature;
         switch ( TN.initialisationParameters.initialTemperatureDeviationMethod )
         {
         case INITIAL_TEMPERATURE_DEVIATION_METHOD::SINGLE_SPH:
            initTemperature = temperatureSingleSPH( *temperatureInitParams, *temperatureReferenceFct );
            break;
         case INITIAL_TEMPERATURE_DEVIATION_METHOD::RANDOM_SUPERPOSITION_SPH:
            initTemperature = temperatureRandomSuperpositioneSPH( *temperatureInitParams, *temperatureReferenceFct );
            break;
         case INITIAL_TEMPERATURE_DEVIATION_METHOD::WHITE_NOISE:
            initTemperature = temperatureWhiteNoise( *temperatureInitParams, *temperatureReferenceFct );
            break;
         default:
            WALBERLA_ABORT( "Unknown initial temperature deviation method" );
         }

         p2ScalarFunctionContainer["TemperatureFE"]->interpolate( initTemperature, l, DirichletBoundary );
      }
      /*############ ENERGY STEP ############*/

      // setupEnergyRHS();
      storage->getTimingTree()->start( "TerraNeo solve Energy equation" );
      solveEnergy();
      storage->getTimingTree()->stop( "TerraNeo solve Energy equation" );

      //######################################################//
      //               STOKES EQUATIONS AGAIN                 //
      //######################################################//

      if ( TN.simulationParameters.simulationType == "CirculationModel" )
      {
         //update plate velocity boundary condition prior to stokes solve
         //update if >=1Myr has passed since Prev update
         if ( ( ( TN.simulationParameters.plateAge - TN.simulationParameters.ageMa ) >= real_c( 1 ) ) )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "" );
            WALBERLA_LOG_INFO_ON_ROOT( "Update plates" );

            //save the age of the current update (rounded to 1Myr intervals)
            TN.simulationParameters.plateAge = std::round( TN.simulationParameters.ageMa );
            updatePlateVelocities( *( p2p1StokesFunctionContainer["VelocityFE"] ) );
         }

         WALBERLA_LOG_INFO_ON_ROOT( "Plate age: " << TN.simulationParameters.plateAge << " Ma" )
      }

      //update ref temp vector based on new temperature field

      temperatureProfiles =
          std::make_shared< RadialProfile >( computeRadialProfile( *( p2ScalarFunctionContainer["TemperatureFE"] ),
                                                                   TN.domainParameters.rMin,
                                                                   TN.domainParameters.rMax,
                                                                   TN.domainParameters.nRad,
                                                                   TN.domainParameters.maxLevel ) );

      TN.physicalParameters.temperatureProfile = temperatureProfiles->mean;

      // Consistency check for unreasonable low min Temperatures of Tmin <= 0 K
      for ( uint_t i = 0; i < temperatureProfiles->min.size(); i++ )
      {
         // Redimensionalise temperature
         real_t dimFactor = TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp;
         if ( ( temperatureProfiles->min[i] ) <= 0 )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Negative Temperature: " << temperatureProfiles->min[i] * dimFactor
                                                                << " detected at shell radii: "
                                                                << temperatureProfiles->shellRadii[i] );
            WALBERLA_LOG_INFO_ON_ROOT( "Dump data" );
            dataOutput();
            WALBERLA_ABORT( "Aborting simulation run" );
         }
      }

      solveStokes();

      // update viscosity Profiles for logging
      if ( TN.simulationParameters.tempDependentViscosity )
      {
         viscosityProfiles =
             std::make_shared< RadialProfile >( computeRadialProfile( *( p2ScalarFunctionContainer["ViscosityFE"] ),
                                                                      TN.domainParameters.rMin,
                                                                      TN.domainParameters.rMax,
                                                                      TN.domainParameters.nRad,
                                                                      TN.domainParameters.maxLevel ) );
      }
   }

   p2ScalarFunctionContainer["TemperaturePrev"]->assign(
       { real_c( 1 ) }, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, TN.domainParameters.maxLevel, All );

   //######################################################//
   //                  DUMP OUTPUT                         //
   //######################################################//

   if ( TN.outputParameters.outputMyr )
   {
      if ( ( ( TN.simulationParameters.modelRunTimeMa - TN.outputParameters.prevOutputTime ) >=
             real_c( TN.outputParameters.outputIntervalMyr ) ) )
      {
         dataOutput();
      }
   }
   else
   {
      dataOutput();
   }

   // Individually decide when checkpoint data is dumped
   if ( TN.outputParameters.ADIOS2StoreCheckpoint )
   {
      outputCheckpoint();
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Finished step: " << TN.simulationParameters.timeStep );

   if ( TN.outputParameters.createTimingDB )
   {
      db->setVariableEntry( "Step_size_Ma", stepSize );
      db->setVariableEntry( "model_runtime_Ma", TN.simulationParameters.modelRunTimeMa );
      db->setVariableEntry( "timestep", TN.simulationParameters.timeStep );
      db->setVariableEntry( "max_magnitude_velocity_cm_a", MaxVelocityMagSI );
      db->setVariableEntry( "Time_solve_Energy", timerSolveEnergy );
      db->writeRowOnRoot();
   }

   ++TN.simulationParameters.timeStep;
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::solveEnergy()
{
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Energy Solve ---------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   const uint_t temperatureMaxLevel = [&] {
      if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
      {
         return TN.domainParameters.maxLevel;
      }
      else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
      {
         return TN.domainParameters.maxLevel + 1;
      }
      else
      {
         WALBERLA_ABORT( "Unknown type" );
      }
   }();

   const std::shared_ptr< TemperatureFunction_T >& temperatureFE = [&] {
      if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
      {
         return p2ScalarFunctionContainer.at( "TemperatureFE" );
      }
      else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
      {
         return p1ScalarFunctionContainer.at( "TemperatureFEP1" );
      }
      else
      {
         WALBERLA_ABORT( "Unknown type" );
      }
   }();

   const std::shared_ptr< TemperatureFunction_T >& temperatureRhs = [&] {
      if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
      {
         return p2ScalarFunctionContainer.at( "EnergyRHSWeak" );
      }
      else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
      {
         return p1ScalarFunctionContainer.at( "EnergyRHSP1" );
      }
      else
      {
         WALBERLA_ABORT( "Unknown type" );
      }
   }();

   temperatureTransportOperator_->setTimestep( TN.simulationParameters.dt );

   p2ScalarFunctionContainer["ShearHeatingTermCoeff"]->interpolate(
       shearHeatingCoeffCalc, { *( p2ScalarFunctionContainer["DensityFE"] ) }, TN.domainParameters.maxLevel, All );

   P2toP1Conversion( *( p2ScalarFunctionContainer["ShearHeatingTermCoeff"] ),
                     *( p1ScalarFunctionContainer["ShearHeatingTermCoeffP1"] ),
                     TN.domainParameters.maxLevel + 1,
                     All );

   // Assemble RHS
   temperatureTransportOperator_->applyRHS( *temperatureRhs, temperatureMaxLevel, All );

   // Solve
   temperatureTransportSolver_->solve( *temperatureTransportOperator_, *temperatureFE, *temperatureRhs, temperatureMaxLevel );

   if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
   {
      P1toP2Conversion( *temperatureFE, *( p2ScalarFunctionContainer["TemperatureFE"] ), TN.domainParameters.maxLevel, All );
   }

   if ( !TN.simulationParameters.compressible && TN.simulationParameters.volAvrgTemperatureDev )
   {
      for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; ++l )
      {
         p2ScalarMassOperator_->apply(
             *( p2ScalarFunctionContainer["TemperatureFE"] ), *( p2ScalarFunctionContainer["TemperatureVolumetric"] ), l, All );
         p2ScalarFunctionContainer["Ones"]->interpolate( real_c( 1 ), l, All );
         p2ScalarMassOperator_->apply( *( p2ScalarFunctionContainer["Ones"] ), *( p2ScalarFunctionContainer["Volume"] ), l, All );
      }
      real_t TemperatureVolumentric =
          p2ScalarFunctionContainer["TemperatureVolumetric"]->sumGlobal( TN.domainParameters.maxLevel );
      real_t Volume                              = p2ScalarFunctionContainer["Volume"]->sumGlobal( TN.domainParameters.maxLevel );
      TN.simulationParameters.avrgTemperatureVol = ( TemperatureVolumentric / Volume );
      WALBERLA_LOG_INFO_ON_ROOT(
          "Average volumentric temperature: " << TN.simulationParameters.avrgTemperatureVol *
                                                     ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) );
   }

   real_t energyResidual = calculateEnergyResidual( TN.domainParameters.maxLevel );

   if ( std::isnan( energyResidual ) )
   {
      WALBERLA_ABORT( "NaN values detected in temperatures after Energy solve" );
   }

   temperatureTransportOperator_->incrementTimestep();
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::setupStokesRHS()
{
   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      ////////////////////
      //    Momentum    //
      ////////////////////
      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > calculateTDev =
          [this]( const Point3D& x, const std::vector< real_t >& vals ) {
             real_t refTemp = referenceTemperatureFct( x );
             return vals[0] - refTemp;
          };

      p2ScalarFunctionContainer["TemperatureDev"]->interpolate(
          calculateTDev, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, l, All );

      p2ScalarMassOperator_->apply(
          *( p2ScalarFunctionContainer["TemperatureDev"] ), p2p1StokesFunctionContainer["StokesRHS"]->uvw()[0], l, All );
      p2ScalarMassOperator_->apply(
          *( p2ScalarFunctionContainer["TemperatureDev"] ), p2p1StokesFunctionContainer["StokesRHS"]->uvw()[1], l, All );
      p2ScalarMassOperator_->apply(
          *( p2ScalarFunctionContainer["TemperatureDev"] ), p2p1StokesFunctionContainer["StokesRHS"]->uvw()[2], l, All );

      // Multiply current RHS with rho and non-dimensionalised numbers
      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > momentumFactors =
          [&]( const Point3D& x, const std::vector< real_t >& deltaT ) {
             if ( TN.simulationParameters.haveThermalExpProfile )
             {
                TN.physicalParameters.thermalExpansivityRadial =
                    terraneo::interpolateDataValues( x,
                                                     TN.physicalParameters.radiusAlpha,
                                                     TN.physicalParameters.thermalExpansivityProfile,
                                                     TN.domainParameters.rMin,
                                                     TN.domainParameters.rMax );
                return ( -( TN.physicalParameters.rayleighNumber / TN.physicalParameters.pecletNumber ) *
                         ( TN.physicalParameters.thermalExpansivityRadial / TN.physicalParameters.thermalExpansivity ) *
                         densityFunction( x ) * deltaT[0] );
             }
             else
             {
                return ( -( TN.physicalParameters.rayleighNumber / TN.physicalParameters.pecletNumber ) * densityFunction( x ) *
                         deltaT[0] );
             }
          };

      // Interpolate functions to RHS
      p2p1StokesFunctionContainer["StokesRHS"]->uvw()[0].interpolate(
          momentumFactors, { p2p1StokesFunctionContainer["StokesRHS"]->uvw()[0] }, l, All );
      p2p1StokesFunctionContainer["StokesRHS"]->uvw()[1].interpolate(
          momentumFactors, { p2p1StokesFunctionContainer["StokesRHS"]->uvw()[1] }, l, All );
      p2p1StokesFunctionContainer["StokesRHS"]->uvw()[2].interpolate(
          momentumFactors, { p2p1StokesFunctionContainer["StokesRHS"]->uvw()[2] }, l, All );

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multiplyWithInwardNormalX =
          []( const Point3D& x, const std::vector< real_t >& vals ) {
             real_t xNorm = x[0] / x.norm();
             return -xNorm * vals[0];
          };

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multiplyWithInwardNormalY =
          []( const Point3D& x, const std::vector< real_t >& vals ) {
             real_t xNorm = x[1] / x.norm();
             return -xNorm * vals[0];
          };

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multiplyWithInwardNormalZ =
          []( const Point3D& x, const std::vector< real_t >& vals ) {
             real_t xNorm = x[2] / x.norm();
             return -xNorm * vals[0];
          };

      p2p1StokesFunctionContainer["StokesTmp1"]->uvw().assign(
          { 1.0 }, { p2p1StokesFunctionContainer["StokesRHS"]->uvw() }, l, All );

      // multply with inward normal (for gravity)
      p2p1StokesFunctionContainer["StokesRHS"]->uvw().component( 0u ).interpolate(
          multiplyWithInwardNormalX, { p2p1StokesFunctionContainer["StokesTmp1"]->uvw().component( 0u ) }, l, All );

      p2p1StokesFunctionContainer["StokesRHS"]->uvw().component( 1u ).interpolate(
          multiplyWithInwardNormalY, { p2p1StokesFunctionContainer["StokesTmp1"]->uvw().component( 1u ) }, l, All );

      p2p1StokesFunctionContainer["StokesRHS"]->uvw().component( 2u ).interpolate(
          multiplyWithInwardNormalZ, { p2p1StokesFunctionContainer["StokesTmp1"]->uvw().component( 2u ) }, l, All );

      /////////////////
      //    Mass    //
      ////////////////

      // Provide the option to run incompressible simulations for test or educational purposes
      if ( TN.simulationParameters.compressible && TN.simulationParameters.frozenVelocity )
      {
         // Update non-dimensional numbers for mass conservation equation
         if ( TN.simulationParameters.haveThermalExpProfile || TN.simulationParameters.haveSpecificHeatCapProfile )
         {
            // Update gradRho/Rho with new non-Dim paramters Di and alpha.
            // grad(rho)/rho = - ( Di / gamma ) * r_hat
            // std::function< real_t( const Point3D& ) > updateDensity = [&]( const Point3D& x ) { return densityFunc( x ); };
            p2ScalarFunctionContainer["DensityFE"]->interpolate( densityFunc, l, All );
         }

         frozenVelocityRHS_->apply(
             p2p1StokesFunctionContainer["VelocityFE"]->uvw(), p2p1StokesFunctionContainer["StokesRHS"]->p(), l, All );
      }
      else
      {
         p2p1StokesFunctionContainer["StokesRHS"]->p().interpolate( real_c( 0 ), l, All );
      }
   }
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::solveStokes()
{
   walberla::WcTimer localTimer;
   if ( TN.simulationParameters.tempDependentViscosity )
   {
      updateViscosity();
      stokesOperator_->getA().computeInverseDiagonalOperatorValues();
      stokesOperatorRotationOpgen_->getA().computeInverseDiagonalOperatorValues();
   }

   localTimer.start();
   setupStokesRHS();
   localTimer.end();
   real_t setupStokesTimer = localTimer.last();

   if ( TN.simulationParameters.timeStep == 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "----------------------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( "------ Initial Stokes solve ------" );
      WALBERLA_LOG_INFO_ON_ROOT( "----------------------------------" );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "--------------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( "------ Stokes solve ------" );
      WALBERLA_LOG_INFO_ON_ROOT( "--------------------------" );
   }

   real_t stokesResidual   = calculateStokesResidual( TN.domainParameters.maxLevel );
   real_t pressureResidual = calculatePressureResidual( TN.domainParameters.maxLevel );

   real_t stokesResidualInitial   = stokesResidual;
   real_t pressureResidualInitial = pressureResidual;

   TN.solverParameters.vCycleResidualUPrev       = stokesResidual;
   TN.solverParameters.initialResidualU          = stokesResidual;
   TN.solverParameters.averageResidualReductionU = real_c( 0 );
   TN.solverParameters.numVCycles                = 0;

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes residual (initial): " << stokesResidual );
   WALBERLA_LOG_INFO_ON_ROOT( "Pressure residual (initial): " << pressureResidual );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   if ( TN.outputParameters.createTimingDB )
   {
      db->setVariableEntry( "Stokes_residual_initial", stokesResidual );
   }

   if ( stokesResidual > TN.solverParameters.stokesKillTolerance && ( TN.simulationParameters.timeStep != 0 ) )
   {
      WALBERLA_ABORT( "Residual " << stokesResidual << " exceeds tolerance of " << TN.solverParameters.stokesKillTolerance
                                  << ". ABORTING SIMULATION. " );
   }

   localTimer.start();
   storage->getTimingTree()->start( "TerraNeo solve Stokes" );

   if ( TN.solverParameters.useRotationWrapper )
   {
      rotationOperator_->rotate( *( p2p1StokesFunctionContainer["VelocityFE"] ), TN.domainParameters.maxLevel, FreeslipBoundary );
      p2p1StokesFunctionContainer["VelocityFERotated"]->assign(
          { 1.0 }, { *( p2p1StokesFunctionContainer["VelocityFE"] ) }, TN.domainParameters.maxLevel, All );

      rotationOperator_->rotate( *( p2p1StokesFunctionContainer["StokesRHS"] ), TN.domainParameters.maxLevel, FreeslipBoundary );
      p2p1StokesFunctionContainer["StokesRHSRotated"]->assign(
          { 1.0 }, { *( p2p1StokesFunctionContainer["StokesRHS"] ) }, TN.domainParameters.maxLevel, All );

      stokesRotationOpgenSolver_->solve( *stokesOperatorRotationOpgen_,
                                         *( p2p1StokesFunctionContainer["VelocityFERotated"] ),
                                         *( p2p1StokesFunctionContainer["StokesRHSRotated"] ),
                                         TN.domainParameters.maxLevel );

      p2p1StokesFunctionContainer["VelocityFE"]->assign(
          { 1.0 }, { *p2p1StokesFunctionContainer["VelocityFERotated"] }, TN.domainParameters.maxLevel, All );
      rotationOperator_->rotate(
          *( p2p1StokesFunctionContainer["VelocityFE"] ), TN.domainParameters.maxLevel, FreeslipBoundary, true );
   }
   else
   {
      projectionOperator_->project(
          *( p2p1StokesFunctionContainer["StokesRHS"] ), TN.domainParameters.maxLevel, FreeslipBoundary );
      stokesSolver_->solve( *stokesOperator_,
                            *( p2p1StokesFunctionContainer["VelocityFE"] ),
                            *( p2p1StokesFunctionContainer["StokesRHS"] ),
                            TN.domainParameters.maxLevel );
   }

   storage->getTimingTree()->stop( "TerraNeo solve Stokes" );
   localTimer.end();

   real_t timeStokesSolve = localTimer.last();
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes time [s]: " << timeStokesSolve );

   if ( TN.solverParameters.useRotationWrapper )
   {
      setupStokesRHS();
   }

   stokesResidual   = calculateStokesResidual( TN.domainParameters.maxLevel );
   pressureResidual = calculatePressureResidual( TN.domainParameters.maxLevel );

   if ( std::isnan( stokesResidual ) )
   {
      WALBERLA_ABORT( "NaN values detected in velocity after Stokes solve" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes residual (final): " << stokesResidual );
   WALBERLA_LOG_INFO_ON_ROOT( "Pressure residual (final): " << pressureResidual );

   real_t relStokesResidual   = stokesResidual / stokesResidualInitial;
   real_t relPressureResidual = pressureResidual / pressureResidualInitial;
   WALBERLA_LOG_INFO_ON_ROOT( "Relative stokes residual (final): " << relStokesResidual );
   WALBERLA_LOG_INFO_ON_ROOT( "Relative pressure residual (final): " << relPressureResidual );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   if ( TN.outputParameters.createTimingDB )
   {
      db->setVariableEntry( "Relative_Stokes_residual", relStokesResidual );
      db->setVariableEntry( "Relative_pressure_residual", relPressureResidual );
      db->setVariableEntry( "Time_solve_Stokes", timeStokesSolve );
   }
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::calculateStokesResidual( uint_t level )
{
   stokesOperator_->apply( *( p2p1StokesFunctionContainer["VelocityFE"] ),
                           *( p2p1StokesFunctionContainer["StokesTmp1"] ),
                           level,
                           Inner | NeumannBoundary | FreeslipBoundary );
   p2p1StokesFunctionContainer["StokesTmp1"]->assign(
       { real_c( 1 ), real_c( -1 ) },
       { *( p2p1StokesFunctionContainer["StokesTmp1"] ), *( p2p1StokesFunctionContainer["StokesRHS"] ) },
       level,
       Inner | NeumannBoundary | FreeslipBoundary );

   p2p1StokesFunctionContainer["StokesResidual"]->assign(
       { real_c( 1 ) }, { *( p2p1StokesFunctionContainer["StokesTmp1"] ) }, level, All );

   return std::sqrt( p2p1StokesFunctionContainer["StokesTmp1"]->dotGlobal(
       *( p2p1StokesFunctionContainer["StokesTmp1"] ), level, Inner | NeumannBoundary | FreeslipBoundary ) );
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::calculatePressureResidual( uint_t level )
{
   stokesOperator_->apply( *( p2p1StokesFunctionContainer["VelocityFE"] ),
                           *( p2p1StokesFunctionContainer["StokesTmp3"] ),
                           level,
                           Inner | NeumannBoundary | FreeslipBoundary );
   p2p1StokesFunctionContainer["StokesTmp3"]->p().assign(
       { real_c( 1 ), real_c( -1 ) },
       { ( p2p1StokesFunctionContainer["StokesTmp3"]->p() ), ( p2p1StokesFunctionContainer["StokesRHS"]->p() ) },
       level,
       Inner | NeumannBoundary | FreeslipBoundary );
   return ( std::sqrt( p2p1StokesFunctionContainer["StokesTmp3"]->p().dotGlobal(
       ( p2p1StokesFunctionContainer["StokesTmp3"]->p() ), level, Inner | NeumannBoundary | FreeslipBoundary ) ) );
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::calculateEnergyResidual( uint_t level )
{
   return 0.0;
   // P2Function< real_t >& tempFunc = p2p1StokesFunctionContainer["StokesTmp1"]->uvw().component( 0u );

   // temperatureTransportOperator_->apply(
   //     *( p2ScalarFunctionContainer["TemperatureFE"] ), tempFunc, level, Inner | NeumannBoundary | FreeslipBoundary );
   // tempFunc.assign( { real_c( 1 ), real_c( -1 ) },
   //                  { tempFunc, *( p2ScalarFunctionContainer["EnergyRHSWeak"] ) },
   //                  level,
   //                  Inner | NeumannBoundary | FreeslipBoundary );
   // return std::sqrt( tempFunc.dotGlobal( tempFunc, level, Inner | NeumannBoundary | FreeslipBoundary ) );
}
} // namespace terraneo