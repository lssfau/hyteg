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

   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureFE     = p2p1StokesFunctionContainer.at( "VelocityFE" );
   std::shared_ptr< P2P1StokesFunction_T >& velocityPressurePrevFE = p2p1StokesFunctionContainer.at( "VelocityFEPrev" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesTmp1             = p2p1StokesFunctionContainer.at( "StokesTmp1" );

   std::shared_ptr< P2ScalarFunction_T >& temperatureP2            = p2ScalarFunctionContainer.at( "TemperatureFE" );
   std::shared_ptr< P2ScalarFunction_T >& temperaturePrevP2        = p2ScalarFunctionContainer.at( "TemperaturePrev" );
   std::shared_ptr< P2ScalarFunction_T >& velocityMagnitudeSquared = p2ScalarFunctionContainer.at( "VelocityMagnitudeSquared" );

   std::shared_ptr< P1VectorFunction_T >& velocityPressureFEP1     = p1VectorFunctionContainer.at( "VelocityFEP1" );
   std::shared_ptr< P1VectorFunction_T >& velocityPressurePrevFEP1 = p1VectorFunctionContainer.at( "VelocityFEPrevP1" );

   std::shared_ptr< P1ScalarFunction_T >& temperatureP1 = p1ScalarFunctionContainer.at( "TemperatureFEP1" );

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
            checkpointImporter->restoreFunction(
                *( temperatureP2 ), TN.domainParameters.minLevel, checkpointDataInfo[0].maxLevel, 0, true );

            WALBERLA_LOG_INFO_ON_ROOT( "Prolongate Temperature field checkpoint data from level: "
                                       << checkpointDataInfo[0].maxLevel
                                       << " to new maxLevel: " << TN.domainParameters.maxLevel );
            for ( uint_t l = checkpointDataInfo[0].maxLevel; l < TN.domainParameters.maxLevel; l++ )
            {
               p2ProlongationOperator_->prolongate( *( temperatureP2 ), l, All );
               WALBERLA_LOG_INFO_ON_ROOT( "Temperatuer field prolongated from level:  " << l << " to level: " << l + 1 );
            }
         }
         else
         {
            WALBERLA_LOG_INFO_ON_ROOT(
                "Restore Temperature field from checkpoint data at level: " << checkpointDataInfo[0].maxLevel );
            checkpointImporter->restoreFunction(
                *( temperatureP2 ), TN.domainParameters.minLevel, TN.domainParameters.maxLevel, 0, true );
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
      velocityPressurePrevFE->assign( { real_c( 1 ) }, { *( velocityPressureFE ) }, TN.domainParameters.maxLevel, All );
   } //end timestep0 stokes

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "-------- Time step: " << TN.simulationParameters.timeStep << " --------" );

   real_t vMax = velocityMaxMagnitude( velocityPressureFE->uvw(),
                                       stokesTmp1->uvw().component( 0u ),
                                       *( velocityMagnitudeSquared ),
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
      if ( TN.simulationParameters.timeStep < TN.simulationParameters.initialNStepsForTimestepLinearIncrease )
      {
         // Initial Step size is 1000 a
         TN.simulationParameters.dt = ( TN.simulationParameters.timeStep * TN.simulationParameters.initialTimestepSize *
                                        TN.physicalParameters.characteristicVelocity *
                                        TN.simulationParameters.plateVelocityScaling * TN.simulationParameters.secondsPerMyr ) /
                                      TN.physicalParameters.mantleThickness;
      }
      else
      {
         TN.simulationParameters.dt = ( TN.simulationParameters.cflMax / vMax ) * TN.simulationParameters.hMin;
      }
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

   temperatureP2->assign( { real_c( 1 ) }, { *( temperaturePrevP2 ) }, TN.domainParameters.maxLevel, All );

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
      temperatureMMOCOperator_->step( *temperatureP2,
                                      velocityPressureFE->uvw(),
                                      velocityPressurePrevFE->uvw(),
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
      P2toP1Conversion( velocityPressurePrevFE->uvw(), *( velocityPressurePrevFEP1 ), TN.domainParameters.maxLevel + 1, All );
      P2toP1Conversion( velocityPressureFE->uvw(), *( velocityPressureFEP1 ), TN.domainParameters.maxLevel + 1, All );
      P2toP1Conversion( *( temperatureP2 ), *temperatureP1, TN.domainParameters.maxLevel + 1, All );

      temperatureMMOCOperator_->step( *temperatureP1,
                                      *( velocityPressureFEP1 ),
                                      *( velocityPressurePrevFEP1 ),
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

   // for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   // {
   temperatureP2->interpolate( temperatureInitialCondition, TN.domainParameters.maxLevel, DirichletBoundary );
   temperatureP1->interpolate( temperatureInitialCondition, TN.domainParameters.maxLevel + 1, DirichletBoundary );
   // }

   /*############ ENERGY STEP ############*/

   localTimerStep.start();
   storage->getTimingTree()->start( "TerraNeo solve Energy equation" );
   solveEnergy();
   storage->getTimingTree()->stop( "TerraNeo solve Energy equation" );
   localTimerStep.end();
   real_t timerSolveEnergy = localTimerStep.last();

   if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
   {
      P1toP2Conversion( *temperatureP1, *temperatureP2, TN.domainParameters.maxLevel, All );
   }

   if ( !TN.simulationParameters.compressible && TN.simulationParameters.volAvrgTemperatureDev )
   {
      P2Function< real_t >& tempP2 = stokesTmp1->uvw().component( 0u );

      for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; ++l )
      {
         p2ScalarMassOperator_->apply( *temperatureP2, tempP2, l, All );
      }
      real_t temperatureVolumetric = tempP2.sumGlobal( TN.domainParameters.maxLevel );

      real_t Volume = real_c( 4.0 / 3.0 ) * walberla::math::pi *
                      ( std::pow( TN.domainParameters.rMax, 3u ) - std::pow( TN.domainParameters.rMin, 3u ) );

      TN.simulationParameters.avrgTemperatureVol = ( temperatureVolumetric / Volume );

      WALBERLA_LOG_INFO_ON_ROOT(
          "Average volumentric temperature: " << TN.simulationParameters.avrgTemperatureVol *
                                                     ( TN.physicalParameters.cmbTemp - TN.physicalParameters.surfaceTemp ) );
   }

   //######################################################//
   //                  STOKES EQUATIONS                    //
   //######################################################//

   // update velocity field storing the velocity field of the Prev timestep
   velocityPressurePrevFE->assign( { real_c( 1 ) }, { *( velocityPressureFE ) }, TN.domainParameters.maxLevel, All );

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
         updatePlateVelocities( *( velocityPressureFE ) );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "Plate age: " << TN.simulationParameters.plateAge << " Ma" )
   }

   temperatureProfiles = std::make_shared< RadialProfile >( computeRadialProfile( *( temperatureP2 ),
                                                                                  TN.domainParameters.rMin,
                                                                                  TN.domainParameters.rMax,
                                                                                  TN.domainParameters.nRad,
                                                                                  TN.domainParameters.maxLevel ) );

   TN.physicalParameters.temperatureProfile = temperatureProfiles->mean;

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
   velocityProfiles                      = std::make_shared< RadialProfile >( computeRadialProfile( velocityPressureFE->uvw(),
                                                                               TN.domainParameters.rMin,
                                                                               TN.domainParameters.rMax,
                                                                               TN.domainParameters.nRad,
                                                                               TN.domainParameters.maxLevel ) );
   TN.physicalParameters.velocityProfile = velocityProfiles->rms;

   // update viscosity Profiles for logging
   if ( TN.simulationParameters.tempDependentViscosity )
   {
      viscosityProfiles = std::make_shared< RadialProfile >( computeRadialProfile( *( velocityPressureFE ),
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
         temperatureMMOCOperator_->step( *temperatureP2,
                                         velocityPressureFE->uvw(),
                                         velocityPressurePrevFE->uvw(),
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
         P2toP1Conversion( velocityPressurePrevFE->uvw(), *( velocityPressurePrevFEP1 ), TN.domainParameters.maxLevel + 1, All );
         P2toP1Conversion( velocityPressureFE->uvw(), *( velocityPressureFEP1 ), TN.domainParameters.maxLevel + 1, All );
         P2toP1Conversion( *( temperatureP2 ), *temperatureP1, TN.domainParameters.maxLevel + 1, All );

         temperatureMMOCOperator_->step( *temperatureP1,
                                         *( velocityPressureFEP1 ),
                                         *( velocityPressurePrevFEP1 ),
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

      // for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
      // {
      temperatureP2->interpolate( temperatureInitialCondition, TN.domainParameters.maxLevel, DirichletBoundary );
      temperatureP1->interpolate( temperatureInitialCondition, TN.domainParameters.maxLevel + 1, DirichletBoundary );
      // }

      /*############ ENERGY STEP ############*/

      // setupEnergyRHS();
      storage->getTimingTree()->start( "TerraNeo solve Energy equation" );
      solveEnergy();
      storage->getTimingTree()->stop( "TerraNeo solve Energy equation" );

      if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
      {
         P1toP2Conversion( *temperatureP1, *temperatureP2, TN.domainParameters.maxLevel, All );
      }

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
            updatePlateVelocities( *( velocityPressureFE ) );
         }

         WALBERLA_LOG_INFO_ON_ROOT( "Plate age: " << TN.simulationParameters.plateAge << " Ma" )
      }

      //update ref temp vector based on new temperature field

      temperatureProfiles = std::make_shared< RadialProfile >( computeRadialProfile( *( temperatureP2 ),
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
         viscosityProfiles = std::make_shared< RadialProfile >( computeRadialProfile( *( velocityPressureFE ),
                                                                                      TN.domainParameters.rMin,
                                                                                      TN.domainParameters.rMax,
                                                                                      TN.domainParameters.nRad,
                                                                                      TN.domainParameters.maxLevel ) );
      }
   }

   temperaturePrevP2->assign( { real_c( 1 ) }, { *( temperatureP2 ) }, TN.domainParameters.maxLevel, All );

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

   std::shared_ptr< P2ScalarFunction_T >& densityFE = p2ScalarFunctionContainer.at( "DensityFE" );

   std::shared_ptr< P2ScalarFunction_T >& shearHeatingCoeffP2 = p2ScalarFunctionContainer["ShearHeatingTermCoeff"];
   std::shared_ptr< P1ScalarFunction_T >& shearHeatingCoeffP1 = p1ScalarFunctionContainer["ShearHeatingTermCoeffP1"];

   temperatureTransportOperator_->setTimestep( TN.simulationParameters.dt );

   shearHeatingCoeffP2->interpolate( shearHeatingCoeffCalc, { *densityFE }, TN.domainParameters.maxLevel, All );

   P2toP1Conversion( *shearHeatingCoeffP2, *shearHeatingCoeffP1, TN.domainParameters.maxLevel + 1, All );

   // Assemble RHS
   temperatureTransportOperator_->applyRHS( *temperatureRhs, temperatureMaxLevel, All );

   // Solve
   temperatureTransportSolver_->solve( *temperatureTransportOperator_, *temperatureFE, *temperatureRhs, temperatureMaxLevel );

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
   // for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   // {

   uint_t l = TN.domainParameters.maxLevel;

   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureFE = p2p1StokesFunctionContainer.at( "VelocityFE" );

   std::shared_ptr< P2P1StokesFunction_T >& stokesRHS  = p2p1StokesFunctionContainer.at( "StokesRHS" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesTmp1 = p2p1StokesFunctionContainer.at( "StokesTmp1" );

   std::shared_ptr< P2ScalarFunction_T >& temperatureP2    = p2ScalarFunctionContainer.at( "TemperatureFE" );
   std::shared_ptr< P2ScalarFunction_T >& temperatureDevP2 = p2ScalarFunctionContainer.at( "TemperatureDev" );

   std::shared_ptr< P2ScalarFunction_T >& densityFE = p2ScalarFunctionContainer.at( "DensityFE" );

   ////////////////////
   //    Momentum    //
   ////////////////////

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > calculateTDev =
       [this]( const Point3D& x, const std::vector< real_t >& vals ) {
          if ( TN.simulationParameters.adaptiveRefTemp )
          {
             real_t radius = x.norm();
             uint_t shell  = nearestShellFromRadius( radius, temperatureProfiles->shellRadii );

             return vals[0] - temperatureProfiles->mean.at( shell );
          }
          else
          {
             real_t refTemp = referenceTemperatureFct( x );
             return vals[0] - refTemp;
          }
       };

   temperatureDevP2->interpolate( calculateTDev, { *temperatureP2 }, l, All );

   p2ScalarMassOperator_->apply( *temperatureDevP2, stokesRHS->uvw()[0], l, All );
   p2ScalarMassOperator_->apply( *temperatureDevP2, stokesRHS->uvw()[1], l, All );
   p2ScalarMassOperator_->apply( *temperatureDevP2, stokesRHS->uvw()[2], l, All );

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
   stokesRHS->uvw()[0].interpolate( momentumFactors, { stokesRHS->uvw()[0] }, l, All );
   stokesRHS->uvw()[1].interpolate( momentumFactors, { stokesRHS->uvw()[1] }, l, All );
   stokesRHS->uvw()[2].interpolate( momentumFactors, { stokesRHS->uvw()[2] }, l, All );

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

   stokesTmp1->uvw().assign( { 1.0 }, { stokesRHS->uvw() }, l, All );

   // multply with inward normal (for gravity)
   stokesRHS->uvw().component( 0u ).interpolate( multiplyWithInwardNormalX, { stokesTmp1->uvw().component( 0u ) }, l, All );
   stokesRHS->uvw().component( 1u ).interpolate( multiplyWithInwardNormalY, { stokesTmp1->uvw().component( 1u ) }, l, All );
   stokesRHS->uvw().component( 2u ).interpolate( multiplyWithInwardNormalZ, { stokesTmp1->uvw().component( 2u ) }, l, All );

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

      frozenVelocityRHS_->apply( velocityPressureFE->uvw(), stokesRHS->p(), l, All );
   }
   else
   {
      stokesRHS->p().interpolate( real_c( 0 ), l, All );
   }
   // }
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::solveStokes()
{
   walberla::WcTimer localTimer;

   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureFE = p2p1StokesFunctionContainer.at( "VelocityFE" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesRHS          = p2p1StokesFunctionContainer.at( "StokesRHS" );

   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureRotatedFE = p2p1StokesFunctionContainer.at( "VelocityFERotated" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesRHSRotated          = p2p1StokesFunctionContainer.at( "StokesRHSRotated" );

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
      rotationOperator_->rotate( *velocityPressureFE, TN.domainParameters.maxLevel, FreeslipBoundary );
      velocityPressureRotatedFE->assign( { 1.0 }, { *velocityPressureFE }, TN.domainParameters.maxLevel, All );

      rotationOperator_->rotate( *stokesRHS, TN.domainParameters.maxLevel, FreeslipBoundary );
      stokesRHSRotated->assign( { 1.0 }, { *stokesRHS }, TN.domainParameters.maxLevel, All );

      stokesRotationOpgenSolver_->solve(
          *stokesOperatorRotationOpgen_, *velocityPressureRotatedFE, *stokesRHSRotated, TN.domainParameters.maxLevel );

      velocityPressureFE->assign( { 1.0 }, { *velocityPressureRotatedFE }, TN.domainParameters.maxLevel, All );
      rotationOperator_->rotate( *( velocityPressureFE ), TN.domainParameters.maxLevel, FreeslipBoundary, true );
   }
   else
   {
      projectionOperator_->project( *stokesRHS, TN.domainParameters.maxLevel, FreeslipBoundary );
      stokesSolver_->solve( *stokesOperator_, *velocityPressureFE, *stokesRHS, TN.domainParameters.maxLevel );
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
   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureFE = p2p1StokesFunctionContainer.at( "VelocityFE" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesRHS          = p2p1StokesFunctionContainer.at( "StokesRHS" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesResidual     = p2p1StokesFunctionContainer.at( "StokesResidual" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesTmp1         = p2p1StokesFunctionContainer.at( "StokesTmp1" );

   stokesOperator_->apply( *velocityPressureFE, *stokesTmp1, level, Inner | NeumannBoundary | FreeslipBoundary );
   stokesTmp1->assign(
       { real_c( 1 ), real_c( -1 ) }, { *stokesTmp1, *stokesRHS }, level, Inner | NeumannBoundary | FreeslipBoundary );

   stokesResidual->assign( { real_c( 1 ) }, { *stokesTmp1 }, level, All );

   return std::sqrt( stokesTmp1->dotGlobal( *stokesTmp1, level, Inner | NeumannBoundary | FreeslipBoundary ) );
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::calculatePressureResidual( uint_t level )
{
   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureFE = p2p1StokesFunctionContainer.at( "VelocityFE" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesRHS          = p2p1StokesFunctionContainer.at( "StokesRHS" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesResidual     = p2p1StokesFunctionContainer.at( "StokesResidual" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesTmp1         = p2p1StokesFunctionContainer.at( "StokesTmp1" );

   stokesOperator_->apply( *velocityPressureFE, *stokesTmp1, level, Inner | NeumannBoundary | FreeslipBoundary );
   stokesTmp1->p().assign( { real_c( 1 ), real_c( -1 ) },
                           { ( stokesTmp1->p() ), ( stokesRHS->p() ) },
                           level,
                           Inner | NeumannBoundary | FreeslipBoundary );
   return ( std::sqrt( stokesTmp1->p().dotGlobal( ( stokesTmp1->p() ), level, Inner | NeumannBoundary | FreeslipBoundary ) ) );
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
real_t ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::calculateEnergyResidual( uint_t level )
{
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

   const TemperatureFunction_T& tempFunc = [&] {
      if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
      {
         return p2p1StokesFunctionContainer.at( "StokesTmp1" )->uvw().component( 0u );
      }
      else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
      {
         return p2p1StokesFunctionContainer.at( "StokesTmp1" )->uvw().component( 0u ).getVertexDoFFunction();
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

   temperatureTransportOperator_->apply( *temperatureFE, tempFunc, level, Inner | NeumannBoundary | FreeslipBoundary );
   tempFunc.assign(
       { real_c( 1 ), real_c( -1 ) }, { tempFunc, *temperatureRhs }, level, Inner | NeumannBoundary | FreeslipBoundary );
   return std::sqrt( tempFunc.dotGlobal( tempFunc, level, Inner | NeumannBoundary | FreeslipBoundary ) );
}
} // namespace terraneo