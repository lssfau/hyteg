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

void ConvectionSimulation::step()
{
   walberla::WcTimer localTimerStep;
   if ( TN.simulationParameters.timeStep == 0 )
   {
      //set up logging
      WALBERLA_ROOT_SECTION()
      {
         walberla::logging::Logging::instance()->includeLoggingToFile( TN.outputParameters.outputDirectory + "/" +
                                                                       TN.outputParameters.outputBaseName + ".out" );
      }
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
               p2ProlongationOperator->prolongate( *( p2ScalarFunctionContainer["TemperatureFE"] ), l, All );
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
      TN.simulationParameters.dt = TN.simulationParameters.dtConstant;
   }

   else
   {
      TN.simulationParameters.dt = ( TN.simulationParameters.cflMax / vMax ) * TN.simulationParameters.hMin;
   }

   TN.simulationParameters.modelTime += TN.simulationParameters.dt;
   real_t stepSize = ( TN.simulationParameters.dt * TN.physicalParameters.mantleThickness ) /
                     ( TN.physicalParameters.characteristicVelocity * TN.simulationParameters.plateVelocityScaling *
                       TN.simulationParameters.secondsPerMyr );
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

   //######################################################//
   //                  ENERGY EQUATION                     //
   //######################################################//

   /*############ ADVECTION STEP ############*/

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Advection Step ------" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   transportOperator->step( *( p2ScalarFunctionContainer["TemperatureFE"] ),
                            p2p1StokesFunctionContainer["VelocityFE"]->uvw(),
                            p2p1StokesFunctionContainer["VelocityFEPrev"]->uvw(),
                            TN.domainParameters.maxLevel,
                            hyteg::Inner,
                            TN.simulationParameters.dt,
                            1,
                            true,
                            true,
                            false );

   // Reset temperature on boundary to initial values

   p2ScalarFunctionContainer["TemperaturePrev"]->assign(
       { real_c( 1 ) }, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, TN.domainParameters.maxLevel, All );

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

   temperatureProfiles = std::make_shared< RadialProfile >( computeRadialProfile( *( p2ScalarFunctionContainer["TemperatureFE"] ),
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
      viscosityProfiles = std::make_shared< RadialProfile >( computeRadialProfile( *( p2ScalarFunctionContainer["ViscosityFE"] ),
                                                                                   TN.domainParameters.rMin,
                                                                                   TN.domainParameters.rMax,
                                                                                   TN.domainParameters.nRad,
                                                                                   TN.domainParameters.maxLevel ) );
   }
   //######################################################//
   //                  DUMP OUTPUT                         //
   //######################################################//

   if ( TN.outputParameters.outputMyr && ( ( TN.simulationParameters.modelRunTimeMa - TN.outputParameters.prevOutputTime ) >=
                                           real_c( TN.outputParameters.outputIntervalMyr ) ) )
   {
      dataOutput();
   }
   // Individually decide when checkpoint data is dumped
   if ( TN.outputParameters.ADIOS2StoreCheckpoint )
   {
      outputCheckpoint();
   }

   if ( !TN.outputParameters.outputMyr )
   {
      dataOutput();
      localTimerStep.end();
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

void ConvectionSimulation::solveEnergy()
{
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Energy Solve ---------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   transportOperatorTALA->setTimestep( TN.simulationParameters.dt );
   // transportOperatorRHS->setTimestep( TN.simulationParameters.dt );

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > shearHeatingCoeffCalc =
       [this]( const Point3D& x, const std::vector< real_t >& density ) {
          real_t radius = x.norm();
          real_t shearHeatingCoeff;

          if ( TN.simulationParameters.lithosphereShearHeatingScaling > 1 )
          {
             WALBERLA_ABORT( "Shear heating scaling factor at Lithosphere > 1 not allowed! --- Abort simulation ---" );
          }

          if ( TN.simulationParameters.haveSpecificHeatCapProfile && TN.simulationParameters.haveDensityProfile )
          {
             TN.physicalParameters.specificHeatCapacityRadial =
                 terraneo::interpolateDataValues( x,
                                                  TN.physicalParameters.radiusCp,
                                                  TN.physicalParameters.specificHeatCapacityProfile,
                                                  TN.domainParameters.rMin,
                                                  TN.domainParameters.rMax );

             shearHeatingCoeff = ( ( TN.physicalParameters.dissipationNumber * TN.physicalParameters.pecletNumber /
                                     TN.physicalParameters.rayleighNumber ) *
                                   ( TN.physicalParameters.specificHeatCapacity /
                                     ( TN.physicalParameters.specificHeatCapacityRadial * densityFunction( x ) ) ) );
          }
          else
          {
             shearHeatingCoeff = ( TN.physicalParameters.dissipationNumber * TN.physicalParameters.pecletNumber /
                                   ( TN.physicalParameters.rayleighNumber * density[0] ) );
          }

          if ( TN.simulationParameters.haveSpecificHeatCapProfile && TN.simulationParameters.haveDensityProfile )
          {
             TN.physicalParameters.specificHeatCapacityRadial =
                 terraneo::interpolateDataValues( x,
                                                  TN.physicalParameters.radiusCp,
                                                  TN.physicalParameters.specificHeatCapacityProfile,
                                                  TN.domainParameters.rMin,
                                                  TN.domainParameters.rMax );

             shearHeatingCoeff = ( ( TN.physicalParameters.dissipationNumber * TN.physicalParameters.pecletNumber /
                                     TN.physicalParameters.rayleighNumber ) *
                                   ( TN.physicalParameters.specificHeatCapacity /
                                     ( TN.physicalParameters.specificHeatCapacityRadial * densityFunction( x ) ) ) );
          }
          else
          {
             shearHeatingCoeff = ( TN.physicalParameters.dissipationNumber * TN.physicalParameters.pecletNumber /
                                   ( TN.physicalParameters.rayleighNumber * density[0] ) );
          }

          if ( radius > TN.domainParameters.rMax -
                            ( TN.simulationParameters.lithosphereThickness * 1000 / TN.physicalParameters.mantleThickness ) )
          {
             return ( shearHeatingCoeff * TN.simulationParameters.lithosphereShearHeatingScaling );
          }
          else
          {
             return shearHeatingCoeff;
          }
       };

   p2ScalarFunctionContainer["ShearHeatingTermCoeff"]->interpolate(
       shearHeatingCoeffCalc, { *( p2ScalarFunctionContainer["DensityFE"] ) }, TN.domainParameters.maxLevel, All );

   // Assemble RHS
   transportOperatorTALA->applyRHS( *( p2ScalarFunctionContainer["EnergyRHSWeak"] ), TN.domainParameters.maxLevel, All );

   // Solve
   transportSolverTALA->solve( *transportOperatorTALA,
                               *( p2ScalarFunctionContainer["TemperatureFE"] ),
                               *( p2ScalarFunctionContainer["EnergyRHSWeak"] ),
                               TN.domainParameters.maxLevel );

   real_t energyResidual = calculateEnergyResidual( TN.domainParameters.maxLevel );

   if ( std::isnan( energyResidual ) )
   {
      WALBERLA_ABORT( "NaN values detected in temperatures after Energy solve" );
   }

   transportOperatorTALA->incrementTimestep();
}

void ConvectionSimulation::setupStokesRHS()
{
   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      ////////////////////
      //    Momentum    //
      ////////////////////
      if ( TN.simulationParameters.adaptiveRefTemp )
      {
         std::function< real_t( const Point3D&, const std::vector< real_t >& ) > adaptiveTemperatureDev =
             [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
                auto radius = x.norm();

                size_t shell = static_cast< size_t >( std::round(
                    real_c( TN.simulationParameters.numLayers ) *
                    ( ( radius - TN.domainParameters.rMin ) / ( TN.domainParameters.rMax - TN.domainParameters.rMin ) ) ) );

                return ( Temperature[0] - temperatureProfiles->mean.at( shell ) );
             };
         p2ScalarFunctionContainer["TemperatureDev"]->interpolate(
             adaptiveTemperatureDev, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, l, All );
      }
      else
      {
         std::function< real_t( const Point3D&, const std::vector< real_t >& ) > calculateTDev =
             [this]( const Point3D& x, const std::vector< real_t >& vals ) {
                real_t refTemp = referenceTemperatureFct( x );
                return vals[0] - refTemp;
             };
         p2ScalarFunctionContainer["TemperatureDev"]->interpolate(
             calculateTDev, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, l, All );

         // p2ScalarFunctionContainer["TemperatureReference"]->interpolate( referenceTemperatureFct, l, All );
         // p2ScalarFunctionContainer["TemperatureDev"]->assign(
         //     { 1.0, -1.0 },
         //     { *( p2ScalarFunctionContainer["TemperatureFE"] ), *( p2ScalarFunctionContainer["TemperatureReference"] ) },
         //     l,
         //     All );
      }

      // Multiply with mass matrix (of velocity space -- P2) to get the weak form

      P2MassOperator->apply(
          *( p2ScalarFunctionContainer["TemperatureDev"] ), p2p1StokesFunctionContainer["StokesRHS"]->uvw()[0], l, All );
      P2MassOperator->apply(
          *( p2ScalarFunctionContainer["TemperatureDev"] ), p2p1StokesFunctionContainer["StokesRHS"]->uvw()[1], l, All );
      P2MassOperator->apply(
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
      if ( TN.simulationParameters.compressible )
      {
         // Update non-dimensional numbers for mass conservation equation
         if ( TN.simulationParameters.haveThermalExpProfile || TN.simulationParameters.haveSpecificHeatCapProfile )
         {
            // Update gradRho/Rho with new non-Dim paramters Di and alpha.
            // grad(rho)/rho = - ( Di / gamma ) * r_hat
            // std::function< real_t( const Point3D& ) > updateDensity = [&]( const Point3D& x ) { return densityFunc( x ); };
            p2ScalarFunctionContainer["DensityFE"]->interpolate( densityFunc, l, All );
         }

         frozenVelocityRHS->apply(
             p2p1StokesFunctionContainer["VelocityFE"]->uvw(), p2p1StokesFunctionContainer["StokesRHS"]->p(), l, All );
             
      }
      else
      {
         p2p1StokesFunctionContainer["StokesRHS"]->p().interpolate( real_c( 0 ), l, All );
      }
   }
}

void ConvectionSimulation::solveStokes()
{
   walberla::WcTimer localTimer;
   if ( TN.simulationParameters.tempDependentViscosity )
   {
      storage->getTimingTree()->start( "TerraNeo update viscosity" );
      updateViscosity();
      storage->getTimingTree()->stop( "TerraNeo update viscosity" );
   }

   localTimer.start();
   storage->getTimingTree()->start( "TerraNeo setup Stokes RHS" );
   setupStokesRHS();
   storage->getTimingTree()->stop( "TerraNeo setup Stokes RHS" );
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

   real_t stokesResidual = calculateStokesResidual( TN.domainParameters.maxLevel );

   TN.solverParameters.vCycleResidualUPrev       = stokesResidual;
   TN.solverParameters.initialResidualU          = stokesResidual;
   TN.solverParameters.averageResidualReductionU = real_c( 0 );
   TN.solverParameters.numVCycles                = 0;

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes residual (initial): " << stokesResidual );
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
   projectionOperator->project( *( p2p1StokesFunctionContainer["StokesRHS"] ), TN.domainParameters.maxLevel, FreeslipBoundary );
   stokesSolverFS->solve( *stokesOperatorFS,
                          *( p2p1StokesFunctionContainer["VelocityFE"] ),
                          *( p2p1StokesFunctionContainer["StokesRHS"] ),
                          TN.domainParameters.maxLevel );
   // stokesSolver->solve( *stokesOperator, *(p2p1StokesFunctionContainer["VelocityFE"]), *(p2p1StokesFunctionContainer["StokesRHS"]), TN.domainParameters.maxLevel );
   storage->getTimingTree()->stop( "TerraNeo solve Stokes" );
   localTimer.end();

   real_t timeStokesSolve = localTimer.last();
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes time [s]: " << timeStokesSolve );

   stokesResidual = calculateStokesResidual( TN.domainParameters.maxLevel );

   if ( std::isnan( stokesResidual ) )
   {
      WALBERLA_ABORT( "NaN values detected in velocity after Stokes solve" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes residual (final): " << stokesResidual );
   if ( TN.outputParameters.createTimingDB )
   {
      db->setVariableEntry( "Stokes_residual_final", stokesResidual );
      db->setVariableEntry( "Time_setup_Stokes_RHS", setupStokesTimer );
      db->setVariableEntry( "Time_solve_Stokes", timeStokesSolve );
   }
}

real_t ConvectionSimulation::calculateStokesResidual( uint_t level )
{
   stokesOperatorFS->apply( *( p2p1StokesFunctionContainer["VelocityFE"] ),
                            *( p2p1StokesFunctionContainer["StokesTmp1"] ),
                            level,
                            Inner | NeumannBoundary | FreeslipBoundary );
   p2p1StokesFunctionContainer["StokesTmp1"]->assign(
       { real_c( 1 ), real_c( -1 ) },
       { *( p2p1StokesFunctionContainer["StokesTmp1"] ), *( p2p1StokesFunctionContainer["StokesRHS"] ) },
       level,
       Inner | NeumannBoundary | FreeslipBoundary );
   return std::sqrt( p2p1StokesFunctionContainer["StokesTmp1"]->dotGlobal(
       *( p2p1StokesFunctionContainer["StokesTmp1"] ), level, Inner | NeumannBoundary | FreeslipBoundary ) );
}

real_t ConvectionSimulation::calculateEnergyResidual( uint_t level )
{
   P2Function< real_t >& tempFunc = p2p1StokesFunctionContainer["StokesTmp1"]->uvw().component( 0u );

   transportOperatorTALA->apply(
       *( p2ScalarFunctionContainer["TemperatureFE"] ), tempFunc, level, Inner | NeumannBoundary | FreeslipBoundary );
   tempFunc.assign( { real_c( 1 ), real_c( -1 ) },
                    { tempFunc, *( p2ScalarFunctionContainer["EnergyRHSWeak"] ) },
                    level,
                    Inner | NeumannBoundary | FreeslipBoundary );
   return std::sqrt( tempFunc.dotGlobal( tempFunc, level, Inner | NeumannBoundary | FreeslipBoundary ) );
}
} // namespace terraneo