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
   if ( TN.simulationParameters.timeStep == 0 )
   {
      //set up logging
      WALBERLA_ROOT_SECTION()
      {
         walberla::logging::Logging::instance()->includeLoggingToFile( TN.outputParameters.outputDirectory + "/" +
                                                                       TN.outputParameters.outputBaseName + ".out" );
      }

      // setupStokesRHS();
      if ( TN.outputParameters.ADIOS2StartFromCheckpoint )
      {
#ifdef HYTEG_BUILD_WITH_ADIOS2
         checkpointImporter->restoreFunction( *( p2ScalarFunctionContainer["TemperatureFE"] ) );
         solveStokes();
         dataOutput();
#else
         WALBERLA_ABORT( " No submodule ADIOS2 enabled! Loading Checkpoint not possible - Abort simulation " );
#endif
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

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT(
       "Step size: " << ( TN.simulationParameters.dt * TN.physicalParameters.mantleThickness ) /
                            ( TN.physicalParameters.characteristicVelocity * TN.simulationParameters.plateVelocityScaling *
                              TN.simulationParameters.secondsPerMyr )
                     << " Ma" );
   WALBERLA_LOG_INFO_ON_ROOT( "Velocity max magnitude: " << vMax << " " );

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
                            true, true, false );

   // Reset temperature on boundary to initial values

   p2ScalarFunctionContainer["TemperaturePrev"]->assign(
       { real_c( 1 ) }, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, TN.domainParameters.maxLevel, All );

   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      if ( TN.initialisationParameters.temperatureNoise )
      {
         p2ScalarFunctionContainer["TemperatureFE"]->interpolate(
             temperatureWhiteNoise( *temperatureInitParams, *temperatureReferenceFct, TN.initialisationParameters.noiseFactor ),
             l,
             DirichletBoundary );
      }
      else
      {
         p2ScalarFunctionContainer["TemperatureFE"]->interpolate(
             temperatureSPH( *temperatureInitParams,
                             *temperatureReferenceFct,
                             TN.initialisationParameters.tempInit,
                             TN.initialisationParameters.deg,
                             TN.initialisationParameters.ord,
                             TN.initialisationParameters.lmax,
                             TN.initialisationParameters.lmin,
                             TN.initialisationParameters.superposition,
                             TN.initialisationParameters.buoyancyFactor,
                             TN.physicalParameters.initialTemperatureSteepness ),
             l,
             DirichletBoundary );
      }
   }

   /*############ ENERGY STEP ############*/

   // setupEnergyRHS();

   storage->getTimingTree()->start( "Solve Energy equation" );
   solveEnergy();
   storage->getTimingTree()->stop( "Solve Energy equation" );

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

   if ( !TN.outputParameters.outputMyr )
   {
      dataOutput();
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Finished step: " << TN.simulationParameters.timeStep );
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
          if ( TN.simulationParameters.radialProfile )
          {
             updateNonDimParameters( x );
          }
          if ( TN.simulationParameters.shearHeatingScaling > 1 )
          {
             WALBERLA_ABORT( "Shear heating scaling factor > 1 not allowed! --- Abort simulation ---" );
          }

          if ( radius > TN.domainParameters.rMax -
                            ( TN.simulationParameters.lithosphereThickness * 1000 / TN.physicalParameters.mantleThickness ) )
          {
             return ( ( TN.physicalParameters.dissipationNumber * TN.physicalParameters.pecletNumber /
                        ( TN.physicalParameters.rayleighNumber * density[0] ) ) *
                      TN.simulationParameters.shearHeatingScaling );
          }
          else
          {
             return TN.physicalParameters.dissipationNumber * TN.physicalParameters.pecletNumber /
                    ( TN.physicalParameters.rayleighNumber * density[0] );
          }
       };

   std::function< real_t( const Point3D& ) > internalHeatingCoeffCalc = [this]( const Point3D& x ) {
      real_t intHeatingFactor = 1.0;
      if ( TN.simulationParameters.radialProfile )
      {
         updateNonDimParameters( x );
      }
      return TN.physicalParameters.hNumber * intHeatingFactor;
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
             if ( TN.simulationParameters.radialProfile )
             {
                updateNonDimParameters( x );
             }
             return ( -( TN.physicalParameters.rayleighNumber / TN.physicalParameters.pecletNumber ) * densityFunction( x ) *
                      deltaT[0] );
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
         if ( TN.simulationParameters.haveThermalExpProfile || TN.simulationParameters.haveSpecificHeatCapProfile )
         {
            // Update gradRho/Rho with new non-Dim paramters Di and alpha.
            // grad(rho)/rho = - ( Di / gamma ) * r_hat
            // std::function< real_t( const Point3D& ) > updateDensity = [&]( const Point3D& x ) { return densityFunc( x ); };
            p2ScalarFunctionContainer["DensityFE"]->interpolate( densityFunc, l, All );
         }

         frozenVelocityRHS->apply(
             p2p1StokesFunctionContainer["VelocityFE"]->uvw(), p2p1StokesFunctionContainer["StokesRHS"]->p(), l, All );

         p2p1StokesFunctionContainer["StokesRHS"]->p().assign(
             { -1.0 }, { p2p1StokesFunctionContainer["StokesRHS"]->p() }, l, All );
      }
      else
      {
         p2p1StokesFunctionContainer["StokesRHS"]->p().interpolate( real_c( 0 ), l, All );
      }
   }
}

void ConvectionSimulation::solveStokes()
{
   if ( TN.simulationParameters.tempDependentViscosity )
   {
      updateViscosity();
   }

   if ( TN.simulationParameters.tempDependentViscosity && TN.simulationParameters.resetSolver &&
        ( TN.simulationParameters.timeStep != 0 ) &&
        ( TN.simulationParameters.timeStep % TN.simulationParameters.resetSolverFrequency == 0 ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( "------ Resetting Solvers and Operators ------" );
      WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------" );

      std::function< real_t( const Point3D& ) > zeros = []( const Point3D& ) { return real_c( 0 ); };

      for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; ++l )
      {
         //save current velocity in temorary, as U and F must be altered for solver setup
         p2p1StokesFunctionContainer["StokesTmp1"]->assign(
             { real_c( 1 ) }, { *( p2p1StokesFunctionContainer["VelocityFE"] ) }, l, All );
         p2p1StokesFunctionContainer["StokesRHS"]->interpolate( { zeros, zeros, zeros }, l, All );
      }

      // for temperature dependent viscosity the spectral radius might change after several solving
      // iterations -> resetting solvers and operators

      setupSolversAndOperators();

      //after setup, reset velocity to values stored in temporary
      for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; ++l )
      {
         p2p1StokesFunctionContainer["VelocityFE"]->assign(
             { real_c( 1 ) }, { *( p2p1StokesFunctionContainer["StokesTmp1"] ) }, l, All );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( "------ Resetting Solvers and Operators: Finished ------" );
      WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------------------------" );
   }

   setupStokesRHS();

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

   walberla::WcTimer localTimer;
   real_t            stokesResidual = calculateStokesResidual( TN.domainParameters.maxLevel );

   TN.solverParameters.vCycleResidualUPrev       = stokesResidual;
   TN.solverParameters.initialResidualU          = stokesResidual;
   TN.solverParameters.averageResidualReductionU = real_c( 0 );
   TN.solverParameters.numVCycles                = 0;

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes residual (initial): " << stokesResidual );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   if ( stokesResidual > TN.solverParameters.stokesKillTolerance && ( TN.simulationParameters.timeStep != 0 ) )
   {
      WALBERLA_ABORT( "Residual " << stokesResidual << " exceeds tolerance of " << TN.solverParameters.stokesKillTolerance
                                  << ". ABORTING SIMULATION. " );
   }

   localTimer.start();
   storage->getTimingTree()->start( "Stokes Solve" );
   projectionOperator->project( *( p2p1StokesFunctionContainer["StokesRHS"] ), TN.domainParameters.maxLevel, FreeslipBoundary );
   stokesSolverFS->solve( *stokesOperatorFS,
                          *( p2p1StokesFunctionContainer["VelocityFE"] ),
                          *( p2p1StokesFunctionContainer["StokesRHS"] ),
                          TN.domainParameters.maxLevel );
   // stokesSolver->solve( *stokesOperator, *(p2p1StokesFunctionContainer["VelocityFE"]), *(p2p1StokesFunctionContainer["StokesRHS"]), TN.domainParameters.maxLevel );
   storage->getTimingTree()->stop( "Stokes Solve" );
   localTimer.end();

   real_t timeStokes = localTimer.last();
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes time [s]: " << timeStokes );

   stokesResidual = calculateStokesResidual( TN.domainParameters.maxLevel );

   if ( std::isnan( stokesResidual ) )
   {
      WALBERLA_ABORT( "NaN values detected in velocity after Stokes solve" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes residual (final): " << stokesResidual );
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