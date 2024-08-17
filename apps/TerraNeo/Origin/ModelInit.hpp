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

namespace terraneo {
////////////////////////
//   Initialisation   //
////////////////////////

ConvectionSimulation::ConvectionSimulation( const walberla::Config::BlockHandle& mainConf )
{
   TN = terraneo::parseConfig( mainConf );

   viscosityFunc  = std::bind( &ConvectionSimulation::viscosityFunction, this, std::placeholders::_1, std::placeholders::_2 );
   densityFunc    = std::bind( &ConvectionSimulation::densityFunction, this, std::placeholders::_1 );
   diffFactorFunc = std::bind( &ConvectionSimulation::diffPreFactorFunction, this, std::placeholders::_1 );
}

void ConvectionSimulation::init()
{
   if ( !std::filesystem::exists( TN.outputParameters.outputDirectory ) )
   {
      std::filesystem::create_directories( TN.outputParameters.outputDirectory );
   }

   if ( !std::filesystem::exists( TN.outputParameters.outputDirectory ) )
   {
      WALBERLA_ABORT( "Failed to create output directory \"" << TN.outputParameters.outputDirectory << "\"" );
   }

   if ( TN.outputParameters.outputProfiles && !std::filesystem::exists( TN.outputParameters.outputDirectory + "/Profiles" ) )
   {
      std::filesystem::create_directories( TN.outputParameters.outputDirectory + "/Profiles" );
   }

   setupDomain();

   TN.simulationParameters.unknownsTemperature = numberOfGlobalDoFs< P2FunctionTag >( *storage, TN.domainParameters.maxLevel );
   TN.simulationParameters.unknownsStokes =
       numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, TN.domainParameters.maxLevel );
   TN.simulationParameters.hMin = MeshQuality::getMinimalEdgeLength( storage, TN.domainParameters.maxLevel );
   TN.simulationParameters.hMax = MeshQuality::getMaximalEdgeLength( storage, TN.domainParameters.maxLevel );

   printConfig( TN );

   setupBoundaryConditions();
   setupFunctions();
   initialiseFunctions();
   setupSolversAndOperators();
   setupOutput();
}

// Setup the domain

void ConvectionSimulation::setupDomain()
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
   meshInfo          = MeshInfo::meshSphericalShell(
       TN.domainParameters.nTan, TN.domainParameters.nRad, TN.domainParameters.rMin, TN.domainParameters.rMax );

   auto setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, walberla::mpi::MPIManager::instance()->numProcesses() );

   loadbalancing::roundRobinVolume( *setupStorage );

   IcosahedralShellMap::setMap( *setupStorage );

   storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );
}

// Setup boundary conditions

void ConvectionSimulation::setupBoundaryConditions()
{
   // For Temperature always use Dirichlet BC on Surface and CMB
   idSurfaceT = bcTemperature.createDirichletBC( "surface", { MeshInfo::hollowFlag::flagOuterBoundary } );
   idCMBT     = bcTemperature.createDirichletBC( "cmb", { MeshInfo::hollowFlag::flagInnerBoundary } );

   // Boundary Conditions for velocity
   // BC 1: No-Slip/No-Slip; 2: Free-Slip/Free-Slip; 3: No-Slip/Free-Slip

   switch ( TN.simulationParameters.boundaryCond )
   {
   case 1:
      idSurface = bcVelocity.createDirichletBC( "surface", { MeshInfo::hollowFlag::flagOuterBoundary } );
      idCMB     = bcVelocity.createDirichletBC( "cmb", { MeshInfo::hollowFlag::flagInnerBoundary } );
      break;

   case 2:
      idSurface = bcVelocity.createFreeslipBC( "surface", { MeshInfo::hollowFlag::flagOuterBoundary } );
      idCMB     = bcVelocity.createFreeslipBC( "cmb", { MeshInfo::hollowFlag::flagInnerBoundary } );
      TN.simulationParameters.boundaryCondFreeSlip = true;
      break;

   case 3:
      idSurface = bcVelocity.createDirichletBC( "surface", { MeshInfo::hollowFlag::flagOuterBoundary } );
      idCMB     = bcVelocity.createFreeslipBC( "cmb", { MeshInfo::hollowFlag::flagInnerBoundary } );
      TN.simulationParameters.boundaryCondFreeSlip = true;
      break;

   default:
      idSurface = bcVelocity.createDirichletBC( "surface", { MeshInfo::hollowFlag::flagOuterBoundary } );
      idCMB     = bcVelocity.createDirichletBC( "cmb", { MeshInfo::hollowFlag::flagInnerBoundary } );
      break;
   }
}

// Setup functions

void ConvectionSimulation::setupFunctions()
{
   if ( !storage )
   {
      setupDomain();
   }

   for ( auto& [p2p1Names, p2p1MinLevel, p2p1MaxLevel, p2p1BcType] : p2p1StokesFunctionDict )
   {
      p2p1MinLevel = TN.domainParameters.minLevel;
      p2p1MaxLevel = TN.domainParameters.maxLevel;
   }

   for ( auto [p2p1Names, p2p1MinLevel, p2p1MaxLevel, p2p1BcType] : p2p1StokesFunctionDict )
   {
      switch ( p2p1BcType )
      {
      case BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION: {
         p2p1StokesFunctionContainer.emplace( std::make_pair(
             p2p1Names, std::make_shared< StokesFunctionP2P1 >( p2p1Names, storage, p2p1MinLevel, p2p1MaxLevel, bcVelocity ) ) );
      }
      break;
      case BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION: {
         WALBERLA_LOG_WARNING_ON_ROOT( "bcTemperature looks odd for P2P1TaylorHoodFunction" );
         p2p1StokesFunctionContainer.emplace( std::make_pair(
             p2p1Names,
             std::make_shared< StokesFunctionP2P1 >( p2p1Names, storage, p2p1MinLevel, p2p1MaxLevel, bcTemperature ) ) );
      }
      break;
      default: {
         p2p1StokesFunctionContainer.emplace( std::make_pair(
             p2p1Names, std::make_shared< StokesFunctionP2P1 >( p2p1Names, storage, p2p1MinLevel, p2p1MaxLevel ) ) );
      }
      }
   }

   for ( auto& [p2Names, p2MinLevel, p2MaxLevel, p2BcType] : p2ScalarFunctionDict )
   {
      p2MinLevel = TN.domainParameters.minLevel;
      p2MaxLevel = TN.domainParameters.maxLevel;
   }

   for ( auto [p2Names, p2MinLevel, p2MaxLevel, p2BcType] : p2ScalarFunctionDict )
   {
      switch ( p2BcType )
      {
      case BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION: {
         WALBERLA_LOG_WARNING_ON_ROOT( "bcVelocity looks odd for P2Function" );
         p2ScalarFunctionContainer.emplace( std::make_pair(
             p2Names, std::make_shared< ScalarFunctionP2 >( p2Names, storage, p2MinLevel, p2MaxLevel, bcVelocity ) ) );
      }
      break;
      case BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION: {
         p2ScalarFunctionContainer.emplace( std::make_pair(
             p2Names, std::make_shared< ScalarFunctionP2 >( p2Names, storage, p2MinLevel, p2MaxLevel, bcTemperature ) ) );
      }
      break;
      default: {
         p2ScalarFunctionContainer.emplace(
             std::make_pair( p2Names, std::make_shared< ScalarFunctionP2 >( p2Names, storage, p2MinLevel, p2MaxLevel ) ) );
      }
      }
   }

   for ( auto& [p2VecNames, p2VecMinLevel, p2VecMaxLevel, p2VecBcType] : p2VectorFunctionDict )
   {
      p2VecMinLevel = TN.domainParameters.minLevel;
      p2VecMaxLevel = TN.domainParameters.maxLevel;
   }

   for ( auto [p2VecNames, p2VecMinLevel, p2VecMaxLevel, p2VecBcType] : p2VectorFunctionDict )
   {
      switch ( p2VecBcType )
      {
      case BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION: {
         p2VectorFunctionContainer.emplace( std::make_pair(
             p2VecNames,
             std::make_shared< VectorFunctionP2 >( p2VecNames, storage, p2VecMinLevel, p2VecMaxLevel, bcVelocity ) ) );
      }
      break;
      case BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION: {
         WALBERLA_LOG_WARNING_ON_ROOT( "bcTemperature looks odd for P2VectorFunction" );
         p2VectorFunctionContainer.emplace( std::make_pair(
             p2VecNames,
             std::make_shared< VectorFunctionP2 >( p2VecNames, storage, p2VecMinLevel, p2VecMaxLevel, bcTemperature ) ) );
      }
      break;
      default: {
         p2VectorFunctionContainer.emplace( std::make_pair(
             p2VecNames, std::make_shared< VectorFunctionP2 >( p2VecNames, storage, p2VecMinLevel, p2VecMaxLevel ) ) );
      }
      }
   }
}

void ConvectionSimulation::initialiseFunctions()
{
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup initial state -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   walberla::math::seedRandomGenerator( 42 );
   std::function< real_t( const Point3D& ) > randFunc = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };
   std::function< real_t( const Point3D& ) > zeros = []( const Point3D& ) { return real_c( 0 ); };

   TemperatureInitializationParameters tempInitParams( TN.physicalParameters.cmbTemp,
                                                       TN.physicalParameters.surfaceTemp,
                                                       TN.physicalParameters.adiabatSurfaceTemp,
                                                       TN.physicalParameters.dissipationNumber,
                                                       TN.domainParameters.rMin,
                                                       TN.domainParameters.rMax );

   temperatureInitParams = std::make_shared< TemperatureInitializationParameters >( tempInitParams );

   temperatureReferenceFct = std::make_shared< std::function< real_t( const Point3D& ) > >(
       terraneo::temperatureReferenceExponential( tempInitParams ) );

   auto referenceTemperature = temperatureReferenceExponential( *temperatureInitParams );

   if ( TN.initialisationParameters.temperatureNoise )
   {
      auto initTemperatureWhiteNoise =
          temperatureWhiteNoise( *temperatureInitParams, *temperatureReferenceFct, TN.initialisationParameters.noiseFactor );

      for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
      {
         p2ScalarFunctionContainer["TemperatureFE"]->interpolate( initTemperatureWhiteNoise, l, All );
      }
   }
   else
   {
      auto initTemperatureSPH = temperatureSPH( *temperatureInitParams,
                                                *temperatureReferenceFct,
                                                TN.initialisationParameters.tempInit,
                                                TN.initialisationParameters.deg,
                                                TN.initialisationParameters.ord,
                                                TN.initialisationParameters.lmax,
                                                TN.initialisationParameters.lmin,
                                                TN.initialisationParameters.superposition,
                                                TN.initialisationParameters.buoyancyFactor,
                                                TN.physicalParameters.initialTemperatureSteepness );

      for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
      {
         p2ScalarFunctionContainer["TemperatureFE"]->interpolate( initTemperatureSPH, l, All );
      }
   }
   if ( TN.simulationParameters.simulationType == "CirculationModel" )
   {
      // initialise plate velocity oracle
      WALBERLA_LOG_INFO_ON_ROOT( "Setup Oracle for Plates" );
      terraneo::oracle = std::make_shared< terraneo::plates::PlateVelocityProvider >(
          TN.simulationParameters.fnameTopologies, TN.simulationParameters.fnameReconstructions );
   }
   // Assign temperature field to temperaturePrev

   for ( uint_t level = TN.domainParameters.minLevel; level <= TN.domainParameters.maxLevel; ++level )
   {
      p2ScalarFunctionContainer["TemperaturePrev"]->assign(
          { real_c( 1 ) }, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, level, All );
      p2ScalarFunctionContainer["TemperatureTmp"]->assign(
          { real_c( 1 ) }, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, level, All );
      p2ScalarFunctionContainer["DensityFE"]->interpolate( densityFunc, level, All );
      p2ScalarFunctionContainer["DiffusionFE"]->interpolate( diffFactorFunc, level, All );

      //set plate velocities for timestep 0 / initialAge
      if ( TN.simulationParameters.simulationType == "CirculationModel" )
      {
         updatePlateVelocities( *( p2p1StokesFunctionContainer["StokesLHS"] ) );
      }
      else
      {
         //currently initialising CMB and surface to zeros for NoSlipNoSlip case
         p2p1StokesFunctionContainer["StokesLHS"]->uvw().interpolate( { zeros, zeros, zeros }, level, idSurface );
         p2p1StokesFunctionContainer["StokesLHSPrev"]->uvw().interpolate( { zeros, zeros, zeros }, level, idSurface );
      }

      p2p1StokesFunctionContainer["StokesLHS"]->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      p2p1StokesFunctionContainer["StokesLHSPrev"]->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      p2p1StokesFunctionContainer["StokesRHS"]->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      p2p1StokesFunctionContainer["StokesTmp"]->uvw().interpolate( { zeros, zeros, zeros }, level, All );

      p2VectorFunctionContainer["InwardNormal"]->interpolate( { normalX, normalY, normalZ }, level, All );

      p2VectorFunctionContainer["OppositeGravityField"]->assign(
          { -1.0 }, { *( p2VectorFunctionContainer["InwardNormal"] ) }, level, All );

      p2ScalarFunctionContainer["OnesFE"]->interpolate( real_c( 1 ), level, All );
      p2VectorFunctionContainer["GradRhoOverRho"]->interpolate( { normalX, normalY, normalZ }, level, All );

      //grad(rho)/rho = - ( Di / gamma ) * r_hat
      p2VectorFunctionContainer["GradRhoOverRho"]->assign(
          { TN.physicalParameters.dissipationNumber / TN.physicalParameters.grueneisenParameter },
          { *( p2VectorFunctionContainer["GradRhoOverRho"] ) },
          level,
          All );
   }

   auto temperatureRadialProfile = computeRadialProfile( *( p2ScalarFunctionContainer["TemperatureFE"] ),
                                                         TN.domainParameters.rMin,
                                                         TN.domainParameters.rMax,
                                                         TN.domainParameters.nRad,
                                                         TN.domainParameters.maxLevel );
   temperatureProfiles           = std::make_shared< RadialProfile >( temperatureRadialProfile );

   if ( TN.outputParameters.outputProfiles && TN.simulationParameters.tempDependentViscosity )
   {
      updateViscosity();
      auto viscosityRadialProfile = computeRadialProfile( *( p2ScalarFunctionContainer["ViscosityFE"] ),
                                                          TN.domainParameters.rMin,
                                                          TN.domainParameters.rMax,
                                                          TN.domainParameters.nRad,
                                                          TN.domainParameters.maxLevel );
      viscosityProfiles           = std::make_shared< RadialProfile >( viscosityRadialProfile );
   }

   referenceTemperatureFct = [this]( const Point3D& x ) {
      real_t radius = x.norm();
      if ( TN.simulationParameters.adaptiveRefTemp )
      {
         uint_t shell = static_cast< uint_t >(
             std::round( real_c( TN.simulationParameters.numLayers ) *
                         ( ( radius - TN.domainParameters.rMin ) / ( TN.domainParameters.rMax - TN.domainParameters.rMin ) ) ) );
         WALBERLA_ASSERT( shell < temperatureProfiles->mean.size() );
         return temperatureProfiles->mean.at( shell );
      }
      else
      {
         return referenceTemperatureFunction( x );
      }
   };

   p2ScalarFunctionContainer["TemperatureReference"]->interpolate( referenceTemperatureFct, TN.domainParameters.maxLevel, All );
}

void ConvectionSimulation::setupSolversAndOperators()
{
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup solvers & operators -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   auto                                      normalFunc_ = [&]( const Point3D& p, Point3D& n ) { normalFunc( p, n ); };
   std::function< real_t( const Point3D& ) > zeros       = []( const Point3D& ) { return real_c( 0 ); };

   walberla::math::seedRandomGenerator( 42 );
   std::function< real_t( const Point3D& ) > randFunc = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };
   //non-dimensionalise viscosity such that minimum value = 1
   updateViscosity();

   projectionOperator = std::make_shared< P2ProjectNormalOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, normalFunc_ );

   stokesOperator = std::make_shared< StokesOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *( p2ScalarFunctionContainer["ViscosityFE"] ) );

   stokesOperatorFS = std::make_shared< StokesOperatorFS >( storage,
                                                            TN.domainParameters.minLevel,
                                                            TN.domainParameters.maxLevel,
                                                            p2ScalarFunctionContainer["ViscosityFE"]->getVertexDoFFunction(),
                                                            p2ScalarFunctionContainer["ViscosityFEInv"]->getVertexDoFFunction(),
                                                            *projectionOperator,
                                                            bcVelocity );

   if ( TN.solverParameters.solverFlag == 0u )
   {
      stokesSolverFS = hyteg::solvertemplates::stokesGMGFSSolver(
          storage,
          TN.domainParameters.minLevel,
          TN.domainParameters.maxLevel,
          stokesOperatorFS,
          projectionOperator,
          bcVelocity,
          TN.solverParameters.estimateUzawaOmega,
          false,
          {
              { solvertemplates::StokesGMGFSSolverParamKey::NUM_POWER_ITERATIONS_SPECTRUM,
                real_c( TN.solverParameters.numPowerIterations ) },
              { solvertemplates::StokesGMGFSSolverParamKey::FGMRES_UZAWA_PRECONDITIONED_OUTER_ITER,
                real_c( TN.solverParameters.FGMRESOuterIterations ) },
              { solvertemplates::StokesGMGFSSolverParamKey::FGMRES_UZAWA_PRECONDITIONED_OUTER_TOLERANCE,
                real_c( TN.solverParameters.FGMRESTolerance ) },
              { solvertemplates::StokesGMGFSSolverParamKey::INEXACT_UZAWA_VELOCITY_ITER,
                real_c( TN.solverParameters.uzawaIterations ) },
              { solvertemplates::StokesGMGFSSolverParamKey::INEXACT_UZAWA_OMEGA, real_c( TN.solverParameters.uzawaOmega ) },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_CG_SOLVER_MG_PRECONDITIONED_ITER,
                real_c( TN.solverParameters.ABlockMGIterations ) },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE,
                real_c( TN.solverParameters.ABlockMGTolerance ) },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_MG_PRESMOOTH,
                real_c( TN.solverParameters.ABlockMGPreSmooth ) },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_MG_POSTSMOOTH,
                real_c( TN.solverParameters.ABlockMGPostSmooth ) },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_COARSE_ITER,
                real_c( TN.solverParameters.ABlockCoarseGridIterations ) },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_COARSE_TOLERANCE,
                real_c( TN.solverParameters.ABlockCoarseGridTolerance ) },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_CG_SOLVER_MG_PRECONDITIONED_ITER,
                real_c( TN.solverParameters.SchurMGIterations ) },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE,
                real_c( TN.solverParameters.SchurMGTolerance ) },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_MG_PRESMOOTH, real_c( TN.solverParameters.SchurMGPreSmooth ) },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_MG_POSTSMOOTH,
                real_c( TN.solverParameters.SchurMGPostSmooth ) },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_COARSE_GRID_CG_ITER,
                real_c( TN.solverParameters.SchurCoarseGridIterations ) },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_COARSE_GRID_CG_TOLERANCE,
                real_c( TN.solverParameters.SchurCoarseGridTolerance ) },
          } );
   }
   else if ( TN.solverParameters.solverFlag == 1u )
   {
      auto stopIterationCallback =
          [&]( const StokesOperatorFS& _A, const StokesFunction& _u, const StokesFunction& _b, const uint_t _level ) {
             WALBERLA_UNUSED( _A );
             WALBERLA_UNUSED( _u );
             WALBERLA_UNUSED( _b );
             real_t stokesResidual;

             stokesResidual = calculateStokesResidual( _level );

             if ( TN.solverParameters.numVCycles == 0 )
             {
                WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
                    "[Uzawa] iter %3d | residual: %10.5e | initial ", 0, TN.solverParameters.vCycleResidualUPrev ) );
             }

             auto reductionRateU = stokesResidual / TN.solverParameters.vCycleResidualUPrev;

             TN.solverParameters.vCycleResidualUPrev = stokesResidual;

             TN.solverParameters.numVCycles++;
             TN.solverParameters.averageResidualReductionU += reductionRateU;

             if ( TN.simulationParameters.verbose )
             {
                WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "[Uzawa] iter %3d | residual: %10.5e | reduction: %10.5e ",
                                                             TN.solverParameters.numVCycles,
                                                             stokesResidual,
                                                             reductionRateU ) );
             }

             if ( stokesResidual / TN.solverParameters.initialResidualU < TN.solverParameters.stokesRelativeResidualUTolerance )
             {
                WALBERLA_LOG_INFO_ON_ROOT( "[Uzawa] reached relative residual threshold" )
                return true;
             }

             if ( stokesResidual < TN.solverParameters.stokesAbsoluteResidualUTolerance )
             {
                WALBERLA_LOG_INFO_ON_ROOT( "[Uzawa] reached absolute residual threshold" )
                return true;
             }
             return false;
          };

      auto multigridSolver = solvertemplates::stokesGMGUzawaFSSolver< StokesOperatorFS, P2ProjectNormalOperator >(
          storage,
          TN.domainParameters.minLevel,
          TN.domainParameters.maxLevel,
          stokesOperatorFS,
          projectionOperator,
          bcVelocity,
          false,
          { { solvertemplates::StokesGMGUzawaFSSolverParamKey::NUM_POWER_ITERATIONS_SPECTRUM,
              real_c( TN.solverParameters.numPowerIterations ) },
            { solvertemplates::StokesGMGUzawaFSSolverParamKey::NUM_COARSE_GRID_ITERATIONS,
              real_c( TN.solverParameters.stokesUzawaCoarseGridIter ) },
            { solvertemplates::StokesGMGUzawaFSSolverParamKey::COARSE_GRID_TOLERANCE,
              real_c( TN.solverParameters.stokesUzawaCoarseGridTol ) },
            { solvertemplates::StokesGMGUzawaFSSolverParamKey::UZAWA_OMEGA, real_c( TN.solverParameters.uzawaOmega ) },
            { solvertemplates::StokesGMGUzawaFSSolverParamKey::MG_PRE_SMOOTH, real_c( TN.solverParameters.ABlockMGPreSmooth ) },
            { solvertemplates::StokesGMGUzawaFSSolverParamKey::MG_POST_SMOOTH, real_c( TN.solverParameters.ABlockMGPostSmooth ) },
            { solvertemplates::StokesGMGUzawaFSSolverParamKey::UZAWA_VELOCITY_ITER,
              real_c( TN.solverParameters.uzawaIterations ) },
            { solvertemplates::StokesGMGUzawaFSSolverParamKey::SMOOTH_INCREMENT_COARSE_GRID,
              real_c( TN.solverParameters.stokesSmoothIncrementCoarseGrid ) } } );

      stokesSolverFS = std::make_shared< SolverLoop< StokesOperatorFS > >(
          multigridSolver, TN.solverParameters.stokesMaxNumIterations, stopIterationCallback );
   }
   else
   {
      WALBERLA_ABORT( "Unknown solver type" );
   }

   P2MassOperator = std::make_shared< P2ElementwiseBlendingMassOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   frozenVelocityRHSX =
       std::make_shared< FrozenVelocityOperator >( storage,
                                                   TN.domainParameters.minLevel,
                                                   TN.domainParameters.maxLevel,
                                                   p2VectorFunctionContainer["GradRhoOverRho"]->component( 0U ) );
   frozenVelocityRHSY =
       std::make_shared< FrozenVelocityOperator >( storage,
                                                   TN.domainParameters.minLevel,
                                                   TN.domainParameters.maxLevel,
                                                   p2VectorFunctionContainer["GradRhoOverRho"]->component( 1U ) );
   frozenVelocityRHSZ =
       std::make_shared< FrozenVelocityOperator >( storage,
                                                   TN.domainParameters.minLevel,
                                                   TN.domainParameters.maxLevel,
                                                   p2VectorFunctionContainer["GradRhoOverRho"]->component( 2U ) );

   /////////////////////////
   // Diffusion Operator //
   ////////////////////////

   transportOperatorTALA = std::make_shared< P2TransportIcosahedralShellMapOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   transportOperatorTALA->setVelocity( p2p1StokesFunctionContainer["StokesLHS"] );
   transportOperatorTALA->setViscosity( p2ScalarFunctionContainer["ViscosityFE"] );
   transportOperatorTALA->setTemperature( p2ScalarFunctionContainer["TemperatureFE"] );

   transportOperatorTALA->setInvGravity( p2VectorFunctionContainer["OppositeGravityField"] );

   transportOperatorTALA->setDiffusivityCoeff( p2ScalarFunctionContainer["DiffusionFE"] );
   transportOperatorTALA->setAdiabaticCoeff( p2ScalarFunctionContainer["AdiabaticTermCoeff"] );
   transportOperatorTALA->setShearHeatingCoeff( p2ScalarFunctionContainer["ShearHeatingTermCoeff"] );
   transportOperatorTALA->setConstEnergyCoeff( p2ScalarFunctionContainer["ConstEnergyCoeff"] );
   transportOperatorTALA->setSurfTempCoeff( p2ScalarFunctionContainer["SurfaceTempCoeff"] );

   transportOperatorTALA->setReferenceTemperature( p2ScalarFunctionContainer["TemperatureReference"] );

   transportOperatorTALA->setTALADict(
       { { TransportOperatorTermKey::ADIABATIC_HEATING_TERM, TN.simulationParameters.adiabaticHeating },
         { TransportOperatorTermKey::SHEAR_HEATING_TERM, TN.simulationParameters.shearHeating },
         { TransportOperatorTermKey::INTERNAL_HEATING_TERM, TN.simulationParameters.internalHeating } } );

   transportOperatorTALA->initializeOperators();

   transportOperatorRHS = std::make_shared< P2TransportRHSIcosahedralShellMapOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   transportOperatorRHS->setVelocity( p2p1StokesFunctionContainer["StokesLHS"] );
   transportOperatorRHS->setViscosity( p2ScalarFunctionContainer["ViscosityFE"] );
   transportOperatorRHS->setTemperature( p2ScalarFunctionContainer["TemperatureFE"] );

   transportOperatorRHS->setInvGravity( p2VectorFunctionContainer["OppositeGravityField"] );

   transportOperatorRHS->setDiffusivityCoeff( p2ScalarFunctionContainer["DiffusionFE"] );
   transportOperatorRHS->setAdiabaticCoeff( p2ScalarFunctionContainer["AdiabaticTermCoeff"] );
   transportOperatorRHS->setShearHeatingCoeff( p2ScalarFunctionContainer["ShearHeatingTermCoeff"] );
   transportOperatorRHS->setConstEnergyCoeff( p2ScalarFunctionContainer["ConstEnergyCoeff"] );

   transportOperatorRHS->setReferenceTemperature( p2ScalarFunctionContainer["TemperatureReference"] );

   transportOperatorRHS->setTALADict(
       { { TransportRHSOperatorTermKey::ADIABATIC_HEATING_TERM, TN.simulationParameters.adiabaticHeating },
         { TransportRHSOperatorTermKey::SHEAR_HEATING_TERM, TN.simulationParameters.shearHeating },
         { TransportRHSOperatorTermKey::INTERNAL_HEATING_TERM, TN.simulationParameters.internalHeating } } );

   transportOperatorRHS->initializeOperators();

   diffusionOperator = std::make_shared< DiffusionOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *( p2ScalarFunctionContainer["DiffusionFE"] ) );

   transportOperator = std::make_shared< MMOCTransport< ScalarFunction > >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, TimeSteppingScheme::RK4 );

   transportSolverTALA = std::make_shared< CGSolver< P2TransportIcosahedralShellMapOperator > >(
       storage,
       TN.domainParameters.minLevel,
       TN.domainParameters.maxLevel,
       TN.solverParameters.diffusionMaxNumIterations,
       TN.solverParameters.diffusionAbsoluteResidualUTolerance );

   transportSolverTALA->setPrintInfo( true );

   diffusionSolver = std::make_shared< CGSolver< DiffusionOperator > >( storage,
                                                                        TN.domainParameters.minLevel,
                                                                        TN.domainParameters.maxLevel,
                                                                        TN.solverParameters.diffusionMaxNumIterations,
                                                                        TN.solverParameters.diffusionAbsoluteResidualUTolerance );

   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup solvers & operators: Finished -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

} // namespace terraneo