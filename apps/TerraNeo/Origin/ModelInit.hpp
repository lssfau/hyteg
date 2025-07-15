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

   // printConfig( TN );

#ifdef HYTEG_BUILD_WITH_ADIOS2
   attributeList_["rCMB"]     = TN.domainParameters.rCMB;
   attributeList_["rSurface"] = TN.domainParameters.rSurface;
   attributeList_["nTan"]     = TN.domainParameters.nTan;
   attributeList_["nRad"]     = TN.domainParameters.nRad;
   attributeList_["minLevel"] = TN.domainParameters.minLevel;
   attributeList_["maxLevel"] = TN.domainParameters.maxLevel;

   attributeList_["simulationType"]             = TN.simulationParameters.simulationType;
   attributeList_["adaptiveRefTemp"]            = TN.simulationParameters.adaptiveRefTemp;
   attributeList_["tempDependentViscosity"]     = TN.simulationParameters.tempDependentViscosity;
   attributeList_["tempDependentViscosityType"] = TN.simulationParameters.tempDependentViscosityType;

   attributeList_["fnameTopologies"]                = TN.simulationParameters.fnameTopologies;
   attributeList_["fnameReconstructions"]           = TN.simulationParameters.fnameReconstructions;
   attributeList_["plateVelocityScaling"]           = TN.simulationParameters.plateVelocityScaling;
   attributeList_["plateSmoothingDistance"]         = TN.simulationParameters.plateSmoothingDistance;
   attributeList_["compressible"]                   = TN.simulationParameters.compressible;
   attributeList_["shearHeating"]                   = TN.simulationParameters.shearHeating;
   attributeList_["adiabaticHeating"]               = TN.simulationParameters.adiabaticHeating;
   attributeList_["internalHeating"]                = TN.simulationParameters.internalHeating;
   attributeList_["boundaryCond"]                   = TN.simulationParameters.boundaryCond;
   attributeList_["boundaryCondFreeSlip"]           = TN.simulationParameters.boundaryCondFreeSlip;
   attributeList_["haveViscosityProfile"]           = TN.simulationParameters.haveViscosityProfile;
   attributeList_["fileViscosityProfile"]           = TN.simulationParameters.fileViscosityProfile;
   attributeList_["lithosphereShearHeatingScaling"] = TN.simulationParameters.lithosphereShearHeatingScaling;
   attributeList_["lithosphereThickness"]           = TN.simulationParameters.lithosphereThickness;

   attributeList_["surfaceTemp"]          = TN.physicalParameters.surfaceTemp;
   attributeList_["cmbTemp"]              = TN.physicalParameters.cmbTemp;
   attributeList_["thermalExpansivity"]   = TN.physicalParameters.thermalExpansivity;
   attributeList_["thermalConductivity"]  = TN.physicalParameters.thermalConductivity;
   attributeList_["specificHeatCapacity"] = TN.physicalParameters.specificHeatCapacity;
   attributeList_["internalHeatingRate"]  = TN.physicalParameters.internalHeatingRate;
   attributeList_["referenceDensity"]     = TN.physicalParameters.referenceDensity;
   attributeList_["surfaceDensity"]       = TN.physicalParameters.surfaceDensity;
   attributeList_["referenceViscosity"]   = TN.physicalParameters.referenceViscosity;
   attributeList_["viscosity"]            = TN.physicalParameters.viscosity;
   attributeList_["grueneisenParameter"]  = TN.physicalParameters.grueneisenParameter;
   attributeList_["adiabatSurfaceTemp"]   = TN.physicalParameters.adiabatSurfaceTemp;
   attributeList_["activationEnergy"]     = TN.physicalParameters.activationEnergy;
   attributeList_["depthViscosityFactor"] = TN.physicalParameters.depthViscosityFactor;
   attributeList_["viscosityLowerBound"]  = TN.physicalParameters.viscosityLowerBound;
   attributeList_["viscosityUpperBound"]  = TN.physicalParameters.viscosityUpperBound;
#endif

   setupBoundaryConditions();
   setupFunctions();
   initialiseFunctions();
   setupSolversAndOperators();
   setupOutput();
   if ( TN.outputParameters.createTimingDB )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Create SQL database for runtime analysis" );
      initTimingDB();
   }

   printConfig( TN, walberla::format( "%s/params.txt", modelPath.c_str() ) );

   std::string logFilename = walberla::format( "%s/logFile.txt", modelPath.c_str() );

   WALBERLA_ROOT_SECTION()
   {
      walberla::logging::Logging::instance()->includeLoggingToFile( logFilename );
   }
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

   storage          = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );
   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------" )
   WALBERLA_LOG_INFO_ON_ROOT( "--- PrimitiveStorage Info ----" )
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------" )
   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );
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
      bcVelocityR.createDirichletBC( "RadialDirichlet",
                                     { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
      bcVelocityThetaPhi.createDirichletBC(
          "TangentialDirichlet", { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );

      idSurface = bcVelocity.createDirichletBC( "surface", { MeshInfo::hollowFlag::flagOuterBoundary } );
      idCMB     = bcVelocity.createDirichletBC( "cmb", { MeshInfo::hollowFlag::flagInnerBoundary } );
      break;

   case 2:
      bcVelocityR.createDirichletBC( "RadialDirichlet",
                                     { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
      bcVelocityThetaPhi.createNeumannBC( "TangentialNeumann",
                                          { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );

      idSurface = bcVelocity.createFreeslipBC( "surface", { MeshInfo::hollowFlag::flagOuterBoundary } );
      idCMB     = bcVelocity.createFreeslipBC( "cmb", { MeshInfo::hollowFlag::flagInnerBoundary } );
      TN.simulationParameters.boundaryCondFreeSlip = true;
      break;

   case 3:
      bcVelocityR.createDirichletBC( "RadialDirichlet",
                                     { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
      bcVelocityThetaPhi.createNeumannBC( "NeumannInner", { MeshInfo::hollowFlag::flagInnerBoundary } );
      bcVelocityThetaPhi.createDirichletBC( "DirichletOuter", { MeshInfo::hollowFlag::flagOuterBoundary } );

      idSurface = bcVelocity.createDirichletBC( "surface", { MeshInfo::hollowFlag::flagOuterBoundary } );
      idCMB     = bcVelocity.createFreeslipBC( "cmb", { MeshInfo::hollowFlag::flagInnerBoundary } );
      TN.simulationParameters.boundaryCondFreeSlip = true;
      break;

   default:
      bcVelocityR.createDirichletBC( "RadialDirichlet",
                                     { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
      bcVelocityThetaPhi.createDirichletBC(
          "TangentialDirichlet", { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );

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

   for ( auto& [p0Names, p0MinLevel, p0MaxLevel, p0BcType] : p0ScalarFunctionDict )
   {
      p0MinLevel = TN.domainParameters.minLevel;
      p0MaxLevel = TN.domainParameters.maxLevel;
   }

   for ( auto [p0Names, p0MinLevel, p0MaxLevel, p0BcType] : p0ScalarFunctionDict )
   {
      switch ( p0BcType )
      {
      case BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION: {
         WALBERLA_LOG_WARNING_ON_ROOT( "bcVelocity looks odd for P2Function" );
         p0ScalarFunctionContainer.emplace( std::make_pair(
             p0Names, std::make_shared< ScalarFunctionP0 >( p0Names, storage, p0MinLevel, p0MaxLevel, bcVelocity ) ) );
      }
      break;
      case BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION: {
         p0ScalarFunctionContainer.emplace( std::make_pair(
             p0Names, std::make_shared< ScalarFunctionP0 >( p0Names, storage, p0MinLevel, p0MaxLevel, bcTemperature ) ) );
      }
      break;
      default: {
         p0ScalarFunctionContainer.emplace(
             std::make_pair( p0Names, std::make_shared< ScalarFunctionP0 >( p0Names, storage, p0MinLevel, p0MaxLevel ) ) );
      }
      }
   }

   p2p1StokesFunctionContainer["VelocityFERotated"]->uvw().component( 0u ).setBoundaryCondition( bcVelocityThetaPhi );
   p2p1StokesFunctionContainer["VelocityFERotated"]->uvw().component( 1u ).setBoundaryCondition( bcVelocityThetaPhi );
   p2p1StokesFunctionContainer["VelocityFERotated"]->uvw().component( 2u ).setBoundaryCondition( bcVelocityR );

   p2p1StokesFunctionContainer["StokesRHSRotated"]->uvw().component( 0u ).setBoundaryCondition( bcVelocityThetaPhi );
   p2p1StokesFunctionContainer["StokesRHSRotated"]->uvw().component( 1u ).setBoundaryCondition( bcVelocityThetaPhi );
   p2p1StokesFunctionContainer["StokesRHSRotated"]->uvw().component( 2u ).setBoundaryCondition( bcVelocityR );
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

   temperatureInitParams = std::make_shared< TemperatureInitializationParameters >( TN.physicalParameters.cmbTemp,
                                                                                    TN.physicalParameters.surfaceTemp,
                                                                                    TN.physicalParameters.adiabatSurfaceTemp,
                                                                                    TN.physicalParameters.dissipationNumber,
                                                                                    TN.domainParameters.rMin,
                                                                                    TN.domainParameters.rMax,
                                                                                    &TN.initialisationParameters );

   if ( TN.simulationParameters.haveTemperatureProfile )
   {
      temperatureReferenceFct =
          std::make_shared< std::function< real_t( const Point3D& ) > >( terraneo::temperatureReferenceProfile(
              *temperatureInitParams, TN.physicalParameters.radiusT, TN.physicalParameters.temperatureInputProfile ) );
   }
   else
   {
      temperatureReferenceFct = std::make_shared< std::function< real_t( const Point3D& ) > >(
          terraneo::temperatureReferenceExponential( *temperatureInitParams ) );
   }

   // Temperature Deivation
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
   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      p2ScalarFunctionContainer["TemperatureFE"]->interpolate( initTemperature, l, All );
   }

   if ( TN.simulationParameters.simulationType == "CirculationModel" )
   {
      // initialise plate velocity oracle
      WALBERLA_LOG_INFO_ON_ROOT( "Setup Oracle for Plates" );
      terraneo::oracle = std::make_shared< terraneo::plates::PlateVelocityProvider >(
          TN.simulationParameters.fnameTopologies, TN.simulationParameters.fnameReconstructions );
      // Plate averaging object
      std::vector< std::pair< real_t, int > > weightPairs = { { 1.0 / 100.0, 6 } };
      terraneo::avgPointProvider = std::make_shared< terraneo::plates::UniformCirclesPointWeightProvider >( weightPairs, 1e-1 );
   }
   // Assign temperature field to temperaturePrev

   for ( uint_t level = TN.domainParameters.minLevel; level <= TN.domainParameters.maxLevel; ++level )
   {
      p2ScalarFunctionContainer["TemperaturePrev"]->assign(
          { real_c( 1 ) }, { *( p2ScalarFunctionContainer["TemperatureFE"] ) }, level, All );
      p2ScalarFunctionContainer["DensityFE"]->interpolate( densityFunc, level, All );

      //set plate velocities for timestep 0 / initialAge
      if ( TN.simulationParameters.simulationType == "CirculationModel" )
      {
         updatePlateVelocities( *( p2p1StokesFunctionContainer["VelocityFE"] ) );
      }
      else
      {
         //currently initialising CMB and surface to zeros for NoSlipNoSlip case
         p2p1StokesFunctionContainer["VelocityFE"]->uvw().interpolate( { zeros, zeros, zeros }, level, idSurface );
         p2p1StokesFunctionContainer["VelocityFEPrev"]->uvw().interpolate( { zeros, zeros, zeros }, level, idSurface );
      }

      p2p1StokesFunctionContainer["VelocityFE"]->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      p2p1StokesFunctionContainer["VelocityFEPrev"]->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      p2p1StokesFunctionContainer["StokesRHS"]->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      p2p1StokesFunctionContainer["StokesTmp1"]->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      // for volumetric average temperature calculation
      p2ScalarFunctionContainer["TemperatureVolumetric"]->interpolate( real_c( 1 ), level, All );
      p2ScalarFunctionContainer["Volume"]->interpolate( real_c( 1 ), level, All );
   }

   auto temperatureRadialProfile            = computeRadialProfile( *( p2ScalarFunctionContainer["TemperatureFE"] ),
                                                         TN.domainParameters.rMin,
                                                         TN.domainParameters.rMax,
                                                         TN.domainParameters.nRad,
                                                         TN.domainParameters.maxLevel );
   temperatureProfiles                      = std::make_shared< RadialProfile >( temperatureRadialProfile );
   TN.physicalParameters.temperatureProfile = temperatureProfiles->mean;

   auto velocityRadialProfile = computeRadialProfile( p2p1StokesFunctionContainer["VelocityFE"]->uvw(),
                                                      TN.domainParameters.rMin,
                                                      TN.domainParameters.rMax,
                                                      TN.domainParameters.nRad,
                                                      TN.domainParameters.maxLevel );
   velocityProfiles           = std::make_shared< RadialProfile >( velocityRadialProfile );

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
   referenceTemperatureFct = [this]( const Point3D& x ) { return referenceTemperatureFunction( x ); };

   std::function< real_t( const Point3D& ) > nX = [&]( const Point3D& x ) {
      Point3D n;
      normalFunc( x, n );
      return n[0];
   };

   std::function< real_t( const Point3D& ) > nY = [&]( const Point3D& x ) {
      Point3D n;
      normalFunc( x, n );
      return n[1];
   };

   std::function< real_t( const Point3D& ) > nZ = [&]( const Point3D& x ) {
      Point3D n;
      normalFunc( x, n );
      return n[2];
   };

   for ( uint_t level = TN.domainParameters.minLevel; level <= TN.domainParameters.maxLevel; level++ )
   {
      p2VectorFunctionContainer["NormalsFS"]->interpolate( { nX, nY, nZ }, level, FreeslipBoundary );
   }
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

   rotationOperator =
       std::make_shared< P2RotationOperator >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, normalFunc_ );

   p2ProlongationOperator = std::make_shared< P2toP2QuadraticProlongation >();

   stokesOperator = std::make_shared< StokesOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *( p2ScalarFunctionContainer["ViscosityFE"] ) );

   stokesOperatorFS = std::make_shared< StokesOperatorFS >( storage,
                                                            TN.domainParameters.minLevel,
                                                            TN.domainParameters.maxLevel,
                                                            p2ScalarFunctionContainer["ViscosityFE"]->getVertexDoFFunction(),
                                                            p2ScalarFunctionContainer["ViscosityFEInv"]->getVertexDoFFunction(),
                                                            *projectionOperator,
                                                            bcVelocity,
                                                            *( p2ScalarFunctionContainer["DensityFE"] ),
                                                            TN.simulationParameters.frozenVelocity );

   if constexpr ( std::is_same_v< ViscosityFunction_T, P1Function< real_t > > )
   {
      stokesOperatorOpgen =
          std::make_shared< StokesOperatorOpgen_T >( storage,
                                                     TN.domainParameters.minLevel,
                                                     TN.domainParameters.maxLevel,
                                                     p2ScalarFunctionContainer["ViscosityFE"]->getVertexDoFFunction(),
                                                     p2ScalarFunctionContainer["ViscosityFEInv"]->getVertexDoFFunction(),
                                                     p2VectorFunctionContainer["NormalsFS"]->component( 0u ),
                                                     p2VectorFunctionContainer["NormalsFS"]->component( 1u ),
                                                     p2VectorFunctionContainer["NormalsFS"]->component( 2u ),
                                                     0.0,
                                                     *rotationOperator,
                                                     bcVelocity );
   }
   else if constexpr ( std::is_same_v< ViscosityFunction_T, P0Function< real_t > > )
   {
      stokesOperatorOpgen =
          std::make_shared< StokesOperatorOpgen_T >( storage,
                                                     TN.domainParameters.minLevel,
                                                     TN.domainParameters.maxLevel,
                                                     *( p0ScalarFunctionContainer["ViscosityFEP0"] ),
                                                     p2ScalarFunctionContainer["ViscosityFEInv"]->getVertexDoFFunction(),
                                                     p2VectorFunctionContainer["NormalsFS"]->component( 0u ),
                                                     p2VectorFunctionContainer["NormalsFS"]->component( 1u ),
                                                     p2VectorFunctionContainer["NormalsFS"]->component( 2u ),
                                                     0.0,
                                                     *rotationOperator,
                                                     bcVelocity );
   }
   else
   {
      WALBERLA_ABORT("Nope");
   }

   if ( TN.solverParameters.solverFlag == 0u )
   {
      stokesSolverClass = std::make_shared< StokesMCFGMRESSolverWithProjection< StokesOperatorFS, P2ProjectNormalOperator > >(
          storage,
          TN.domainParameters.minLevel,
          TN.domainParameters.maxLevel,
          stokesOperatorFS,
          projectionOperator,
          p2p1StokesFunctionContainer["StokesTmpProlongation"],
          p2p1StokesFunctionContainer["StokesTmp1"],
          p2p1StokesFunctionContainer["StokesTmp2"],
          TN );
   }
   else if ( TN.solverParameters.solverFlag == 1u )
   {
      stokesSolverClass = std::make_shared< StokesMCUzawaSolver< StokesOperatorFS, P2ProjectNormalOperator > >(
          storage,
          TN.domainParameters.minLevel,
          TN.domainParameters.maxLevel,
          stokesOperatorFS,
          projectionOperator,
          p2p1StokesFunctionContainer["StokesTmpProlongation"],
          p2p1StokesFunctionContainer["StokesTmp1"],
          TN );
   }
   else
   {
      WALBERLA_ABORT( "Unknown solver type" );
   }

   stokesSolverFS = stokesSolverClass->getSolver();

   stokesSolverOpgenClass =
       std::make_shared< StokesMCFGMRESSolver< StokesOperatorOpgen_T, P2ProjectNormalOperator > >(
           storage,
           TN.domainParameters.minLevel,
           TN.domainParameters.maxLevel,
           stokesOperatorOpgen,
           p2p1StokesFunctionContainer["StokesTmpProlongation"],
           p2p1StokesFunctionContainer["StokesTmp1"],
           p2p1StokesFunctionContainer["StokesTmp2"],
           TN );

   stokesSolverOpgen = stokesSolverOpgenClass->getSolver();

   P2MassOperator = std::make_shared< P2ElementwiseBlendingMassOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   frozenVelocityRHS = std::make_shared< FrozenVelocityFullOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *( p2ScalarFunctionContainer["DensityFE"] ) );

   /////////////////////////
   // Diffusion Operator //
   ////////////////////////

   transportOperatorTALA = std::make_shared< P2TransportIcosahedralShellMapOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   transportOperatorP1Coefficients = std::make_shared< P2TransportP1CoefficientsIcosahedralShellMapOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   transportOperatorTALA->setVelocity( p2p1StokesFunctionContainer["VelocityFE"] );
   transportOperatorTALA->setVelocityPrev( p2p1StokesFunctionContainer["VelocityFEPrev"] );
   transportOperatorTALA->setViscosity( p2ScalarFunctionContainer["ViscosityFE"] );
   transportOperatorTALA->setTemperature( p2ScalarFunctionContainer["TemperatureFE"] );

   transportOperatorTALA->setReferenceTemperature(
       std::make_shared< std::function< real_t( const Point3D& ) > >( referenceTemperatureFct ) );

   oppositeGravityX = [this]( const Point3D& x ) {
      Point3D n;
      oppositeGravity( x, n );

      return n[0];
   };

   oppositeGravityY = [this]( const Point3D& x ) {
      Point3D n;
      oppositeGravity( x, n );

      return n[1];
   };

   oppositeGravityZ = [this]( const Point3D& x ) {
      Point3D n;
      oppositeGravity( x, n );

      return n[2];
   };

   adiabaticCoeffFunc = [this]( const Point3D& x ) { return adiabaticCoefficientFunction( x ); };

   constEnergyCoeffFunc = [this]( const Point3D& x ) { return constantEnergyCoefficientFunction( x ); };

   surfTempCoeffFunc = [this]( const Point3D& x ) { return surfaceTempCoefficientFunction( x ); };

   transportOperatorTALA->setInvGravity( { std::make_shared< std::function< real_t( const Point3D& ) > >( oppositeGravityX ),
                                           std::make_shared< std::function< real_t( const Point3D& ) > >( oppositeGravityY ),
                                           std::make_shared< std::function< real_t( const Point3D& ) > >( oppositeGravityZ ) } );

   transportOperatorTALA->setDiffusivityCoeff( std::make_shared< std::function< real_t( const Point3D& ) > >( diffFactorFunc ) );
   transportOperatorTALA->setAdiabaticCoeff(
       std::make_shared< std::function< real_t( const Point3D& ) > >( adiabaticCoeffFunc ) );
   transportOperatorTALA->setConstEnergyCoeff(
       std::make_shared< std::function< real_t( const Point3D& ) > >( constEnergyCoeffFunc ) );
   transportOperatorTALA->setSurfTempCoeff( std::make_shared< std::function< real_t( const Point3D& ) > >( surfTempCoeffFunc ) );
   transportOperatorTALA->setShearHeatingCoeff( p2ScalarFunctionContainer["ShearHeatingTermCoeff"] );

   if ( TN.simulationParameters.compressible && !TN.simulationParameters.adiabaticHeating )
   {
      TN.simulationParameters.adiabaticHeating = true;
   }

   transportOperatorTALA->setTALADict(
       { { TransportOperatorTermKey::ADIABATIC_HEATING_TERM, TN.simulationParameters.adiabaticHeating },
         { TransportOperatorTermKey::SHEAR_HEATING_TERM, TN.simulationParameters.shearHeating },
         { TransportOperatorTermKey::INTERNAL_HEATING_TERM, TN.simulationParameters.internalHeating } } );

   transportOperatorTALA->initializeOperators();

   transportOperator = std::make_shared< MMOCTransport< ScalarFunction > >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, TimeSteppingScheme::RK4 );

   transportOperator->setP1Evaluate( false );
   transportOperator->setParticleLocalRadiusTolerance( 1e-3 );

   std::function< void( const hyteg::Point3D&, hyteg::Point3D& ) > projectPointsBack = [&]( const hyteg::Point3D& xOld,
                                                                                            hyteg::Point3D&       xNew ) {
      xNew = xOld;

      real_t r = xOld.norm();

      real_t eps = 1e-8;

      if ( r < TN.domainParameters.rMax - eps && r > TN.domainParameters.rMin + eps )
      {
         WALBERLA_ABORT( "Particle is inside the domain, but seems like neighbour search failed due to tighter tolerance" );
      }
      else if ( r > TN.domainParameters.rMax - eps )
      {
         real_t dr = r - TN.domainParameters.rMax;
         xNew      = xOld - dr * xOld / r;
      }
      else if ( r < TN.domainParameters.rMin + eps )
      {
         real_t dr = TN.domainParameters.rMin - r;
         xNew      = xOld + dr * xOld / r;
      }
      else
      {
         WALBERLA_ABORT( "Cannot be here" );
      }
   };

   transportOperator->setProjectPointsBackOutsideDomainFunction( projectPointsBack );

   transportSolverTALA = std::make_shared< CGSolver< P2TransportIcosahedralShellMapOperator > >(
       storage,
       TN.domainParameters.minLevel,
       TN.domainParameters.maxLevel,
       TN.solverParameters.diffusionMaxNumIterations,
       real_c( 0 ),
       TN.solverParameters.diffusionAbsoluteResidualUTolerance );

   transportSolverTALA->setPrintInfo( true );

   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup solvers & operators: Finished -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

} // namespace terraneo