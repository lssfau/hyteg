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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::ConvectionSimulation(
    const walberla::Config::BlockHandle& mainConf )
{
   TN = terraneo::parseConfig( mainConf );

   densityFunc    = std::bind( &ConvectionSimulation::densityFunction, this, std::placeholders::_1 );
   diffFactorFunc = std::bind( &ConvectionSimulation::diffPreFactorFunction, this, std::placeholders::_1 );
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::init()
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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::setupDomain()
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
   meshInfo          = MeshInfo::meshSphericalShell( TN.domainParameters.nTan, TN.domainParameters.macroLayers );

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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::setupBoundaryConditions()
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

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::setupFunctions()
{
   if ( !storage )
   {
      setupDomain();
   }

   auto initialiseFunctionContainers = [this]< typename FunctionType >(
                                           const std::vector< std::tuple< std::string, Level_T, Level_T, BCs_T > >& functionDict,
                                           std::map< std::string, std::shared_ptr< FunctionType > >& functionContainer ) {
      for ( auto [functionName, functionMinLevelType, functionMaxLevelType, functionBcType] : functionDict )
      {
         uint_t functionMinLevel = 0u;
         uint_t functionMaxLevel = 0u;

         if ( functionMinLevelType == Level_T::MINLEVEL )
         {
            functionMinLevel = TN.domainParameters.minLevel;
         }
         else if ( functionMinLevelType == Level_T::MAXLEVEL )
         {
            functionMinLevel = TN.domainParameters.maxLevel;
         }
         else if ( functionMinLevelType == Level_T::MAXLEVELPLUSONE )
         {
            functionMinLevel = TN.domainParameters.maxLevel + 1;
         }
         else
         {
            WALBERLA_ABORT( "Unknown Level_T" );
         }

         if ( functionMaxLevelType == Level_T::MINLEVEL )
         {
            functionMaxLevel = TN.domainParameters.minLevel;
         }
         else if ( functionMaxLevelType == Level_T::MAXLEVEL )
         {
            functionMaxLevel = TN.domainParameters.maxLevel;
         }
         else if ( functionMaxLevelType == Level_T::MAXLEVELPLUSONE )
         {
            functionMaxLevel = TN.domainParameters.maxLevel + 1;
         }
         else
         {
            WALBERLA_ABORT( "Unknown Level_T" );
         }

         switch ( functionBcType )
         {
         case BCs_T::VELOCITY_BC: {
            functionContainer.emplace( std::make_pair(
                functionName,
                std::make_shared< FunctionType >( functionName, storage, functionMinLevel, functionMaxLevel, bcVelocity ) ) );
         }
         break;
         case BCs_T::VELOCITYROTATION_BC: {
            functionContainer.emplace( std::make_pair(
                functionName,
                std::make_shared< FunctionType >( functionName, storage, functionMinLevel, functionMaxLevel, bcVelocity ) ) );

            if constexpr ( std::is_same_v< FunctionType, P2P1StokesFunction_T > )
            {
               functionContainer.at( functionName )->uvw().component( 0u ).setBoundaryCondition( bcVelocityThetaPhi );
               functionContainer.at( functionName )->uvw().component( 1u ).setBoundaryCondition( bcVelocityThetaPhi );
               functionContainer.at( functionName )->uvw().component( 2u ).setBoundaryCondition( bcVelocityR );
            }
            else if constexpr ( std::is_same_v< FunctionType, P2VectorFunction_T > )
            {
               functionContainer.at( functionName )->component( 0u ).setBoundaryCondition( bcVelocityThetaPhi );
               functionContainer.at( functionName )->component( 1u ).setBoundaryCondition( bcVelocityThetaPhi );
               functionContainer.at( functionName )->component( 2u ).setBoundaryCondition( bcVelocityR );
            }
            else
            {
               WALBERLA_ABORT( "VELOCITYROTATION BC does not apply for this function type" );
            }
         }
         break;
         case BCs_T::TEMPERATURE_BC: {
            functionContainer.emplace( std::make_pair(
                functionName,
                std::make_shared< FunctionType >( functionName, storage, functionMinLevel, functionMaxLevel, bcTemperature ) ) );
         }
         break;
         default: {
            functionContainer.emplace( std::make_pair(
                functionName, std::make_shared< FunctionType >( functionName, storage, functionMinLevel, functionMaxLevel ) ) );
         }
         }
      }
   };

   initialiseFunctionContainers( p2p1StokesFunctionDict, p2p1StokesFunctionContainer );
   initialiseFunctionContainers( p2ScalarFunctionDict, p2ScalarFunctionContainer );
   initialiseFunctionContainers( p2VectorFunctionDict, p2VectorFunctionContainer );
   initialiseFunctionContainers( p1ScalarFunctionDict, p1ScalarFunctionContainer );
   initialiseFunctionContainers( p1VectorFunctionDict, p1VectorFunctionContainer );
   initialiseFunctionContainers( p0ScalarFunctionDict, p0ScalarFunctionContainer );
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::initialiseFunctions()
{
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup initial state -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   walberla::math::seedRandomGenerator( 42 );
   randFunc = []( const Point3D& ) { return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) ); };

   temperatureInitParams = std::make_shared< TemperatureInitializationParameters >( TN.physicalParameters.cmbTemp,
                                                                                    TN.physicalParameters.surfaceTemp,
                                                                                    TN.physicalParameters.adiabatSurfaceTemp,
                                                                                    TN.physicalParameters.dissipationNumber,
                                                                                    TN.domainParameters.rMin,
                                                                                    TN.domainParameters.rMax,
                                                                                    &TN.initialisationParameters );

   referenceTemperatureFct = [this]( const Point3D& x ) { return referenceTemperatureFunction( x ); };

   // Temperature Deivation
   switch ( TN.initialisationParameters.initialTemperatureDeviationMethod )
   {
   case INITIAL_TEMPERATURE_DEVIATION_METHOD::SINGLE_SPH:
      temperatureInitialCondition = temperatureSingleSPH( *temperatureInitParams, referenceTemperatureFct );
      break;
   case INITIAL_TEMPERATURE_DEVIATION_METHOD::RANDOM_SUPERPOSITION_SPH:
      temperatureInitialCondition = temperatureRandomSuperpositioneSPH( *temperatureInitParams, referenceTemperatureFct );
      break;
   case INITIAL_TEMPERATURE_DEVIATION_METHOD::WHITE_NOISE:
      temperatureInitialCondition = temperatureWhiteNoise( *temperatureInitParams, referenceTemperatureFct );
      break;
   default:
      WALBERLA_ABORT( "Unknown initial temperature deviation method" );
   }

   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureRHS    = p2p1StokesFunctionContainer.at( "StokesRHS" );
   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureTemp   = p2p1StokesFunctionContainer.at( "StokesTmp1" );
   std::shared_ptr< P2P1StokesFunction_T >& velocityPressureFE     = p2p1StokesFunctionContainer.at( "VelocityFE" );
   std::shared_ptr< P2P1StokesFunction_T >& velocityPressurePrevFE = p2p1StokesFunctionContainer.at( "VelocityFEPrev" );

   std::shared_ptr< P2VectorFunction_T >& normalsFE = p2VectorFunctionContainer.at( "NormalsFS" );

   std::shared_ptr< P2ScalarFunction_T >& viscosityP2       = p2ScalarFunctionContainer.at( "ViscosityFE" );
   std::shared_ptr< P2ScalarFunction_T >& temperatureP2     = p2ScalarFunctionContainer.at( "TemperatureFE" );
   std::shared_ptr< P2ScalarFunction_T >& temperaturePrevP2 = p2ScalarFunctionContainer.at( "TemperaturePrev" );
   std::shared_ptr< P2ScalarFunction_T >& densityFE         = p2ScalarFunctionContainer.at( "DensityFE" );

   std::shared_ptr< P1ScalarFunction_T >& temperatureP1     = p1ScalarFunctionContainer.at( "TemperatureFEP1" );
   std::shared_ptr< P1ScalarFunction_T >& temperaturePrevP1 = p1ScalarFunctionContainer.at( "TemperatureFEPrevP1" );

   temperatureP2->interpolate( temperatureInitialCondition, TN.domainParameters.maxLevel, All );
   temperatureP1->interpolate( temperatureInitialCondition, TN.domainParameters.maxLevel + 1, All );

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

   // for ( uint_t level = TN.domainParameters.minLevel; level <= TN.domainParameters.maxLevel; ++level )
   // {
   temperaturePrevP2->assign( { real_c( 1 ) }, { *( temperatureP2 ) }, TN.domainParameters.maxLevel, All );
   temperaturePrevP1->assign( { real_c( 1 ) }, { *( temperatureP1 ) }, TN.domainParameters.maxLevel, All );

   if ( TN.simulationParameters.simulationType == "CirculationModel" )
   {
      updatePlateVelocities( *( velocityPressureFE ) );
   }
   else
   {
      velocityPressureFE->uvw().interpolate( 0.0, TN.domainParameters.maxLevel, idSurface );
      velocityPressurePrevFE->uvw().interpolate( 0.0, TN.domainParameters.maxLevel, idSurface );
   }

   velocityPressureFE->uvw().interpolate( 0.0, TN.domainParameters.maxLevel, All );
   velocityPressurePrevFE->uvw().interpolate( 0.0, TN.domainParameters.maxLevel, All );
   velocityPressureRHS->uvw().interpolate( 0.0, TN.domainParameters.maxLevel, All );
   velocityPressureTemp->uvw().interpolate( 0.0, TN.domainParameters.maxLevel, All );
   // }

   auto temperatureRadialProfile            = computeRadialProfile( *( temperatureP2 ),
                                                         TN.domainParameters.rMin,
                                                         TN.domainParameters.rMax,
                                                         TN.domainParameters.macroLayers,
                                                         TN.domainParameters.maxLevel );
   temperatureProfiles                      = std::make_shared< RadialProfile >( temperatureRadialProfile );
   TN.physicalParameters.temperatureProfile = temperatureProfiles->mean;

   auto velocityRadialProfile = computeRadialProfile( velocityPressureFE->uvw(),
                                                      TN.domainParameters.rMin,
                                                      TN.domainParameters.rMax,
                                                      TN.domainParameters.macroLayers,
                                                      TN.domainParameters.maxLevel );
   velocityProfiles           = std::make_shared< RadialProfile >( velocityRadialProfile );

   if ( TN.outputParameters.outputProfiles && TN.simulationParameters.tempDependentViscosity )
   {
      updateViscosity();
      auto viscosityRadialProfile = computeRadialProfile( *( viscosityP2 ),
                                                          TN.domainParameters.rMin,
                                                          TN.domainParameters.rMax,
                                                          TN.domainParameters.macroLayers,
                                                          TN.domainParameters.maxLevel );
      viscosityProfiles           = std::make_shared< RadialProfile >( viscosityRadialProfile );
   }

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
      densityFE->interpolate( densityFunc, level, All );
      normalsFE->interpolate( { nX, nY, nZ }, level, FreeslipBoundary );
   }
}

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
void ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T >::setupSolversAndOperators()
{
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup solvers & operators -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   std::shared_ptr< P2P1StokesFunction_T >& stokesTmpProlongation = p2p1StokesFunctionContainer.at( "StokesTmpProlongation" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesTmp1            = p2p1StokesFunctionContainer.at( "StokesTmp1" );
   std::shared_ptr< P2P1StokesFunction_T >& stokesTmp2            = p2p1StokesFunctionContainer.at( "StokesTmp2" );

   std::shared_ptr< P2VectorFunction_T >& normalsFE = p2VectorFunctionContainer.at( "NormalsFS" );

   std::shared_ptr< P2ScalarFunction_T >& densityP2      = p2ScalarFunctionContainer.at( "DensityFE" );
   std::shared_ptr< P2ScalarFunction_T >& viscosityP2    = p2ScalarFunctionContainer.at( "ViscosityFE" );
   std::shared_ptr< P2ScalarFunction_T >& viscosityP2Inv = p2ScalarFunctionContainer.at( "ViscosityFEInv" );

   std::shared_ptr< P0ScalarFunction_T >& viscosityP0 = p0ScalarFunctionContainer.at( "ViscosityFEP0" );

   auto normalFunc_ = [&]( const Point3D& p, Point3D& n ) { normalFunc( p, n ); };

   updateViscosity();

   projectionOperator_ = std::make_shared< P2ProjectNormalOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, normalFunc_ );

   rotationOperator_ =
       std::make_shared< P2RotationOperator >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, normalFunc_ );

   p2ProlongationOperator_ = std::make_shared< P2toP2QuadraticProlongation >();

   if constexpr ( std::is_same_v< ViscosityFunction_T, P1Function< real_t > > )
   {
      stokesOperator_ = std::make_shared< StokesOperator_T >( storage,
                                                              TN.domainParameters.minLevel,
                                                              TN.domainParameters.maxLevel,
                                                              viscosityP2->getVertexDoFFunction(),
                                                              viscosityP2Inv->getVertexDoFFunction(),
                                                              *projectionOperator_,
                                                              bcVelocity,
                                                              *( densityP2 ),
                                                              false,
                                                              nullptr,
                                                              nullptr,
                                                              nullptr );

      stokesOperatorRotation_ = std::make_shared< StokesOperatorRotation_T >( storage,
                                                                              TN.domainParameters.minLevel,
                                                                              TN.domainParameters.maxLevel,
                                                                              viscosityP2->getVertexDoFFunction(),
                                                                              viscosityP2Inv->getVertexDoFFunction(),
                                                                              *rotationOperator_,
                                                                              bcVelocity,
                                                                              *( densityP2 ),
                                                                              false,
                                                                              nullptr,
                                                                              nullptr,
                                                                              nullptr );
   }
   else if constexpr ( std::is_same_v< ViscosityFunction_T, P0Function< real_t > > )
   {
      stokesOperator_ = std::make_shared< StokesOperator_T >( storage,
                                                              TN.domainParameters.minLevel,
                                                              TN.domainParameters.maxLevel,
                                                              *( viscosityP0 ),
                                                              viscosityP2Inv->getVertexDoFFunction(),
                                                              *projectionOperator_,
                                                              bcVelocity,
                                                              *( densityP2 ),
                                                              false,
                                                              nullptr,
                                                              nullptr,
                                                              nullptr );

      stokesOperatorRotation_ = std::make_shared< StokesOperatorRotation_T >( storage,
                                                                              TN.domainParameters.minLevel,
                                                                              TN.domainParameters.maxLevel,
                                                                              *( viscosityP0 ),
                                                                              viscosityP2Inv->getVertexDoFFunction(),
                                                                              *rotationOperator_,
                                                                              bcVelocity,
                                                                              *( densityP2 ),
                                                                              false,
                                                                              nullptr,
                                                                              nullptr,
                                                                              nullptr );
   }
   else
   {
      WALBERLA_ABORT( "Nope" );
   }

   if ( TN.solverParameters.solverFlag == 0u )
   {
      stokesMCSolver_ = std::make_shared< StokesMCFGMRESSolver< StokesOperator_T > >( storage,
                                                                                      TN.domainParameters.minLevel,
                                                                                      TN.domainParameters.maxLevel,
                                                                                      stokesOperator_,
                                                                                      stokesTmpProlongation,
                                                                                      stokesTmp1,
                                                                                      stokesTmp2,
                                                                                      TN );

      stokesRotationMCSolver_ =
          std::make_shared< StokesMCFGMRESSolver< StokesOperatorRotation_T > >( storage,
                                                                                TN.domainParameters.minLevel,
                                                                                TN.domainParameters.maxLevel,
                                                                                stokesOperatorRotation_,
                                                                                stokesTmpProlongation,
                                                                                stokesTmp1,
                                                                                stokesTmp2,
                                                                                TN );
   }
   else if ( TN.solverParameters.solverFlag == 1u )
   {
      stokesMCSolver_ =
          std::make_shared< StokesMCUzawaSolver< StokesOperator_T, ProjectionOperator_T > >( storage,
                                                                                             TN.domainParameters.minLevel,
                                                                                             TN.domainParameters.maxLevel,
                                                                                             stokesOperator_,
                                                                                             projectionOperator_,
                                                                                             stokesTmpProlongation,
                                                                                             stokesTmp1,
                                                                                             TN );

      WALBERLA_LOG_WARNING_ON_ROOT( "Uzawa is not tested/initialized for all operators, so the application can crash" );
   }
   else if ( TN.solverParameters.solverFlag == 2u )
   {
#ifdef HYTEG_BUILD_WITH_PETSC
      stokesMCSolver_ = std::make_shared< MCSolverBase< StokesOperator_T > >(
          storage,
          TN.domainParameters.minLevel,
          TN.domainParameters.maxLevel,
          std::make_shared< PETScLUSolver< StokesOperator_T > >( storage, TN.domainParameters.maxLevel ) );

      stokesRotationMCSolver_ = std::make_shared< MCSolverBase< StokesOperatorRotation_T > >(
          storage,
          TN.domainParameters.minLevel,
          TN.domainParameters.maxLevel,
          std::make_shared< PETScLUSolver< StokesOperatorRotation_T > >( storage, TN.domainParameters.maxLevel ) );
#else
      WALBERLA_ABORT( "PETSc solver requested but not built" );
#endif
   }
   else
   {
      WALBERLA_ABORT( "Unknown solver type" );
   }

   stokesSolver_         = stokesMCSolver_->getSolver();
   stokesRotationSolver_ = stokesRotationMCSolver_->getSolver();

   p2ScalarMassOperator_ =
       std::make_shared< P2ScalarMassOperator_T >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   frozenVelocityRHS_ = std::make_shared< FrozenVelocityFullOperator_T >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *( densityP2 ) );

   ////////////////////////
   // Diffusion Operator //
   ////////////////////////

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

   adiabaticCoeffFunc   = [this]( const Point3D& x ) { return adiabaticCoefficientFunction( x ); };
   constEnergyCoeffFunc = [this]( const Point3D& x ) { return constantEnergyCoefficientFunction( x ); };
   surfTempCoeffFunc    = [this]( const Point3D& x ) { return surfaceTempCoefficientFunction( x ); };

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

      const std::shared_ptr< TemperatureFunction_T >& viscosityFE = [&] {
         if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
         {
            return p2ScalarFunctionContainer.at( "ViscosityFE" );
         }
         else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
         {
            return p1ScalarFunctionContainer.at( "ViscosityFEP1" );
         }
         else
         {
            WALBERLA_ABORT( "Unknown type" );
         }
      }();

      const std::shared_ptr< TemperatureFunction_T >& shearHeatingCoeff = [&] {
         if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
         {
            return p2ScalarFunctionContainer.at( "ShearHeatingTermCoeff" );
         }
         else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
         {
            return p1ScalarFunctionContainer.at( "ShearHeatingTermCoeffP1" );
         }
         else
         {
            WALBERLA_ABORT( "Unknown type" );
         }
      }();

      temperatureTransportOperator_ =
          std::make_shared< TransportOperator_T >( storage, TN.domainParameters.minLevel, temperatureMaxLevel );

      temperatureTransportOperator_->setVelocity( p2p1StokesFunctionContainer["VelocityFE"] );
      temperatureTransportOperator_->setVelocityPrev( p2p1StokesFunctionContainer["VelocityFEPrev"] );
      temperatureTransportOperator_->setViscosity( viscosityFE );
      temperatureTransportOperator_->setTemperature( temperatureFE );

      temperatureTransportOperator_->setReferenceTemperature(
          std::make_shared< std::function< real_t( const Point3D& ) > >( referenceTemperatureFct ) );

      temperatureTransportOperator_->setInvGravity(
          { std::make_shared< std::function< real_t( const Point3D& ) > >( oppositeGravityX ),
            std::make_shared< std::function< real_t( const Point3D& ) > >( oppositeGravityY ),
            std::make_shared< std::function< real_t( const Point3D& ) > >( oppositeGravityZ ) } );

      temperatureTransportOperator_->setDiffusivityCoeff(
          std::make_shared< std::function< real_t( const Point3D& ) > >( diffFactorFunc ) );
      temperatureTransportOperator_->setAdiabaticCoeff(
          std::make_shared< std::function< real_t( const Point3D& ) > >( adiabaticCoeffFunc ) );
      temperatureTransportOperator_->setConstEnergyCoeff(
          std::make_shared< std::function< real_t( const Point3D& ) > >( constEnergyCoeffFunc ) );
      temperatureTransportOperator_->setSurfTempCoeff(
          std::make_shared< std::function< real_t( const Point3D& ) > >( surfTempCoeffFunc ) );
      temperatureTransportOperator_->setShearHeatingCoeff( shearHeatingCoeff );
      temperatureTransportOperator_->setTALADict(
          { { TransportOperatorTermKey::ADIABATIC_HEATING_TERM, TN.simulationParameters.compressible },
            { TransportOperatorTermKey::SHEAR_HEATING_TERM, TN.simulationParameters.shearHeating },
            { TransportOperatorTermKey::INTERNAL_HEATING_TERM, TN.simulationParameters.internalHeating } } );

      temperatureTransportOperator_->initializeOperators();

      temperatureMMOCOperator_ = std::make_shared< MMOCTransport< TemperatureFunction_T > >(
          storage, TN.domainParameters.minLevel, temperatureMaxLevel, TimeSteppingScheme::RK4 );
      temperatureMMOCOperator_->setCautionedEvaluate( TN.simulationParameters.cautionedEvaluate );

      temperatureTransportSolver_ =
          std::make_shared< CGSolver< TransportOperator_T > >( storage,
                                                               TN.domainParameters.minLevel,
                                                               temperatureMaxLevel,
                                                               TN.solverParameters.diffusionMaxNumIterations,
                                                               TN.solverParameters.diffusionAbsoluteResidualUTolerance );

      temperatureTransportSolver_->setPrintInfo( true );
   }

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

   temperatureMMOCOperator_->setParticleLocalRadiusTolerance( TN.simulationParameters.particleLocationRadius );

   if ( TN.simulationParameters.projectPointsBackToDomain )
   {
      temperatureMMOCOperator_->setProjectPointsBackOutsideDomainFunction( projectPointsBack );
   }

   shearHeatingCoeffCalc = [this]( const Point3D& x, const std::vector< real_t >& density ) {
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

   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup solvers & operators: Finished -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

} // namespace terraneo
