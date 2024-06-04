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
   // referenceTemperatureFct = std::bind( &ConvectionSimulation::referenceTemperatureFunction, this, std::placeholders::_1 );
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

   temperature = std::make_shared< ScalarFunction >(
       "temperature", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   temperaturePrev = std::make_shared< ScalarFunction >(
       "temperature_Prev_tstep", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   temperatureExtrapolated = std::make_shared< ScalarFunction >(
       "temperature_extrapolated", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   temperatureDev = std::make_shared< ScalarFunction >(
       "temperature_dev", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   temperatureTmp = std::make_shared< ScalarFunction >(
       "temperatureTmp", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   temperatureReference = std::make_shared< ScalarFunction >(
       "temperatureReference", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   viscosityFE = std::make_shared< ScalarFunction >(
       "viscosity", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   viscosityFEInv = std::make_shared< ScalarFunction >(
       "viscosityInv", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   energyRHS = std::make_shared< ScalarFunction >(
       "q", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   energyRHSWeak = std::make_shared< ScalarFunction >(
       "qWeak", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   onesFE = std::make_shared< ScalarFunction >(
       "onesFE", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   densityFE = std::make_shared< ScalarFunction >(
       "rho", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   diffusionFE = std::make_shared< ScalarFunction >(
       "k", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   adiabaticTermCoeff = std::make_shared< ScalarFunction >(
       "adiabaticTermCoeff", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   shearHeatingTermCoeff = std::make_shared< ScalarFunction >(
       "shearHeatingTermCoeff", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   constEnergyCoeff = std::make_shared< ScalarFunction >(
       "constEnergyCoeff", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   surfTempCoeff = std::make_shared< ScalarFunction >(
       "surfTempCoeff", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcTemperature );

   stokesLHS =
       std::make_shared< StokesFunction >( "u", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   stokesLHSP1 = std::make_shared< P1VectorFunction< real_t > >(
       "uP1", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   stokesLHSP1Weak = std::make_shared< P1VectorFunction< real_t > >(
       "uP1Weak", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   stokesLHSExtrapolated = std::make_shared< StokesFunction >(
       "u_extrapolated", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   stokesLHSPrev = std::make_shared< StokesFunction >(
       "uPrev", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   stokesRHS =
       std::make_shared< StokesFunction >( "f", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   stokesTmp = std::make_shared< StokesFunction >(
       "fTmp", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   gradRhoOverRho = std::make_shared< P2VectorFunction< real_t > >(
       "gradRho", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   inwardNormal = std::make_shared< P2VectorFunction< real_t > >(
       "outwardNormal", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   oppositeGravityField = std::make_shared< P2VectorFunction< real_t > >(
       "oppositeGravityField", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   velocityMagnitudeSquared = std::make_shared< ScalarFunction >(
       "uMagnitudeSquared", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   scalarTmp = std::make_shared< ScalarFunction >(
       "uTmp", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );
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
         temperature->interpolate( initTemperatureWhiteNoise, l, All );
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
         temperature->interpolate( initTemperatureSPH, l, All );
      }
   }

   // Assign temperature field to temperaturePrev

   for ( uint_t level = TN.domainParameters.minLevel; level <= TN.domainParameters.maxLevel; ++level )
   {
      temperaturePrev->assign( { real_c( 1 ) }, { *temperature }, level, All );
      temperatureTmp->assign( { real_c( 1 ) }, { *temperature }, level, All );
      densityFE->interpolate( densityFunc, level, All );
      diffusionFE->interpolate( diffFactorFunc, level, All );

      //set plate velocities for timestep 0 / initialAge
      if ( TN.simulationParameters.simulationType == "CirculationModel" )
      {
         updatePlateVelocities( *stokesLHS );
      }
      else
      {
         //currently initialising CMB and surface to zeros for NoSlipNoSlip case
         stokesLHS->uvw().interpolate( { zeros, zeros, zeros }, level, idSurface );
         stokesLHSPrev->uvw().interpolate( { zeros, zeros, zeros }, level, idSurface );
      }

      stokesLHS->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      stokesLHSPrev->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      stokesRHS->uvw().interpolate( { zeros, zeros, zeros }, level, All );
      stokesTmp->uvw().interpolate( { zeros, zeros, zeros }, level, All );

      inwardNormal->interpolate( { normalX, normalY, normalZ }, level, All );

      oppositeGravityField->assign( { -1.0 }, { *inwardNormal }, level, All );

      onesFE->interpolate( real_c( 1 ), level, All );
      gradRhoOverRho->interpolate( { normalX, normalY, normalZ }, level, All );

      //grad(rho)/rho = - ( Di / gamma ) * r_hat
      gradRhoOverRho->assign( { TN.physicalParameters.dissipationNumber / TN.physicalParameters.grueneisenParameter },
                              { *gradRhoOverRho },
                              level,
                              All );
   }

   auto temperatureRadialProfile = computeRadialProfile(
       *temperature, TN.domainParameters.rMin, TN.domainParameters.rMax, TN.domainParameters.nRad, TN.domainParameters.maxLevel );
   temperatureProfiles = std::make_shared< RadialProfile >( temperatureRadialProfile );

   if ( TN.outputParameters.outputProfiles & TN.simulationParameters.tempDependentViscosity )
   {
      auto viscosityRadialProfile = computeRadialProfile( *viscosityFE,
                                                          TN.domainParameters.rMin,
                                                          TN.domainParameters.rMax,
                                                          TN.domainParameters.nRad,
                                                          TN.domainParameters.maxLevel );
      viscosityProfiles           = std::make_shared< RadialProfile >( viscosityRadialProfile );
   }

   referenceTemperatureFct = [=]( const Point3D& x ) {
      real_t radius = x.norm();
      if ( TN.simulationParameters.adaptiveRefTemp )
      {
         uint_t shell = static_cast< uint_t >(
             std::round( real_c( TN.simulationParameters.numLayers ) *
                         ( ( radius - TN.domainParameters.rMin ) / ( TN.domainParameters.rMax - TN.domainParameters.rMin ) ) ) );
         WALBERLA_ASSERT( shell < T.size() );
         return temperatureProfiles->mean.at( shell );
      }
      else
      {
         return referenceTemperatureFunction( x );
      }
   };

   temperatureReference->interpolate( referenceTemperatureFct, TN.domainParameters.maxLevel, All );
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
   updateRefViscosity();

   P1Function< real_t > muInv( "muInv", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   for ( uint_t l = TN.domainParameters.minLevel; l <= TN.domainParameters.maxLevel; l++ )
   {
      muInv.interpolate( 1.0 / TN.physicalParameters.referenceViscosity, l, All );
   }

   auto viscInvTest = viscosityFEInv->getVertexDoFFunction();

   projectionOperator = std::make_shared< P2ProjectNormalOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, normalFunc_ );

   stokesOperator =
       std::make_shared< StokesOperator >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *viscosityFE );

   stokesOperatorFS = std::make_shared< StokesOperatorFS >( storage,
                                                            TN.domainParameters.minLevel,
                                                            TN.domainParameters.maxLevel,
                                                            *viscosityFE,
                                                            viscosityFEInv->getVertexDoFFunction(),
                                                            *projectionOperator,
                                                            bcVelocity );

   stokesSolverFS = temporary::stokesGMGFSSolver< StokesOperatorFS, P2ProjectNormalOperator >(
       storage,
       TN.domainParameters.minLevel,
       TN.domainParameters.maxLevel,
       stokesOperatorFS,
       projectionOperator,
       TN.solverParameters.stokesMaxNumIterations,
       TN.solverParameters.stokesRelativeResidualUTolerance,
       0.3,
       bcVelocity );

   P2MassOperator = std::make_shared< P2ElementwiseBlendingMassOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
   MassOperatorVelocityP1 =
       std::make_shared< P1MassOperatorVelocity >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   frozenVelocityRHSX = std::make_shared< FrozenVelocityOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, gradRhoOverRho->component( 0U ) );
   frozenVelocityRHSY = std::make_shared< FrozenVelocityOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, gradRhoOverRho->component( 1U ) );
   frozenVelocityRHSZ = std::make_shared< FrozenVelocityOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, gradRhoOverRho->component( 2U ) );

   /////////////////////////
   // Diffusion Operator //
   ////////////////////////

   transportOperatorTALA =
       std::make_shared< P2TransportOperatorTALA >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   transportOperatorTALA->setVelocity( stokesLHS );
   transportOperatorTALA->setViscosity( viscosityFE );
   transportOperatorTALA->setTemperature( temperature );

   transportOperatorTALA->setInvGravity( oppositeGravityField );

   transportOperatorTALA->setDiffusivityCoeff( diffusionFE );
   transportOperatorTALA->setAdiabaticCoeff( adiabaticTermCoeff );
   transportOperatorTALA->setShearHeatingCoeff( shearHeatingTermCoeff );
   transportOperatorTALA->setConstEnergyCoeff( constEnergyCoeff );
   transportOperatorTALA->setSurfTempCoeff( surfTempCoeff );

   transportOperatorTALA->setReferenceTemperature( temperatureReference );

   transportOperatorTALA->setTALADict( { { OperatorTermKey::ADIABATIC_HEATING_TERM, false },
                                         { OperatorTermKey::SHEAR_HEATING_TERM, false },
                                         { OperatorTermKey::INTERNAL_HEATING_TERM, false } } );

   transportOperatorTALA->initializeOperators();

   diffusionOperator =
       std::make_shared< DiffusionOperator >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *diffusionFE );

   transportOperator = std::make_shared< MMOCTransport< ScalarFunction > >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, TimeSteppingScheme::RK4 );

   transportSolverTALA =
       std::make_shared< CGSolver< P2TransportOperatorTALA > >( storage,
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