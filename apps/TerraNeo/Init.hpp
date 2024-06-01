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

   gradRhoOverRho = std::make_shared< P1VectorFunction< real_t > >(
       "gradRho", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

   inwardNormal = std::make_shared< P2VectorFunction< real_t > >(
       "outwardNormal", storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, bcVelocity );

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

   using SubstAType = BlockLaplaceOperator;

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

   auto prolongationOperator = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto restrictionOperator  = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );

   stokesOperator =
       std::make_shared< StokesOperator >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *viscosityFE );

   stokesOperatorFS = std::make_shared< StokesOperatorFS >( storage,
                                                            TN.domainParameters.minLevel,
                                                            TN.domainParameters.maxLevel,
                                                            *viscosityFE,
                                                            viscosityFEInv->getVertexDoFFunction(),
                                                            *projectionOperator,
                                                            bcVelocity );

   stokesSolverFS = stokesGMGFSSolver( storage,
                                       TN.domainParameters.minLevel,
                                       TN.domainParameters.maxLevel,
                                       stokesOperatorFS,
                                       projectionOperator,
                                       25,
                                       0.3,
                                       bcVelocity );

   schurOperator =
       std::make_shared< SchurOperator >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, viscInvTest );
   auto APrecOperator = std::make_shared< SubstAType >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
   APrecOperator->computeInverseDiagonalOperatorValues();

   auto ABlockProlongationOperator = std::make_shared< P2toP2QuadraticVectorProlongation >();
   auto ABlockRestrictionOperator  = std::make_shared< P2toP2QuadraticVectorRestriction >();
   auto ABlockCoarseGridSolver     = std::make_shared< CGSolver< SubstAType > >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, 10, 1e-8 );
   auto ABlockSmoother =
       std::make_shared< ChebyshevSmoother< SubstAType > >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   //temporarily alter U so so it is not in the kernel of stokesOperator, required for spectral radius estimation
   stokesTmp->uvw().interpolate( randFunc, TN.domainParameters.maxLevel, All );

   auto spectralRadiusA = chebyshev::estimateRadius( *APrecOperator,
                                                     TN.domainParameters.maxLevel,
                                                     TN.solverParameters.chebyshevIterations,
                                                     storage,
                                                     stokesTmp->uvw(),
                                                     stokesRHS->uvw() );

   WALBERLA_LOG_INFO_ON_ROOT( "------ Chebyshev spectral radius: " << spectralRadiusA << " ------" );
   WALBERLA_LOG_INFO_ON_ROOT( " " );

   ABlockSmoother->setupCoefficients( 3, spectralRadiusA );

   //reset plates

   if ( TN.simulationParameters.simulationType == "CirculationModel" )
   {
      updatePlateVelocities( *stokesLHS );
   }

   auto ABlockMultigridSolver = std::make_shared< GeometricMultigridSolver< SubstAType > >( storage,
                                                                                            ABlockSmoother,
                                                                                            ABlockCoarseGridSolver,
                                                                                            ABlockRestrictionOperator,
                                                                                            ABlockProlongationOperator,
                                                                                            TN.domainParameters.minLevel,
                                                                                            TN.domainParameters.maxLevel,
                                                                                            TN.solverParameters.uzawaPreSmooth,
                                                                                            TN.solverParameters.uzawaPostSmooth,
                                                                                            0,
                                                                                            CycleType::VCYCLE );

   auto SubstABlockSolver =
       std::make_shared< SubstitutePreconditioner< typename StokesOperator::ViscousOperator_T, SubstAType > >(
           ABlockMultigridSolver, APrecOperator );

   auto ABlockSolver = std::make_shared< CGSolver< typename StokesOperator::ViscousOperator_T > >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, 3, 1e-2, SubstABlockSolver );

   auto SchurSolver = std::make_shared< CGSolver< SchurOperator > >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, 50, 1e-8 );

   uint_t chebyshevIterationsPerUzawa = 1;
   uzawaSmoother =
       std::make_shared< InexactUzawaPreconditioner< StokesOperator, StokesOperator::ViscousOperator_T, SchurOperator > >(
           storage,
           TN.domainParameters.minLevel,
           TN.domainParameters.maxLevel,
           *schurOperator,
           ABlockSolver,
           SchurSolver,
           1.0,
           1.0,
           chebyshevIterationsPerUzawa );

   if ( TN.solverParameters.solverType == 1 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: MINRES PETSc" );
#ifdef HYTEG_BUILD_WITH_PETSC
      PETScManager petscManager;
      coarseGridSolver =
          std::make_shared< PETScMinResSolver< StokesOperator > >( storage,
                                                                   TN.domainParameters.minLevel,
                                                                   TN.solverParameters.coarseGridRelativeResidualTolerance,
                                                                   TN.solverParameters.coarseGridAbsoluteResidualTolerance,
                                                                   TN.solverParameters.stokesMaxNumIterations );
#else
      WALBERLA_ABORT( "PETSc module not found or enabled." );
#endif
   }
   else
   {
      coarseGridSolver =
          solvertemplates::stokesMinResSolver< StokesOperator >( storage,
                                                                 TN.domainParameters.minLevel,
                                                                 TN.solverParameters.coarseGridAbsoluteResidualTolerance,
                                                                 TN.solverParameters.stokesMaxNumIterations,
                                                                 TN.simulationParameters.verbose );
   }

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >( storage,
                                                                                          uzawaSmoother,
                                                                                          coarseGridSolver,
                                                                                          restrictionOperator,
                                                                                          prolongationOperator,
                                                                                          TN.domainParameters.minLevel,
                                                                                          TN.domainParameters.maxLevel,
                                                                                          2,
                                                                                          2,
                                                                                          2,
                                                                                          CycleType::VCYCLE );

   stokesSolver = std::make_shared< FGMRESSolver< StokesOperator > >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, 100, 50, 1e-6, 1e-6, 0, multigridSolver );
   stokesSolver->setPrintInfo( true );

   P2MassOperator = std::make_shared< P2ElementwiseBlendingMassOperator >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );
   MassOperatorVelocityP1 =
       std::make_shared< P1MassOperatorVelocity >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel );

   /////////////////////////
   // Diffusion Operator //
   ////////////////////////

   diffusionOperator =
       std::make_shared< DiffusionOperator >( storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, *diffusionFE );

   transportOperator = std::make_shared< MMOCTransport< ScalarFunction > >(
       storage, TN.domainParameters.minLevel, TN.domainParameters.maxLevel, TimeSteppingScheme::RK4 );

   diffusionSolver = std::make_shared< CGSolver< DiffusionOperator > >( storage,
                                                                        TN.domainParameters.minLevel,
                                                                        TN.domainParameters.maxLevel,
                                                                        TN.solverParameters.diffusionMaxNumIterations,
                                                                        TN.solverParameters.diffusionAbsoluteResidualUTolerance );

   diffusionSolverTest =
       std::make_shared< CGSolver< P2DiffusionOperatorWrapper > >( storage,
                                                                   TN.domainParameters.minLevel,
                                                                   TN.domainParameters.maxLevel,
                                                                   TN.solverParameters.diffusionMaxNumIterations,
                                                                   TN.solverParameters.diffusionAbsoluteResidualUTolerance );

   diffusionSolverTest->setPrintInfo( true );

   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup solvers & operators: Finished -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

} // namespace terraneo