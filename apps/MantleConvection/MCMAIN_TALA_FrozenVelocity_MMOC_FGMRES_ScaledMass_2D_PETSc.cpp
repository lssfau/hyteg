/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include "hyteg/functions/FunctionHistory.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2InjectionRestriction.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/numerictools/BDFScheme.hpp"
#include "hyteg/operators/GEMV.hpp"
#include "hyteg/p2functionspace/P2FullViscousTDependentOperator.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

#include "MantleConvection.hpp"
#include "Utility/Data/DataLoader.hpp"
#include "Utility/Density/ConstantDensity.hpp"
#include "Utility/Density/ExponentialDensity.hpp"
#include "Utility/Density/LinearDensityInterpolation.hpp"
#include "Utility/Density/PressureProfileDensityBilinearInterpolation.hpp"
#include "Utility/LHS/AdvectionDiffusionOperator.hpp"
#include "Utility/LHS/SaddlePointOperator.hpp"
#include "Utility/OperatorTools/CombinedABlockOperator.hpp"
#include "Utility/OperatorTools/OperatorTypedefs.hpp"
#include "Utility/Parameters/NondimensionalisationParameters.hpp"
#include "Utility/Pressure/ConstantPressure.hpp"
#include "Utility/Pressure/LinearPressureInterpolation.hpp"
#include "Utility/RHS/SaddlePointOperatorRHS.hpp"
#include "Utility/Solver/ABlock/ABlockAgglomeratedCGCoarseGridSolver.hpp"
#include "Utility/Solver/ABlock/ABlockCGCoarseGridSolver.hpp"
#include "Utility/Solver/ABlock/ABlockCGOuterLoopSolver.hpp"
#include "Utility/Solver/ABlock/ABlockCGVectorMassSolver.hpp"
#include "Utility/Solver/ABlock/ABlockChebyshevSmoother.hpp"
#include "Utility/Solver/ABlock/ABlockInverseDiagSolver.hpp"
#include "Utility/Solver/ABlock/ABlockMinResCoarseGridSolver.hpp"
#include "Utility/Solver/ABlock/ABlockMinResOuterLoopSolver.hpp"
#include "Utility/Solver/ABlock/ABlockMultigridSolver.hpp"
#include "Utility/Solver/ABlock/ABlockPETScCoarseGridSolver.hpp"
#include "Utility/Solver/ABlock/ABlockWeightedJacobiSolver.hpp"
#include "Utility/Solver/AdvectionDiffusion/AdvectionDiffusionCGSPDPreconditioner.hpp"
#include "Utility/Solver/AdvectionDiffusion/AdvectionDiffusionFGMRESSolver.hpp"
#include "Utility/Solver/SaddlePoint/SaddlePointAdjointInexactUzawaSmoother.hpp"
#include "Utility/Solver/SaddlePoint/SaddlePointBlockApproximationFactorisationSmoother.hpp"
#include "Utility/Solver/SaddlePoint/SaddlePointFGMRESSolver.hpp"
#include "Utility/Solver/SaddlePoint/SaddlePointInexactUzawaSmoother.hpp"
#include "Utility/Solver/SaddlePoint/SaddlePointSubstituteSolver.hpp"
#include "Utility/Solver/SaddlePoint/SaddlePointSymmetricUzawaSmoother.hpp"
#include "Utility/Solver/Schur/SchurBFBTSolver.hpp"
#include "Utility/Solver/Schur/SchurCGCoarseGridSolver.hpp"
#include "Utility/Solver/Schur/SchurCGMassSolver.hpp"
#include "Utility/Solver/Schur/SchurCGOuterLoopSolver.hpp"
#include "Utility/Solver/Schur/SchurCGPoissonNeumannSolver.hpp"
#include "Utility/Solver/Schur/SchurChebyshevSmoother.hpp"
#include "Utility/Solver/Schur/SchurInverseDiagMassSolver.hpp"
#include "Utility/Solver/Schur/SchurInverseLumpedMassSolver.hpp"
#include "Utility/Solver/Schur/SchurMultigridSolver.hpp"
#include "Utility/Solver/Schur/SchurPoissonNeumannChebyshevSmoother.hpp"
#include "Utility/Solver/Schur/SchurPoissonNeumannMultigridSolver.hpp"
#include "Utility/Temperature/ConstantTemperature.hpp"
#include "Utility/Temperature/ExponentialTemperature.hpp"
#include "Utility/Temperature/LinearTemperature.hpp"
#include "Utility/Temperature/LinearTemperatureInterpolation.hpp"
#include "Utility/Temperature/RelativeRandomTemperature.hpp"
#include "Utility/Temperature/SphericalHarmonicsTemperature.hpp"
#include "Utility/Viscosity/ConstantViscosity.hpp"
#include "Utility/Viscosity/ExponentialViscosity.hpp"
#include "Utility/Viscosity/LinearViscosityInterpolation.hpp"
#include "Utility/Viscosity/SimpleSpaceDependentViscosityProfileWithJumps.hpp"
#include "terraneo/helpers/ConvectionToolbox.hpp"

using namespace hyteg::convectionToolbox;
using namespace MantleConvection;

// Define Combined A Block type
typedef MC_ABlock_Vec_AnnulusMap                                                 NewAType;
typedef P2ElementwiseBlendingFullViscousTDependentOperator_Centroid_ScaledDivDiv OldAType;

typedef CombinedABlockOperator< NewAType, OldAType, hyteg::P2VectorFunction< real_t > > CombinedAType;

// Define Mantle Convection operator types
typedef MantleConvection::
    SaddlePointOperator< CombinedAType, MC_BTBlock_AnnulusMap, MC_BBlock_AnnulusMap, MC_Projection, MC_NoOp >
        stokesLHSType;

typedef MantleConvection::SaddlePointOperatorRHS< MC_NoOp,
                                                  MC_NoOp,
                                                  MC_TemperatureToVelocityRHS_AnnulusMap,
                                                  MC_VelocityToPressureRHS_AnnulusMap,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_Projection >
    stokesRHSType;

typedef MantleConvection::AdvectionDiffusionOperator< MC_P2Mass_AnnulusMap,
                                                      MC_NoOp,
                                                      MC_DivKGrad_AnnulusMap,
                                                      MC_DiffusionAdditional_AnnulusMap,
                                                      MC_AdiabaticHeating_AnnulusMap,

                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp >
    transportLHSType;

typedef MantleConvection::AdvectionDiffusionOperatorRHS< MC_P2Mass_AnnulusMap,
                                                         MC_P2Mass_AnnulusMap,
                                                         MC_ShearHeating_NoSurface_AnnulusMap,
                                                         MC_NoOp,

                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp >
    transportRHSType;

// scaled mass type
typedef MC_KMass_AnnulusMap invKMassType;

typedef ScaledOperator< invKMassType, typename invKMassType::srcType, typename invKMassType::dstType > scaledInvKMassType;
typedef StabilisationProjectionWrapper< invKMassType >                                                 WrappedMassType;

// define model type
using SaddlePointOperatorType_           = stokesLHSType;
using SaddlePointRHSOperatorType_        = stokesRHSType;
using AdvectionDiffusionOperatorType_    = transportLHSType;
using AdvectionDiffusionRHSOperatorType_ = transportRHSType;

typedef MantleConvectionModel< SaddlePointOperatorType_,
                               SaddlePointRHSOperatorType_,
                               AdvectionDiffusionOperatorType_,
                               AdvectionDiffusionRHSOperatorType_,
                               MC_P2Mass_AnnulusMap,
                               hyteg::P2toP2QuadraticRestriction >
    MCModel;

typedef MCModel::SaddlePointOperatorType stokesType;

// define custom operator updater
template < class SaddlePointOperatorType_         = hyteg::NoOperator,
           class InvKMassOperatorType_            = hyteg::NoOperator,
           class ABlockChebyshevSmootherType_     = hyteg::NoOperator,
           class ABlockPETScCoarseGridSolverType_ = hyteg::NoOperator >
class CustomOperatorUpdater : public OperatorUpdater
{
 public:
   CustomOperatorUpdater( const std::shared_ptr< SaddlePointOperatorType_ >&         saddlePointOp,
                          const std::shared_ptr< InvKMassOperatorType_ >&            invKMass,
                          const std::shared_ptr< ABlockChebyshevSmootherType_ >&     ABlockChebyshevSmoother,
                          const std::shared_ptr< ABlockPETScCoarseGridSolverType_ >& ABlockPETScCoarseGridSolver,
                          const std::shared_ptr< ABlockChebyshevSmootherType_ >&     ABlockCoarseCorrectedChebyshevSmoother,
                          const std::shared_ptr< MCModel >&                          _MCModel,
                          real_t                                                     _chebyshevUpdateRateMyrs,
                          std::function< void() >                                    _auxFct )
   : saddlePointOp_( saddlePointOp )
   , invKMass_( invKMass )
   , ABlockChebyshevSmoother_( ABlockChebyshevSmoother )
   , ABlockPETScCoarseGridSolver_( ABlockPETScCoarseGridSolver )
   , ABlockCoarseCorrectedChebyshevSmoother_( ABlockCoarseCorrectedChebyshevSmoother )
   , MCModel_( _MCModel )
   , chebyshevUpdateRateMyrs_( _chebyshevUpdateRateMyrs )
   , alwaysUpdateChebyshev_( _chebyshevUpdateRateMyrs <= real_c( 0 ) )
   , auxFct_( _auxFct )
   {}

   void updateBeforeSaddlePointSolve() override {}

   void updateAfterCheckPointLoad() override
   {
      if ( saddlePointOp_ != nullptr )
      {
         saddlePointOp_->computeInverseDiagonalOperatorValues();
         saddlePointOp_->getA().getOperatorPtr()->getOperatorPtr()->getCoarseOperator()->computeAndStoreLocalElementMatrices();
      }

      if ( invKMass_ != nullptr )
      {
         invKMass_->computeInverseDiagonalOperatorValues();
      }

      if ( ABlockChebyshevSmoother_ != nullptr )
      {
         ABlockChebyshevSmoother_->regenerate();
      }

      if ( ABlockPETScCoarseGridSolver_ != nullptr )
      {
         ABlockPETScCoarseGridSolver_->reassembleMatrix( saddlePointOp_->getA() );
      }

      if ( ABlockCoarseCorrectedChebyshevSmoother_ != nullptr )
      {
         ABlockCoarseCorrectedChebyshevSmoother_->regenerate();
      }

      auxFct_();

      lastAge_ = MCModel_->getCurrentAge();
   }

   void updateAfterSaddlePointSolve() override {}

   void updateAfterAdvectionDiffusionSolve() override {}

   void updateAfterViscosityRecalculation() override
   {
      currentAge_ = MCModel_->getCurrentAge();

      if ( saddlePointOp_ != nullptr )
      {
         saddlePointOp_->computeInverseDiagonalOperatorValues();
         saddlePointOp_->getA().getOperatorPtr()->getOperatorPtr()->getCoarseOperator()->computeAndStoreLocalElementMatrices();
      }

      if ( invKMass_ != nullptr )
      {
         invKMass_->computeInverseDiagonalOperatorValues();
      }

      if ( ABlockChebyshevSmoother_ != nullptr )
      {
         if ( ( alwaysUpdateChebyshev_ ) ||
              ( std::floor( currentAge_ / chebyshevUpdateRateMyrs_ ) > std::floor( lastAge_ / chebyshevUpdateRateMyrs_ ) ) )
         {
            ABlockChebyshevSmoother_->regenerate();
         }
      }

      if ( ABlockPETScCoarseGridSolver_ != nullptr )
      {
         ABlockPETScCoarseGridSolver_->reassembleMatrix( saddlePointOp_->getA() );
      }

      if ( ABlockCoarseCorrectedChebyshevSmoother_ != nullptr )
      {
         if ( ( alwaysUpdateChebyshev_ ) ||
              ( std::floor( currentAge_ / chebyshevUpdateRateMyrs_ ) > std::floor( lastAge_ / chebyshevUpdateRateMyrs_ ) ) )
         {
            ABlockCoarseCorrectedChebyshevSmoother_->regenerate();
         }
      }

      auxFct_();

      lastAge_ = currentAge_;
   }

   void updateAfterDensityRecalculation() override {}

   void updateAfterViscosityExtrapolationRecalculation() override {}

   void updateAfterDensityExtrapolationRecalculation() override {}

   void updateAfterUpExtrapolationRecalculation() override {}

   void updateAfterTemperatureExtrapolationRecalculation() override {}

 private:
   std::shared_ptr< SaddlePointOperatorType_ >         saddlePointOp_;
   std::shared_ptr< InvKMassOperatorType_ >            invKMass_;
   std::shared_ptr< ABlockChebyshevSmootherType_ >     ABlockChebyshevSmoother_;
   std::shared_ptr< ABlockPETScCoarseGridSolverType_ > ABlockPETScCoarseGridSolver_;
   std::shared_ptr< ABlockChebyshevSmootherType_ >     ABlockCoarseCorrectedChebyshevSmoother_;

   std::shared_ptr< MCModel > MCModel_;
   real_t                     chebyshevUpdateRateMyrs_;
   bool                       alwaysUpdateChebyshev_;

   real_t currentAge_;
   real_t lastAge_;

   std::function< void() > auxFct_;
};

typedef CustomOperatorUpdater< SaddlePointOperatorType_,
                               WrappedMassType,
                               ABlockChebyshevSmoother< SaddlePointOperatorType_ >,
                               ABlockPETScCoarseGridSolver< stokesType::AOperatorType > >
    OperatorUpdaterType;

// ##################################
// ############## Main ##############
// ##################################

int main( int argc, char** argv )
{
   // Create walberla & MPI environment
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   // PETSc
   PETScManager petscManager( &argc, &argv );
#endif

   // load scenario file
   std::string parameterfile;
   std::string vtkfile;
   if ( argc < 2 )
   {
      std::shared_ptr< walberla::config::Config > cfgScenario_ = std::make_shared< walberla::config::Config >();
      cfgScenario_->readParameterFile( "./ScenarioSelector.prm" );
      walberla::config::Config::BlockHandle parametersScenario_ = cfgScenario_->getOneBlock( "Parameters" );
      std::string                           s                   = parametersScenario_.getParameter< std::string >( "Scenario" );

      parameterfile = s + std::string( ".prm" );
      vtkfile       = s;
   }
   else
   {
      parameterfile = std::string( argv[1] ) + std::string( ".prm" );
      vtkfile       = argv[1];
   }

   WALBERLA_LOG_INFO_ON_ROOT( "parameterfile: " << parameterfile );

   std::shared_ptr< MCModel > Model = std::make_shared< MCModel >( 2, parameterfile );

   auto  storage_       = Model->getStorage();
   auto  minLevel_      = Model->getMinLevel();
   auto  maxLevel_      = Model->getMaxLevel();
   auto  lowMemoryMode_ = Model->getLowMemoryMode();
   auto& ND_            = Model->getNondimensionalisation();
   auto& parameters_    = Model->getParameters();
   auto& surfaceFct_    = Model->getSurfaceFct();
   auto& CMBFct_        = Model->getCMBFct();
   auto  projection_    = Model->getVelocityProjection();
   auto& eta_           = Model->getEta();
   auto& rho_           = Model->getRho();
   auto& THistory       = Model->getTHistory();
   auto& up_            = Model->getUp();
   auto& up_extra_      = Model->getUpExtra();
   auto& eta_extra_     = Model->getEtaExtra();
   auto& rho_extra_     = Model->getRhoExtra();
   auto& BC_            = Model->getBC();
   //    auto  SUPG_scaling_ = Model->getSUPG_scaling();
   auto const_H_     = Model->getConst_H();
   auto const_alpha_ = Model->getConst_alpha();
   auto const_k_     = Model->getConst_k();
   auto const_C_p_   = Model->getConst_C_p();

   auto inv_eta_       = Model->getInvEta();
   auto inv_rho_       = Model->getInvRho();
   auto inv_eta_extra_ = Model->getInvEtaExtra();
   auto inv_rho_extra_ = Model->getInvRhoExtra();

   auto BDFOrder_ = Model->getBDFOrder();

   // ###############################
   // ### Model specific settings ###
   // ###############################

   WALBERLA_LOG_INFO_ON_ROOT( "--- Model specific settings ---" );

   uint_t blockpreconditionerType_ = parameters_.getParameter< uint_t >( "blockpreconditionerType" );
   WALBERLA_LOG_INFO_ON_ROOT( "Blockpreconditioner type: " << blockpreconditionerType_ );

   const int  OldAParam   = parameters_.getParameter< int >( "OldAMaxLevel" );
   const bool disableOldA = parameters_.getParameter< bool >( "disableOldA" );
   WALBERLA_LOG_INFO_ON_ROOT( "disableOldA: " << disableOldA );
   uint_t OldAMaxLevel;
   if ( OldAParam < 0 )
   {
      OldAMaxLevel = std::min( minLevel_ + walberla::numeric_cast< uint_t >( -OldAParam ), maxLevel_ );
   }
   else
   {
      OldAMaxLevel = walberla::numeric_cast< uint_t >( OldAParam );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "OldA MaxLevel: " << OldAMaxLevel );

   real_t chebyshevUpdateRateMyrs = parameters_.getParameter< real_t >( "chebyshevUpdateRateMyrs" );
   WALBERLA_LOG_INFO_ON_ROOT( "chebyshevUpdateRateMyrs: " << chebyshevUpdateRateMyrs );

   const real_t shearHeatingCutoffDimensional = parameters_.getParameter< real_t >( "shearHeatingCutoff" );
   const real_t shearHeatingCutoff            = shearHeatingCutoffDimensional / ND_.d_;
   WALBERLA_LOG_INFO_ON_ROOT( "shearHeatingCutoff: " << std::scientific << shearHeatingCutoffDimensional << " meters" );

   // #################################################
   // ### Temperature, Density and Viscosity models ###
   // #################################################

   auto ExpTemp     = std::make_shared< ExponentialTemperature >( ND_, parameters_, surfaceFct_, CMBFct_ );
   auto RandomTemp  = std::make_shared< RelativeRandomTemperature >( ND_, parameters_, ExpTemp );
   auto ExpDensity  = std::make_shared< ExponentialDensity >( ND_, parameters_ );
   auto ViscProfile = std::make_shared< SimpleSpaceDependentViscosityProfileWithJumps >( ND_, parameters_ );
   auto ExpVisc     = std::make_shared< ExponentialViscosity >( ND_, parameters_, ViscProfile );

   // some extra options:
   //    auto data = MantleConvection::loadCSV( "./Utility/Data/TemperatureProfiles/SimpleQualitativeProfile.csv" );
   //    MantleConvection::nondimensionaliseCSV( parameters_, data );

   //    auto LinearTempInter = std::make_shared< MantleConvection::LinearTemperatureInterpolation >(
   //        ND_, parameters_, surfaceFct_, CMBFct_, data.range.front(), data.values.front(), false );
   //    auto SphericalHarmonicsTemp = std::make_shared< SphericalHarmonicsTemperature >( ND_, parameters_, ExpTemp );

   // ########################
   // ### Initialise model ###
   // ########################

   auto& initTempModel  = RandomTemp;
   auto& refTempModel   = ExpTemp;
   auto& densityModel   = ExpDensity;
   auto& viscosityModel = ExpVisc;

   Model->setModels( initTempModel, refTempModel, densityModel, viscosityModel );
   Model->init();

   // ############################
   // ###### Define A Block ######
   // ############################

   const real_t divdivScaling = real_c( 1 );

   auto                                              viscosityFct             = viscosityModel->getViscosityFct();
   std::function< real_t( const Point3D&, real_t ) > viscosityFunctionWrapper = [=]( const Point3D& x, real_t temp ) {
      return viscosityFct( x, { temp } );
   };

   std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > visc_A =
       [&]( const Point3D& x, real_t temp, real_t temp_cent, real_t centX, real_t centY, real_t centZ ) {
          WALBERLA_UNUSED( centX );
          WALBERLA_UNUSED( centY );
          WALBERLA_UNUSED( centZ );
          WALBERLA_UNUSED( temp_cent );
          return viscosityFunctionWrapper( x, temp );
       };

   std::function< real_t( const Point3D& ) > k_ = [=]( const Point3D& x ) {
      WALBERLA_UNUSED( x );
      return divdivScaling;
   };

   std::shared_ptr< stokesType::AOperatorTypeInternal > ABlock;
   if ( disableOldA )
   {
      auto NewABlock = std::make_shared< NewAType >( storage_, minLevel_, maxLevel_, eta_, divdivScaling );

      ABlock = std::make_shared< stokesType::AOperatorTypeInternal >( storage_, minLevel_, maxLevel_, NewABlock, nullptr );
   }
   else
   {
      auto OldABlock =
          std::make_shared< OldAType >( storage_, minLevel_, OldAMaxLevel, THistory->getState( 0 ), visc_A, k_, false, nullptr );
      OldABlock->computeAndStoreLocalElementMatrices();

      auto NewABlock = std::make_shared< NewAType >( storage_, minLevel_, maxLevel_, eta_, divdivScaling );

      ABlock = std::make_shared< stokesType::AOperatorTypeInternal >( storage_, minLevel_, maxLevel_, NewABlock, OldABlock );
   }

   // ############################
   // ###### Define B Block ######
   // ############################

   auto BBlock = std::make_shared< stokesType::BOperatorTypeInternal >( storage_, minLevel_, maxLevel_ );

   // #############################
   // ###### Define BT Block ######
   // #############################

   auto BTBlock = std::make_shared< stokesType::BTOperatorTypeInternal >( storage_, minLevel_, maxLevel_ );

   real_t AScaling  = real_c( 1 );
   real_t BScaling  = real_c( 1 );
   real_t BTScaling = real_c( 1 );

   // ####################################
   // ### Define saddle point operator ###
   // ####################################

   auto stokesOperator_ = std::make_shared< stokesLHSType >( storage_,
                                                             minLevel_,
                                                             maxLevel_,
                                                             ABlock,
                                                             BBlock,
                                                             BTBlock,
                                                             nullptr,
                                                             AScaling,
                                                             BScaling,
                                                             BTScaling,
                                                             real_c( 0 ),
                                                             lowMemoryMode_,
                                                             projection_,
                                                             hyteg::FreeslipBoundary );
   stokesOperator_->computeInverseDiagonalOperatorValues();

   // ########################################
   // ### Define saddle point RHS operator ###
   // ########################################

   typedef MCModel::SaddlePointRHSOperatorType stokesTypeRHS;
   auto                                        TALA_RHS =
       std::make_shared< stokesTypeRHS::TemperatureToVelocityRHSOperatorTypeInternal >( storage_, minLevel_, maxLevel_, rho_ );
   auto gradRhoRho = std::make_shared< stokesTypeRHS::VelocityToPressureRHSOperatorTypeInternal >(
       storage_, minLevel_, maxLevel_, *inv_rho_, rho_ );

   real_t TALAScaling       = ND_.Ra_ / ND_.Pe_ * const_alpha_;
   real_t gradRhoRhoScaling = real_c( 1 );

   auto stokesOperatorRHS_ = std::make_shared< stokesTypeRHS >( storage_,
                                                                minLevel_,
                                                                maxLevel_,
                                                                nullptr,
                                                                nullptr,
                                                                TALA_RHS,
                                                                gradRhoRho,
                                                                nullptr,
                                                                nullptr,
                                                                nullptr,
                                                                real_c( 0 ),
                                                                real_c( 0 ),
                                                                TALAScaling,
                                                                gradRhoRhoScaling,
                                                                real_c( 0 ),
                                                                real_c( 0 ),
                                                                real_c( 0 ),
                                                                projection_,
                                                                hyteg::FreeslipBoundary,
                                                                lowMemoryMode_ );

   // ###########################################
   // ### Define advection diffusion operator ###
   // ###########################################

   typedef MCModel::AdvectionDiffusionOperatorType transportType;

   auto timeScheme = std::make_shared< BDFScheme< P2Function< real_t >, transportType::AdditiveMassType > >( BDFOrder_ );

   auto MassOperator = std::make_shared< transportType::MassOperatorTypeInternal >( storage_, minLevel_, maxLevel_ );
   auto DiffusionOperator =
       std::make_shared< transportType::DiffusionOperatorTypeInternal >( storage_, minLevel_, maxLevel_, *inv_rho_extra_ );
   auto DiffusionAdditionalOperator = std::make_shared< transportType::DiffusionAdditionalOperatorTypeInternal >(
       storage_, minLevel_, maxLevel_, rho_extra_ );
   auto AdiabaticHeatingOperator = std::make_shared< transportType::AdiabaticHeatingOperatorTypeInternal >(
       storage_, minLevel_, maxLevel_, up_extra_.uvw()[0], up_extra_.uvw()[1] );

   real_t AdvectionScaling           = real_c( 1 );
   real_t DiffusionScaling           = const_k_ / ND_.Pe_ / const_C_p_;
   real_t DiffusionAdditionalScaling = const_k_ / ND_.Pe_ / const_C_p_;
   real_t AdiabaticHeatingScaling    = ND_.Di_ * const_alpha_ / const_C_p_;

   auto transportOperator_ = std::make_shared< transportType >( storage_,
                                                                minLevel_,
                                                                maxLevel_,
                                                                THistory,
                                                                MassOperator,
                                                                nullptr,
                                                                DiffusionOperator,
                                                                DiffusionAdditionalOperator,
                                                                AdiabaticHeatingOperator,
                                                                AdvectionScaling,
                                                                DiffusionScaling,
                                                                DiffusionAdditionalScaling,
                                                                AdiabaticHeatingScaling,
                                                                timeScheme );

   // ###############################################
   // ### Define advection diffusion RHS operator ###
   // ###############################################

   typedef MCModel::AdvectionDiffusionRHSOperatorType transportType_RHS;

   auto ShearHeatingOperator = std::make_shared< transportType_RHS::ShearHeatingOperatorTypeInternal >( storage_,
                                                                                                        minLevel_,
                                                                                                        maxLevel_,
                                                                                                        eta_extra_,
                                                                                                        *inv_rho_extra_,
                                                                                                        up_extra_.uvw()[0],
                                                                                                        up_extra_.uvw()[1],
                                                                                                        ND_.radiusSurface_,
                                                                                                        shearHeatingCutoff );

   real_t InternalHeatingScaling = const_H_ / const_C_p_;
   real_t ShearHeatingScaling    = ND_.Pe_ * ND_.Di_ / const_C_p_ / ND_.Ra_;

   auto transportOperator_RHS_ = std::make_shared< transportType_RHS >( storage_,
                                                                        minLevel_,
                                                                        maxLevel_,
                                                                        THistory,
                                                                        MassOperator,
                                                                        MassOperator,
                                                                        ShearHeatingOperator,
                                                                        nullptr,
                                                                        InternalHeatingScaling,
                                                                        ShearHeatingScaling,
                                                                        AdiabaticHeatingScaling,
                                                                        timeScheme );

   // ###################################
   // ###### Define A Block solver ######
   // ###################################

   // PETSC Solver template:
   PETScSolverOptions solverOptions;
   solverOptions.setKspType( "cg" );
   solverOptions.setPcType( "gamg" );

   auto ABlockPETScCoarseSolver = std::make_shared< ABlockPETScCoarseGridSolver< stokesType::AOperatorType > >(
       parameters_, storage_, minLevel_, solverOptions );

   ABlockPETScCoarseSolver->getSolver()->getEnumerator().copyBoundaryConditionFromFunction( up_.uvw() );

   PETScSparseMatrix< P2ProjectNormalOperator > PMat;
   PMat.createMatrixFromOperator(
       *projection_, minLevel_, ABlockPETScCoarseSolver->getSolver()->getEnumerator(), FreeslipBoundary );

   auto PVec = std::make_shared< PETScVector< real_t, P2VectorFunction > >();
   PVec->createVectorFromFunction( up_.uvw(), ABlockPETScCoarseSolver->getSolver()->getEnumerator(), minLevel_, All );

   std::shared_ptr< std::function< void( Vec ) > > petscFreeslipFunction =
       std::make_shared< std::function< void( Vec v ) > >( [&]( Vec v ) {
          VecCopy( v, PVec->get() );
          MatMult( PMat.get(), PVec->get(), v );
       } );

   ABlockPETScCoarseSolver->getSolver()->setNullSpaceFunction( petscFreeslipFunction );

   std::shared_ptr< ABlockSolver< stokesType::AOperatorType > > ABlockCoarseSolver;

   ABlockCoarseSolver = ABlockPETScCoarseSolver;

   auto ABlockRestrictionOperator =
       std::make_shared< P2toP2QuadraticVectorRestrictionWithProjection< stokesType::VelocityProjectionOperatorType > >(
           storage_, minLevel_, maxLevel_, stokesOperator_->getProjPtr(), stokesOperator_->getProjFlag(), lowMemoryMode_ );
   auto ABlockProlongationOperator =
       std::make_shared< P2toP2QuadraticVectorProlongationWithProjection< stokesType::VelocityProjectionOperatorType > >(
           storage_, minLevel_, maxLevel_, stokesOperator_->getProjPtr(), stokesOperator_->getProjFlag(), lowMemoryMode_ );

   auto ABlockSmoother = std::make_shared< ABlockChebyshevSmoother< stokesType > >(
       parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, true, lowMemoryMode_ );

   auto ABlockMGSolver = std::make_shared<
       ABlockMultigridSolver< stokesType::AOperatorType,
                              P2toP2QuadraticVectorRestrictionWithProjection< stokesType::VelocityProjectionOperatorType >,
                              P2toP2QuadraticVectorProlongationWithProjection< stokesType::VelocityProjectionOperatorType > > >(
       parameters_,
       storage_,
       minLevel_,
       maxLevel_,
       ABlockSmoother,
       ABlockCoarseSolver,
       ABlockRestrictionOperator,
       ABlockProlongationOperator,
       lowMemoryMode_ );

   auto ABlockOuterLoop = std::make_shared< ABlockCGOuterLoopSolver< stokesType::AOperatorType > >(
       parameters_, storage_, minLevel_, maxLevel_, lowMemoryMode_, ABlockMGSolver );

   // #####################################################
   // ### Define coarse corrected A Block approximation ###
   // #####################################################

   auto CoarseCorrectedABlockSmoother = std::make_shared< ABlockChebyshevSmoother< stokesType > >(
       parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, true, lowMemoryMode_, true, 35466, "CoarseCorrected" );

   auto CoarseCorrectedABlockSolver = std::make_shared<
       ABlockMultigridSolver< stokesType::AOperatorType,
                              P2toP2QuadraticVectorRestrictionWithProjection< stokesType::VelocityProjectionOperatorType >,
                              P2toP2QuadraticVectorProlongationWithProjection< stokesType::VelocityProjectionOperatorType > > >(
       parameters_,
       storage_,
       minLevel_,
       maxLevel_,
       CoarseCorrectedABlockSmoother,
       ABlockCoarseSolver,
       ABlockRestrictionOperator,
       ABlockProlongationOperator,
       lowMemoryMode_,
       "CoarseCorrected" );

   // #########################################################
   // ### Define scaled mass Schur complement approximation ###
   // #########################################################

   auto InvKMass = std::make_shared< invKMassType >( storage_, minLevel_, maxLevel_, *inv_eta_ );
   InvKMass->computeInverseDiagonalOperatorValues();
   auto WrappedMass = std::make_shared< WrappedMassType >( InvKMass, lowMemoryMode_ );

   auto SchurInverseDiagonalSolver = std::make_shared< SchurInverseDiagMassSolver< WrappedMassType > >(
       parameters_, storage_, minLevel_, maxLevel_, *WrappedMass );

   auto SchurMassSolver = std::make_shared< SchurCGMassSolver< WrappedMassType > >(
       parameters_, storage_, minLevel_, maxLevel_, WrappedMass, lowMemoryMode_, SchurInverseDiagonalSolver );

   // ###########################
   // ### Define main solvers ###
   // ###########################

   // Stokes solver
   std::shared_ptr< SaddlePointSolver< stokesType > > SaddlePointPreconditioner;
   switch ( blockpreconditionerType_ )
   {
   case 0: {
      SaddlePointPreconditioner = std::make_shared< SaddlePointInexactUzawaSmoother< stokesType > >(
          parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, ABlockOuterLoop, SchurMassSolver, lowMemoryMode_ );
   }
   break;
   case 1: {
      SaddlePointPreconditioner = std::make_shared< SaddlePointAdjointInexactUzawaSmoother< stokesType > >(
          parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, ABlockOuterLoop, SchurMassSolver, lowMemoryMode_ );
   }
   break;
   case 2: {
      SaddlePointPreconditioner = std::make_shared< SaddlePointBlockApproximationFactorisationSmoother< stokesType > >(
          parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, ABlockOuterLoop, SchurMassSolver, lowMemoryMode_ );
   }
   break;
   case 3: {
      SaddlePointPreconditioner = std::make_shared< SaddlePointSymmetricUzawaSmoother< stokesType > >(
          parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, ABlockOuterLoop, SchurMassSolver, lowMemoryMode_ );
   }
   break;
   default: {
      WALBERLA_ABORT( "Unknown block preconditioner type " << blockpreconditionerType_ << "!" );
   }
   break;
   }

   auto FGMRES = std::make_shared< SaddlePointFGMRESSolver< stokesType > >(
       parameters_, storage_, minLevel_, maxLevel_, lowMemoryMode_, SaddlePointPreconditioner );

   // Advection Diffusion solver
   auto AdvectionDiffusionPreconditioner = std::make_shared< AdvectionDiffusionCGSPDPreconditioner< transportType > >(
       parameters_, storage_, minLevel_, maxLevel_, transportOperator_ );

   auto AdvectionDiffusionFGMRES = std::make_shared< AdvectionDiffusionFGMRESSolver< transportType > >(
       parameters_, storage_, minLevel_, maxLevel_, lowMemoryMode_, AdvectionDiffusionPreconditioner );

   // ###################
   // ### Final setup ###
   // ###################
   Model->setOperators( stokesOperator_, stokesOperatorRHS_, transportOperator_, transportOperator_RHS_ );
   Model->setSolvers( FGMRES, AdvectionDiffusionFGMRES );

   std::function< void() > auxFct           = [&]() {};
   auto                    operatorUpdater_ = std::make_shared< OperatorUpdaterType >( stokesOperator_,
                                                                    WrappedMass,
                                                                    ABlockSmoother,
                                                                    ABlockPETScCoarseSolver,
                                                                    CoarseCorrectedABlockSmoother,
                                                                    Model,
                                                                    chebyshevUpdateRateMyrs,
                                                                    auxFct );
   Model->setOperatorUpdater( operatorUpdater_ );

   // ###############################
   // ### Print Model Information ###
   // ###############################
   WALBERLA_LOG_INFO_ON_ROOT( *Model );

   // #################
   // ### Run Model ###
   // #################
   Model->runSimulation();

   return EXIT_SUCCESS;
}