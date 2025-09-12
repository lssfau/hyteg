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
#include <map>

#include "hyteg/functions/FunctionHistory.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2InjectionRestriction.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/numerictools/BDFScheme.hpp"
#include "hyteg/operators/GEMV.hpp"
#include "hyteg/p2functionspace/P2FullViscousTDependentOperator.hpp"
#include "hyteg/p2functionspace/P2VectorMassOperator.hpp"
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
#include "Utility/Viscosity/ArrheniusViscosity.hpp"
#include "Utility/Viscosity/BoundedViscosity.hpp"
#include "Utility/Viscosity/ConstantViscosity.hpp"
#include "Utility/Viscosity/ExponentialViscosity.hpp"
#include "Utility/Viscosity/LinearViscosityInterpolation.hpp"
#include "Utility/Viscosity/SimpleSpaceDependentViscosityProfileWithJumps.hpp"
#include "terraneo/helpers/ConvectionToolbox.hpp"

using namespace hyteg::convectionToolbox;
using namespace MantleConvection;

// Define Combined A Block type
typedef MC_ABlock_Vec_IcosahedralShellMap                                        NewAType;
typedef P2ElementwiseBlendingFullViscousTDependentOperator_Centroid_ScaledDivDiv OldAType;

typedef CombinedABlockOperator< NewAType, OldAType, hyteg::P2VectorFunction< real_t > > CombinedAType;

// Define Mantle Convection operator types
typedef MantleConvection::SaddlePointOperator< CombinedAType,
                                               MC_BTBlock_IcosahedralShellMap,
                                               MC_GradRhoRhoDivergence_IcosahedralShellMap,
                                               MC_Projection,
                                               MC_NoOp >
    stokesLHSType;

typedef MantleConvection::
    SaddlePointOperator< CombinedAType, MC_BTBlock_IcosahedralShellMap, MC_BBlock_IcosahedralShellMap, MC_Projection, MC_NoOp >
        stokesLHSTypeSym;

typedef MantleConvection::SaddlePointOperatorRHS< MC_NoOp,
                                                  MC_NoOp,
                                                  MC_TemperatureToVelocityRHS_IcosahedralShellMap,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_Projection >
    stokesRHSType;

typedef MantleConvection::AdvectionDiffusionOperator< MC_P2Mass_IcosahedralShellMap,
                                                      MC_NoOp,
                                                      MC_DivKGrad_IcosahedralShellMap,
                                                      MC_DiffusionAdditional_IcosahedralShellMap,
                                                      MC_AdiabaticHeating_IcosahedralShellMap,

                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp >
    transportLHSType;

typedef MantleConvection::AdvectionDiffusionOperatorRHS< MC_P2Mass_IcosahedralShellMap,
                                                         MC_P2Mass_IcosahedralShellMap,
                                                         MC_ShearHeating_NoSurface_IcosahedralShellMap,
                                                         MC_NoOp,

                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp >
    transportRHSType;

// scaled mass type
typedef MC_KMass_IcosahedralShellMap invKMassType;

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
                               MC_P2Mass_IcosahedralShellMap,
                               hyteg::P2toP2QuadraticRestriction >
    MCModel;

typedef MCModel::SaddlePointOperatorType stokesType;

// WBFBT typedefs
typedef MC_P2KP1Mass_IcosahedralShellMap                                                              ViscSqrtMassType;
typedef P2VectorMassOperator< ViscSqrtMassType >                                                      BaseVectorSqrtMassType;
typedef ABlockProjectionWrapper< BaseVectorSqrtMassType, stokesType::VelocityProjectionOperatorType > VectorSqrtMassType;

typedef MC_P1DivKGrad_IcosahedralShellMap                     InvSqrtViscMassType;
typedef StabilisationProjectionWrapper< InvSqrtViscMassType > WrappedDivKGradType;

// define custom operator updater
template < class SaddlePointOperatorType_         = hyteg::NoOperator,
           class InvKMassOperatorType_            = hyteg::NoOperator,
           class ABlockChebyshevSmootherType_     = hyteg::NoOperator,
           class ABlockPETScCoarseGridSolverType_ = hyteg::NoOperator >
class CustomOperatorUpdater : public OperatorUpdater
{
 public:
   CustomOperatorUpdater(
       const std::shared_ptr< SaddlePointOperatorType_ >&         saddlePointOp,
       const std::shared_ptr< InvKMassOperatorType_ >&            invKMass,
       const std::shared_ptr< ABlockChebyshevSmootherType_ >&     ABlockChebyshevSmoother,
       const std::shared_ptr< ABlockPETScCoarseGridSolverType_ >& ABlockPETScCoarseGridSolver,
       const std::shared_ptr< ABlockChebyshevSmootherType_ >&     ABlockCoarseCorrectedChebyshevSmoother,
       const std::shared_ptr< MCModel >&                          _MCModel,
       real_t                                                     _chebyshevUpdateRateMyrs,
       std::function< void() >                                    _auxFct,

       uint_t                                                                                   minLevel,
       uint_t                                                                                   maxLevel,
       const std::shared_ptr< hyteg::FunctionHistory< hyteg::P2Function< real_t > > >&          THist,
       const std::shared_ptr< MantleConvection::TemperatureDependentViscosityModel< real_t > >& viscModel,

       const std::shared_ptr< MCModel::ViscosityFunctionType >&                                  sqrtEtaRight,
       const std::shared_ptr< MCModel::ViscosityFunctionType >&                                  invSqrtEtaRight,
       const std::shared_ptr< ViscSqrtMassType >&                                                SqrtViscMassOperatorRight,
       const std::shared_ptr< InvSqrtViscMassType >&                                             invSqrtViscMassOperatorRight,
       const std::shared_ptr< MantleConvection::ABlockChebyshevSmoother< VectorSqrtMassType > >& VectorMassSmootherRight,
       const std::shared_ptr< MantleConvection::SchurPoissonNeumannChebyshevSmoother< WrappedDivKGradType > >&
           SchurPoissonNeumannChebyshevRight,

       const std::shared_ptr< MCModel::ViscosityFunctionType >&                                  sqrtEtaLeft,
       const std::shared_ptr< MCModel::ViscosityFunctionType >&                                  invSqrtEtaLeft,
       const std::shared_ptr< ViscSqrtMassType >&                                                SqrtViscMassOperatorLeft,
       const std::shared_ptr< InvSqrtViscMassType >&                                             invSqrtViscMassOperatorLeft,
       const std::shared_ptr< MantleConvection::ABlockChebyshevSmoother< VectorSqrtMassType > >& VectorMassSmootherLeft,
       const std::shared_ptr< MantleConvection::SchurPoissonNeumannChebyshevSmoother< WrappedDivKGradType > >&
           SchurPoissonNeumannChebyshevLeft,

       const std::shared_ptr< hyteg::P1Function< real_t > >& maskRight,
       const std::shared_ptr< hyteg::P1Function< real_t > >& maskLeft )
   : saddlePointOp_( saddlePointOp )
   , invKMass_( invKMass )
   , ABlockChebyshevSmoother_( ABlockChebyshevSmoother )
   , ABlockPETScCoarseGridSolver_( ABlockPETScCoarseGridSolver )
   , ABlockCoarseCorrectedChebyshevSmoother_( ABlockCoarseCorrectedChebyshevSmoother )
   , MCModel_( _MCModel )
   , chebyshevUpdateRateMyrs_( _chebyshevUpdateRateMyrs )
   , alwaysUpdateChebyshev_( _chebyshevUpdateRateMyrs <= real_c( 0 ) )
   , auxFct_( _auxFct )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , THist_( THist )
   , viscModel_( viscModel )

   , sqrtEtaRight_( sqrtEtaRight )
   , invSqrtEtaRight_( invSqrtEtaRight )
   , SqrtViscMassOperatorRight_( SqrtViscMassOperatorRight )
   , invSqrtViscMassOperatorRight_( invSqrtViscMassOperatorRight )
   , VectorMassSmootherRight_( VectorMassSmootherRight )
   , SchurPoissonNeumannChebyshevRight_( SchurPoissonNeumannChebyshevRight )

   , sqrtEtaLeft_( sqrtEtaLeft )
   , invSqrtEtaLeft_( invSqrtEtaLeft )
   , SqrtViscMassOperatorLeft_( SqrtViscMassOperatorLeft )
   , invSqrtViscMassOperatorLeft_( invSqrtViscMassOperatorLeft )
   , VectorMassSmootherLeft_( VectorMassSmootherLeft )
   , SchurPoissonNeumannChebyshevLeft_( SchurPoissonNeumannChebyshevLeft )

   , maskRight_( maskRight )
   , maskLeft_( maskLeft )
   {}

   void updateBeforeSaddlePointSolve() override {}

   void updateAfterCheckPointLoad() override
   {
      auto                                              viscFct             = viscModel_->getViscosityFct();
      std::function< real_t( const Point3D&, real_t ) > viscFunctionWrapper = [=]( const Point3D& x, real_t temp ) {
         return viscFct( x, { temp } );
      };

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > sqrtEtaFct =
          [&]( const Point3D& x, const std::vector< real_t >& temperature ) {
             return std::sqrt( viscFunctionWrapper( x, temperature[0] ) );
          };

      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         sqrtEtaRight_->interpolate( sqrtEtaFct, { THist_->getState( 0 ).getVertexDoFFunction() }, level, All );
         sqrtEtaRight_->multElementwise( { *sqrtEtaRight_, *maskRight_ }, level, hyteg::All );
         sqrtEtaLeft_->interpolate( sqrtEtaFct, { THist_->getState( 0 ).getVertexDoFFunction() }, level, All );
         sqrtEtaLeft_->multElementwise( { *sqrtEtaLeft_, *maskLeft_ }, level, hyteg::All );

         invSqrtEtaRight_->assign( { real_c( 1 ) }, { *sqrtEtaRight_ }, level, hyteg::All );
         invSqrtEtaRight_->invertElementwise( level, hyteg::All );
         invSqrtEtaLeft_->assign( { real_c( 1 ) }, { *sqrtEtaLeft_ }, level, hyteg::All );
         invSqrtEtaLeft_->invertElementwise( level, hyteg::All );
      }

      SqrtViscMassOperatorRight_->computeInverseDiagonalOperatorValues();
      invSqrtViscMassOperatorRight_->computeInverseDiagonalOperatorValues();
      SqrtViscMassOperatorLeft_->computeInverseDiagonalOperatorValues();
      invSqrtViscMassOperatorLeft_->computeInverseDiagonalOperatorValues();

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

      if ( VectorMassSmootherRight_ != nullptr )
      {
         VectorMassSmootherRight_->regenerate();
      }
      if ( VectorMassSmootherLeft_ != nullptr )
      {
         VectorMassSmootherLeft_->regenerate();
      }

      if ( SchurPoissonNeumannChebyshevRight_ != nullptr )
      {
         SchurPoissonNeumannChebyshevRight_->regenerate();
      }
      if ( SchurPoissonNeumannChebyshevLeft_ != nullptr )
      {
         SchurPoissonNeumannChebyshevLeft_->regenerate();
      }

      auxFct_();

      lastAge_ = MCModel_->getCurrentAge();
   }

   void updateAfterSaddlePointSolve() override {}

   void updateAfterAdvectionDiffusionSolve() override {}

   void updateAfterViscosityRecalculation() override
   {
      currentAge_ = MCModel_->getCurrentAge();

      auto                                              viscFct             = viscModel_->getViscosityFct();
      std::function< real_t( const Point3D&, real_t ) > viscFunctionWrapper = [=]( const Point3D& x, real_t temp ) {
         return viscFct( x, { temp } );
      };

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > sqrtEtaFct =
          [&]( const Point3D& x, const std::vector< real_t >& temperature ) {
             return std::sqrt( viscFunctionWrapper( x, temperature[0] ) );
          };

      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         sqrtEtaRight_->interpolate( sqrtEtaFct, { THist_->getState( 0 ).getVertexDoFFunction() }, level, All );
         sqrtEtaRight_->multElementwise( { *sqrtEtaRight_, *maskRight_ }, level, hyteg::All );
         sqrtEtaLeft_->interpolate( sqrtEtaFct, { THist_->getState( 0 ).getVertexDoFFunction() }, level, All );
         sqrtEtaLeft_->multElementwise( { *sqrtEtaLeft_, *maskLeft_ }, level, hyteg::All );

         invSqrtEtaRight_->assign( { real_c( 1 ) }, { *sqrtEtaRight_ }, level, hyteg::All );
         invSqrtEtaRight_->invertElementwise( level, hyteg::All );
         invSqrtEtaLeft_->assign( { real_c( 1 ) }, { *sqrtEtaLeft_ }, level, hyteg::All );
         invSqrtEtaLeft_->invertElementwise( level, hyteg::All );
      }

      SqrtViscMassOperatorRight_->computeInverseDiagonalOperatorValues();
      invSqrtViscMassOperatorRight_->computeInverseDiagonalOperatorValues();
      SqrtViscMassOperatorLeft_->computeInverseDiagonalOperatorValues();
      invSqrtViscMassOperatorLeft_->computeInverseDiagonalOperatorValues();

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

      if ( VectorMassSmootherRight_ != nullptr )
      {
         VectorMassSmootherRight_->regenerate();
      }
      if ( VectorMassSmootherLeft_ != nullptr )
      {
         VectorMassSmootherLeft_->regenerate();
      }

      if ( SchurPoissonNeumannChebyshevRight_ != nullptr )
      {
         SchurPoissonNeumannChebyshevRight_->regenerate();
      }
      if ( SchurPoissonNeumannChebyshevLeft_ != nullptr )
      {
         SchurPoissonNeumannChebyshevLeft_->regenerate();
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

   uint_t                                                                            minLevel_;
   uint_t                                                                            maxLevel_;
   std::shared_ptr< hyteg::FunctionHistory< hyteg::P2Function< real_t > > >          THist_;
   std::shared_ptr< MantleConvection::TemperatureDependentViscosityModel< real_t > > viscModel_;

   std::shared_ptr< MCModel::ViscosityFunctionType >                                  sqrtEtaRight_;
   std::shared_ptr< MCModel::ViscosityFunctionType >                                  invSqrtEtaRight_;
   std::shared_ptr< ViscSqrtMassType >                                                SqrtViscMassOperatorRight_;
   std::shared_ptr< InvSqrtViscMassType >                                             invSqrtViscMassOperatorRight_;
   std::shared_ptr< MantleConvection::ABlockChebyshevSmoother< VectorSqrtMassType > > VectorMassSmootherRight_;
   std::shared_ptr< MantleConvection::SchurPoissonNeumannChebyshevSmoother< WrappedDivKGradType > >
       SchurPoissonNeumannChebyshevRight_;

   std::shared_ptr< MCModel::ViscosityFunctionType >                                  sqrtEtaLeft_;
   std::shared_ptr< MCModel::ViscosityFunctionType >                                  invSqrtEtaLeft_;
   std::shared_ptr< ViscSqrtMassType >                                                SqrtViscMassOperatorLeft_;
   std::shared_ptr< InvSqrtViscMassType >                                             invSqrtViscMassOperatorLeft_;
   std::shared_ptr< MantleConvection::ABlockChebyshevSmoother< VectorSqrtMassType > > VectorMassSmootherLeft_;
   std::shared_ptr< MantleConvection::SchurPoissonNeumannChebyshevSmoother< WrappedDivKGradType > >
       SchurPoissonNeumannChebyshevLeft_;

   std::shared_ptr< hyteg::P1Function< real_t > > maskRight_;
   std::shared_ptr< hyteg::P1Function< real_t > > maskLeft_;
};

typedef CustomOperatorUpdater< SaddlePointOperatorType_,
                               WrappedMassType,
                               ABlockChebyshevSmoother< SaddlePointOperatorType_ >,
                               hyteg::NoOperator >
    OperatorUpdaterType;

// Weighted BFBT Tools

// Using a std::map for a sparse vector might not the most efficient idea but since we
// only use this once this avoids using PETSc vectors and should be good enough to store
// which DoFs are flagged as belonging to a boundary touching element.
class CustomVector : public VectorProxy
{
 public:
   CustomVector( real_t a = real_c( 10 ) )
   : a_( a )
   {}
   virtual ~CustomVector() {}

   void setValue( uint_t idx, real_t value ) override
   {
      if ( value > 0 )
      {
         valueMap[idx] = true;
      }
   }

   real_t getValue( uint_t idx ) const override
   {
      if ( valueMap.find( idx ) == valueMap.end() )
      {
         return real_c( 1 );
      }
      else
      {
         return a_;
      }
   }

 private:
   real_t                           a_;
   mutable std::map< uint_t, bool > valueMap;
};

void createBoundaryScalingMask( const std::shared_ptr< P1Function< real_t > >& mask,
                                uint_t                                         level,
                                uint_t                                         depth,
                                real_t                                         a,
                                hyteg::DoFType                                 flag = hyteg::DirichletBoundary )
{
   auto storage = mask->getStorage();
   auto bc      = mask->getBoundaryCondition();

   // mark initial dirichlet boundary
   mask->interpolate( real_c( 1.0 ), level, flag );

   P1Function< real_t > temp( "mask tmp", storage, level, level );

   if ( storage->hasGlobalCells() )
   {
      MC_P1Mass_IcosahedralShellMap massOp( storage, level, level );

      for ( uint_t k = 0; k < depth; k++ )
      {
         temp.assign( { real_c( 1 ) }, { *mask }, level, hyteg::All );
         massOp.apply( temp, *mask, level, hyteg::All ^ flag );
      }
   }
   else
   {
      MC_P1Mass_AnnulusMap massOp( storage, level, level );

      for ( uint_t k = 0; k < depth; k++ )
      {
         temp.assign( { real_c( 1 ) }, { *mask }, level, hyteg::All );
         massOp.apply( temp, *mask, level, hyteg::All ^ flag );
      }
   }

   // use custom vector to create mask
   auto                vec = std::make_shared< CustomVector >( a );
   P1Function< idx_t > enumerator( "enumerator", storage, level, level, bc );

   enumerator.enumerate( level );

   mask->toVector( enumerator, vec, level, All );
   mask->fromVector( enumerator, vec, level, All );
}

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

   std::shared_ptr< MCModel > Model = std::make_shared< MCModel >( 3, parameterfile );

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

   const bool agglomerateACoarseSolver = parameters_.getParameter< bool >( "agglomerateACoarseSolver" );
   WALBERLA_LOG_INFO_ON_ROOT( "agglomerateACoarseSolver: " << agglomerateACoarseSolver );

   real_t chebyshevUpdateRateMyrs = parameters_.getParameter< real_t >( "chebyshevUpdateRateMyrs" );
   WALBERLA_LOG_INFO_ON_ROOT( "chebyshevUpdateRateMyrs: " << chebyshevUpdateRateMyrs );

   // Asymmetric Preconditioner mode (if used)
   // 0 = Asymmetric system is solved with the preconditioner for the symmetric system,
   // unmodified (divergence) B Block gets used for the block preconditioner and BFBT Schur Complement approx.
   // 1 = Asymmetric system is solved with the preconditioner for the asymmetric system,
   // modified B Block gets used for the block preconditioner and BFBT Schur Complement approx.
   // Note: For mode 1 the BFBT suboperator should be solved via a GMRES solver
   uint_t AsymmetricPreconditionerMode = parameters_.getParameter< uint_t >( "AsymmetricPreconditionerMode" );
   WALBERLA_LOG_INFO_ON_ROOT( "AsymmetricPreconditionerMode: " << AsymmetricPreconditionerMode );

   const real_t shearHeatingCutoffDimensional = parameters_.getParameter< real_t >( "shearHeatingCutoff" );
   const real_t shearHeatingCutoff            = shearHeatingCutoffDimensional / ND_.d_;
   WALBERLA_LOG_INFO_ON_ROOT( "shearHeatingCutoff: " << std::scientific << shearHeatingCutoffDimensional << " meters" );

   uint_t WBFBTType_ = parameters_.getParameter< uint_t >( "WBFBTType" );

   WALBERLA_LOG_INFO_ON_ROOT( "WBFBTType: " << WBFBTType_ );

   real_t ar = parameters_.getParameter< real_t >( "WBFBT_ar" );
   WALBERLA_LOG_INFO_ON_ROOT( "WBFBT_ar: " << ar );
   real_t al = parameters_.getParameter< real_t >( "WBFBT_al" );
   WALBERLA_LOG_INFO_ON_ROOT( "WBFBT_al: " << al );

   // #################################################
   // ### Temperature, Density and Viscosity models ###
   // #################################################

   auto dataTemp = MantleConvection::loadCSV( "./Utility/Data/TemperatureProfiles/InitialTemperatureProfileForArrhenius.csv" );
   MantleConvection::nondimensionaliseCSV( parameters_, dataTemp );
   auto dataVisc = MantleConvection::loadCSV( "./Utility/Data/ViscosityProfiles/viscosityStotz2018Cutoff_Arrhenius.csv" );
   MantleConvection::nondimensionaliseCSV( parameters_, dataVisc );

   auto ExpTemp  = std::make_shared< ExponentialTemperature >( ND_, parameters_, surfaceFct_, CMBFct_ );
   auto InitTemp = std::make_shared< LinearTemperatureInterpolation >(
       ND_, parameters_, surfaceFct_, CMBFct_, dataTemp.range.front(), dataTemp.values.front(), false );
   auto RandomTemp = std::make_shared< RelativeRandomTemperature >( ND_, parameters_, InitTemp );
   auto ExpDensity = std::make_shared< ExponentialDensity >( ND_, parameters_ );

   auto ViscProfile = std::make_shared< MantleConvection::LinearViscosityInterpolation >(
       ND_, parameters_, dataVisc.range.front(), dataVisc.values.front() );

   auto ArrVisc     = std::make_shared< ArrheniusViscosity >( ND_, parameters_, ViscProfile );
   auto BoundedVisc = std::make_shared< BoundedViscosity >( ND_, parameters_, ArrVisc );

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
   auto& viscosityModel = BoundedVisc;

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

   std::shared_ptr< stokesLHSTypeSym::BOperatorTypeInternal > BBlockSym;
   if ( AsymmetricPreconditionerMode == 0 )
   {
      BBlockSym = std::make_shared< stokesLHSTypeSym::BOperatorTypeInternal >( storage_, minLevel_, maxLevel_ );
   }

   auto BBlock = std::make_shared< stokesType::BOperatorTypeInternal >( storage_, minLevel_, maxLevel_, *inv_rho_, rho_ );

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

   std::shared_ptr< stokesLHSTypeSym > stokesOperatorSym_;
   if ( AsymmetricPreconditionerMode == 0 )
   {
      stokesOperatorSym_ = std::make_shared< stokesLHSTypeSym >( storage_,
                                                                 minLevel_,
                                                                 maxLevel_,
                                                                 ABlock,
                                                                 BBlockSym,
                                                                 BTBlock,
                                                                 nullptr,
                                                                 AScaling,
                                                                 BScaling,
                                                                 BTScaling,
                                                                 real_c( 0 ),
                                                                 lowMemoryMode_,
                                                                 projection_,
                                                                 hyteg::FreeslipBoundary );
      stokesOperatorSym_->computeInverseDiagonalOperatorValues();
   }

   auto stokesOperator_ = std::make_shared< stokesType >( storage_,
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

   real_t TALAScaling       = ND_.Ra_ / ND_.Pe_ * const_alpha_;
   real_t gradRhoRhoScaling = real_c( 0 );

   auto stokesOperatorRHS_ = std::make_shared< stokesTypeRHS >( storage_,
                                                                minLevel_,
                                                                maxLevel_,
                                                                nullptr,
                                                                nullptr,
                                                                TALA_RHS,
                                                                nullptr,
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
       storage_, minLevel_, maxLevel_, up_extra_.uvw()[0], up_extra_.uvw()[1], up_extra_.uvw()[2] );

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
                                                                                                        up_extra_.uvw()[2],
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

   std::shared_ptr< ABlockSolver< stokesType::AOperatorType > > ABlockCoarseSolver;
   std::shared_ptr< stokesType::AOperatorTypeInternal >         AuxABlock;
   std::function< void() >                                      auxFct = [&]() {};
   std::shared_ptr< P2Function< real_t > >                      AuxTemp;
   std::shared_ptr< P1Function< real_t > >                      AuxEta;
   if ( agglomerateACoarseSolver )
   {
      auto ABlockAgglomeratedCoarseSolver = std::make_shared< ABlockAgglomeratedCGCoarseGridSolver< stokesType::AOperatorType > >(
          parameters_, storage_, minLevel_, maxLevel_, lowMemoryMode_ );

      AuxTemp = std::make_shared< P2Function< real_t > >(
          "AuxTemp", ABlockAgglomeratedCoarseSolver->getAgglomerationStorage(), minLevel_, maxLevel_ );

      AuxEta = std::make_shared< P1Function< real_t > >(
          "AuxEta", ABlockAgglomeratedCoarseSolver->getAgglomerationStorage(), minLevel_, maxLevel_ );

      if ( disableOldA )
      {
         auto AuxNewABlock = std::make_shared< NewAType >(
             ABlockAgglomeratedCoarseSolver->getAgglomerationStorage(), minLevel_, maxLevel_, *AuxEta, divdivScaling );

         AuxABlock = std::make_shared< stokesType::AOperatorTypeInternal >(
             ABlockAgglomeratedCoarseSolver->getAgglomerationStorage(), minLevel_, maxLevel_, AuxNewABlock, nullptr );
      }
      else
      {
         auto AuxOldABlock = std::make_shared< OldAType >( ABlockAgglomeratedCoarseSolver->getAgglomerationStorage(),
                                                           minLevel_,
                                                           OldAMaxLevel,
                                                           *AuxTemp,
                                                           visc_A,
                                                           k_,
                                                           false,
                                                           nullptr );
         AuxOldABlock->computeAndStoreLocalElementMatrices();

         auto AuxNewABlock = std::make_shared< NewAType >(
             ABlockAgglomeratedCoarseSolver->getAgglomerationStorage(), minLevel_, maxLevel_, *AuxEta, divdivScaling );

         AuxABlock = std::make_shared< stokesType::AOperatorTypeInternal >(
             ABlockAgglomeratedCoarseSolver->getAgglomerationStorage(), minLevel_, maxLevel_, AuxNewABlock, AuxOldABlock );
      }

      typedef MantleConvection::SaddlePointOperator< CombinedAType, MC_NoOp, MC_NoOp, MC_Projection, MC_NoOp > AuxStokesLHSType;

      auto AuxProjection_ = std::make_shared< P2ProjectNormalOperator >(
          ABlockAgglomeratedCoarseSolver->getAgglomerationStorage(), minLevel_, maxLevel_, Model->getSurfaceNormal() );
      auto AuxStokesOperator_ = std::make_shared< AuxStokesLHSType >( ABlockAgglomeratedCoarseSolver->getAgglomerationStorage(),
                                                                      minLevel_,
                                                                      minLevel_,
                                                                      AuxABlock,
                                                                      nullptr,
                                                                      nullptr,
                                                                      nullptr,
                                                                      AScaling,
                                                                      real_c( 0 ),
                                                                      real_c( 0 ),
                                                                      real_c( 0 ),
                                                                      lowMemoryMode_,
                                                                      AuxProjection_,
                                                                      hyteg::FreeslipBoundary );

      auxFct = [&]() {
         ABlockAgglomeratedCoarseSolver->copyFunctionToAgglomerationStorage< P2Function< real_t > >(
             THistory->getState( 0 ), *AuxTemp, minLevel_ );
         ABlockAgglomeratedCoarseSolver->copyFunctionToAgglomerationStorage< P1Function< real_t > >( eta_, *AuxEta, minLevel_ );

         if ( AuxABlock->getCoarseOperator() != nullptr )
         {
            AuxABlock->getCoarseOperator()->computeAndStoreLocalElementMatrices();
         }
         AuxABlock->computeInverseDiagonalOperatorValues();
      };

      ABlockAgglomeratedCoarseSolver->setAgglomerationOperator( AuxStokesOperator_->getAPtr() );
      //   ABlockAgglomeratedCoarseSolver->setAuxiliaryUpdateFunction( auxFct );

      auxFct();
      ABlockCoarseSolver = ABlockAgglomeratedCoarseSolver;
   }
   else
   {
      ABlockCoarseSolver = std::make_shared< ABlockCGCoarseGridSolver< stokesType::AOperatorType > >(
          parameters_, storage_, minLevel_, maxLevel_, lowMemoryMode_ );
   }

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

   // ######################################
   // ### Define sqrt(eta) velocity mass ###
   // ######################################

   // create scaling mask
   auto maskRight = std::make_shared< P1Function< real_t > >( "maskRight", storage_, minLevel_, maxLevel_, BC_.bcVelocityX_ );
   auto maskLeft  = std::make_shared< P1Function< real_t > >( "maskLeft", storage_, minLevel_, maxLevel_, BC_.bcVelocityX_ );

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      createBoundaryScalingMask( maskRight, level, 1, ar );
      createBoundaryScalingMask( maskLeft, level, 1, al );
   }

   // sqrt vector mass
   auto sqrtEtaRight_ = std::make_shared< MCModel::ViscosityFunctionType >( "sqrt eta right", storage_, minLevel_, maxLevel_ );
   auto sqrtEtaLeft_  = std::make_shared< MCModel::ViscosityFunctionType >( "sqrt eta left", storage_, minLevel_, maxLevel_ );

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > sqrtEtaFct =
       [&]( const Point3D& x, const std::vector< real_t >& temperature ) {
          return std::sqrt( viscosityFunctionWrapper( x, temperature[0] ) );
       };

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      sqrtEtaRight_->interpolate( sqrtEtaFct, { THistory->getState( 0 ).getVertexDoFFunction() }, level, hyteg::All );
      sqrtEtaRight_->multElementwise( { *sqrtEtaRight_, *maskRight }, level, hyteg::All );
      sqrtEtaLeft_->interpolate( sqrtEtaFct, { THistory->getState( 0 ).getVertexDoFFunction() }, level, hyteg::All );
      sqrtEtaLeft_->multElementwise( { *sqrtEtaRight_, *maskLeft }, level, hyteg::All );
   }

   auto SqrtMassRight = std::make_shared< ViscSqrtMassType >( storage_, minLevel_, maxLevel_, *sqrtEtaRight_ );
   SqrtMassRight->computeInverseDiagonalOperatorValues();
   auto SqrtMassLeft = std::make_shared< ViscSqrtMassType >( storage_, minLevel_, maxLevel_, *sqrtEtaLeft_ );
   SqrtMassLeft->computeInverseDiagonalOperatorValues();

   auto BaseVectorMassRight = std::make_shared< BaseVectorSqrtMassType >( storage_, minLevel_, maxLevel_, SqrtMassRight );
   auto BaseVectorMassLeft  = std::make_shared< BaseVectorSqrtMassType >( storage_, minLevel_, maxLevel_, SqrtMassLeft );

   auto VectorMassRight = std::make_shared< VectorSqrtMassType >(
       BaseVectorMassRight, stokesOperator_->getProjPtr(), stokesOperator_->getProjFlag(), lowMemoryMode_ );
   auto VectorMassLeft = std::make_shared< VectorSqrtMassType >(
       BaseVectorMassLeft, stokesOperator_->getProjPtr(), stokesOperator_->getProjFlag(), lowMemoryMode_ );

   auto VectorMassCoarseSolver = std::make_shared< ABlockCGCoarseGridSolver< VectorSqrtMassType > >(
       parameters_,
       storage_,
       minLevel_,
       maxLevel_,
       lowMemoryMode_,
       std::make_shared< ABlockIdentityPreconditioner< VectorSqrtMassType > >(),
       "VectorMass" );

   auto VectorMassRestrictionOperator =
       std::make_shared< P2toP2QuadraticVectorRestrictionWithProjection< stokesType::VelocityProjectionOperatorType > >(
           storage_, minLevel_, maxLevel_, stokesOperator_->getProjPtr(), stokesOperator_->getProjFlag(), lowMemoryMode_ );
   auto VectorMassProlongationOperator =
       std::make_shared< P2toP2QuadraticVectorProlongationWithProjection< stokesType::VelocityProjectionOperatorType > >(
           storage_, minLevel_, maxLevel_, stokesOperator_->getProjPtr(), stokesOperator_->getProjFlag(), lowMemoryMode_ );

   auto VectorMassSmootherRight = std::make_shared< ABlockChebyshevSmoother< VectorSqrtMassType > >(
       parameters_, storage_, minLevel_, maxLevel_, VectorMassRight, true, lowMemoryMode_, true, 45489, "VectorMass" );
   auto VectorMassSmootherLeft = std::make_shared< ABlockChebyshevSmoother< VectorSqrtMassType > >(
       parameters_, storage_, minLevel_, maxLevel_, VectorMassLeft, true, lowMemoryMode_, true, 45489, "VectorMass" );

   auto VectorMassMGSolverRight = std::make_shared<
       ABlockMultigridSolver< VectorSqrtMassType,
                              P2toP2QuadraticVectorRestrictionWithProjection< stokesType::VelocityProjectionOperatorType >,
                              P2toP2QuadraticVectorProlongationWithProjection< stokesType::VelocityProjectionOperatorType > > >(
       parameters_,
       storage_,
       minLevel_,
       maxLevel_,
       VectorMassSmootherRight,
       VectorMassCoarseSolver,
       VectorMassRestrictionOperator,
       VectorMassProlongationOperator,
       lowMemoryMode_,
       "VectorMass" );
   auto VectorMassMGSolverLeft = std::make_shared<
       ABlockMultigridSolver< VectorSqrtMassType,
                              P2toP2QuadraticVectorRestrictionWithProjection< stokesType::VelocityProjectionOperatorType >,
                              P2toP2QuadraticVectorProlongationWithProjection< stokesType::VelocityProjectionOperatorType > > >(
       parameters_,
       storage_,
       minLevel_,
       maxLevel_,
       VectorMassSmootherLeft,
       VectorMassCoarseSolver,
       VectorMassRestrictionOperator,
       VectorMassProlongationOperator,
       lowMemoryMode_,
       "VectorMass" );

   auto ABlockCGVectorMassSolverRight_ =
       std::make_shared< ABlockCGVectorMassSolver< stokesType::AOperatorType, VectorSqrtMassType > >(
           parameters_, storage_, minLevel_, maxLevel_, VectorMassRight, VectorMassMGSolverRight, lowMemoryMode_ );
   auto ABlockCGVectorMassSolverLeft_ =
       std::make_shared< ABlockCGVectorMassSolver< stokesType::AOperatorType, VectorSqrtMassType > >(
           parameters_, storage_, minLevel_, maxLevel_, VectorMassLeft, VectorMassMGSolverLeft, lowMemoryMode_ );

   // Poisson Neumann prec
   auto invSqrtEtaRight_ =
       std::make_shared< MCModel::ViscosityFunctionType >( "inv sqrt eta right", storage_, minLevel_, maxLevel_ );
   auto invSqrtEtaLeft_ =
       std::make_shared< MCModel::ViscosityFunctionType >( "inv sqrt eta left", storage_, minLevel_, maxLevel_ );

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      invSqrtEtaRight_->assign( { real_c( 1 ) }, { *sqrtEtaRight_ }, level, hyteg::All );
      invSqrtEtaRight_->invertElementwise( level, hyteg::All );

      invSqrtEtaLeft_->assign( { real_c( 1 ) }, { *sqrtEtaLeft_ }, level, hyteg::All );
      invSqrtEtaLeft_->invertElementwise( level, hyteg::All );
   }

   auto InvSqrtDiffusionRight = std::make_shared< InvSqrtViscMassType >( storage_, minLevel_, maxLevel_, *invSqrtEtaRight_ );
   InvSqrtDiffusionRight->computeInverseDiagonalOperatorValues();
   auto InvSqrtDiffusionLeft = std::make_shared< InvSqrtViscMassType >( storage_, minLevel_, maxLevel_, *invSqrtEtaLeft_ );
   InvSqrtDiffusionLeft->computeInverseDiagonalOperatorValues();

   auto WrappedDivKGradRight = std::make_shared< WrappedDivKGradType >( InvSqrtDiffusionRight, lowMemoryMode_ );
   auto WrappedDivKGradLeft  = std::make_shared< WrappedDivKGradType >( InvSqrtDiffusionLeft, lowMemoryMode_ );

   // Poisson Neumann Multigrid
   auto SchurPoissonNeumannChebyshevRight = std::make_shared< SchurPoissonNeumannChebyshevSmoother< WrappedDivKGradType > >(
       parameters_, storage_, minLevel_, maxLevel_, WrappedDivKGradRight, true, lowMemoryMode_ );
   auto SchurPoissonNeumannChebyshevLeft = std::make_shared< SchurPoissonNeumannChebyshevSmoother< WrappedDivKGradType > >(
       parameters_, storage_, minLevel_, maxLevel_, WrappedDivKGradLeft, true, lowMemoryMode_ );

   auto SchurPoissonNeumannRestrictionOperator =
       std::make_shared< P1toP1LinearRestrictionWithProjection< real_t > >( storage_, minLevel_, maxLevel_, lowMemoryMode_ );
   auto SchurPoissonNeumannProlongationOperator =
       std::make_shared< P1toP1LinearProlongationWithProjection< real_t > >( storage_, minLevel_, maxLevel_, lowMemoryMode_ );

   auto SchurCGPoissonNeumannCoarseGridSolverRight = std::make_shared< SchurCGPoissonNeumannSolver< WrappedDivKGradType > >(
       parameters_,
       storage_,
       minLevel_,
       minLevel_,
       WrappedDivKGradRight,
       lowMemoryMode_,
       std::make_shared< SchurIdentityPreconditioner< P1Function< real_t > > >(),
       "Coarse" );
   auto SchurCGPoissonNeumannCoarseGridSolverLeft = std::make_shared< SchurCGPoissonNeumannSolver< WrappedDivKGradType > >(
       parameters_,
       storage_,
       minLevel_,
       minLevel_,
       WrappedDivKGradLeft,
       lowMemoryMode_,
       std::make_shared< SchurIdentityPreconditioner< P1Function< real_t > > >(),
       "Coarse" );

   auto SchurCGPoissonNeumannMGSolverRight =
       std::make_shared< SchurPoissonNeumannMultigridSolver< WrappedDivKGradType,
                                                             P1toP1LinearRestrictionWithProjection< real_t >,
                                                             P1toP1LinearProlongationWithProjection< real_t > > >(
           parameters_,
           storage_,
           minLevel_,
           maxLevel_,
           WrappedDivKGradRight,
           SchurPoissonNeumannChebyshevRight,
           SchurCGPoissonNeumannCoarseGridSolverRight,
           SchurPoissonNeumannRestrictionOperator,
           SchurPoissonNeumannProlongationOperator,
           lowMemoryMode_ );
   auto SchurCGPoissonNeumannMGSolverLeft =
       std::make_shared< SchurPoissonNeumannMultigridSolver< WrappedDivKGradType,
                                                             P1toP1LinearRestrictionWithProjection< real_t >,
                                                             P1toP1LinearProlongationWithProjection< real_t > > >(
           parameters_,
           storage_,
           minLevel_,
           maxLevel_,
           WrappedDivKGradLeft,
           SchurPoissonNeumannChebyshevLeft,
           SchurCGPoissonNeumannCoarseGridSolverLeft,
           SchurPoissonNeumannRestrictionOperator,
           SchurPoissonNeumannProlongationOperator,
           lowMemoryMode_ );

   auto SchurCGPoissonNeumannSolverRight_ = std::make_shared< SchurCGPoissonNeumannSolver< WrappedDivKGradType > >(
       parameters_, storage_, minLevel_, maxLevel_, WrappedDivKGradRight, lowMemoryMode_ );
   auto SchurCGPoissonNeumannSolverLeft_ = std::make_shared< SchurCGPoissonNeumannSolver< WrappedDivKGradType > >(
       parameters_, storage_, minLevel_, maxLevel_, WrappedDivKGradLeft, lowMemoryMode_ );

   // ##################################################
   // ### Define main Schur complement approximation ###
   // ##################################################

   // asymmetrical BFBT type 1
   typedef BFBTOperator< stokesType::AOperatorType,
                         VectorSqrtMassType,
                         VectorSqrtMassType,
                         stokesType::BTOperatorType,
                         stokesType::BOperatorType >
       BFBTType1Right;
   typedef BFBTOperator< stokesType::AOperatorType,
                         VectorSqrtMassType,
                         VectorSqrtMassType,
                         stokesType::BTOperatorType,
                         stokesType::BOperatorType >
       BFBTType1Left;
   typedef SchurBFBTSolver< stokesType,
                            VectorSqrtMassType,
                            VectorSqrtMassType,
                            BFBTType1Right,
                            BFBTType1Left,
                            hyteg::FGMRESSolver< BFBTType1Right >,
                            hyteg::FGMRESSolver< BFBTType1Left > >
       BFBTSolverType1;

   // symmetrical BFBT type 1
   typedef BFBTOperator< stokesLHSTypeSym::AOperatorType,
                         VectorSqrtMassType,
                         VectorSqrtMassType,
                         stokesLHSTypeSym::BTOperatorType,
                         stokesLHSTypeSym::BOperatorType >
       BFBTType1RightSym;
   typedef BFBTOperator< stokesLHSTypeSym::AOperatorType,
                         VectorSqrtMassType,
                         VectorSqrtMassType,
                         stokesLHSTypeSym::BTOperatorType,
                         stokesLHSTypeSym::BOperatorType >
       BFBTType1LeftSym;
   typedef SchurBFBTSolver< stokesLHSTypeSym,
                            VectorSqrtMassType,
                            VectorSqrtMassType,
                            BFBTType1RightSym,
                            BFBTType1LeftSym,
                            hyteg::CGSolver< BFBTType1RightSym >,
                            hyteg::CGSolver< BFBTType1LeftSym > >
       BFBTSolverType1Sym;

   // asymmetrical BFBT type 2
   typedef SchurBFBTSolver< stokesType,
                            VectorSqrtMassType,
                            VectorSqrtMassType,
                            WrappedDivKGradType,
                            WrappedDivKGradType,
                            hyteg::CGSolver< WrappedDivKGradType >,
                            hyteg::CGSolver< WrappedDivKGradType > >
       BFBTSolverType2;

   // asymmetrical BFBT type 2
   typedef SchurBFBTSolver< stokesLHSTypeSym,
                            VectorSqrtMassType,
                            VectorSqrtMassType,
                            WrappedDivKGradType,
                            WrappedDivKGradType,
                            hyteg::CGSolver< WrappedDivKGradType >,
                            hyteg::CGSolver< WrappedDivKGradType > >
       BFBTSolverType2Sym;

   std::shared_ptr< SchurSolver< SchurOperator< P1Function< real_t > > > > SchurBFBTSolver_;

   if ( WBFBTType_ == 1 )
   {
      if ( AsymmetricPreconditionerMode == 0 )
      {
         SchurBFBTSolver_ = std::make_shared< BFBTSolverType1Sym >( parameters_,
                                                                    storage_,
                                                                    minLevel_,
                                                                    maxLevel_,
                                                                    stokesOperatorSym_,
                                                                    VectorMassRight,
                                                                    VectorMassLeft,
                                                                    ABlockCGVectorMassSolverRight_,
                                                                    ABlockCGVectorMassSolverLeft_,
                                                                    BC_.bcVelocityX_,
                                                                    BC_.bcVelocityY_,
                                                                    BC_.bcVelocityZ_,
                                                                    lowMemoryMode_,
                                                                    SchurCGPoissonNeumannMGSolverRight,
                                                                    nullptr,
                                                                    SchurCGPoissonNeumannMGSolverLeft,
                                                                    nullptr,
                                                                    "WBFBT_Type1" );
      }
      else
      {
         SchurBFBTSolver_ = std::make_shared< BFBTSolverType1 >( parameters_,
                                                                 storage_,
                                                                 minLevel_,
                                                                 maxLevel_,
                                                                 stokesOperator_,
                                                                 VectorMassRight,
                                                                 VectorMassLeft,
                                                                 ABlockCGVectorMassSolverRight_,
                                                                 ABlockCGVectorMassSolverLeft_,
                                                                 BC_.bcVelocityX_,
                                                                 BC_.bcVelocityY_,
                                                                 BC_.bcVelocityZ_,
                                                                 lowMemoryMode_,
                                                                 SchurCGPoissonNeumannMGSolverRight,
                                                                 nullptr,
                                                                 SchurCGPoissonNeumannMGSolverLeft,
                                                                 nullptr,
                                                                 "WBFBT_Type1" );
      }
   }
   else
   {
      if ( AsymmetricPreconditionerMode == 0 )
      {
         SchurBFBTSolver_ = std::make_shared< BFBTSolverType2Sym >( parameters_,
                                                                    storage_,
                                                                    minLevel_,
                                                                    maxLevel_,
                                                                    stokesOperatorSym_,
                                                                    VectorMassRight,
                                                                    VectorMassLeft,
                                                                    ABlockCGVectorMassSolverRight_,
                                                                    ABlockCGVectorMassSolverLeft_,
                                                                    BC_.bcVelocityX_,
                                                                    BC_.bcVelocityY_,
                                                                    BC_.bcVelocityZ_,
                                                                    lowMemoryMode_,
                                                                    SchurCGPoissonNeumannMGSolverRight,
                                                                    WrappedDivKGradRight,
                                                                    SchurCGPoissonNeumannMGSolverLeft,
                                                                    WrappedDivKGradLeft,
                                                                    "WBFBT_Type2" );
      }
      else
      {
         SchurBFBTSolver_ = std::make_shared< BFBTSolverType2 >( parameters_,
                                                                 storage_,
                                                                 minLevel_,
                                                                 maxLevel_,
                                                                 stokesOperator_,
                                                                 VectorMassRight,
                                                                 VectorMassLeft,
                                                                 ABlockCGVectorMassSolverRight_,
                                                                 ABlockCGVectorMassSolverLeft_,
                                                                 BC_.bcVelocityX_,
                                                                 BC_.bcVelocityY_,
                                                                 BC_.bcVelocityZ_,
                                                                 lowMemoryMode_,
                                                                 SchurCGPoissonNeumannMGSolverRight,
                                                                 WrappedDivKGradRight,
                                                                 SchurCGPoissonNeumannMGSolverLeft,
                                                                 WrappedDivKGradLeft,
                                                                 "WBFBT_Type2" );
      }
   }

   // ###########################
   // ### Define main solvers ###
   // ###########################

   // Stokes solver
   std::shared_ptr< SaddlePointSolver< stokesType > > SaddlePointPreconditioner;

   if ( AsymmetricPreconditionerMode == 0 )
   {
      std::shared_ptr< SaddlePointSolver< stokesLHSTypeSym > > SaddlePointPreconditionerSym;

      switch ( blockpreconditionerType_ )
      {
      case 0: {
         SaddlePointPreconditionerSym = std::make_shared< SaddlePointInexactUzawaSmoother< stokesLHSTypeSym > >(
             parameters_, storage_, minLevel_, maxLevel_, stokesOperatorSym_, ABlockOuterLoop, SchurBFBTSolver_, lowMemoryMode_ );
      }
      break;
      case 1: {
         SaddlePointPreconditionerSym = std::make_shared< SaddlePointAdjointInexactUzawaSmoother< stokesLHSTypeSym > >(
             parameters_, storage_, minLevel_, maxLevel_, stokesOperatorSym_, ABlockOuterLoop, SchurBFBTSolver_, lowMemoryMode_ );
      }
      break;
      case 2: {
         SaddlePointPreconditionerSym =
             std::make_shared< SaddlePointBlockApproximationFactorisationSmoother< stokesLHSTypeSym > >( parameters_,
                                                                                                         storage_,
                                                                                                         minLevel_,
                                                                                                         maxLevel_,
                                                                                                         stokesOperatorSym_,
                                                                                                         ABlockOuterLoop,
                                                                                                         SchurBFBTSolver_,
                                                                                                         lowMemoryMode_ );
      }
      break;
      case 3: {
         SaddlePointPreconditionerSym = std::make_shared< SaddlePointSymmetricUzawaSmoother< stokesLHSTypeSym > >(
             parameters_, storage_, minLevel_, maxLevel_, stokesOperatorSym_, ABlockOuterLoop, SchurBFBTSolver_, lowMemoryMode_ );
      }
      break;
      default: {
         WALBERLA_ABORT( "Unknown block preconditioner type " << blockpreconditionerType_ << "!" );
      }
      break;
      }

      SaddlePointPreconditioner = std::make_shared< SaddlePointSubstituteSolver< stokesLHSTypeSym, stokesType > >(
          stokesOperatorSym_, SaddlePointPreconditionerSym );
   }
   else
   {
      switch ( blockpreconditionerType_ )
      {
      case 0: {
         SaddlePointPreconditioner = std::make_shared< SaddlePointInexactUzawaSmoother< stokesType > >(
             parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, ABlockOuterLoop, SchurBFBTSolver_, lowMemoryMode_ );
      }
      break;
      case 1: {
         SaddlePointPreconditioner = std::make_shared< SaddlePointAdjointInexactUzawaSmoother< stokesType > >(
             parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, ABlockOuterLoop, SchurBFBTSolver_, lowMemoryMode_ );
      }
      break;
      case 2: {
         SaddlePointPreconditioner = std::make_shared< SaddlePointBlockApproximationFactorisationSmoother< stokesType > >(
             parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, ABlockOuterLoop, SchurBFBTSolver_, lowMemoryMode_ );
      }
      break;
      case 3: {
         SaddlePointPreconditioner = std::make_shared< SaddlePointSymmetricUzawaSmoother< stokesType > >(
             parameters_, storage_, minLevel_, maxLevel_, stokesOperator_, ABlockOuterLoop, SchurBFBTSolver_, lowMemoryMode_ );
      }
      break;
      default: {
         WALBERLA_ABORT( "Unknown block preconditioner type " << blockpreconditionerType_ << "!" );
      }
      break;
      }
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

   auto operatorUpdater_ = std::make_shared< OperatorUpdaterType >( stokesOperator_,
                                                                    nullptr,
                                                                    ABlockSmoother,
                                                                    nullptr,
                                                                    nullptr,
                                                                    Model,
                                                                    chebyshevUpdateRateMyrs,
                                                                    auxFct,
                                                                    minLevel_,
                                                                    maxLevel_,
                                                                    THistory,
                                                                    viscosityModel,
                                                                    sqrtEtaRight_,
                                                                    invSqrtEtaRight_,
                                                                    SqrtMassRight,
                                                                    InvSqrtDiffusionRight,
                                                                    VectorMassSmootherRight,
                                                                    SchurPoissonNeumannChebyshevRight,
                                                                    sqrtEtaLeft_,
                                                                    invSqrtEtaLeft_,
                                                                    SqrtMassLeft,
                                                                    InvSqrtDiffusionLeft,
                                                                    VectorMassSmootherLeft,
                                                                    SchurPoissonNeumannChebyshevLeft,
                                                                    maskRight,
                                                                    maskLeft );
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