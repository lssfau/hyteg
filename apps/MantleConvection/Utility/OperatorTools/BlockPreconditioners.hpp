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
#pragma once

#include "core/math/Random.h"

#include "hyteg/functions/PressureMeanProjection.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/SubstitutePreconditioner.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

#include "../Solver/ABlock/ABlockSolver.hpp"
#include "../Solver/Schur/SchurSolver.hpp"
#include "SchurOperator.hpp"

namespace MantleConvection {

// -------------------------------------------------------------
// -------------------------- Preface --------------------------
// -------------------------------------------------------------

// The block preconditioners in this header expect a MantleConvecion::SaddlePointOperator
//    {{A, B^T},{B, C}}
//
// NOTE: A is assumed to be positive semidefinite and C is assumed to be negative semidefinite.
// A possible stabilisation for the P1P1 case is the pspg stabilisation for which the form already
// contains a negative sign hence it is negative semidefinite.
//
// The operator in question is expected to operate on functions that provide
// access to their respective parts via *.uvw() and *.p() with typedefs
// VelocityFunction_T and PressureFunction_T for the respective types.
// Alos the source and destination function types of your operator should be the identical.
//
// Additionally the operator is expected to define the typedefs
//    AOperatorType             for the ABlock                 Operator,
//    BOperatorType             for the BBlock                 Operator,
//    BTOperatorType            for the BTBlock                Operator and
//    StabilisationOperatorType for the optional stabilisation Operator
//
// and provide access to its inner operators via the functions
//    AOperatorType&              getA()     const
//    BOperatorType&              getB()     const
//    BTOperatorType&             getBT()    const
//    StabilisationOperatorType&  getStab()  const
//
// If one of the blocks is 0, you can substitute the operator
// with the
//    NoOperator
// class found in src/hyteg/operators/.
//
// For more information about the block preconditioners see
// Daniel Drzisga, Lorenz John, Ulrich Rüde, Barbara Wohlmuth, and Walter Zulehner
// On the Analysis of Block Smoothers for Saddle Point Problems
// SIAM Journal on Matrix Analysis and Applications, Volume 39, Number 2, Pages 932-960, 2018
// https://doi.org/10.1137/16M1106304

// --------------------------------------------------------------
// ---------------------------- Note ----------------------------
// --------------------------------------------------------------

// The loops applying multiple velocity/pressure updated could also be achieved via loops of the form
//
// A.getBT().apply( x.p(), residual_.uvw(), level, flag_, Replace );
// residual_.uvw().assign( { 1.0, -1.0 }, { b.uvw(), residual_.uvw() }, level, flag_ );
//
// for ( uint_t i = 0; i < 2; i++ )
// {
//    ABlockSolver_->solve( A.getA(), x.uvw(), residual_.uvw(), level );
// }
//
// but here we want to specifically enforce e.g. the usage of A.getA() for the residual calculation
// (some solver wrappers just ignore the given operator and use something else).

// -------------------------------------------------------------
// ---------------------- Preconditioners ----------------------
// -------------------------------------------------------------

template < class OperatorType,
           bool projectASolverRHS_     = true,
           bool projectSchurSolverRHS_ = true,
           bool postProjectPressure_   = true,
           bool postProjectVelocity_   = true >
class InexactUzawaPreconditioner : public hyteg::Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType                     SrcFunctionType;
   typedef typename OperatorType::AOperatorType               AOperatorType;
   typedef typename OperatorType::srcType::PressureFunction_T PressureFunctionType;

   InexactUzawaPreconditioner(
       const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
       const uint_t                                                                   minLevel,
       const uint_t                                                                   maxLevel,
       const std::shared_ptr< OperatorType >&                                         SaddleOp,
       const std::shared_ptr< ABlockSolver< AOperatorType > >&                        ABlockSolver,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& SchurComplementSolver,
       real_t                                                                         relaxParamA,
       real_t                                                                         relaxParamSchur,
       bool                                                                           lowMemoryMode      = false,
       uint_t                                                                         VelocityIterations = 1,
       hyteg::DoFType flag = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , schurOp_( storage, minLevel, maxLevel )
   , ABlockSolver_( ABlockSolver )
   , SchurComplementSolver_( SchurComplementSolver )
   , flag_( flag )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParamA_( relaxParamA )
   , relaxParamSchur_( relaxParamSchur )
   , VelocityIterations_( VelocityIterations )
   , lowMemoryMode_( lowMemoryMode )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      if ( SaddleOp == nullptr )
      {
         WALBERLA_ABORT( "Operator is nullptr!" );
      }

      if ( ABlockSolver == nullptr )
      {
         WALBERLA_ABORT( "ABlockSolver is nullptr!" );
      }

      if ( SchurComplementSolver == nullptr )
      {
         WALBERLA_ABORT( "SchurComplementSolver is nullptr!" );
      }

      // A Block
      if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getAPtr() == nullptr )
         {
            WALBERLA_ABORT( "A Block set to nullptr but type is not NoOperator!" );
         }
      }

      // BT Block
      if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getBTPtr() == nullptr )
         {
            WALBERLA_ABORT( "BT Block set to nullptr but type is not NoOperator!" );
         }
      }

      // B Block
      if constexpr ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getBPtr() == nullptr )
         {
            WALBERLA_ABORT( "B Block set to nullptr but type is not NoOperator!" );
         }
      }

      // Stabilisation Block
      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getStabPtr() == nullptr )
         {
            WALBERLA_ABORT( "Stabilisation Block set to nullptr but type is not NoOperator!" );
         }
      }

      // Projection
      if constexpr ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getProjPtr() == nullptr )
         {
            WALBERLA_ABORT( "Velocity projection set to nullptr but type is not NoOperator!" );
         }
         // set projection flag
         projectionFlag_ = SaddleOp->getProjFlag();
      }

      // init functions
      if ( !lowMemoryMode_ )
      {
         residual_ = std::make_shared< SrcFunctionType >( "InexactUzawaPreconditioner residual", storage, minLevel, maxLevel );
         tmp_      = std::make_shared< SrcFunctionType >( "InexactUzawaPreconditioner tmp", storage, minLevel, maxLevel );
      }
   }

   virtual void solve( const OperatorType&                   A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const uint_t                          level ) override
   {
      std::shared_ptr< SrcFunctionType > residualSolve;
      std::shared_ptr< SrcFunctionType > tmpSolve;

      if ( !lowMemoryMode_ )
      {
         residualSolve = residual_;
         tmpSolve      = tmp_;
      }
      else
      {
         residualSolve = hyteg::getTemporaryFunction< SrcFunctionType >( storage_, minLevel_, maxLevel_ );
         tmpSolve      = hyteg::getTemporaryFunction< SrcFunctionType >( storage_, minLevel_, maxLevel_ );
      }

      residualSolve->copyBoundaryConditionFromFunction( x );
      tmpSolve->copyBoundaryConditionFromFunction( x );

      // Pre projections are expected to be handled via the projection wrappers

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual

         if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getBT().apply( x.p(), residualSolve->uvw(), level, flag_, hyteg::Replace );
         }
         if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getA().apply( x.uvw(),
                            residualSolve->uvw(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
         }

         residualSolve->uvw().assign( { real_c( 1 ), real_c( -1 ) }, { b.uvw(), residualSolve->uvw() }, level, flag_ );

         // Usually we can assume that b is already projected if the user correctly applied the RHS operator
         // and applied the projection after prolongation / restriction. However if you do not have any previous
         // knowledge about the nature of b you can set projectASolverRHS_ to true.
         if constexpr ( ( projectASolverRHS_ ) &&
                        ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
         {
            A.getProj().project( residualSolve->uvw(), level, projectionFlag_ );
         }

         // apply A_hat
         tmpSolve->uvw().setToZero( level );
         ABlockSolver_->solve( A.getA(), tmpSolve->uvw(), residualSolve->uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1 ), relaxParamA_ }, { x.uvw(), tmpSolve->uvw() }, level, flag_ );
      }

      // Projecting x.uvw() before applying the B Block is expected to be handled by the projection wrappers

      // --------------------------Pressure--------------------------

      // calculate residual
      if constexpr ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         A.getB().apply( x.uvw(), residualSolve->p(), level, flag_, hyteg::Replace );
      }
      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         A.getStab().apply( x.p(),
                            residualSolve->p(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
      }

      residualSolve->p().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.p(), residualSolve->p() }, level, flag_ );

      // Usually we can assume that b is already projected if the user correctly applied the RHS operator
      // and applied the projection after prolongation / restriction. However if you do not have any previous
      // knowledge about the nature of b you can set projectSchurSolverRHS_ to true.
      if constexpr ( projectSchurSolverRHS_ )
      {
         hyteg::projectPressureMean( residualSolve->p(), level );
      }

      // apply S_hat
      tmpSolve->p().setToZero( level );
      SchurComplementSolver_->solve( schurOp_, tmpSolve->p(), residualSolve->p(), level );

      // update x with the scaled result
      x.p().assign( { real_c( 1 ), -relaxParamSchur_ }, { x.p(), tmpSolve->p() }, level, flag_ );

      // ------------------------Post Projection------------------------

      // Usually we can assume that x is already projected if the user correctly uses a projected initial guess for x.
      // However if you do not have any previous knowledge about the nature of x you can set postProjectPressure_ and postProjectVelocity_ to true.
      if constexpr ( postProjectPressure_ )
      {
         hyteg::projectPressureMean( x.p(), level );
      }

      if constexpr ( ( postProjectVelocity_ ) &&
                     ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
      {
         A.getProj().project( x.uvw(), level, projectionFlag_ );
      }
   }

 protected:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   SchurOperator< PressureFunctionType > schurOp_;

   std::shared_ptr< ABlockSolver< AOperatorType > >                        ABlockSolver_;
   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > SchurComplementSolver_;

   hyteg::DoFType flag_;
   hyteg::DoFType projectionFlag_;
   bool           hasGlobalCells_;

   real_t                             relaxParamA_;
   real_t                             relaxParamSchur_;
   uint_t                             VelocityIterations_;
   std::shared_ptr< SrcFunctionType > residual_;
   std::shared_ptr< SrcFunctionType > tmp_;

   bool lowMemoryMode_;
};

template < class OperatorType,
           bool projectASolverRHS_     = true,
           bool projectSchurSolverRHS_ = true,
           bool postProjectPressure_   = true,
           bool postProjectVelocity_   = true >
class AdjointInexactUzawaPreconditioner : public hyteg::Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType                     SrcFunctionType;
   typedef typename OperatorType::AOperatorType               AOperatorType;
   typedef typename OperatorType::srcType::PressureFunction_T PressureFunctionType;

   AdjointInexactUzawaPreconditioner(
       const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
       const uint_t                                                                   minLevel,
       const uint_t                                                                   maxLevel,
       const std::shared_ptr< OperatorType >&                                         SaddleOp,
       const std::shared_ptr< ABlockSolver< AOperatorType > >&                        ABlockSolver,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& SchurComplementSolver,
       real_t                                                                         relaxParamA,
       real_t                                                                         relaxParamSchur,
       bool                                                                           lowMemoryMode      = false,
       uint_t                                                                         VelocityIterations = 1,
       hyteg::DoFType flag = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , schurOp_( storage, minLevel, maxLevel )
   , ABlockSolver_( ABlockSolver )
   , SchurComplementSolver_( SchurComplementSolver )
   , flag_( flag )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParamA_( relaxParamA )
   , relaxParamSchur_( relaxParamSchur )
   , VelocityIterations_( VelocityIterations )
   , lowMemoryMode_( lowMemoryMode )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      if ( SaddleOp == nullptr )
      {
         WALBERLA_ABORT( "Operator is nullptr!" );
      }

      if ( ABlockSolver == nullptr )
      {
         WALBERLA_ABORT( "ABlockSolver is nullptr!" );
      }

      if ( SchurComplementSolver == nullptr )
      {
         WALBERLA_ABORT( "SchurComplementSolver is nullptr!" );
      }

      // A Block
      if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getAPtr() == nullptr )
         {
            WALBERLA_ABORT( "A Block set to nullptr but type is not NoOperator!" );
         }
      }

      // BT Block
      if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getBTPtr() == nullptr )
         {
            WALBERLA_ABORT( "BT Block set to nullptr but type is not NoOperator!" );
         }
      }

      // B Block
      if constexpr ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getBPtr() == nullptr )
         {
            WALBERLA_ABORT( "B Block set to nullptr but type is not NoOperator!" );
         }
      }

      // Stabilisation Block
      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getStabPtr() == nullptr )
         {
            WALBERLA_ABORT( "Stabilisation Block set to nullptr but type is not NoOperator!" );
         }
      }

      // Projection
      if constexpr ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getProjPtr() == nullptr )
         {
            WALBERLA_ABORT( "Velocity projection set to nullptr but type is not NoOperator!" );
         }
         // set projection flag
         projectionFlag_ = SaddleOp->getProjFlag();
      }

      // init functions
      if ( !lowMemoryMode_ )
      {
         residual_ =
             std::make_shared< SrcFunctionType >( "AdjointInexactUzawaPreconditioner residual", storage, minLevel, maxLevel );
         tmp_ = std::make_shared< SrcFunctionType >( "AdjointInexactUzawaPreconditioner tmp", storage, minLevel, maxLevel );
      }
   }

   virtual void solve( const OperatorType&                   A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const uint_t                          level ) override
   {
      std::shared_ptr< SrcFunctionType > residualSolve;
      std::shared_ptr< SrcFunctionType > tmpSolve;

      if ( !lowMemoryMode_ )
      {
         residualSolve = residual_;
         tmpSolve      = tmp_;
      }
      else
      {
         residualSolve = hyteg::getTemporaryFunction< SrcFunctionType >( storage_, minLevel_, maxLevel_ );
         tmpSolve      = hyteg::getTemporaryFunction< SrcFunctionType >( storage_, minLevel_, maxLevel_ );
      }

      residualSolve->copyBoundaryConditionFromFunction( x );
      tmpSolve->copyBoundaryConditionFromFunction( x );

      // Pre projections are expected to be handled via the projection wrappers

      // --------------------------Pressure--------------------------

      // calculate residual
      if constexpr ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         A.getB().apply( x.uvw(), residualSolve->p(), level, flag_, hyteg::Replace );
      }
      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         A.getStab().apply( x.p(),
                            residualSolve->p(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
      }

      residualSolve->p().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.p(), residualSolve->p() }, level, flag_ );

      // Usually we can assume that b is already projected if the user correctly applied the RHS operator
      // and applied the projection after prolongation / restriction. However if you do not have any previous
      // knowledge about the nature of b you can set projectSchurSolverRHS_ to true.
      if constexpr ( projectSchurSolverRHS_ )
      {
         hyteg::projectPressureMean( residualSolve->p(), level );
      }

      // apply S_hat
      tmpSolve->p().setToZero( level );
      SchurComplementSolver_->solve( schurOp_, tmpSolve->p(), residualSolve->p(), level );

      // update x with the scaled result
      x.p().assign( { real_c( 1 ), -relaxParamSchur_ }, { x.p(), tmpSolve->p() }, level, flag_ );

      // Projecting x.p() before applying the TB Block is expected to be handled by the projection wrappers

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual

         if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getBT().apply( x.p(), residualSolve->uvw(), level, flag_, hyteg::Replace );
         }
         if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getA().apply( x.uvw(),
                            residualSolve->uvw(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
         }

         residualSolve->uvw().assign( { real_c( 1 ), real_c( -1 ) }, { b.uvw(), residualSolve->uvw() }, level, flag_ );

         // Usually we can assume that b is already projected if the user correctly applied the RHS operator
         // and applied the projection after prolongation / restriction. However if you do not have any previous
         // knowledge about the nature of b you can set projectASolverRHS_ to true.
         if constexpr ( ( projectASolverRHS_ ) &&
                        ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
         {
            A.getProj().project( residualSolve->uvw(), level, projectionFlag_ );
         }

         // apply A_hat
         tmpSolve->uvw().setToZero( level );
         ABlockSolver_->solve( A.getA(), tmpSolve->uvw(), residualSolve->uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1 ), relaxParamA_ }, { x.uvw(), tmpSolve->uvw() }, level, flag_ );
      }

      // ------------------------Post Projection------------------------

      // Usually we can assume that x is already projected if the user correctly uses a projected initial guess for x.
      // However if you do not have any previous knowledge about the nature of x you can set postProjectPressure_ and postProjectVelocity_ to true.
      if constexpr ( postProjectPressure_ )
      {
         hyteg::projectPressureMean( x.p(), level );
      }

      if constexpr ( ( postProjectVelocity_ ) &&
                     ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
      {
         A.getProj().project( x.uvw(), level, projectionFlag_ );
      }
   }

 protected:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   SchurOperator< PressureFunctionType > schurOp_;

   std::shared_ptr< ABlockSolver< AOperatorType > >                        ABlockSolver_;
   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > SchurComplementSolver_;

   hyteg::DoFType flag_;
   hyteg::DoFType projectionFlag_;
   bool           hasGlobalCells_;

   real_t                             relaxParamA_;
   real_t                             relaxParamSchur_;
   uint_t                             VelocityIterations_;
   std::shared_ptr< SrcFunctionType > residual_;
   std::shared_ptr< SrcFunctionType > tmp_;

   bool lowMemoryMode_;
};

template < class OperatorType,
           bool projectASolverRHS_     = true,
           bool projectSchurSolverRHS_ = true,
           bool postProjectPressure_   = true,
           bool postProjectVelocity_   = true >
class BlockApproximateFactorisationPreconditioner : public hyteg::Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType                     SrcFunctionType;
   typedef typename OperatorType::AOperatorType               AOperatorType;
   typedef typename OperatorType::srcType::VelocityFunction_T VelocityFunctionType;
   typedef typename OperatorType::srcType::PressureFunction_T PressureFunctionType;

   BlockApproximateFactorisationPreconditioner(
       const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
       const uint_t                                                                   minLevel,
       const uint_t                                                                   maxLevel,
       const std::shared_ptr< OperatorType >&                                         SaddleOp,
       const std::shared_ptr< ABlockSolver< AOperatorType > >&                        ABlockSolver,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& SchurComplementSolver,
       real_t                                                                         relaxParamA,
       real_t                                                                         relaxParamSchur,
       bool                                                                           lowMemoryMode      = false,
       uint_t                                                                         VelocityIterations = 1,
       hyteg::DoFType flag = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , schurOp_( storage, minLevel, maxLevel )
   , ABlockSolver_( ABlockSolver )
   , SchurComplementSolver_( SchurComplementSolver )
   , flag_( flag )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParamA_( relaxParamA )
   , relaxParamSchur_( relaxParamSchur )
   , VelocityIterations_( VelocityIterations )
   , lowMemoryMode_( lowMemoryMode )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      if ( SaddleOp == nullptr )
      {
         WALBERLA_ABORT( "Operator is nullptr!" );
      }

      if ( ABlockSolver == nullptr )
      {
         WALBERLA_ABORT( "ABlockSolver is nullptr!" );
      }

      if ( SchurComplementSolver == nullptr )
      {
         WALBERLA_ABORT( "SchurComplementSolver is nullptr!" );
      }

      // A Block
      if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getAPtr() == nullptr )
         {
            WALBERLA_ABORT( "A Block set to nullptr but type is not NoOperator!" );
         }
      }

      // BT Block
      if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getBTPtr() == nullptr )
         {
            WALBERLA_ABORT( "BT Block set to nullptr but type is not NoOperator!" );
         }
      }

      // B Block
      if constexpr ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getBPtr() == nullptr )
         {
            WALBERLA_ABORT( "B Block set to nullptr but type is not NoOperator!" );
         }
      }

      // Stabilisation Block
      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getStabPtr() == nullptr )
         {
            WALBERLA_ABORT( "Stabilisation Block set to nullptr but type is not NoOperator!" );
         }
      }

      // Projection
      if constexpr ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getProjPtr() == nullptr )
         {
            WALBERLA_ABORT( "Velocity projection set to nullptr but type is not NoOperator!" );
         }
         // set projection flag
         projectionFlag_ = SaddleOp->getProjFlag();
      }

      // init functions
      if ( !lowMemoryMode_ )
      {
         residual_ = std::make_shared< SrcFunctionType >(
             "BlockApproximateFactorisationPreconditioner residual", storage, minLevel, maxLevel );
         tmp_ = std::make_shared< SrcFunctionType >(
             "BlockApproximateFactorisationPreconditioner tmp", storage, minLevel, maxLevel );
         store_ = std::make_shared< VelocityFunctionType >(
             "BlockApproximateFactorisationPreconditioner store", storage, minLevel, maxLevel );
      }
   }

   virtual void solve( const OperatorType&                   A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const uint_t                          level ) override
   {
      std::shared_ptr< SrcFunctionType >      residualSolve;
      std::shared_ptr< SrcFunctionType >      tmpSolve;
      std::shared_ptr< VelocityFunctionType > storeSolve;

      if ( !lowMemoryMode_ )
      {
         residualSolve = residual_;
         tmpSolve      = tmp_;
         storeSolve    = store_;
      }
      else
      {
         residualSolve = hyteg::getTemporaryFunction< SrcFunctionType >( storage_, minLevel_, maxLevel_ );
         tmpSolve      = hyteg::getTemporaryFunction< SrcFunctionType >( storage_, minLevel_, maxLevel_ );
         storeSolve    = hyteg::getTemporaryFunction< VelocityFunctionType >( storage_, minLevel_, maxLevel_ );
      }

      residualSolve->copyBoundaryConditionFromFunction( x );
      tmpSolve->copyBoundaryConditionFromFunction( x );
      storeSolve->copyBoundaryConditionFromFunction( x.uvw() );

      // Pre projections are expected to be handled via the projection wrappers

      // --------------------------Store Velocity--------------------------

      // we store the velocity for later
      storeSolve->assign( { real_c( 1 ) }, { x.uvw() }, level, flag_ );

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual

         if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getBT().apply( x.p(), residualSolve->uvw(), level, flag_, hyteg::Replace );
         }
         if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getA().apply( x.uvw(),
                            residualSolve->uvw(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
         }

         residualSolve->uvw().assign( { real_c( 1 ), real_c( -1 ) }, { b.uvw(), residualSolve->uvw() }, level, flag_ );

         // Usually we can assume that b is already projected if the user correctly applied the RHS operator
         // and applied the projection after prolongation / restriction. However if you do not have any previous
         // knowledge about the nature of b you can set projectASolverRHS_ to true.
         if constexpr ( ( projectASolverRHS_ ) &&
                        ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
         {
            A.getProj().project( residualSolve->uvw(), level, projectionFlag_ );
         }

         // apply A_hat
         tmpSolve->uvw().setToZero( level );
         ABlockSolver_->solve( A.getA(), tmpSolve->uvw(), residualSolve->uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1 ), relaxParamA_ }, { x.uvw(), tmpSolve->uvw() }, level, flag_ );
      }

      // Projecting x.uvw() before applying the B Block is expected to be handled by the projection wrappers

      // --------------------------Pressure--------------------------

      // calculate residual
      if constexpr ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         A.getB().apply( x.uvw(), residualSolve->p(), level, flag_, hyteg::Replace );
      }
      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         A.getStab().apply( x.p(),
                            residualSolve->p(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
      }

      residualSolve->p().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.p(), residualSolve->p() }, level, flag_ );

      // Usually we can assume that b is already projected if the user correctly applied the RHS operator
      // and applied the projection after prolongation / restriction. However if you do not have any previous
      // knowledge about the nature of b you can set projectSchurSolverRHS_ to true.
      if constexpr ( projectSchurSolverRHS_ )
      {
         hyteg::projectPressureMean( residualSolve->p(), level );
      }

      // apply S_hat
      tmpSolve->p().setToZero( level );
      SchurComplementSolver_->solve( schurOp_, tmpSolve->p(), residualSolve->p(), level );

      // update x with the scaled result
      x.p().assign( { real_c( 1 ), -relaxParamSchur_ }, { x.p(), tmpSolve->p() }, level, flag_ );

      // --------------------------Restore Velocity--------------------------

      // we restore the velocity from earlier
      x.uvw().assign( { real_c( 1 ) }, { *storeSolve }, level, flag_ );

      // Projecting x.p() and x.uvw() before applying the BT Block and A Block is expected to be handled by the projection wrappers

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual

         if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getBT().apply( x.p(), residualSolve->uvw(), level, flag_, hyteg::Replace );
         }
         if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getA().apply( x.uvw(),
                            residualSolve->uvw(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
         }

         residualSolve->uvw().assign( { real_c( 1 ), real_c( -1 ) }, { b.uvw(), residualSolve->uvw() }, level, flag_ );

         // Usually we can assume that b is already projected if the user correctly applied the RHS operator
         // and applied the projection after prolongation / restriction. However if you do not have any previous
         // knowledge about the nature of b you can set projectASolverRHS_ to true.
         if constexpr ( ( projectASolverRHS_ ) &&
                        ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
         {
            A.getProj().project( residualSolve->uvw(), level, projectionFlag_ );
         }

         // apply A_hat
         tmpSolve->uvw().setToZero( level );
         ABlockSolver_->solve( A.getA(), tmpSolve->uvw(), residualSolve->uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1 ), relaxParamA_ }, { x.uvw(), tmpSolve->uvw() }, level, flag_ );
      }

      // ------------------------Post Projection------------------------

      // Usually we can assume that x is already projected if the user correctly uses a projected initial guess for x.
      // However if you do not have any previous knowledge about the nature of x you can set postProjectPressure_ and postProjectVelocity_ to true.
      if constexpr ( postProjectPressure_ )
      {
         hyteg::projectPressureMean( x.p(), level );
      }

      if constexpr ( ( postProjectVelocity_ ) &&
                     ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
      {
         A.getProj().project( x.uvw(), level, projectionFlag_ );
      }
   }

 protected:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   SchurOperator< PressureFunctionType > schurOp_;

   std::shared_ptr< ABlockSolver< AOperatorType > >                        ABlockSolver_;
   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > SchurComplementSolver_;

   hyteg::DoFType flag_;
   hyteg::DoFType projectionFlag_;
   bool           hasGlobalCells_;

   real_t                                  relaxParamA_;
   real_t                                  relaxParamSchur_;
   uint_t                                  VelocityIterations_;
   std::shared_ptr< SrcFunctionType >      residual_;
   std::shared_ptr< SrcFunctionType >      tmp_;
   std::shared_ptr< VelocityFunctionType > store_;

   bool lowMemoryMode_;
};

template < class OperatorType,
           bool projectASolverRHS_     = true,
           bool projectSchurSolverRHS_ = true,
           bool postProjectPressure_   = true,
           bool postProjectVelocity_   = true >
class SymmetricUzawaPreconditioner : public hyteg::Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType                     SrcFunctionType;
   typedef typename OperatorType::AOperatorType               AOperatorType;
   typedef typename OperatorType::srcType::PressureFunction_T PressureFunctionType;

   SymmetricUzawaPreconditioner(
       const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
       const uint_t                                                                   minLevel,
       const uint_t                                                                   maxLevel,
       const std::shared_ptr< OperatorType >&                                         SaddleOp,
       const std::shared_ptr< ABlockSolver< AOperatorType > >&                        ABlockSolver,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& SchurComplementSolver,
       real_t                                                                         relaxParamA,
       real_t                                                                         relaxParamSchur,
       bool                                                                           lowMemoryMode      = false,
       uint_t                                                                         VelocityIterations = 1,
       hyteg::DoFType flag = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , schurOp_( storage, minLevel, maxLevel )
   , ABlockSolver_( ABlockSolver )
   , SchurComplementSolver_( SchurComplementSolver )
   , flag_( flag )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParamA_( relaxParamA )
   , relaxParamSchur_( relaxParamSchur )
   , VelocityIterations_( VelocityIterations )
   , lowMemoryMode_( lowMemoryMode )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      if ( SaddleOp == nullptr )
      {
         WALBERLA_ABORT( "Operator is nullptr!" );
      }

      if ( ABlockSolver == nullptr )
      {
         WALBERLA_ABORT( "ABlockSolver is nullptr!" );
      }

      if ( SchurComplementSolver == nullptr )
      {
         WALBERLA_ABORT( "SchurComplementSolver is nullptr!" );
      }

      // A Block
      if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getAPtr() == nullptr )
         {
            WALBERLA_ABORT( "A Block set to nullptr but type is not NoOperator!" );
         }
      }

      // BT Block
      if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getBTPtr() == nullptr )
         {
            WALBERLA_ABORT( "BT Block set to nullptr but type is not NoOperator!" );
         }
      }

      // B Block
      if constexpr ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getBPtr() == nullptr )
         {
            WALBERLA_ABORT( "B Block set to nullptr but type is not NoOperator!" );
         }
      }

      // Stabilisation Block
      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getStabPtr() == nullptr )
         {
            WALBERLA_ABORT( "Stabilisation Block set to nullptr but type is not NoOperator!" );
         }
      }

      // Projection
      if constexpr ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getProjPtr() == nullptr )
         {
            WALBERLA_ABORT( "Velocity projection set to nullptr but type is not NoOperator!" );
         }
         // set projection flag
         projectionFlag_ = SaddleOp->getProjFlag();
      }

      // init functions
      if ( !lowMemoryMode_ )
      {
         residual_ = std::make_shared< SrcFunctionType >( "SymmetricUzawaPreconditioner residual", storage, minLevel, maxLevel );
         tmp_      = std::make_shared< SrcFunctionType >( "SymmetricUzawaPreconditioner tmp", storage, minLevel, maxLevel );
      }
   }

   virtual void solve( const OperatorType&                   A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const uint_t                          level ) override
   {
      std::shared_ptr< SrcFunctionType > residualSolve;
      std::shared_ptr< SrcFunctionType > tmpSolve;

      if ( !lowMemoryMode_ )
      {
         residualSolve = residual_;
         tmpSolve      = tmp_;
      }
      else
      {
         residualSolve = hyteg::getTemporaryFunction< SrcFunctionType >( storage_, minLevel_, maxLevel_ );
         tmpSolve      = hyteg::getTemporaryFunction< SrcFunctionType >( storage_, minLevel_, maxLevel_ );
      }

      residualSolve->copyBoundaryConditionFromFunction( x );
      tmpSolve->copyBoundaryConditionFromFunction( x );

      // Pre projections are expected to be handled via the projection wrappers

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual

         if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getBT().apply( x.p(), residualSolve->uvw(), level, flag_, hyteg::Replace );
         }
         if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getA().apply( x.uvw(),
                            residualSolve->uvw(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
         }

         residualSolve->uvw().assign( { real_c( 1 ), real_c( -1 ) }, { b.uvw(), residualSolve->uvw() }, level, flag_ );

         // Usually we can assume that b is already projected if the user correctly applied the RHS operator
         // and applied the projection after prolongation / restriction. However if you do not have any previous
         // knowledge about the nature of b you can set projectASolverRHS_ to true.
         if constexpr ( ( projectASolverRHS_ ) &&
                        ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
         {
            A.getProj().project( residualSolve->uvw(), level, projectionFlag_ );
         }

         // apply A_hat
         tmpSolve->uvw().setToZero( level );
         ABlockSolver_->solve( A.getA(), tmpSolve->uvw(), residualSolve->uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1 ), relaxParamA_ }, { x.uvw(), tmpSolve->uvw() }, level, flag_ );
      }

      // Projecting x.uvw() before applying the B Block is expected to be handled by the projection wrappers

      // --------------------------Pressure--------------------------

      // calculate residual
      if constexpr ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         A.getB().apply( x.uvw(), residualSolve->p(), level, flag_, hyteg::Replace );
      }
      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         A.getStab().apply( x.p(),
                            residualSolve->p(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
      }

      residualSolve->p().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.p(), residualSolve->p() }, level, flag_ );

      // Usually we can assume that b is already projected if the user correctly applied the RHS operator
      // and applied the projection after prolongation / restriction. However if you do not have any previous
      // knowledge about the nature of b you can set projectSchurSolverRHS_ to true.
      if constexpr ( projectSchurSolverRHS_ )
      {
         hyteg::projectPressureMean( residualSolve->p(), level );
      }

      // apply S_hat
      tmpSolve->p().setToZero( level );
      SchurComplementSolver_->solve( schurOp_, tmpSolve->p(), residualSolve->p(), level );

      // update x with the scaled result
      x.p().assign( { real_c( 1 ), -relaxParamSchur_ }, { x.p(), tmpSolve->p() }, level, flag_ );

      // Projecting x.p() and x.uvw() before applying the BT Block and A Block is expected to be handled by the projection wrappers

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual

         if constexpr ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getBT().apply( x.p(), residualSolve->uvw(), level, flag_, hyteg::Replace );
         }
         if constexpr ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value )
         {
            A.getA().apply( x.uvw(),
                            residualSolve->uvw(),
                            level,
                            flag_,
                            ( ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value ) ?
                                  hyteg::UpdateType::Add :
                                  hyteg::UpdateType::Replace ) );
         }

         residualSolve->uvw().assign( { real_c( 1 ), real_c( -1 ) }, { b.uvw(), residualSolve->uvw() }, level, flag_ );

         // Usually we can assume that b is already projected if the user correctly applied the RHS operator
         // and applied the projection after prolongation / restriction. However if you do not have any previous
         // knowledge about the nature of b you can set projectASolverRHS_ to true.
         if constexpr ( ( projectASolverRHS_ ) &&
                        ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
         {
            A.getProj().project( residualSolve->uvw(), level, projectionFlag_ );
         }

         // apply A_hat
         tmpSolve->uvw().setToZero( level );
         ABlockSolver_->solve( A.getA(), tmpSolve->uvw(), residualSolve->uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1 ), relaxParamA_ }, { x.uvw(), tmpSolve->uvw() }, level, flag_ );
      }

      // ------------------------Post Projection------------------------

      // Usually we can assume that x is already projected if the user correctly uses a projected initial guess for x.
      // However if you do not have any previous knowledge about the nature of x you can set postProjectPressure_ and postProjectVelocity_ to true.
      if constexpr ( postProjectPressure_ )
      {
         hyteg::projectPressureMean( x.p(), level );
      }

      if constexpr ( ( postProjectVelocity_ ) &&
                     ( !std::is_same< typename OperatorType::VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
      {
         A.getProj().project( x.uvw(), level, projectionFlag_ );
      }
   }

 protected:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   SchurOperator< PressureFunctionType > schurOp_;

   std::shared_ptr< ABlockSolver< AOperatorType > >                        ABlockSolver_;
   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > SchurComplementSolver_;

   hyteg::DoFType flag_;
   hyteg::DoFType projectionFlag_;
   bool           hasGlobalCells_;

   real_t                             relaxParamA_;
   real_t                             relaxParamSchur_;
   uint_t                             VelocityIterations_;
   std::shared_ptr< SrcFunctionType > residual_;
   std::shared_ptr< SrcFunctionType > tmp_;

   bool lowMemoryMode_;
};

} // namespace MantleConvection
