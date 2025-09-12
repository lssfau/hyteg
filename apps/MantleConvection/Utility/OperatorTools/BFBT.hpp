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

#include "hyteg/functions/PressureMeanProjection.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/SubstitutePreconditioner.hpp"

namespace MantleConvection {

template < class AOperatorType_,
           class AApproximationOperatorTypeRight_,
           class AApproximationOperatorTypeLeft_,
           class BTOperatorType_,
           class BOperatorType_ >
class BFBTOperator : public hyteg::Operator< typename BTOperatorType_::srcType, typename BOperatorType_::dstType >
{
 public:
   BFBTOperator( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                 uint_t                                            minLevel,
                 uint_t                                            maxLevel,

                 const std::shared_ptr< AOperatorType_ >&                   AOperator,
                 const std::shared_ptr< AApproximationOperatorTypeRight_ >& AApproximationOperatorRight,
                 const std::shared_ptr< AApproximationOperatorTypeLeft_ >&  AApproximationOperatorLeft,
                 const std::shared_ptr< BTOperatorType_ >&                  BTOperator,
                 const std::shared_ptr< BOperatorType_ >&                   BOperator,

                 const std::shared_ptr< hyteg::Solver< AApproximationOperatorTypeRight_ > >& AApproximationSolverRight,
                 const std::shared_ptr< hyteg::Solver< AApproximationOperatorTypeLeft_ > >&  AApproximationSolverLeft,

                 hyteg::BoundaryCondition VelocityBC,

                 bool             lowMemoryMode = false,
                 hyteg::BFBTState state         = hyteg::BFBTState::InnerBFBT )
   : BFBTOperator( storage,
                   minLevel,
                   maxLevel,
                   AOperator,
                   AApproximationOperatorRight,
                   AApproximationOperatorLeft,
                   BTOperator,
                   BOperator,
                   AApproximationSolverRight,
                   AApproximationSolverLeft,
                   VelocityBC,
                   VelocityBC,
                   VelocityBC,
                   lowMemoryMode,
                   state )
   {}

   BFBTOperator( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                 uint_t                                            minLevel,
                 uint_t                                            maxLevel,

                 const std::shared_ptr< AOperatorType_ >&                   AOperator,
                 const std::shared_ptr< AApproximationOperatorTypeRight_ >& AApproximationOperatorRight,
                 const std::shared_ptr< AApproximationOperatorTypeLeft_ >&  AApproximationOperatorLeft,
                 const std::shared_ptr< BTOperatorType_ >&                  BTOperator,
                 const std::shared_ptr< BOperatorType_ >&                   BOperator,

                 const std::shared_ptr< hyteg::Solver< AApproximationOperatorTypeRight_ > >& AApproximationSolverRight,
                 const std::shared_ptr< hyteg::Solver< AApproximationOperatorTypeLeft_ > >&  AApproximationSolverLeft,

                 hyteg::BoundaryCondition VelocityBCx,
                 hyteg::BoundaryCondition VelocityBCy,
                 hyteg::BoundaryCondition VelocityBCz,

                 bool             lowMemoryMode = false,
                 hyteg::BFBTState state         = hyteg::BFBTState::InnerBFBT )
   : hyteg::Operator< typename BTOperatorType_::srcType, typename BOperatorType_::dstType >( storage, minLevel, maxLevel )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , AOperator_( AOperator )
   , AApproximationOperatorRight_( AApproximationOperatorRight )
   , AApproximationOperatorLeft_( AApproximationOperatorLeft )
   , BTOperator_( BTOperator )
   , BOperator_( BOperator )
   , AApproximationSolverRight_( AApproximationSolverRight )
   , AApproximationSolverLeft_( AApproximationSolverLeft )
   , lowMemoryMode_( lowMemoryMode )
   , state_( state )
   , VelocityBCx_( VelocityBCx )
   , VelocityBCy_( VelocityBCy )
   , VelocityBCz_( VelocityBCz )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      if ( AOperator_ == nullptr )
      {
         WALBERLA_ABORT( "AOperator_ is nullptr!" );
      }
      if ( AApproximationOperatorRight_ == nullptr )
      {
         WALBERLA_ABORT( "AApproximationOperatorRight_ is nullptr!" );
      }
      if ( AApproximationOperatorLeft_ == nullptr )
      {
         WALBERLA_ABORT( "AApproximationOperatorLeft_ is nullptr!" );
      }

      if ( BTOperator == nullptr )
      {
         WALBERLA_ABORT( "BTOperator is nullptr!" );
      }

      if ( BOperator == nullptr )
      {
         WALBERLA_ABORT( "BOperator is nullptr!" );
      }

      if ( AApproximationSolverRight_ == nullptr )
      {
         WALBERLA_ABORT( "AApproximationSolverRight_ is nullptr!" );
      }
      if ( AApproximationSolverLeft_ == nullptr )
      {
         WALBERLA_ABORT( "AApproximationSolverLeft_ is nullptr!" );
      }

      // #############################
      // #### Boundary conditions ####
      // #############################

      if ( !lowMemoryMode_ )
      {
         temporary_ = std::make_shared< typename AApproximationOperatorTypeRight_::srcType >(
             "schur temporary", storage, minLevel, maxLevel );
         temporarySolution_ = std::make_shared< typename AApproximationOperatorTypeRight_::srcType >(
             "schur temporary solution", storage, minLevel, maxLevel );

         temporary_->component( 0 ).setBoundaryCondition( VelocityBCx_ );
         temporary_->component( 1 ).setBoundaryCondition( VelocityBCy_ );
         if ( storage_->hasGlobalCells() )
         {
            temporary_->component( 2 ).setBoundaryCondition( VelocityBCz_ );
         }

         temporarySolution_->component( 0 ).setBoundaryCondition( VelocityBCx_ );
         temporarySolution_->component( 1 ).setBoundaryCondition( VelocityBCy_ );
         if ( storage_->hasGlobalCells() )
         {
            temporarySolution_->component( 2 ).setBoundaryCondition( VelocityBCz_ );
         }
      }
   }

   void apply( const typename BTOperatorType_::srcType& src,
               const typename BOperatorType_::dstType&  dst,
               const uint_t                             level,
               const hyteg::DoFType                     flag,
               const hyteg::UpdateType                  updateType = hyteg::Replace ) const
   {
      std::shared_ptr< typename AApproximationOperatorTypeRight_::srcType > temporaryApply;
      std::shared_ptr< typename AApproximationOperatorTypeRight_::srcType > temporarySolutionApply;

      if ( !lowMemoryMode_ )
      {
         temporaryApply         = temporary_;
         temporarySolutionApply = temporarySolution_;
      }
      else
      {
         temporaryApply =
             hyteg::getTemporaryFunction< typename AApproximationOperatorTypeRight_::srcType >( storage_, minLevel_, maxLevel_ );
         temporarySolutionApply =
             hyteg::getTemporaryFunction< typename AApproximationOperatorTypeRight_::srcType >( storage_, minLevel_, maxLevel_ );

         temporaryApply->component( 0 ).setBoundaryCondition( VelocityBCx_ );
         temporaryApply->component( 1 ).setBoundaryCondition( VelocityBCy_ );
         if ( storage_->hasGlobalCells() )
         {
            temporaryApply->component( 2 ).setBoundaryCondition( VelocityBCz_ );
         }

         temporarySolutionApply->component( 0 ).setBoundaryCondition( VelocityBCx_ );
         temporarySolutionApply->component( 1 ).setBoundaryCondition( VelocityBCy_ );
         if ( storage_->hasGlobalCells() )
         {
            temporarySolutionApply->component( 2 ).setBoundaryCondition( VelocityBCz_ );
         }
      }

      switch ( state_ )
      {
      case hyteg::BFBTState::RightBFBT: {
         // apply BT Block
         BTOperator_->apply( src, *temporaryApply, level, hyteg::All, hyteg::Replace );

         // apply A^-1 Block
         temporarySolutionApply->setToZero( level );
         AApproximationSolverRight_->solve( *AApproximationOperatorRight_, *temporarySolutionApply, *temporaryApply, level );

         // apply B Block
         BOperator_->apply( *temporarySolutionApply, dst, level, flag, updateType );
      }
      break;
      case hyteg::BFBTState::InnerBFBT: {
         // apply BT Block
         BTOperator_->apply( src, *temporaryApply, level, hyteg::All, hyteg::Replace );

         // apply A^-1 Block
         temporarySolutionApply->setToZero( level );
         AApproximationSolverRight_->solve( *AApproximationOperatorRight_, *temporarySolutionApply, *temporaryApply, level );

         // apply A Block
         AOperator_->apply( *temporarySolutionApply,
                            *temporaryApply,
                            level,
                            hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
                            hyteg::Replace );
         temporaryApply->assign( { real_c( 1 ) }, { *temporarySolutionApply }, level, hyteg::DirichletBoundary );

         // apply A^-1 Block
         temporarySolutionApply->setToZero( level );
         AApproximationSolverLeft_->solve( *AApproximationOperatorLeft_, *temporarySolutionApply, *temporaryApply, level );

         // apply B Block
         BOperator_->apply( *temporarySolutionApply, dst, level, flag, updateType );
      }
      break;
      case hyteg::BFBTState::LeftBFBT: {
         // apply BT Block
         BTOperator_->apply( src, *temporaryApply, level, hyteg::All, hyteg::Replace );

         // apply A^-1 Block
         temporarySolutionApply->setToZero( level );
         AApproximationSolverLeft_->solve( *AApproximationOperatorLeft_, *temporarySolutionApply, *temporaryApply, level );

         // apply B Block
         BOperator_->apply( *temporarySolutionApply, dst, level, flag, updateType );
      }
      break;

      default:
         break;
      }
   }

   /// 0 = Inner, 1 = Right, 2 = Left
   void setState( hyteg::BFBTState state ) { state_ = state; }
   /// 0 = Inner, 1 = Right, 2 = Left
   hyteg::BFBTState getInner() { return state_; }

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< AOperatorType_ >                   AOperator_;
   std::shared_ptr< AApproximationOperatorTypeRight_ > AApproximationOperatorRight_;
   std::shared_ptr< AApproximationOperatorTypeLeft_ >  AApproximationOperatorLeft_;
   std::shared_ptr< BTOperatorType_ >                  BTOperator_;
   std::shared_ptr< BOperatorType_ >                   BOperator_;

   std::shared_ptr< hyteg::Solver< AApproximationOperatorTypeRight_ > > AApproximationSolverRight_;
   std::shared_ptr< hyteg::Solver< AApproximationOperatorTypeLeft_ > >  AApproximationSolverLeft_;

   std::shared_ptr< typename AApproximationOperatorTypeRight_::srcType > temporary_;
   std::shared_ptr< typename AApproximationOperatorTypeRight_::srcType > temporarySolution_;

   bool             lowMemoryMode_;
   hyteg::BFBTState state_;

   hyteg::BoundaryCondition VelocityBCx_;
   hyteg::BoundaryCondition VelocityBCy_;
   hyteg::BoundaryCondition VelocityBCz_;
};

} // namespace MantleConvection
