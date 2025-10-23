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

template < typename OperatorType >
class ApproximateSchurComplementOperator
: public hyteg::Operator< typename OperatorType::srcType::PressureFunction_T,
                          typename OperatorType::srcType::PressureFunction_T >,
  public hyteg::OperatorWithInverseDiagonal< typename OperatorType::srcType::PressureFunction_T >
{
 public:
   using hyteg::Operator< typename OperatorType::srcType::PressureFunction_T,
                          typename OperatorType::srcType::PressureFunction_T >::storage_;
   using hyteg::Operator< typename OperatorType::srcType::PressureFunction_T,
                          typename OperatorType::srcType::PressureFunction_T >::minLevel_;
   using hyteg::Operator< typename OperatorType::srcType::PressureFunction_T,
                          typename OperatorType::srcType::PressureFunction_T >::maxLevel_;
   typedef typename OperatorType::AOperatorType                  AOperatorType;
   typedef typename OperatorType::VelocityProjectionOperatorType VelocityProjectionOperatorType_;

   ApproximateSchurComplementOperator( const std::shared_ptr< hyteg::PrimitiveStorage >&        storage,
                                       uint_t                                                   minLevel,
                                       uint_t                                                   maxLevel,
                                       const std::shared_ptr< OperatorType >&                   SaddleOp,
                                       const std::shared_ptr< hyteg::Solver< AOperatorType > >& ABlockSolver,
                                       hyteg::BoundaryCondition                                 VelocityBC,
                                       bool                                                     lowMemoryMode = false )
   : ApproximateSchurComplementOperator( storage,
                                         minLevel,
                                         maxLevel,
                                         SaddleOp,
                                         ABlockSolver,
                                         VelocityBC,
                                         VelocityBC,
                                         VelocityBC,
                                         lowMemoryMode )
   {}

   ApproximateSchurComplementOperator( const std::shared_ptr< hyteg::PrimitiveStorage >&        storage,
                                       uint_t                                                   minLevel,
                                       uint_t                                                   maxLevel,
                                       const std::shared_ptr< OperatorType >&                   SaddleOp,
                                       const std::shared_ptr< hyteg::Solver< AOperatorType > >& ABlockSolver,
                                       hyteg::BoundaryCondition                                 VelocityBCx,
                                       hyteg::BoundaryCondition                                 VelocityBCy,
                                       hyteg::BoundaryCondition                                 VelocityBCz,
                                       bool                                                     lowMemoryMode = false )
   : hyteg::Operator< typename OperatorType::srcType::PressureFunction_T, typename OperatorType::srcType::PressureFunction_T >(
         storage,
         minLevel,
         maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , SaddleOp_( SaddleOp )
   , ABlockSolver_( ABlockSolver )
   , lowMemoryMode_( lowMemoryMode )
   , VelocityBCx_( VelocityBCx )
   , VelocityBCy_( VelocityBCy )
   , VelocityBCz_( VelocityBCz )
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

      if constexpr ( BTInvABBlock )
      {
         if ( ( SaddleOp->getBTPtr() == nullptr ) || ( SaddleOp->getAPtr() == nullptr ) || ( SaddleOp->getBPtr() == nullptr ) )
         {
            WALBERLA_ABORT( "BT, A and B must not be nullptr if all three operators are unequal to NoOperator in type!" );
         }
      }

      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( SaddleOp->getStabPtr() == nullptr )
         {
            WALBERLA_ABORT( "Stabilisation set to nullptr but type is not NoOperator!" );
         }
      }

      // #############################
      // #### Boundary conditions ####
      // #############################

      if ( !lowMemoryMode_ )
      {
         temporary_ = std::make_shared< typename OperatorType::srcType::VelocityFunction_T >(
             "schur temporary", storage, minLevel, maxLevel );
         temporarySolution_ = std::make_shared< typename OperatorType::srcType::VelocityFunction_T >(
             "schur temporary solution", storage, minLevel, maxLevel );

         temporary_->component( 0 ).setBoundaryCondition( VelocityBCx_ );
         temporary_->component( 1 ).setBoundaryCondition( VelocityBCy_ );
         if ( hasGlobalCells_ )
         {
            temporary_->component( 2 ).setBoundaryCondition( VelocityBCz_ );
         }

         temporarySolution_->component( 0 ).setBoundaryCondition( VelocityBCx_ );
         temporarySolution_->component( 1 ).setBoundaryCondition( VelocityBCy_ );
         if ( hasGlobalCells_ )
         {
            temporarySolution_->component( 2 ).setBoundaryCondition( VelocityBCz_ );
         }
      }
   }

   void apply( const typename OperatorType::srcType::PressureFunction_T& src,
               const typename OperatorType::srcType::PressureFunction_T& dst,
               const uint_t                                              level,
               const hyteg::DoFType                                      flag,
               const hyteg::UpdateType                                   updateType = hyteg::Replace ) const
   {
      std::shared_ptr< typename OperatorType::srcType::VelocityFunction_T > temporaryApply;
      std::shared_ptr< typename OperatorType::srcType::VelocityFunction_T > temporarySolutionApply;

      if ( !lowMemoryMode_ )
      {
         temporaryApply         = temporary_;
         temporarySolutionApply = temporarySolution_;
      }
      else
      {
         temporaryApply =
             hyteg::getTemporaryFunction< typename OperatorType::srcType::VelocityFunction_T >( storage_, minLevel_, maxLevel_ );
         temporarySolutionApply =
             hyteg::getTemporaryFunction< typename OperatorType::srcType::VelocityFunction_T >( storage_, minLevel_, maxLevel_ );

         temporaryApply->component( 0 ).setBoundaryCondition( VelocityBCx_ );
         temporaryApply->component( 1 ).setBoundaryCondition( VelocityBCy_ );
         if ( hasGlobalCells_ )
         {
            temporaryApply->component( 2 ).setBoundaryCondition( VelocityBCz_ );
         }

         temporarySolutionApply->component( 0 ).setBoundaryCondition( VelocityBCx_ );
         temporarySolutionApply->component( 1 ).setBoundaryCondition( VelocityBCy_ );
         if ( hasGlobalCells_ )
         {
            temporarySolutionApply->component( 2 ).setBoundaryCondition( VelocityBCz_ );
         }
      }

      if constexpr ( BTInvABBlock )
      {
         // apply BT Block
         SaddleOp_->getBT().apply( src, *temporaryApply, level, hyteg::All, hyteg::Replace );

         // apply A^-1 Block
         //temporaryApply->interpolate( real_c( 0.0 ), level, hyteg::DirichletBoundary );
         temporarySolutionApply->setToZero( level );
         ABlockSolver_->solve( SaddleOp_->getA(), *temporarySolutionApply, *temporaryApply, level );
         //temporarySolutionApply->interpolate( real_c( 0.0 ), level, hyteg::DirichletBoundary );

         // apply B Block
         SaddleOp_->getB().apply( *temporarySolutionApply, dst, level, flag, updateType );
      }

      // apply Stabilisation Block
      if constexpr ( !std::is_same< typename OperatorType::StabilisationOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         SaddleOp_->getStab().apply( src, dst, level, flag, ( BTInvABBlock ? hyteg::UpdateType::Add : updateType ) );
      }
   }

   std::shared_ptr< typename OperatorType::srcType::PressureFunction_T > getInverseDiagonalValues() const { return invDiag_; }

   void computeInverseDiagonalOperatorValues()
   {
      if ( invDiag_ == nullptr )
      {
         invDiag_ = std::make_shared< typename OperatorType::srcType::PressureFunction_T >(
             "schur temporary", storage_, minLevel_, maxLevel_ );
      }
   }
   void setInverseDiagonal( typename OperatorType::srcType::PressureFunction_T diag )
   {
      computeInverseDiagonalOperatorValues();
      for ( uint_t l = minLevel_; l <= maxLevel_; l++ )
      {
         invDiag_->assign( { real_c( 1.0 ) }, { diag }, l, hyteg::All );
      }
   }
   void setInverseDiagonal( std::shared_ptr< typename OperatorType::srcType::PressureFunction_T > diagPtr )
   {
      invDiag_ = diagPtr;
   }

   static constexpr bool BTInvABBlock =
       ( !std::is_same< typename OperatorType::BTOperatorTypeInternal, hyteg::NoOperator >::value ) &&
       ( !std::is_same< typename OperatorType::AOperatorTypeInternal, hyteg::NoOperator >::value ) &&
       ( !std::is_same< typename OperatorType::BOperatorTypeInternal, hyteg::NoOperator >::value );

 private:
   bool                                                                  hasGlobalCells_;
   std::shared_ptr< typename OperatorType::srcType::VelocityFunction_T > temporary_;
   std::shared_ptr< typename OperatorType::srcType::VelocityFunction_T > temporarySolution_;
   std::shared_ptr< OperatorType >                                       SaddleOp_;
   std::shared_ptr< hyteg::Solver< AOperatorType > >                     ABlockSolver_;

   std::shared_ptr< typename OperatorType::srcType::PressureFunction_T > invDiag_;

   bool lowMemoryMode_;

   hyteg::BoundaryCondition VelocityBCx_;
   hyteg::BoundaryCondition VelocityBCy_;
   hyteg::BoundaryCondition VelocityBCz_;
};

} // namespace MantleConvection
