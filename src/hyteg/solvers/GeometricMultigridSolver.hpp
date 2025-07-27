/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Nils Kohl, Andreas Burkhart.
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

#include "core/DataTypes.h"
#include "core/timing/TimingTree.h"

#include "hyteg/functions/FunctionTools.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

template < class OperatorType >
class GeometricMultigridSolver : public Solver< OperatorType >
{
 public:
   using FunctionType = typename OperatorType::srcType;
   using ValueType    = typename FunctionTrait< FunctionType >::ValueType;

   /// \brief A generic geometric multigrid solver.
   ///
   /// \param storage                       A PrimitiveStorage instance.
   /// \param tmpFunction                   The geoemetric multigrid solver requires a temporary function. This can either be
   ///                                      constructed internally (using a different constructor), or can be passed explicitly
   ///                                      with this constructor.
   /// \param smoother                      A Solver instance that is employed as smoother.
   /// \param smoothIncrementOnCoarserGrids For each coarser level than the invoked max level, the number of pre- and post-smoothing
   ///                                      steps is increased by this value.
   /// \param constantRHS                   In most cases, this parameter can be ignored and set to false. If true, the rhs
   ///                                      function passed into the solve call does not need to be allocated on the finest level
   ///                                      IF the rhs is a constant function.
   ///                                      The sole purpose is to reduce the total allocated memory of the application, so no
   ///                                      performance advantage should be expected.
   /// \param constantRHSScalar             The constant RHS, if constantRHS is true.
   ///
   GeometricMultigridSolver( const std::shared_ptr< PrimitiveStorage >&              storage,
                             const std::shared_ptr< FunctionType >&                  tmpFunction,
                             std::shared_ptr< Solver< OperatorType > >               smoother,
                             std::shared_ptr< Solver< OperatorType > >               coarseSolver,
                             std::shared_ptr< RestrictionOperator< FunctionType > >  restrictionOperator,
                             std::shared_ptr< ProlongationOperator< FunctionType > > prolongationOperator,
                             uint_t                                                  minLevel,
                             uint_t                                                  maxLevel,
                             uint_t                                                  preSmoothSteps                = 3,
                             uint_t                                                  postSmoothSteps               = 3,
                             uint_t                                                  smoothIncrementOnCoarserGrids = 0,
                             CycleType                                               cycleType         = CycleType::VCYCLE,
                             bool                                                    lowMemoryMode     = false,
                             bool                                                    constantRHS       = false,
                             real_t                                                  constantRHSScalar = real_c( 0 ) )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , preSmoothSteps_( preSmoothSteps )
   , postSmoothSteps_( postSmoothSteps )   
   , smoothIncrement_( smoothIncrementOnCoarserGrids )   
   , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )   
   , cycleType_( cycleType )   
   , smoother_( smoother )
   , coarseSolver_( coarseSolver )
   , restrictionOperator_( restrictionOperator )
   , prolongationOperator_( prolongationOperator )
   , timingTree_( storage->getTimingTree() )
   , constantRHS_( constantRHS )
   , constantRHSScalar_( constantRHSScalar )
   , lowMemoryMode_( lowMemoryMode )
   {
      if ( !lowMemoryMode_ )
      {
         tmp_ = tmpFunction;
      }
   }


   //  private:
   // std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   // uint_t                                     minLevel_;
   // uint_t                                     maxLevel_;
   // uint_t                                     preSmoothSteps_;
   // uint_t                                     postSmoothSteps_;
   // uint_t                                     smoothIncrement_;
   // uint_t                                     invokedLevel_;

   // hyteg::DoFType flag_;
   // CycleType      cycleType_;

   // std::shared_ptr< hyteg::Solver< OperatorType > >               smoother_;
   // std::shared_ptr< hyteg::Solver< OperatorType > >               coarseSolver_;
   // std::shared_ptr< hyteg::RestrictionOperator< FunctionType > >  restrictionOperator_;
   // std::shared_ptr< hyteg::ProlongationOperator< FunctionType > > prolongationOperator_;

   // std::shared_ptr< FunctionType > tmp_;

   // std::shared_ptr< walberla::WcTimingTree > timingTree_;

   // bool   constantRHS_;
   // real_t constantRHSScalar_;

   // bool lowMemoryMode_;

   GeometricMultigridSolver( const std::shared_ptr< PrimitiveStorage >&              storage,
                             std::shared_ptr< Solver< OperatorType > >               smoother,
                             std::shared_ptr< Solver< OperatorType > >               coarseSolver,
                             std::shared_ptr< RestrictionOperator< FunctionType > >  restrictionOperator,
                             std::shared_ptr< ProlongationOperator< FunctionType > > prolongationOperator,
                             uint_t                                                  minLevel,
                             uint_t                                                  maxLevel,
                             uint_t                                                  preSmoothSteps                = 3,
                             uint_t                                                  postSmoothSteps               = 3,
                             uint_t                                                  smoothIncrementOnCoarserGrids = 0,
                             CycleType                                               cycleType         = CycleType::VCYCLE,
                             bool                                                    lowMemoryMode     = false,
                             bool                                                    constantRHS       = false,
                             real_t                                                  constantRHSScalar = real_c( 0 ) )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , preSmoothSteps_( preSmoothSteps )
   , postSmoothSteps_( postSmoothSteps )
   , smoothIncrement_( smoothIncrementOnCoarserGrids )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   , cycleType_( cycleType )
   , smoother_( smoother )
   , coarseSolver_( coarseSolver )
   , restrictionOperator_( restrictionOperator )
   , prolongationOperator_( prolongationOperator )
   , timingTree_( storage->getTimingTree() )
   , constantRHS_( constantRHS )
   , constantRHSScalar_( constantRHSScalar )
   , lowMemoryMode_( lowMemoryMode )
   {
      if ( !lowMemoryMode_ )
      {
         tmp_ = std::make_shared< FunctionType >( "gmg_tmp", storage, minLevel, maxLevel );
      }
   }

   ~GeometricMultigridSolver() = default;

   void setSmoothingSteps( const uint_t& preSmoothingSteps, const uint_t& postSmoothingSteps, const uint_t smoothIncrement = 0 )
   {
      preSmoothSteps_  = preSmoothingSteps;
      postSmoothSteps_ = postSmoothingSteps;
      smoothIncrement_ = smoothIncrement;
   }

   /// \brief applies the generic geometric multigrid solver to the LSE Ax=b.
   ///     Even though, the solution is computed for the finest level 'level', all parameters must be allocated between the levels {'minLevel_', ..., 'maxLevel'}.
   ///
   /// \param A        Operator matrix A.
   /// \param x        Solution vector x.
   /// \param b        Right-hand-side vector b. After execution of solve, only the finest level 'level' still contains the passed values, all courser values are overwritten with intermediate values.
   /// \param level    Finest level for which 'x' is approximated.
   ///
   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      std::shared_ptr< FunctionType > tmpSolve;

      if ( !lowMemoryMode_ )
      {
         tmpSolve = tmp_;
      }
      else
      {
         tmpSolve = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );

         // // Since the temporary function only serves as an apply target or RHS
         // // it is unnecessary to zero the tmpSolve first.
         // for ( uint_t l = minLevel_; l <= maxLevel_; l++ )
         // {
         //    tmpSolve->setToZero( l );
         // }
      }

      timingTree_->start( "Geometric Multigrid Solver" );
      copyBCs( x, *tmpSolve );
      invokedLevel_ = level;
      solveRecursively( A, x, b, level, tmpSolve );
      timingTree_->stop( "Geometric Multigrid Solver" );
   }

   void solveRecursively( const OperatorType&                    A,
                          const FunctionType&                    x,
                          const FunctionType&                    b,
                          const uint_t                           level,
                          const std::shared_ptr< FunctionType >& tmpSolve ) const
   {
      if ( level == minLevel_ )
      {
         timingTree_->start( "Level " + std::to_string( level ) );
         timingTree_->start( "Coarse Grid Solver" );
         coarseSolver_->solve( A, x, b, minLevel_ );
         timingTree_->stop( "Coarse Grid Solver" );
         timingTree_->stop( "Level " + std::to_string( level ) );
      }
      else
      {
         timingTree_->start( "Level " + std::to_string( level ) );

         if ( constantRHS_ && level == invokedLevel_ )
         {
            tmpSolve->interpolate( constantRHSScalar_, level, flag_ );
         }

         // pre-smooth
         const uint_t preSmoothingSteps = preSmoothSteps_ + smoothIncrement_ * ( invokedLevel_ - level );
         for ( uint_t i = 0; i < preSmoothingSteps; ++i )
         {
            timingTree_->start( "Smoother" );
            if ( constantRHS_ && level == invokedLevel_ )
            {
               smoother_->solve( A, x, *tmpSolve, level );
            }
            else
            {
               smoother_->solve( A, x, b, level );
            }

            timingTree_->stop( "Smoother" );
         }

         if ( constantRHS_ && level == invokedLevel_ )
         {
            A.apply( x, *tmpSolve, level, flag_ );
            tmpSolve->assign( { walberla::numeric_cast< ValueType >( -1.0 ) }, { *tmpSolve }, level, flag_ );
            tmpSolve->add( constantRHSScalar_, level, flag_ );
         }
         else
         {
            A.apply( x, *tmpSolve, level, flag_ );
            tmpSolve->assign( { walberla::numeric_cast< ValueType >( 1.0 ), walberla::numeric_cast< ValueType >( -1.0 ) },
                              { b, *tmpSolve },
                              level,
                              flag_ );
         }

         // restrict
         timingTree_->start( "Restriction" );
         restrictionOperator_->restrict( *tmpSolve, level, flag_ );
         timingTree_->stop( "Restriction" );

         b.assign( { walberla::numeric_cast< ValueType >( 1.0 ) }, { *tmpSolve }, level - 1, flag_ );

         x.interpolate( 0, level - 1 );

         timingTree_->stop( "Level " + std::to_string( level ) );
         solveRecursively( A, x, b, level - 1, tmpSolve );

         if ( cycleType_ == CycleType::WCYCLE )
         {
            solveRecursively( A, x, b, level - 1, tmpSolve );
         }
         timingTree_->start( "Level " + std::to_string( level ) );

         // prolongate
         timingTree_->start( "Prolongation" );
         prolongationOperator_->prolongateAndAdd( x, level - 1, flag_ );
         timingTree_->stop( "Prolongation" );

         if ( constantRHS_ && level == invokedLevel_ )
         {
            tmpSolve->interpolate( constantRHSScalar_, level, flag_ );
         }

         // post-smooth
         const uint_t postSmoothingSteps = postSmoothSteps_ + smoothIncrement_ * ( invokedLevel_ - level );
         for ( uint_t i = 0; i < postSmoothingSteps; ++i )
         {
            timingTree_->start( "Smoother" );
            if ( constantRHS_ && level == invokedLevel_ )
            {
               smoother_->solve( A, x, *tmpSolve, level );
            }
            else
            {
               smoother_->solve( A, x, b, level );
            }
            timingTree_->stop( "Smoother" );
         }

         timingTree_->stop( "Level " + std::to_string( level ) );
      }
   }

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;
   uint_t                                     preSmoothSteps_;
   uint_t                                     postSmoothSteps_;
   uint_t                                     smoothIncrement_;
   uint_t                                     invokedLevel_;

   hyteg::DoFType flag_;
   CycleType      cycleType_;

   std::shared_ptr< hyteg::Solver< OperatorType > >               smoother_;
   std::shared_ptr< hyteg::Solver< OperatorType > >               coarseSolver_;
   std::shared_ptr< hyteg::RestrictionOperator< FunctionType > >  restrictionOperator_;
   std::shared_ptr< hyteg::ProlongationOperator< FunctionType > > prolongationOperator_;

   std::shared_ptr< FunctionType > tmp_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

   bool   constantRHS_;
   real_t constantRHSScalar_;

   bool lowMemoryMode_;
};

} // namespace hyteg
