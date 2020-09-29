/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/types/pointnd.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

template < class OperatorType >
class GeometricMultigridSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   //  static_assert( !std::is_same< FunctionType, typename OperatorType::dstType >::value,
   //                 "CGSolver does not work for Operator with different src and dst FunctionTypes" );

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
                             CycleType                                               cycleType = CycleType::VCYCLE )
   : GeometricMultigridSolver( storage,
                               FunctionType( "gmg_tmp", storage, minLevel, maxLevel ),
                               smoother,
                               coarseSolver,
                               restrictionOperator,
                               prolongationOperator,
                               minLevel,
                               maxLevel,
                               preSmoothSteps,
                               postSmoothSteps,
                               smoothIncrementOnCoarserGrids,
                               cycleType )
   {}

   GeometricMultigridSolver( const std::shared_ptr< PrimitiveStorage >&              storage,
                             const FunctionType&                                     tmpFunction,
                             std::shared_ptr< Solver< OperatorType > >               smoother,
                             std::shared_ptr< Solver< OperatorType > >               coarseSolver,
                             std::shared_ptr< RestrictionOperator< FunctionType > >  restrictionOperator,
                             std::shared_ptr< ProlongationOperator< FunctionType > > prolongationOperator,
                             uint_t                                                  minLevel,
                             uint_t                                                  maxLevel,
                             uint_t                                                  preSmoothSteps                = 3,
                             uint_t                                                  postSmoothSteps               = 3,
                             uint_t                                                  smoothIncrementOnCoarserGrids = 0,
                             CycleType                                               cycleType = CycleType::VCYCLE )
   : minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , smoother_( smoother )
   , coarseSolver_( coarseSolver )
   , restrictionOperator_( restrictionOperator )
   , prolongationOperator_( prolongationOperator )
   , tmp_( tmpFunction )
   , preSmoothSteps_( preSmoothSteps )
   , postSmoothSteps_( postSmoothSteps )
   , smoothIncrement_( smoothIncrementOnCoarserGrids )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary )
   , cycleType_( cycleType )
   , timingTree_( storage->getTimingTree() )
   {}

   ~GeometricMultigridSolver() = default;

   void setSmoothingSteps( const uint_t& preSmoothingSteps, const uint_t& postSmoothingSteps, const uint_t smoothIncrement = 0 )
   {
      preSmoothSteps_  = preSmoothingSteps;
      postSmoothSteps_ = postSmoothingSteps;
      smoothIncrement_ = smoothIncrement;
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      timingTree_->start( "Geometric Multigrid Solver" );
      tmp_.copyBoundaryConditionFromFunction( x );
      invokedLevel_ = level;
      solveRecursively( A, x, b, level );
      timingTree_->stop( "Geometric Multigrid Solver" );
   }

   void solveRecursively( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) const
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

         // pre-smooth
         const uint_t preSmoothingSteps = preSmoothSteps_ + smoothIncrement_ * ( invokedLevel_ - level );
         for ( uint_t i = 0; i < preSmoothingSteps; ++i )
         {
            timingTree_->start( "Smoother" );
            smoother_->solve( A, x, b, level );
            timingTree_->stop( "Smoother" );
         }

         A.apply( x, tmp_, level, flag_ );
         tmp_.assign( {1.0, -1.0}, {b, tmp_}, level, flag_ );

         // restrict
         timingTree_->start( "Restriction" );
         restrictionOperator_->restrict( tmp_, level, flag_ );
         timingTree_->stop( "Restriction" );

         b.assign( {1.0}, {tmp_}, level - 1, flag_ );

         x.interpolate( 0, level - 1 );

         timingTree_->stop( "Level " + std::to_string( level ) );
         solveRecursively( A, x, b, level - 1 );

         if ( cycleType_ == CycleType::WCYCLE )
         {
            solveRecursively( A, x, b, level - 1 );
         }
         timingTree_->start( "Level " + std::to_string( level ) );

         // prolongate
         tmp_.assign( {1.0}, {x}, level, flag_ );
         timingTree_->start( "Prolongation" );
         prolongationOperator_->prolongate( x, level - 1, flag_ );
         timingTree_->stop( "Prolongation" );
         x.add( {1.0}, {tmp_}, level, flag_ );

         // post-smooth
         const uint_t postSmoothingSteps = postSmoothSteps_ + smoothIncrement_ * ( invokedLevel_ - level );
         for ( uint_t i = 0; i < postSmoothingSteps; ++i )
         {
            timingTree_->start( "Smoother" );
            smoother_->solve( A, x, b, level );
            timingTree_->stop( "Smoother" );
         }

         timingTree_->stop( "Level " + std::to_string( level ) );
      }
   }

 private:
   uint_t minLevel_;
   uint_t maxLevel_;
   uint_t preSmoothSteps_;
   uint_t postSmoothSteps_;
   uint_t smoothIncrement_;
   uint_t invokedLevel_;

   hyteg::DoFType flag_;
   CycleType    cycleType_;

   std::shared_ptr< hyteg::Solver< OperatorType > >               smoother_;
   std::shared_ptr< hyteg::Solver< OperatorType > >               coarseSolver_;
   std::shared_ptr< hyteg::RestrictionOperator< FunctionType > >  restrictionOperator_;
   std::shared_ptr< hyteg::ProlongationOperator< FunctionType > > prolongationOperator_;

   FunctionType tmp_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;
};

} // namespace hyteg
