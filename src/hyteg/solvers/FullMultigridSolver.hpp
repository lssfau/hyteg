/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Nils Kohl, Daniel Bauer.
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

#include <vector>

#include "core/Abort.h"
#include "core/DataTypes.h"

#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/Solver.hpp"

namespace hyteg {

template < class OperatorType,
           class MultigridSolverType = GeometricMultigridSolver<
               OperatorType > > // MultigridSolverType: chose GeometricMultigridSolver or FlexibleMultigridSolver
class FullMultigridSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   FullMultigridSolver(
       const std::shared_ptr< PrimitiveStorage >&                     storage,
       const std::shared_ptr< MultigridSolverType >&                  gmgSolver,
       const std::shared_ptr< ProlongationOperator< FunctionType > >& fmgProlongation,
       const uint_t&                                                  minLevel,
       const uint_t&                                                  maxLevel,
       const uint_t&                                                  cyclesPerLevel         = 1,
       const std::function< void( uint_t currentLevel ) >&            postCycleCallback      = []( uint_t ) {},
       const std::function< void( uint_t currentLevel ) >&            postProlongateCallback = []( uint_t ) {} )
   : FullMultigridSolver( storage,
                          gmgSolver,
                          fmgProlongation,
                          minLevel,
                          maxLevel,
                          std::vector< uint_t >( maxLevel + 1, cyclesPerLevel ),
                          postCycleCallback )
   {}

   FullMultigridSolver(
       const std::shared_ptr< PrimitiveStorage >&                     storage,
       const std::shared_ptr< MultigridSolverType >&                  gmgSolver,
       const std::shared_ptr< ProlongationOperator< FunctionType > >& fmgProlongation,
       const uint_t&                                                  minLevel,
       const uint_t&                                                  maxLevel,
       const std::vector< uint_t >&                                   cyclesPerLevel,
       const std::function< void( uint_t currentLevel ) >&            postCycleCallback      = []( uint_t ) {},
       const std::function< void( uint_t currentLevel ) >&            postProlongateCallback = []( uint_t ) {} )
   : gmgSolver_( gmgSolver )
   , fmgProlongation_( fmgProlongation )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , cyclesPerLevel_( cyclesPerLevel )
   , flag_( Inner | NeumannBoundary )
   , postCycleCallback_( postCycleCallback )
   , timingTree_( storage->getTimingTree() )
   {
      if ( cyclesPerLevel.size() != maxLevel + 1 )
      {
         WALBERLA_ABORT(
             "Specify the number of V-Cycles for every level in {0, ..., maxLevel}. Values below minLevel must be given but are ignored." )
      }
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      timingTree_->start( "FMG Solver" );
      for ( uint_t currentLevel = minLevel_; currentLevel <= level; currentLevel++ )
      {
         timingTree_->start( "GMG Solver" );
         for ( uint_t cycle = 0; cycle < cyclesPerLevel_[currentLevel]; cycle++ )
         {
            gmgSolver_->solve( A, x, b, currentLevel );
         }
         timingTree_->stop( "GMG Solver" );

         timingTree_->start( "Post-cycle callback" );
         postCycleCallback_( currentLevel );
         timingTree_->stop( "Post-cycle callback" );

         timingTree_->start( "FMG Prolongation" );
         if ( currentLevel < maxLevel_ )
         {
            fmgProlongation_->prolongate( x, currentLevel, flag_ );
         }
         timingTree_->stop( "FMG Prolongation" );

         timingTree_->start( "Post-prolongate callback" );
         postCycleCallback_( currentLevel );
         timingTree_->stop( "Post-prolongate callback" );
      }
      timingTree_->stop( "FMG Solver" );
   }

 private:
   std::shared_ptr< MultigridSolverType >                  gmgSolver_;
   std::shared_ptr< ProlongationOperator< FunctionType > > fmgProlongation_;

   uint_t minLevel_;
   uint_t maxLevel_;

   std::vector< uint_t > cyclesPerLevel_;

   DoFType flag_;

   std::function< void( uint_t currentLevel ) > postCycleCallback_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;
};

} // namespace hyteg