/*
 * Copyright (c) 2024 Dinesh Parthasarathy.
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
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::int8_t;
using walberla::real_t;
using walberla::uint_t;

template < class OperatorType >
class FlexibleMultigridSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   FlexibleMultigridSolver( const std::shared_ptr< PrimitiveStorage >&               storage,
                            std::vector< std::shared_ptr< Solver< OperatorType > > > smootherList,
                            std::shared_ptr< Solver< OperatorType > >                coarseSolver,
                            std::shared_ptr< RestrictionOperator< FunctionType > >   restrictionOperator,
                            std::shared_ptr< ProlongationOperator< FunctionType > >  prolongationOperator,
                            uint_t                                                   minLevel,
                            uint_t                                                   maxLevel,
                            std::vector< int8_t >                                    cyclingStructure,
                            bool                                                     constantRHS       = false,
                            real_t                                                   constantRHSScalar = real_c( 0 ) )
   : FlexibleMultigridSolver( storage,
                              FunctionType( "flexmg_tmp", storage, minLevel, maxLevel ),
                              smootherList,
                              coarseSolver,
                              restrictionOperator,
                              prolongationOperator,
                              minLevel,
                              maxLevel,
                              cyclingStructure,
                              constantRHS,
                              constantRHSScalar )
   {}

   /// \brief A flexible multigrid solver i.e. the cycle structure can be chosen freely (outside of regular recursive  V-, W-, F-, kappa-cycles)
   /// as a graph with nodal (smoothing) and internodal (restriction, prolongation) operations.
   ///
   /// \param storage                       A PrimitiveStorage instance.
   /// \param tmpFunction                   The multigrid solver requires a temporary function. This can either be
   ///                                      constructed internally (using a different constructor), or can be passed explicitly
   ///                                      with this constructor.
   /// \param smootherList                  List of Solver instances that is employed as smoother in every node of the MG cycle graph.
   /// \param cyclingStructure              List of edge information in the MG cycle graph. -1 implies restriction, 1 implies prolongation and 0 implies staying on the same level.
   /// \param constantRHS                   In most cases, this parameter can be ignored and set to false. If true, the rhs
   ///                                      function passed into the solve call does not need to be allocated on the finest level
   ///                                      IF the rhs is a constant function.
   ///                                      The sole purpose is to reduce the total allocated memory of the application, so no
   ///                                      performance advantage should be expected.
   /// \param constantRHSScalar             The constant RHS, if constantRHS is true.
   ///
   FlexibleMultigridSolver( const std::shared_ptr< PrimitiveStorage >&               storage,
                            const FunctionType&                                      tmpFunction,
                            std::vector< std::shared_ptr< Solver< OperatorType > > > smootherList,
                            std::shared_ptr< Solver< OperatorType > >                coarseSolver,
                            std::shared_ptr< RestrictionOperator< FunctionType > >   restrictionOperator,
                            std::shared_ptr< ProlongationOperator< FunctionType > >  prolongationOperator,
                            uint_t                                                   minLevel,
                            uint_t                                                   maxLevel,
                            std::vector< int8_t >                                    cyclingStructure,
                            bool                                                     constantRHS       = false,
                            real_t                                                   constantRHSScalar = real_c( 0 ) )
   : minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , smootherList_( smootherList )
   , cyclingStructure_( cyclingStructure )
   , coarseSolver_( coarseSolver )
   , restrictionOperator_( restrictionOperator )
   , prolongationOperator_( prolongationOperator )
   , tmp_( tmpFunction )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary )
   , timingTree_( storage->getTimingTree() )
   , constantRHS_( constantRHS )
   , constantRHSScalar_( constantRHSScalar )
   {}

   ~FlexibleMultigridSolver() = default;

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      uint_t invokedLevel_;
      uint_t level_;

      timingTree_->start( "Flexible Multigrid Solver" );
      copyBCs( x, tmp_ );
      invokedLevel_ = level;
      level_        = maxLevel_;

      // iterate  over the nodes of the mg cycle from maxLevel_ to minLevel_,
      for ( uint_t i = 0; i < smootherList_.size(); ++i )
      {
         if ( level_ <= invokedLevel_ ) // solve the system on the level it is invoked
         {
            if ( level_ == minLevel_ )
            {
               timingTree_->start( "Level " + std::to_string( level_ ) );
               timingTree_->start( "Coarse Grid Solver" );
               coarseSolver_->solve( A, x, b, minLevel_ );
               timingTree_->stop( "Coarse Grid Solver" );
               timingTree_->stop( "Level " + std::to_string( level_ ) );

               level_++;
               timingTree_->start( "Prolongation" );
               prolongationOperator_->prolongateAndAdd( x, level_ - 1, flag_ );
               timingTree_->stop( "Prolongation" );
            }
            else
            {
               if ( i == 0 || cyclingStructure_[i - 1] != 0 )
               {
                  timingTree_->start( "Level " + std::to_string( level_ ) );
               }

               if ( constantRHS_ && level == invokedLevel_ )
               {
                  tmp_.interpolate( constantRHSScalar_, level, flag_ );
               }

               // I. Smoothing
               timingTree_->start( "Smoother" );
               if ( constantRHS_ && level_ == invokedLevel_ )
               {
                  smootherList_[i]->solve( A, x, tmp_, level_ );
               }
               else
               {
                  smootherList_[i]->solve( A, x, b, level_ );
               }
               timingTree_->stop( "Smoother" );

               // II. Internodal operations (Restriction, Prolongation)
               if ( i < cyclingStructure_.size() )
               {
                  if ( cyclingStructure_[i] == 1 ) // prolongation
                  {
                     timingTree_->stop( "Level " + std::to_string( level_ ) );
                     level_++;
                     timingTree_->start( "Prolongation" );
                     prolongationOperator_->prolongateAndAdd( x, level_ - 1, flag_ );
                     timingTree_->stop( "Prolongation" );
                  }
                  else if ( cyclingStructure_[i] == -1 ) // restriction
                  {
                     A.apply( x, tmp_, level_, flag_ );
                     if ( constantRHS_ && level_ == invokedLevel_ )
                     {
                        tmp_.assign( { -1.0 }, { tmp_ }, level_, flag_ );
                        tmp_.add( constantRHSScalar_, level_, flag_ );
                     }
                     else
                     {
                        tmp_.assign( { 1.0, -1.0 }, { b, tmp_ }, level_, flag_ );
                     }
                     // restrict
                     timingTree_->start( "Restriction" );
                     restrictionOperator_->restrict( tmp_, level_, flag_ );
                     timingTree_->stop( "Restriction" );

                     b.assign( { 1.0 }, { tmp_ }, level_ - 1, flag_ );
                     x.interpolate( 0, level_ - 1 );
                     timingTree_->stop( "Level " + std::to_string( level_ ) );
                     level_--;
                  }
               }
               else
               {
                  timingTree_->stop( "Level " + std::to_string( level_ ) );
               }
            }
         }
         else if ( i < cyclingStructure_.size() )
         {
            level_ += cyclingStructure_[i];
         }
      }
      timingTree_->stop( "Flexible Multigrid Solver" );
   }

 private:
   uint_t minLevel_;
   uint_t maxLevel_;

   hyteg::DoFType flag_;

   std::vector< std::shared_ptr< hyteg::Solver< OperatorType > > > smootherList_;
   std::vector< int8_t >                                           cyclingStructure_;
   std::shared_ptr< hyteg::Solver< OperatorType > >                coarseSolver_;
   std::shared_ptr< hyteg::RestrictionOperator< FunctionType > >   restrictionOperator_;
   std::shared_ptr< hyteg::ProlongationOperator< FunctionType > >  prolongationOperator_;

   FunctionType tmp_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

   bool   constantRHS_;
   real_t constantRHSScalar_;
};

} // namespace hyteg
