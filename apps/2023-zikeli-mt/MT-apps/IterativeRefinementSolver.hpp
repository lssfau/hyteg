/*
 * Copyright (c) 2024 Michael Zikeli
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
#include "core/debug/Debug.h"
#include "core/timing/TimingTree.h"

#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/types/PointND.hpp"

#ifndef PERFORMANCE_RUN
#include "Tree.h"
#endif

namespace hyteg {

using walberla::uint_t;

template < class SmootherOperator_t, class ResidualOperator_t = SmootherOperator_t >
class IterativeRefinementSolver : public Solver< ResidualOperator_t >
{
 public:
   using BarFunctionType = typename ResidualOperator_t::srcType;
   using DotFunctionType = typename SmootherOperator_t::srcType;
   using BarValueType    = typename FunctionTrait< BarFunctionType >::ValueType;
   using DotValueType    = typename FunctionTrait< DotFunctionType >::ValueType;

   /// \brief A generic iterative refinement solver.
   ///
   ///     Note: This solver allocates 4 temporary functions for the levels {minLevel, ..., maxLevel}
   ///     TODO add paper describing IR.
   ///
   /// \param storage                       A PrimitiveStorage instance.
   /// \param smoother                      A Solver instance that is employed as smoother.
   ///
   IterativeRefinementSolver( const std::shared_ptr< PrimitiveStorage >&      storage,
                              const SmootherOperator_t&                       dotA,
                              std::shared_ptr< Solver< SmootherOperator_t > > smoother,
                              uint_t                                          minLevel,
                              uint_t maxLevel, // FIXME this is actually not necessary since it is never used.
                              uint_t smoothingSteps )
   : minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , correctionSteps_( smoothingSteps )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary ) // FIXME I have no idea, I just copied this from GMG
   , smoother_( smoother )
   , dotA_( dotA )
   , storage_( storage )
   , dotCorrection_( "DotCorrection", storage, minLevel, maxLevel )
   , dotResidual_( "DotResidual", storage, minLevel, maxLevel )
   , barCorrection_( "BarCorrection", storage, minLevel, maxLevel ) // FIXME only necessary for fine level
   , barResidual_( "BarResidual", storage, minLevel, maxLevel )     // FIXME only necessary for fine level
   , residualIsPrecalculated_( false )
   , levelPrecalculatedResidual_()
   , timingTree_( storage->getTimingTree() )
#ifndef PERFORMANCE_RUN
   , residualTree_( nullptr )
   , residualNode_( nullptr )
#endif
   {}

   ~IterativeRefinementSolver() = default;

   /// \brief Helper function to reduce the runtime, if the residual in Bar-Precision of Ax=b is already precomputed.
   ///
   ///     If this function is called before solve, the given residual will be used as RHS of the residual function,
   ///         i.e. the RHS function for the inner solve.
   ///
   /// \param residual      The residual in Bar-Precision.
   /// \param level         The level for which the residual was precomputed.
   void setPrecalculatedResidual( const BarFunctionType& residual, const uint_t level )
   {
      levelPrecalculatedResidual_ = level;
      barResidual_.assign( { walberla::numeric_cast< BarValueType >( 1 ) }, { residual }, level );
      residualIsPrecalculated_ = true;
   }

#ifndef PERFORMANCE_RUN
   void setResidualNode( TreeNode< BarValueType >* node, Tree< BarValueType >* tree )
   {
      residualNode_ = node;
      residualTree_ = tree;
   }
#endif

   void solve( const ResidualOperator_t& barA, const BarFunctionType& x, const BarFunctionType& b, const uint_t level ) override
   {
      timingTree_->start( "Iterative Refinement Solver" );

      // TODO allocate variables
      //    res = calculation // @Nils the rhs is overwritten before restriction, so it doesn't matter what the RHS is.
      //    c = 0             // @Nils same argument

      if ( residualIsPrecalculated_ )
      {
         // Check if the precomputed residual was given for the right level
         WALBERLA_CHECK_EQUAL( level, levelPrecalculatedResidual_ );
         // Reset the precomputed flag, since the precomputed residual is not correct anymore, once the approximation was updated.
         residualIsPrecalculated_ = false;
      }
      else
      {
         timingTree_->start( "Inner Residual Computation" );
         // r = Ax - b
         barA.apply( x, barResidual_, level, All, Replace );
         barResidual_.assign( { walberla::numeric_cast< BarValueType >( 1 ), walberla::numeric_cast< BarValueType >( -1 ) },
                              { barResidual_, b },
                              level,
                              All );
         timingTree_->stop( "Inner Residual Computation" );
      }

      timingTree_->start( "Copying-Casting" );
      dotCorrection_.setToZero( level );
      dotResidual_.copyFrom( barResidual_, level );
      timingTree_->stop( "Copying-Casting" );

      //      // FIXME how to fix underflow or do we even want to fix underflow?
      //      {
      //         const auto avgResBef = dotResidual_.sumGlobal( level ) / walberla::numeric_cast<DotValueType>( dotResidual_.getNumberOfGlobalDoFs( level ) );
      //         if ( l2normResidualTooSmall(avgResBef) )
      //         {
      //            timingTree_->stop( "Iterative Refinement Solver" );
      //            return;
      //         }
      //      } // FIXME can I remove this?
      timingTree_->start( "InnerLoop" );
      for ( uint_t i = 0; i < correctionSteps_; ++i )
      {
         // c <- solve(Ac, r)
         //         timingTree_->start( "Inner Solve i=" + std::to_string(i) );
         //         timingTree_->start( "Inner Solve" );
         smoother_->solve( dotA_, dotCorrection_, dotResidual_, level );
         //         timingTree_->stop( "Inner Solve" );
         //         timingTree_->stop( "Inner Solve i=" + std::to_string(i) );

#ifndef PERFORMANCE_RUN
         // +++ Calculate Correction Norm +++
         if ( residualNode_ != nullptr && residualTree_ != nullptr )
         {
            DotValueType correctionL2Norm =
                std::sqrt( dotCorrection_.dotGlobal( dotCorrection_, level, DoFType::Inner ) /
                           walberla::numeric_cast< DotValueType >( dotCorrection_.getNumberOfGlobalInnerDoFs( level ) ) );
            WALBERLA_ROOT_SECTION()
            {
               [[maybe_unused]] auto CorrectionNodePtr =
                   residualTree_->insert( walberla::numeric_cast< BarValueType >( correctionL2Norm ),
                                          "correction_" + std::to_string( i ),
                                          "correction_lev" + std::to_string( level ),
                                          residualNode_ );
            }
         }
#endif
      }
      timingTree_->stop( "InnerLoop" );

      timingTree_->start( "Copying-Casting" );
      barCorrection_.copyFrom( dotCorrection_, level );
      timingTree_->stop( "Copying-Casting" );

      // x = x - c
      timingTree_->start( "Correct Approximation" );
      x.assign( { walberla::numeric_cast< BarValueType >( 1 ), walberla::numeric_cast< BarValueType >( -1 ) },
                { x, barCorrection_ },
                level,
                flag_ );
      timingTree_->stop( "Correct Approximation" );

      timingTree_->stop( "Iterative Refinement Solver" );
   }

 private:
   template < typename ValueType >
   bool l2normResidualTooSmall( const ValueType l2DotResidual )
   {
      real_t smallestPossibleFloatValue = 0;
      // FIXME use numerics_limits library instead of own values

#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
      if constexpr ( std::is_same_v< ValueType, walberla::float16 > )
      {
         smallestPossibleFloatValue = 1e-14;
      }
#endif // def WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
      if constexpr ( std::is_same_v< ValueType, walberla::float32 > )
      {
         smallestPossibleFloatValue = 1e-127;
      }
      if constexpr ( std::is_same_v< ValueType, walberla::float64 > )
      {
         smallestPossibleFloatValue = 1e-308;
      }

      if ( std::abs( l2DotResidual ) <= smallestPossibleFloatValue )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "\nResidual " << l2DotResidual << " in dot precision (" << typeid( DotValueType ).name()
                                                  << ") is too close to zero to continue.\n" );
         return true;
      }
      else
      {
         return false;
      }
   }

   uint_t minLevel_;
   uint_t maxLevel_;
   uint_t correctionSteps_;

   hyteg::DoFType flag_;

   std::shared_ptr< hyteg::Solver< SmootherOperator_t > > smoother_;

   const SmootherOperator_t& dotA_;

   const std::shared_ptr< PrimitiveStorage > storage_;
   DotFunctionType                           dotCorrection_;
   DotFunctionType                           dotResidual_;
   BarFunctionType                           barCorrection_;
   BarFunctionType                           barResidual_;
   bool                                      residualIsPrecalculated_;
   uint_t                                    levelPrecalculatedResidual_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

#ifndef PERFORMANCE_RUN
   Tree< BarValueType >*     residualTree_;
   TreeNode< BarValueType >* residualNode_;
#endif
};

} // namespace hyteg
