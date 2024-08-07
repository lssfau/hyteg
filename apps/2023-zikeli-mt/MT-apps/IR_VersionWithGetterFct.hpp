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

namespace hyteg {

using walberla::uint_t;

template < class SmootherOperator_t, class ResidualOperator_t = SmootherOperator_t >
class IterativeRefinementSolver : public Solver< ResidualOperator_t >
{
 public:
   using DashFunctionType = typename ResidualOperator_t::srcType;
   using DotFunctionType  = typename SmootherOperator_t::srcType;
   using DashValueType    = typename FunctionTrait< DashFunctionType >::ValueType;
   using DotValueType     = typename FunctionTrait< DotFunctionType >::ValueType;

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
                              uint_t              maxLevel, // FIXME this is actually not necessary since it is never used.
                              uint_t              smoothingSteps,
                              const bool          computeAbsoluteL2ResidualNorm = true,
                              const DashValueType absoluteL2ResidualThreshold   = 2 *
                                                                                std::numeric_limits< DashValueType >::epsilon(),
                              const real_t absoluteL2ResidualRateTerminationThreshold = 1.09 )
   : minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , smoothingSteps_( smoothingSteps )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary ) // FIXME I have no idea, I just copied this from GMG
   , smoother_( smoother )
   , dotA_( dotA )
   , storage_( storage )
   , dotCorrection_( "DotCorrection", storage, minLevel, maxLevel )
   , dotResidual_( "DotResidual", storage, minLevel, maxLevel )
   , dashCorrection_( "DashCorrection", storage, minLevel, maxLevel )
   , dashResidual_( "DashResidual", storage, minLevel, maxLevel )
   , computeAbsoluteL2ResidualNorm_( computeAbsoluteL2ResidualNorm )
   , absoluteL2ResidualNorm_( std::numeric_limits< DashValueType >::max() )
   , oldAbsoluteL2ResidualNorm_( std::numeric_limits< DashValueType >::max() )
   , absoluteL2ResidualRate_( std::numeric_limits< real_t >::max() )
   , levelOfAbsoluteL2ResidualNorm_()
   , absoluteL2ResidualThreshold_( absoluteL2ResidualThreshold )
   , absoluteL2ResidualRateTerminationThreshold_( absoluteL2ResidualRateTerminationThreshold )
   , timingTree_( storage->getTimingTree() )
   {}

   ~IterativeRefinementSolver() = default;

   void
       solve( const ResidualOperator_t& dashA, const DashFunctionType& x, const DashFunctionType& b, const uint_t level ) override
   {
      timingTree_->start( "Iterative Refinement Solver" );

      // TODO allocate variables
      //    res = calculation // @Nils the rhs is overwritten before restriction, so it doesn't matter what the RHS is.
      //    c = 0             // @Nils same argument

      // r = Ax - b
      dashA.apply( x, dashResidual_, level, All, Replace );
      dashResidual_.assign( { walberla::numeric_cast< DashValueType >( 1 ), walberla::numeric_cast< DashValueType >( -1 ) },
                            { dashResidual_, b },
                            level,
                            All );

      if ( computeAbsoluteL2ResidualNorm_ )
      {
         absoluteL2ResidualNorm_ =
             std::sqrt( dashResidual_.dotGlobal( dashResidual_, level, DoFType::All ) /
                        walberla::numeric_cast< DashValueType >()( dashResidual_.getNumberOfGlobalDoFs( level ) ) );
         levelOfAbsoluteL2ResidualNorm_ = level;
         absoluteL2ResidualRate_        = oldAbsoluteL2ResidualNorm_ / absoluteL2ResidualNorm_;
         oldAbsoluteL2ResidualNorm_     = absoluteL2ResidualNorm_;
         if ( absoluteL2ResidualNorm_ < absoluteL2ResidualThreshold_ |
              absoluteL2ResidualRateTerminationThreshold_ < absoluteL2ResidualRate_ )
         {
            return;
         }
      }

      dotCorrection_.setToZero( level );

      // FIXME how to fix underflow or do we even want to fix underflow?
      {
         dotResidual_.copyFrom( dashResidual_, level );
         const auto avgResBef = dotResidual_.sumGlobal( level ) /
                                walberla::numeric_cast< DotValueType >( dotResidual_.getNumberOfGlobalDoFs( level ) );
         WALBERLA_LOG_INFO_ON_ROOT( "Residual before in dot : " << avgResBef );
         if ( l2normResidualTooSmall( avgResBef ) )
         {
            return;
         }
      } // FIXME can I remove this?
      for ( uint_t i = 0; i < smoothingSteps_; ++i )
      {
         timingTree_->start( "Inner Solve" );
         // c <- solve(Ac, r)
         // TODO look up to throw an exception once an nan ocurres.
         smoother_->solve( dotA_, dotCorrection_, dotResidual_, level );
         timingTree_->stop( "Inner Solve" );
      }
      dashCorrection_.copyFrom( dotCorrection_, level );

      // x = x - c
      x.assign( { walberla::numeric_cast< DashValueType >( 1 ), walberla::numeric_cast< DashValueType >( -1 ) },
                { x, dashCorrection_ },
                level,
                flag_ );
      timingTree_->stop( "Iterative Refinement Solver" );
   }

   /// \brief Getter function for the residual in all levels.
   ///
   ///     FIXME: This function is horrible...
   ///     NOTE: This gives you only the current residual, if the inner solve step has not been executed
   ///             otherwise it will return the residual before execution.
   ///
   const DashFunctionType& getResidual() const { return dashResidual_; }

   /// \brief Getter function for the residual.
   ///
   ///     If the required level does not match the previously computed one, the program terminates.
   ///
   /// \param level                       level for which the residual norm is required.
   ///
   const DashValueType& getAbsoluteL2ResidualNorm( const uint_t level ) const
   {
      WALBERLA_CHECK( computeAbsoluteL2ResidualNorm_ );
      WALBERLA_CHECK_EQUAL( level, levelOfAbsoluteL2ResidualNorm_ );
      return absoluteL2ResidualNorm_;
   }

   /// \brief Getter function for the residual rate.
   const real_t& getAbsoluteL2ResidualRate() const
   {
      WALBERLA_CHECK( computeAbsoluteL2ResidualNorm_ );
      return absoluteL2ResidualRate_;
   }

 private:
   uint_t minLevel_;
   uint_t maxLevel_;
   uint_t smoothingSteps_;

   hyteg::DoFType flag_;

   std::shared_ptr< hyteg::Solver< SmootherOperator_t > > smoother_;

   const SmootherOperator_t& dotA_;

   const std::shared_ptr< PrimitiveStorage > storage_;
   DotFunctionType                           dotCorrection_;
   DotFunctionType                           dotResidual_;
   DashFunctionType                          dashCorrection_;
   DashFunctionType                          dashResidual_;
   const bool                                computeAbsoluteL2ResidualNorm_;
   DashValueType                             absoluteL2ResidualNorm_;
   DashValueType                             oldAbsoluteL2ResidualNorm_;
   real_t                                    absoluteL2ResidualRate_;
   uint_t                                    levelOfAbsoluteL2ResidualNorm_;
   const DashValueType                       absoluteL2ResidualThreshold_;
   const DashValueType                       absoluteL2ResidualRateTerminationThreshold_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

   template < typename ValueType >
   bool l2normResidualTooSmall( const ValueType l2DotResidual )
   {
      real_t smallestPossibleFloatValue = 0;
      // FIXME use numerics_limits library instead of own values
      if constexpr ( std::is_same_v< ValueType, walberla::float16 > )
      {
         smallestPossibleFloatValue = std::numeric_limits< ValueType >::min();
      }

      if ( std::abs( l2DotResidual ) <= smallestPossibleFloatValue )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "\nResidual in dot precision (" << typeid( DotValueType ).name()
                                                                    << ") is too close to zero to continue.\n" );
         return true;
      }
      else
      {
         return false;
      }
   }
};

} // namespace hyteg
