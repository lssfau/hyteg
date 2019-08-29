
#pragma once

#include "core/DataTypes.h"
#include "core/timing/TimingTree.h"

#include "tinyhhg_core/gridtransferoperators/ProlongationOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/RestrictionOperator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/solvers/Solver.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {

using walberla::real_t;
using walberla::uint_t;

template < class OperatorType >
class FASSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   FASSolver( const std::shared_ptr< PrimitiveStorage >&              storage,
              std::shared_ptr< Solver< OperatorType > >               smoother,
              std::shared_ptr< Solver< OperatorType > >               coarseSolver,
              std::shared_ptr< RestrictionOperator< FunctionType > >  restrictionOperator,
              std::shared_ptr< RestrictionOperator< FunctionType > >  solutionRestrictionOperator,
              std::shared_ptr< ProlongationOperator< FunctionType > > prolongationOperator,
              uint_t                                                  minLevel,
              uint_t                                                  maxLevel,
              uint_t                                                  preSmoothSteps                = 3,
              uint_t                                                  postSmoothSteps               = 3,
              uint_t                                                  smoothIncrementOnCoarserGrids = 0,
              CycleType                                               cycleType                     = CycleType::VCYCLE )
   : minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , smoother_( smoother )
   , coarseSolver_( coarseSolver )
   , restrictionOperator_( restrictionOperator )
   , solutionRestrictionOperator_( solutionRestrictionOperator )
   , prolongationOperator_( prolongationOperator )
   , tmp_( "fas_tmp", storage, minLevel, maxLevel )
   , d_( "fas_d", storage, minLevel, maxLevel )
   , w_( "fas_w", storage, minLevel, maxLevel )
   , preSmoothSteps_( preSmoothSteps )
   , postSmoothSteps_( postSmoothSteps )
   , smoothIncrement_( smoothIncrementOnCoarserGrids )
   , flag_( hhg::Inner | hhg::NeumannBoundary )
   , cycleType_( cycleType )
   , timingTree_( storage->getTimingTree() )
   {}

   ~FASSolver() = default;

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      timingTree_->start( "FAS Multigrid Solver" );
      invokedLevel_ = level;
      solveRecursively( A, x, b, level );
      timingTree_->stop( "FAS Multigrid Solver" );
   }

   void solveRecursively( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) const
   {
      if ( level == minLevel_ )
      {
         timingTree_->start( "Coarse Grid Solver" );
         coarseSolver_->solve( A, x, b, minLevel_ );
         timingTree_->stop( "Coarse Grid Solver" );
      }
      else
      {
         // pre-smooth
         const uint_t preSmoothingSteps = preSmoothSteps_ + smoothIncrement_ * ( invokedLevel_ - level );
         for ( uint_t i = 0; i < preSmoothingSteps; ++i )
         {
            timingTree_->start( "Smoother" );
            smoother_->solve( A, x, b, level );
            timingTree_->stop( "Smoother" );
         }

         A.apply( x, tmp_, level, flag_ );
         d_.assign( {1.0, -1.0}, {b, tmp_}, level, flag_ );

         // restrict
         timingTree_->start( "Restriction" );
         restrictionOperator_->restrict( d_, level, flag_ );
         solutionRestrictionOperator_->restrict( x, level, flag_ );
         timingTree_->stop( "Restriction" );

         A.apply( x, tmp_, level - 1, flag_ );
         b.assign( {1.0, 1.0}, {d_, tmp_}, level - 1, flag_ );

         w_.assign( {1.0}, {x}, level - 1, All );

         solveRecursively( A, x, b, level - 1 );

         if ( cycleType_ == CycleType::WCYCLE )
         {
            solveRecursively( A, x, b, level - 1 );
         }

         // coarse grid correction
         tmp_.assign( {1.0, -1.0}, {x, w_}, level - 1, flag_ );
         timingTree_->start( "Prolongation" );
         prolongationOperator_->prolongate( tmp_, level - 1, flag_ );
         timingTree_->stop( "Prolongation" );
         x.add( {1.0}, {tmp_}, level, flag_ );

         // post-smooth
         const uint_t postSmoothingSteps = postSmoothSteps_ + smoothIncrement_ * ( invokedLevel_ - level );
         for ( size_t i = 0; i < postSmoothingSteps; ++i )
         {
            timingTree_->start( "Smoother" );
            smoother_->solve( A, x, b, level );
            timingTree_->stop( "Smoother" );
         }
      }
   }

 private:
   uint_t minLevel_;
   uint_t maxLevel_;
   uint_t preSmoothSteps_;
   uint_t postSmoothSteps_;
   uint_t smoothIncrement_;
   uint_t invokedLevel_;

   hhg::DoFType flag_;
   CycleType    cycleType_;

   std::shared_ptr< hhg::Solver< OperatorType > >               smoother_;
   std::shared_ptr< hhg::Solver< OperatorType > >               coarseSolver_;
   std::shared_ptr< hhg::RestrictionOperator< FunctionType > >  restrictionOperator_;
   std::shared_ptr< hhg::RestrictionOperator< FunctionType > >  solutionRestrictionOperator_;
   std::shared_ptr< hhg::ProlongationOperator< FunctionType > > prolongationOperator_;

   FunctionType tmp_;
   FunctionType d_;
   FunctionType w_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;
};

} // namespace hhg
