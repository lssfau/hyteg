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
   : minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , smoother_( smoother )
   , coarseSolver_( coarseSolver )
   , restrictionOperator_( restrictionOperator )
   , prolongationOperator_( prolongationOperator )
   , ax_( "gmg_ax", storage, minLevel, maxLevel )
   , tmp_( "gmg_tmp", storage, minLevel, maxLevel )
   , preSmoothSteps_( preSmoothSteps )
   , postSmoothSteps_( postSmoothSteps )
   , smoothIncrement_( smoothIncrementOnCoarserGrids )
   , flag_( hhg::Inner | hhg::NeumannBoundary )
   , cycleType_( cycleType )
   , timingTree_( storage->getTimingTree() )
   {
      uint_t inc = 0;
      for ( uint_t level = maxLevel; level >= minLevel; level-- )
      {
         preSmoothingPerLevel_[level]  = preSmoothSteps + inc * smoothIncrementOnCoarserGrids;
         postSmoothingPerLevel_[level] = postSmoothSteps + inc * smoothIncrementOnCoarserGrids;
         inc++;
      }
   }

   ~GeometricMultigridSolver() = default;

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      timingTree_->start( "Geometric Multigrid Solver" );
      solveRecursively( A, x, b, level );
      timingTree_->stop( "Geometric Multigrid Solver" );
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
         for ( size_t i = 0; i < preSmoothingPerLevel_.at( level ); ++i )
         {
            timingTree_->start( "Smoother" );
            smoother_->solve( A, x, b, level );
            timingTree_->stop( "Smoother" );
         }

         A.apply( x, ax_, level, flag_ );
         tmp_.assign( {1.0, -1.0}, {b, ax_}, level, flag_ );

         // restrict
         timingTree_->start( "Restriction" );
         restrictionOperator_->restrict( tmp_, level, flag_ );
         timingTree_->stop( "Restriction" );

         b.assign( {1.0}, {tmp_}, level - 1, flag_ );

         x.interpolate( 0, level - 1 );

         solveRecursively( A, x, b, level - 1 );

         if ( cycleType_ == CycleType::WCYCLE )
         {
            solveRecursively( A, x, b, level - 1 );
         }

         // prolongate
         tmp_.assign( {1.0}, {x}, level, flag_ );
         timingTree_->start( "Prolongation" );
         prolongationOperator_->prolongate( x, level - 1, flag_ );
         timingTree_->stop( "Prolongation" );
         x.add( {1.0}, {tmp_}, level, flag_ );

         // post-smooth
         for ( size_t i = 0; i < postSmoothingPerLevel_.at( level ); ++i )
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

   hhg::DoFType flag_;
   CycleType    cycleType_;

   std::shared_ptr< hhg::Solver< OperatorType > >               smoother_;
   std::shared_ptr< hhg::Solver< OperatorType > >               coarseSolver_;
   std::shared_ptr< hhg::RestrictionOperator< FunctionType > >  restrictionOperator_;
   std::shared_ptr< hhg::ProlongationOperator< FunctionType > > prolongationOperator_;

   FunctionType ax_;
   FunctionType tmp_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

   std::map< uint_t, uint_t > preSmoothingPerLevel_;
   std::map< uint_t, uint_t > postSmoothingPerLevel_;
};

} // namespace hhg
