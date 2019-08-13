
#pragma once

#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/solvers/Solver.hpp"

namespace hhg {

template < class OperatorType >
class FullMultigridSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   FullMultigridSolver( const std::shared_ptr< PrimitiveStorage >&                         storage,
                        const std::shared_ptr< GeometricMultigridSolver< OperatorType > >& gmgSolver,
                        const std::shared_ptr< ProlongationOperator< FunctionType > >&     fmgProlongation,
                        const uint_t&                                                      minLevel,
                        const uint_t&                                                      maxLevel,
                        const uint_t&                                                      cyclesPerLevel = 1,
                        const std::function< void( uint_t currentLevel ) >&                postCycleCallback = []( uint_t ){} )
   : gmgSolver_( gmgSolver )
   , fmgProlongation_( fmgProlongation )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , cyclesPerLevel_( cyclesPerLevel )
   , flag_( Inner | NeumannBoundary )
   , postCycleCallback_( postCycleCallback )
   , timingTree_( storage->getTimingTree() )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      timingTree_->start( "FMG Solver" );
      for ( uint_t currentLevel = minLevel_; currentLevel <= level; currentLevel++ )
      {
         timingTree_->start( "GMG Solver" );
         for ( uint_t cycle = 0; cycle < cyclesPerLevel_; cycle++ )
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
      }
      timingTree_->stop( "FMG Solver" );
   }

 private:
   std::shared_ptr< GeometricMultigridSolver< OperatorType > > gmgSolver_;
   std::shared_ptr< ProlongationOperator< FunctionType > >     fmgProlongation_;

   uint_t minLevel_;
   uint_t maxLevel_;

   uint_t cyclesPerLevel_;

   DoFType      flag_;

   std::function< void( uint_t currentLevel ) > postCycleCallback_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;
};

} // namespace hhg