
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
                        const uint_t&                                                      cyclesPerLevel = 1 )
   : gmgSolver_( gmgSolver )
   , fmgProlongation_( fmgProlongation )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , cyclesPerLevel_( cyclesPerLevel )
   , Ax_( "Ax", storage, minLevel, maxLevel )
   , flag_( Inner | NeumannBoundary )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         for ( uint_t cycle = 0; cycle < cyclesPerLevel_; cycle++ )
         {
            gmgSolver_->solve( A, x, b, level );
         }

         if ( level < maxLevel_ )
         {
            fmgProlongation_->prolongate( x, level, flag_ );
         }
      }
   }

 private:
   std::shared_ptr< GeometricMultigridSolver< OperatorType > > gmgSolver_;
   std::shared_ptr< ProlongationOperator< FunctionType > >     fmgProlongation_;

   uint_t minLevel_;
   uint_t maxLevel_;

   uint_t cyclesPerLevel_;

   FunctionType Ax_;
   DoFType      flag_;
};

} // namespace hhg