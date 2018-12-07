#pragma once

#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"

namespace hhg {

template < class OperatorType, class pressureBlockPreconditionerType >
class StokesBlockDiagonalPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   StokesBlockDiagonalPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage,
                                      uint_t                                     minLevel,
                                      uint_t                                     maxLevel,
                                      uint_t                                     velocityPreconditionSteps,
                                      std::shared_ptr< hhg::Solver< typename OperatorType::VelocityOperator_T > >
                                          velocityBlockPreconditioner = std::make_shared< hhg::EmptySolver< typename OperatorType::VelocityOperator_T > >() )
   : velocityPreconditionSteps_( velocityPreconditionSteps )
   , flag_( hhg::Inner | hhg::NeumannBoundary )
   , velocityBlockPreconditioner_( velocityBlockPreconditioner )
   , pressureBlockPreconditioner_( std::make_shared< pressureBlockPreconditionerType >( storage, minLevel, maxLevel ) )
   {}

   // y = M^{-1} * x
   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, uint_t level ) const override
   {
      b.assign( {1.0}, {x}, level, flag_ );

      for( uint_t steps = 0; steps < velocityPreconditionSteps_; steps++ )
      {
         velocityBlockPreconditioner_->solve( A.A, b.u, x.u, level );
         velocityBlockPreconditioner_->solve( A.A, b.v, x.v, level );
         velocityBlockPreconditioner_->solve( A.A, b.w, x.w, level );
      }

      pressureBlockPreconditioner_->apply( x.p, b.p, level, flag_, Replace );
   }

 private:
   uint_t                                                                      velocityPreconditionSteps_;
   hhg::DoFType                                                                flag_;
   std::shared_ptr< hhg::Solver< typename OperatorType::VelocityOperator_T > > velocityBlockPreconditioner_;
   std::shared_ptr< pressureBlockPreconditionerType >                          pressureBlockPreconditioner_;
};

} // namespace hhg