#pragma once

#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/solvers/Solver.hpp"

namespace hhg {

template < class OperatorType, class pressureBlockPreconditionerType >
class StokesPressureBlockPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   StokesPressureBlockPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : pressureBlockPreconditioner_( std::make_shared< pressureBlockPreconditionerType >( storage, minLevel, maxLevel ) )
   , flag_( hhg::Inner | hhg::NeumannBoundary )
   {}

   void solve( const OperatorType&, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      b.assign( {1.0}, {x}, level, flag_ );
      pressureBlockPreconditioner_->apply( x.p, b.p, level, flag_, Replace );
   }

 private:
   std::shared_ptr< pressureBlockPreconditionerType > pressureBlockPreconditioner_;
   hhg::DoFType                                       flag_;
};

} // namespace hhg