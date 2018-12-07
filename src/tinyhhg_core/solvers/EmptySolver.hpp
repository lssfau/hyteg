
#pragma once

#include "tinyhhg_core/solvers/Solver.hpp"

namespace hhg {

template < class OperatorType >
class EmptySolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   void solve( const OperatorType&, const FunctionType&, const FunctionType&, uint_t ) const override
   {
      // does nothing
   }
};

} // namespace hhg