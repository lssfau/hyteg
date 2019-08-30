#pragma once
#include <vector>

#include "tinyhhg_core/solvers/Solver.hpp"

namespace hyteg {

template < class OperatorType >
class IdentityPreconditioner : public Solver< OperatorType >
{
 public:
   IdentityPreconditioner()
   : updateType_( Replace )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary )
   {}

   void solve( const OperatorType&,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const uint_t                          level ) override
   {
      b.assign( {1.0}, {x}, level, flag_ );
   }

 private:
   UpdateType updateType_;
   DoFType    flag_;
};

} // namespace hyteg