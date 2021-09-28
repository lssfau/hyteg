#pragma once

#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"

namespace hyteg {

class P1StokesBlockPreconditioner : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   typedef P1ConstantLaplaceOperator VelocityOperator_T;

   P1StokesBlockPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , P( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   P1ConstantLaplaceOperator A;
   P1LumpedMassOperator      P;

   bool hasGlobalCells_;
};

} // namespace hyteg
