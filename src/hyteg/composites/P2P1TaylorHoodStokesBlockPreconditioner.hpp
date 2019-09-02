#pragma once

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

class P2P1TaylorHoodStokesBlockPreconditioner : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:

   typedef P2ConstantLaplaceOperator VelocityOperator_T;

    P2P1TaylorHoodStokesBlockPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , P( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   P2ConstantLaplaceOperator A;
   P1LumpedMassOperator      P;

   bool hasGlobalCells_;
};

} // namespace hyteg
