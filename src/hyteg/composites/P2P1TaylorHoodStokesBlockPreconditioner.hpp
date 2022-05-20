#pragma once

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

class P2P1TaylorHoodStokesBlockPreconditioner
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ConstantVectorLaplaceOperator VelocityOperator_T;

   P2P1TaylorHoodStokesBlockPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , P( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
     A.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
     P.toMatrix( mat, src.p(), dst.p(), level, flag );
   }

   P2ConstantVectorLaplaceOperator A;
   P1LumpedMassOperator      P;

   bool hasGlobalCells_;
};

} // namespace hyteg
