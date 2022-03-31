#pragma once

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

class P1StokesBlockPreconditioner : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   typedef P1ConstantVectorLaplaceOperator VelocityOperator_T;

   P1StokesBlockPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , P( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1StokesFunction< idx_t >&            src,
                  const P1StokesFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
     A.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
     P.toMatrix( mat, src.p(), dst.p(), level, flag );
   }

   P1ConstantVectorLaplaceOperator A;
   P1LumpedMassOperator      P;

   bool hasGlobalCells_;
};

} // namespace hyteg
