#pragma once

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

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

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1StokesFunction< idx_t >&            src,
                  const P1StokesFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      for ( uint_t k = 0; k < src.uvw().getDimension(); ++k )
      {
         A.toMatrix( mat, src.uvw()[0], dst.uvw()[0], level, flag );
      }

      P.toMatrix( mat, src.p(), dst.p(), level, flag );
   }

   P1ConstantLaplaceOperator A;
   P1LumpedMassOperator      P;

   bool hasGlobalCells_;
};

} // namespace hyteg
