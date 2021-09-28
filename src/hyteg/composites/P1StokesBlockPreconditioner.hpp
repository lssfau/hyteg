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
                  const P1StokesFunction< matIdx_t >&         src,
                  const P1StokesFunction< matIdx_t >&         dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      A.toMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
      A.toMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );

      if ( src.uvw[0].getStorage()->hasGlobalCells() )
      {
         A.toMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
      }

      P.toMatrix( mat, src.p, dst.p, level, flag );
   }

   P1ConstantLaplaceOperator A;
   P1LumpedMassOperator      P;

   bool hasGlobalCells_;
};

} // namespace hyteg
