#pragma once

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/elementwiseoperators/DiagonalNonConstantOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp"

namespace hyteg {

class P2P1ElementwiseBlendingStokesBlockPreconditioner
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ElementwiseBlendingLaplaceOperator VelocityOperator_T;

   P2P1ElementwiseBlendingStokesBlockPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage,
                                                     size_t                                     minLevel,
                                                     size_t                                     maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , P( storage,
        minLevel,
        maxLevel,
        // a lower quadrature order might well be sufficient for a preconditioner
        std::make_shared< P1RowSumForm >( std::make_shared< forms::p1_mass_blending_q4 >() ) )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      for ( uint_t dim = 0; dim < src.uvw().getDimension(); dim++ )
      {
         A.toMatrix( mat, src.uvw()[dim], dst.uvw()[dim], level, flag );
      }
      P.toMatrix( mat, src.p(), dst.p(), level, flag );
   }

   P2ElementwiseBlendingLaplaceOperator A;
   P1BlendingLumpedDiagonalOperator     P;

   bool hasGlobalCells_;
};

} // namespace hyteg
