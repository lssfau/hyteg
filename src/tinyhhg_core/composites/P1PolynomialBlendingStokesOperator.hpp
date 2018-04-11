#pragma once

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/p1functionspace/P1PolynomialBlendingOperatorNew.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"

namespace hhg {

class P1PolynomialBlendingStokesOperator
{
 public:
   P1PolynomialBlendingStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                       uint_t                                     minLevel,
                                       uint_t                                     maxLevel,
                                       uint_t                                     interpolationLevel )
   : A_uu( storage, minLevel, maxLevel, interpolationLevel )
   , A_uv( storage, minLevel, maxLevel, interpolationLevel )
   , A_vu( storage, minLevel, maxLevel, interpolationLevel )
   , A_vv( storage, minLevel, maxLevel, interpolationLevel )
   , div_x( storage, minLevel, maxLevel, interpolationLevel )
   , div_y( storage, minLevel, maxLevel, interpolationLevel )
   , divT_x( storage, minLevel, maxLevel, interpolationLevel )
   , divT_y( storage, minLevel, maxLevel, interpolationLevel )
   , pspg( storage, minLevel, maxLevel )
   {}

   void interpolateStencils( uint_t polyDegree )
   {
      A_uu.interpolateStencils( polyDegree );
      A_uv.interpolateStencils( polyDegree );
      A_vu.interpolateStencils( polyDegree );
      A_vv.interpolateStencils( polyDegree );
      div_x.interpolateStencils( polyDegree );
      div_y.interpolateStencils( polyDegree );
      divT_x.interpolateStencils( polyDegree );
      divT_y.interpolateStencils( polyDegree );
   }

   void useDegree( uint_t degree )
   {
      A_uu.useDegree( degree );
      A_uv.useDegree( degree );
      A_vu.useDegree( degree );
      A_vv.useDegree( degree );
      div_x.useDegree( degree );
      div_y.useDegree( degree );
      divT_x.useDegree( degree );
      divT_y.useDegree( degree );
   }

   void apply( P1StokesFunction< real_t >& src, P1StokesFunction< real_t >& dst, size_t level, DoFType flag )
   {
      A_uu.apply( src.u, dst.u, level, flag, Replace );
      A_uv.apply( src.v, dst.u, level, flag, Add );
      divT_x.apply( src.p, dst.u, level, flag, Add );

      A_vu.apply( src.u, dst.v, level, flag, Replace );
      A_vv.apply( src.v, dst.v, level, flag, Add );
      divT_y.apply( src.p, dst.v, level, flag, Add );

      div_x.apply( src.u, dst.p, level, flag | DirichletBoundary, Replace );
      div_y.apply( src.v, dst.p, level, flag | DirichletBoundary, Add );
      pspg.apply( src.p, dst.p, level, flag | DirichletBoundary, Add );
   }

   P1PolynomialBlendingEpsilonOperator_11 A_uu;
   P1PolynomialBlendingEpsilonOperator_12 A_uv;
   P1PolynomialBlendingEpsilonOperator_21 A_vu;
   P1PolynomialBlendingEpsilonOperator_22 A_vv;
   P1PolynomialBlendingDivOperator_1      div_x;
   P1PolynomialBlendingDivOperator_2      div_y;
   P1PolynomialBlendingDivTOperator_1     divT_x;
   P1PolynomialBlendingDivTOperator_2     divT_y;
   P1PSPGOperator                         pspg;
};

} // namespace hhg