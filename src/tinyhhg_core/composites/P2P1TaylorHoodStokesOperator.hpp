#pragma once

#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"

namespace hhg
{

class P2P1TaylorHoodStokesOperator
{
public:

  P2P1TaylorHoodStokesOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : A(storage, minLevel, maxLevel),
      div_x(storage, minLevel, maxLevel),
      div_y(storage, minLevel, maxLevel),
      divT_x(storage, minLevel, maxLevel),
      divT_y(storage, minLevel, maxLevel)
  {
  }

  void apply(P2P1TaylorHoodFunction<real_t>& src, P2P1TaylorHoodFunction<real_t>& dst, size_t level, DoFType flag)
  {
    A.apply(src.u, dst.u, level, flag, Replace);
    divT_x.getVertexToVertexOpr().apply(src.p, *dst.u.getVertexDoFFunction(), level, flag, Add);
    divT_x.getVertexToEdgeOpr().apply(src.p, *dst.u.getEdgeDoFFunction(), level, flag, Add);

    A.apply(src.v, dst.v, level, flag, Replace);
    divT_y.getVertexToVertexOpr().apply(src.p, *dst.v.getVertexDoFFunction(), level, flag, Add);
    divT_y.getVertexToEdgeOpr().apply(src.p, *dst.v.getEdgeDoFFunction(), level, flag, Add);

    div_x.getVertexToVertexOpr().apply(*src.u.getVertexDoFFunction(), dst.p, level, flag | DirichletBoundary, Replace);
    div_x.getEdgeToVertexOpr().apply(*src.u.getEdgeDoFFunction(), dst.p, level, flag | DirichletBoundary, Add);
    div_y.getVertexToVertexOpr().apply(*src.v.getVertexDoFFunction(), dst.p, level, flag | DirichletBoundary, Add);
    div_y.getEdgeToVertexOpr().apply(*src.v.getEdgeDoFFunction(), dst.p, level, flag | DirichletBoundary, Add);
  }

  P2ConstantLaplaceOperator A;
  P2ConstantDivxOperator div_x;
  P2ConstantDivyOperator div_y;
  P2ConstantDivTxOperator divT_x;
  P2ConstantDivTyOperator divT_y;
};

}
