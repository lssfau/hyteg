#pragma once

#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/mixedoperators/P2ToP1Operator.hpp"
#include "tinyhhg_core/mixedoperators/P1ToP2Operator.hpp"

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
    divT_x.apply(src.p, dst.u, level, flag, Add);

    A.apply(src.v, dst.v, level, flag, Replace);
    divT_y.apply(src.p, dst.v, level, flag, Add);

    div_x.apply(src.u, dst.p, level, flag | DirichletBoundary, Replace);
    div_y.apply(src.v, dst.p, level, flag | DirichletBoundary, Add);
  }

  P2ConstantLaplaceOperator A;
  P2ToP1ConstantDivxOperator div_x;
  P2ToP1ConstantDivyOperator div_y;
  P1ToP2ConstantDivTxOperator divT_x;
  P1ToP2ConstantDivTyOperator divT_y;
};

}
