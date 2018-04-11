#pragma once

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1BlendingOperatorNew.hpp"

namespace hhg
{

class P1BlendingStokesOperator
{
public:

  P1BlendingStokesOperator(const std::shared_ptr< PrimitiveStorage > & storage, uint_t minLevel, uint_t maxLevel)
    : A_uu(storage, minLevel, maxLevel),
      A_uv(storage, minLevel, maxLevel),
      A_vu(storage, minLevel, maxLevel),
      A_vv(storage, minLevel, maxLevel),
      div_x(storage, minLevel, maxLevel),
      div_y(storage, minLevel, maxLevel),
      divT_x(storage, minLevel, maxLevel),
      divT_y(storage, minLevel, maxLevel),
      pspg(storage, minLevel, maxLevel)
  {
  }

  void apply(P1StokesFunction<real_t>& src, P1StokesFunction<real_t>& dst, size_t level, DoFType flag)
  {
    A_uu.apply(src.u, dst.u, level, flag, Replace);
    A_uv.apply(src.v, dst.u, level, flag, Add);
    divT_x.apply(src.p, dst.u, level, flag, Add);

    A_vu.apply(src.u, dst.v, level, flag, Replace);
    A_vv.apply(src.v, dst.v, level, flag, Add);
    divT_y.apply(src.p, dst.v, level, flag, Add);

    div_x.apply(src.u, dst.p, level, flag | DirichletBoundary, Replace);
    div_y.apply(src.v, dst.p, level, flag | DirichletBoundary, Add);
    pspg.apply(src.p, dst.p, level, flag | DirichletBoundary, Add);
  }

  P1BlendingEpsilonOperator_11 A_uu;
  P1BlendingEpsilonOperator_12 A_uv;
  P1BlendingEpsilonOperator_21 A_vu;
  P1BlendingEpsilonOperator_22 A_vv;
  P1BlendingDivOperator_1 div_x;
  P1BlendingDivOperator_2 div_y;
  P1BlendingDivTOperator_1 divT_x;
  P1BlendingDivTOperator_2 divT_y;
  P1PSPGOperator pspg;
};

}
