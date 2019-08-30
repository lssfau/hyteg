#pragma once

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/StokesOperatorTraits.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1CoefficientOperator.hpp"

namespace hyteg {

class P1CoefficientStokesOperator
{
public:

  P1CoefficientStokesOperator(const std::shared_ptr< PrimitiveStorage > & storage, const std::shared_ptr<P1Function< real_t >>& coefficient, size_t minLevel, size_t maxLevel)
    : A_uu(storage, coefficient, minLevel, maxLevel),
      A_uv(storage, coefficient, minLevel, maxLevel),
      A_vu(storage, coefficient, minLevel, maxLevel),
      A_vv(storage, coefficient, minLevel, maxLevel),
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

  P1CoefficientEpsilonOperator_uu A_uu;
  P1CoefficientEpsilonOperator_uv A_uv;
  P1CoefficientEpsilonOperator_vu A_vu;
  P1CoefficientEpsilonOperator_vv A_vv;
  P1DivxOperator div_x;
  P1DivyOperator div_y;
  P1DivTxOperator divT_x;
  P1DivTyOperator divT_y;
  P1PSPGOperator pspg;
};

template<>
struct has_pspg_block< P1CoefficientStokesOperator > {
    static const bool value = true;
};

template<>
struct tensor_variant< P1CoefficientStokesOperator > {
  static const bool value = true;
};

}