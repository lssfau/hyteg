#pragma once

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

class P1EpsilonStokesOperator
{
public:

  P1EpsilonStokesOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
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

  P1ConstantEpsilonOperator_11 A_uu;
  P1ConstantEpsilonOperator_12 A_uv;
  P1ConstantEpsilonOperator_21 A_vu;
  P1ConstantEpsilonOperator_22 A_vv;
  P1DivxOperator div_x;
  P1DivyOperator div_y;
  P1DivTxOperator divT_x;
  P1DivTyOperator divT_y;
  P1PSPGOperator pspg;
};

template<>
struct has_pspg_block< P1EpsilonStokesOperator > {
    static const bool value = true;
};

struct tensor_variant< P1EpsilonStokesOperator > {
  static const bool value = true;
};

}
