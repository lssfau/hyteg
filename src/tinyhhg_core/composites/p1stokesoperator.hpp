#ifndef TINYHHG_P1STOKESOPERATOR_HPP
#define TINYHHG_P1STOKESOPERATOR_HPP

#include "tinyhhg_core/composites/p1stokesfunction.hpp"
#include "tinyhhg_core/p1functionspace/p1operator.hpp"

namespace hhg
{

class P1StokesOperator
{
public:

  P1StokesOperator(Mesh& mesh, size_t minLevel, size_t maxLevel)
    : A(mesh, minLevel, maxLevel),
      div_x(mesh, minLevel, maxLevel),
      div_y(mesh, minLevel, maxLevel),
      divT_x(mesh, minLevel, maxLevel),
      divT_y(mesh, minLevel, maxLevel),
      pspg(mesh, minLevel, maxLevel)
  {
  }

  void apply(P1StokesFunction& src, P1StokesFunction& dst, size_t level, DoFType flag)
  {
    A.apply(src.u, dst.u, level, flag, Replace);
    divT_x.apply(src.p, dst.u, level, flag, Add);

    A.apply(src.v, dst.v, level, flag, Replace);
    divT_y.apply(src.p, dst.v, level, flag, Add);

    div_x.apply(src.u, dst.p, level, flag | DirichletBoundary, Replace);
    div_y.apply(src.v, dst.p, level, flag | DirichletBoundary, Add);
    pspg.apply(src.p, dst.p, level, flag | DirichletBoundary, Add);
  }

  P1LaplaceOperator A;
  P1DivxOperator div_x;
  P1DivyOperator div_y;
  P1DivTxOperator divT_x;
  P1DivTyOperator divT_y;
  P1PSPGOperator pspg;
};

}

#endif //TINYHHG_P1STOKESOPERATOR_HPP
