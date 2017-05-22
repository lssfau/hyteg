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

  template<size_t Level>
  void apply(P1StokesFunction& src, P1StokesFunction& dst, DoFType flag)
  {
    A.apply<Level>(src.u, dst.u, flag, Replace);
    divT_x.apply<Level>(src.p, dst.u, flag, Add);

    A.apply<Level>(src.v, dst.v, flag, Replace);
    divT_y.apply<Level>(src.p, dst.v, flag, Add);

    div_x.apply<Level>(src.u, dst.p, flag | DirichletBoundary, Replace);
    div_y.apply<Level>(src.v, dst.p, flag | DirichletBoundary, Add);
    pspg.apply<Level>(src.p, dst.p, flag | DirichletBoundary, Add);
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
