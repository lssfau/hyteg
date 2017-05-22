#ifndef TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP
#define TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP

#include "tinyhhg_core/composites/p1stokesfunction.hpp"
#include "tinyhhg_core/p1functionspace/p1operator.hpp"

namespace hhg
{

class P1BlockLaplaceOperator
{
public:

  P1BlockLaplaceOperator(Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
    : A(_mesh, _minLevel, _maxLevel)
  {
  }

  template<size_t Level>
  void apply(P1StokesFunction& src, P1StokesFunction& dst, DoFType flag)
  {
    A.apply<Level>(src.u, dst.u, flag);
    A.apply<Level>(src.v, dst.v, flag);
    A.apply<Level>(src.p, dst.p, flag);
  }

  P1LaplaceOperator A;
};

}

#endif //TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP
