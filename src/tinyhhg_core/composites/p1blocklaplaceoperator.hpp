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

  void apply(P1StokesFunction& src, P1StokesFunction& dst, size_t level, DoFType flag)
  {
    A.apply(src.u, dst.u, level, flag);
    A.apply(src.v, dst.v, level, flag);
    A.apply(src.p, dst.p, level, flag);
  }

  P1LaplaceOperator A;
};

}

#endif //TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP
