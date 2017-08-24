#ifndef TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP
#define TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP

#include "tinyhhg_core/composites/p1stokesfunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"

namespace hhg
{

class P1BlockLaplaceOperator
{
public:

  P1BlockLaplaceOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : A(storage, minLevel, maxLevel)
  {
  }

  void apply(P1StokesFunction& src, P1StokesFunction& dst, size_t level, DoFType flag, UpdateType updateType)
  {
    A.apply(src.u, dst.u, level, flag, updateType);
    A.apply(src.v, dst.v, level, flag, updateType);
    A.apply(src.p, dst.p, level, flag, updateType);
  }

  void save(P1StokesFunction& src, P1StokesFunction& dst, size_t level, DoFType flag)
  {
    WALBERLA_UNUSED(src);
    WALBERLA_UNUSED(dst);
    WALBERLA_UNUSED(level);
    WALBERLA_UNUSED(flag);
    WALBERLA_ABORT("not implemented");
  }

  P1LaplaceOperator A;
};

}

#endif //TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP
