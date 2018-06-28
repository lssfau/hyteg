#ifndef TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP
#define TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"

namespace hhg
{

class P1BlockLaplaceOperator
{
public:

  P1BlockLaplaceOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : A(storage, minLevel, maxLevel)
  {
  }

  void apply(P1StokesFunction<real_t>& src, P1StokesFunction<real_t>& dst, size_t level, DoFType flag, UpdateType updateType)
  {
    A.apply(src.u, dst.u, level, flag, updateType);
    A.apply(src.v, dst.v, level, flag, updateType);
    A.apply(src.p, dst.p, level, flag, updateType);
  }

  P1ConstantLaplaceOperator A;
};

}

#endif //TINYHHG_P1BLOCKLAPLACEOPERATOR_HPP
