#pragma once

#include "tinyhhg_core/solvers/GeometricMultiGrid.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"

namespace hhg {

template< class F, class PressureBlockPreconditioner_T >
class StokesPressureBlockPreconditioner
{
public:

    StokesPressureBlockPreconditioner(PressureBlockPreconditioner_T pressureBlockPreconditioner,
                                    const std::shared_ptr<PrimitiveStorage> & storage,
                                    uint_t minLevel, uint_t maxLevel)
      : pressureBlockPreconditioner_( pressureBlockPreconditioner ),
        r("r", storage, minLevel, maxLevel)
  {}

  // y = M^{-1} * x
  void apply(F &x, F &y, uint_t level, DoFType flag) {

    y.assign({1.0}, {&x}, level, flag);
    pressureBlockPreconditioner_.apply( x.p, y.p, level, flag, Replace );
  }

private:

  F r;

  PressureBlockPreconditioner_T pressureBlockPreconditioner_;
};

}