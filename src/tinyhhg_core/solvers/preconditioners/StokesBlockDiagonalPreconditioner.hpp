#pragma once

#include "tinyhhg_core/solvers/cgsolver.hpp"

namespace hhg {

template<class F, class O1, class O2>
class StokesBlockDiagonalPreconditioner {
public:

  StokesBlockDiagonalPreconditioner(O1& velocityOpr, O2& pressureOpr, const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel)
      : velocityOpr_(velocityOpr), pressureOpr_(pressureOpr),
        velocityCG(storage, minLevel, maxLevel),
        pressureCG(storage, minLevel, maxLevel),
        r("r", storage, minLevel, maxLevel)
  {}

  // y = M^{-1} * x
  void apply(F &x, F &y, uint_t level, DoFType flag) {

    y.assign({1.0}, {&x}, level, flag);

    velocityCG.solve(velocityOpr_.A, y.u, x.u, r.u, level, 1e-3, 100, flag, false);
    velocityCG.solve(velocityOpr_.A, y.v, x.v, r.v, level, 1e-3, 100, flag, false);
    pressureCG.solve(pressureOpr_, y.p, x.p, r.p, level, 1e-3, 100, flag | hhg::DirichletBoundary, false);


  }

private:
  O1& velocityOpr_;
  O2& pressureOpr_;

  F r;

  CGSolver<decltype(F::u), decltype(O1::A)> velocityCG;
  CGSolver<decltype(F::u), O2> pressureCG;
};

}