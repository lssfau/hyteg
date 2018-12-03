#pragma once

#include "tinyhhg_core/solvers/GeometricMultigrid.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"

namespace hhg {

template< class F, class O1, class VelocityBlockSolver_T, class PressureBlockPreconditioner_T >
class StokesBlockDiagonalPreconditioner
{
public:

  StokesBlockDiagonalPreconditioner(O1& velocityOpr,
                                    VelocityBlockSolver_T velocityBlockSolver,
                                    PressureBlockPreconditioner_T pressureBlockPreconditioner,
                                    const std::shared_ptr<PrimitiveStorage> & storage,
                                    uint_t minLevel, uint_t maxLevel, uint_t numVcycles)
      : velocityOpr_(velocityOpr),
        velocityBlockSolver_( velocityBlockSolver ),
        pressureBlockPreconditioner_( pressureBlockPreconditioner ),
        r("r", storage, minLevel, maxLevel),
        numVcycles_( numVcycles )
  {}

  // y = M^{-1} * x
  void apply(F &x, F &y, uint_t level, DoFType flag) {

    y.assign({1.0}, {&x}, level, flag);

    for ( uint_t vcycles = 0; vcycles < numVcycles_; vcycles++ )
    {
      velocityBlockSolver_.solve(velocityOpr_, y.u, x.u, r.u, level, 1e-16, 10000, flag);
      velocityBlockSolver_.solve(velocityOpr_, y.v, x.v, r.v, level, 1e-16, 10000, flag);
      velocityBlockSolver_.solve(velocityOpr_, y.w, x.w, r.w, level, 1e-16, 10000, flag);
    }

    pressureBlockPreconditioner_.apply( x.p, y.p, level, flag | DirichletBoundary, Replace );
  }

private:

  O1& velocityOpr_;

  F r;
  uint_t numVcycles_;

  VelocityBlockSolver_T         velocityBlockSolver_;
  PressureBlockPreconditioner_T pressureBlockPreconditioner_;
};

}