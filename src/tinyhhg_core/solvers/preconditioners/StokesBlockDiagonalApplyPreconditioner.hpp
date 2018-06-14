#pragma once

namespace hhg {

template<class F, class O1, class O2 >
class StokesBlockDiagonalApplyPreconditioner {
public:

    StokesBlockDiagonalApplyPreconditioner(O1 velocityOpr, O2 pressureOpr, const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel)
      : velocityOpr_(velocityOpr), pressureOpr_(pressureOpr),
        r("r", storage, minLevel, maxLevel)
  {}

  // y = M^{-1} * x
  void apply(F &x, F &y, uint_t level, DoFType flag)
  {

    y.assign({1.0}, {&x}, level, flag);

    velocityOpr_.apply( x.u, y.u, level, flag, Replace );
    velocityOpr_.apply( x.v, y.v, level, flag, Replace );
    pressureOpr_.apply( x.p, y.p, level, flag | DirichletBoundary, Replace );


  }

private:
  O1 velocityOpr_;
  O2 pressureOpr_;

  F r;
};

}