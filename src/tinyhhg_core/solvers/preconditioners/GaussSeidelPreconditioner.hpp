#pragma once

namespace hhg {

template<class F, class O>
class GaussSeidelPreconditioner {
public:

  GaussSeidelPreconditioner(O& opr, uint_t iterations) : opr_(opr), iterations_(iterations) {}

  // y = M^{-1} * x
  void apply(F &x, F &y, uint_t level, DoFType flag) {
    y.assign({1.0}, {&x}, level, flag);
    for (uint_t i = 0; i < iterations_; ++i) {
      opr_.smooth_gs(y, x, level, flag);
    }
  }

private:
  O& opr_;
  uint_t iterations_;
};

}