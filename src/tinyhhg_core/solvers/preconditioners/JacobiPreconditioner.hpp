#pragma once

namespace hhg {

template<class F, class O>
class JacobiPreconditioner {
public:

  JacobiPreconditioner(const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel, O& opr, uint_t iterations) : opr_(opr), iterations_(iterations),
                                                                                                                                         tmp("jac_tmp", storage, minLevel, maxLevel)
  {}

  // y = M^{-1} * x
  void apply(F &x, F &y, uint_t level, DoFType flag) {
    y.assign({1.0}, {&x}, level, flag);

    for (uint_t i = 0; i < iterations_; ++i) {
      tmp.assign({1.0}, {&y}, level, flag);
      opr_.smooth_jac(y, x, tmp, level, flag);
    }
  }

private:
  O& opr_;
  uint_t iterations_;
  F tmp;
};

}