#pragma once

#include "PETScWrapper.hpp"

#include "PETScLUSolver.hpp"

#ifdef HHG_BUILD_WITH_PETSC

namespace hhg {

template<class F, class O>
class PETScPreconditioner {
public:

  PETScPreconditioner(O& opr, std::shared_ptr<F> &numerator, uint_t localSize, uint_t globalSize)
   : opr_(opr), petscSolver(numerator, localSize, globalSize) {}

  // y = M^{-1} * x
  void apply(F &x, F &y, uint_t level, DoFType flag) {
    petscSolver.solve(opr_,y,x,x,level,0,0,flag);
  }

private:
  O opr_;
  PETScLUSolver<F, O> petscSolver;
};

}

#endif