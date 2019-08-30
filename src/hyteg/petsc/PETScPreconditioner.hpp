#pragma once

#include "PETScWrapper.hpp"

#include "PETScLUSolver.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

template<typename ValueType, template <class> class F, class O>
class PETScPreconditioner {
public:


  PETScPreconditioner(O& opr, std::shared_ptr<F<PetscInt>> &numerator, uint_t localSize, uint_t globalSize)
   : opr_(opr), petscSolver(numerator, localSize, globalSize) {}

  // y = M^{-1} * x
  void apply(F<ValueType> &x, F<ValueType> &y, uint_t level, DoFType flag) {
    petscSolver.solve(opr_,y,x,x,level,0,0,flag);
  }

private:
  O opr_;
  PETScLUSolver<ValueType, F, O> petscSolver;
};

}

#endif