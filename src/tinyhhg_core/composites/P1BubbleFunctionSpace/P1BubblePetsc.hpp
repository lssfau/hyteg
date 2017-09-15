#pragma once

namespace hhg {
namespace petsc {

inline void createVectorFromFunction(P1BubbleFunction<PetscScalar> &function,
                                     P1BubbleFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createVectorFromFunction(function.p1, numerator.p1, vec, level, flag);
  createVectorFromFunction(function.b, numerator.b, vec, level, flag);
}

inline void createFunctionFromVector(P1BubbleFunction<PetscScalar> &function,
                                     P1BubbleFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createFunctionFromVector(function.p1, numerator.p1, vec, level, flag);
  createFunctionFromVector(function.b, numerator.b, vec, level, flag);
}

inline void applyDirichletBC(P1BubbleFunction<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {
  applyDirichletBC(numerator.p1, mat, level);
//  applyDirichletBC(numerator.b, mat, level);
}

template<class OperatorType>
inline void createMatrix(OperatorType& opr, P1BubbleFunction< PetscInt > & src, P1BubbleFunction< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  createMatrix(opr.A_p1, src.p1, dst.p1, mat, level, flag);
  createMatrix(opr.A_b, src.b, dst.b, mat, level, flag);
}

}
}