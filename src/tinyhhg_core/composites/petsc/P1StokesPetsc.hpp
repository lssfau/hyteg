#pragma once

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Petsc.hpp"

namespace hhg {
namespace petsc {

inline void createVectorFromFunction(P1StokesFunction<PetscScalar> &function,
                                     P1StokesFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createVectorFromFunction(function.u, numerator.u, vec, level, flag);
  createVectorFromFunction(function.v, numerator.v, vec, level, flag);
  createVectorFromFunction(function.p, numerator.p, vec, level, flag);
}

inline void createFunctionFromVector(P1StokesFunction<PetscScalar> &function,
                                     P1StokesFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createFunctionFromVector(function.u, numerator.u, vec, level, flag);
  createFunctionFromVector(function.v, numerator.v, vec, level, flag);
  createFunctionFromVector(function.p, numerator.p, vec, level, flag);
}

inline void applyDirichletBC(P1StokesFunction<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {
  applyDirichletBC(numerator.u, mat, level);
  applyDirichletBC(numerator.v, mat, level);
//  applyDirichletBC(numerator.p, mat, level);
}

template<class OperatorType>
inline void createMatrix(OperatorType& opr, P1StokesFunction< PetscInt > & src, P1StokesFunction< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  createMatrix(opr.A, src.u, dst.u, mat, level, flag);
  createMatrix(opr.divT_x, src.p, dst.u, mat, level, flag);

  createMatrix(opr.A, src.v, dst.v, mat, level, flag);
  createMatrix(opr.divT_y, src.p, dst.v, mat, level, flag);

  createMatrix(opr.div_x, src.u, dst.p, mat, level, flag | DirichletBoundary);
  createMatrix(opr.div_y, src.v, dst.p, mat, level, flag | DirichletBoundary);
  createMatrix(opr.pspg, src.p, dst.p, mat, level, flag | DirichletBoundary);

}

}
}