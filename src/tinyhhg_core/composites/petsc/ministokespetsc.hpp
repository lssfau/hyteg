#pragma once

namespace hhg {
namespace petsc {

inline void createVectorFromFunction(MiniStokesFunction<PetscScalar> &function,
                                     MiniStokesFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createVectorFromFunction(function.u, numerator.u, vec, level, flag);
  createVectorFromFunction(function.v, numerator.v, vec, level, flag);
  createVectorFromFunction(function.p, numerator.p, vec, level, flag);
}

inline void createFunctionFromVector(MiniStokesFunction<PetscScalar> &function,
                                     MiniStokesFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createFunctionFromVector(function.u, numerator.u, vec, level, flag);
  createFunctionFromVector(function.v, numerator.v, vec, level, flag);
  createFunctionFromVector(function.p, numerator.p, vec, level, flag);
}

inline void applyDirichletBC(MiniStokesFunction<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {
  applyDirichletBC(numerator.u, mat, level);
  applyDirichletBC(numerator.v, mat, level);
//  applyDirichletBC(numerator.p, mat, level);
}

template<class OperatorType>
inline void createMatrix(OperatorType& opr, MiniStokesFunction< PetscInt > & src, MiniStokesFunction< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  createMatrix(opr.A, src.u, dst.u, mat, level, flag);
  createMatrix(opr.divT_x_p1, src.p, dst.u.p1, mat, level, flag);
  createMatrix(opr.divT_x_b, src.p, dst.u.b, mat, level, flag);

  createMatrix(opr.A, src.v, dst.v, mat, level, flag);
  createMatrix(opr.divT_y_p1, src.p, dst.v.p1, mat, level, flag);
  createMatrix(opr.divT_y_b, src.p, dst.v.b, mat, level, flag);

  createMatrix(opr.div_x_p1, src.u.p1, dst.p, mat, level, flag | DirichletBoundary);
  createMatrix(opr.div_x_b, src.u.b, dst.p, mat, level, flag | DirichletBoundary);
  createMatrix(opr.div_y_p1, src.v.p1, dst.p, mat, level, flag | DirichletBoundary);
  createMatrix(opr.div_y_b, src.v.b, dst.p, mat, level, flag | DirichletBoundary);

}

}
}