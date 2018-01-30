#pragma once

namespace hhg {
namespace petsc {

inline void createVectorFromFunction(P2P1TaylorHoodFunction<PetscScalar> &function,
                                     P2P1TaylorHoodFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createVectorFromFunction(function.u, numerator.u, vec, level, flag);
  createVectorFromFunction(function.v, numerator.v, vec, level, flag);
  createVectorFromFunction(function.p, numerator.p, vec, level, flag);
}

inline void createFunctionFromVector(P2P1TaylorHoodFunction<PetscScalar> &function,
                                     P2P1TaylorHoodFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createFunctionFromVector(function.u, numerator.u, vec, level, flag);
  createFunctionFromVector(function.v, numerator.v, vec, level, flag);
  createFunctionFromVector(function.p, numerator.p, vec, level, flag);
}

inline void applyDirichletBC(P2P1TaylorHoodFunction<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {
  applyDirichletBC(numerator.u, mat, level);
  applyDirichletBC(numerator.v, mat, level);
  applyDirichletBC(numerator.p, mat, level);
}

template<class OperatorType>
inline void createMatrix(OperatorType& opr, P2P1TaylorHoodFunction< PetscInt > & src, P2P1TaylorHoodFunction< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  createMatrix(opr.A, src.u, dst.u, mat, level, flag);
  createMatrix(opr.divT_x_vertexToVertex, src.p, *dst.u.getVertexDoFFunction(), mat, level, flag);
  VertexDoFToEdgeDoF::createMatrix(opr.divT_x_vertexToEdge, src.p, *dst.u.getEdgeDoFFunction(), mat, level, flag);

  createMatrix(opr.A, src.v, dst.v, mat, level, flag);
  createMatrix(opr.divT_y_vertexToVertex, src.p, *dst.v.getVertexDoFFunction(), mat, level, flag);
  VertexDoFToEdgeDoF::createMatrix(opr.divT_y_vertexToEdge, src.p, *dst.v.getEdgeDoFFunction(), mat, level, flag);

  createMatrix(opr.div_x_vertexToVertex, *src.u.getVertexDoFFunction(), dst.p, mat, level, flag | DirichletBoundary);
  EdgeDoFToVertexDoF::createMatrix(opr.div_x_edgeToVertex, *src.u.getEdgeDoFFunction(), dst.p, mat, level, flag | DirichletBoundary);
  createMatrix(opr.div_y_vertexToVertex, *src.v.getVertexDoFFunction(), dst.p, mat, level, flag | DirichletBoundary);
  EdgeDoFToVertexDoF::createMatrix(opr.div_y_edgeToVertex, *src.v.getEdgeDoFFunction(), dst.p, mat, level, flag | DirichletBoundary);
}

}
}