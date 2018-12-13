#pragma once

#include <tinyhhg_core/p2functionspace/P2Function.hpp>

#include <tinyhhg_core/p1functionspace/P1Petsc.hpp>
#include <tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFPetsc.hpp>
#include <tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFPetsc.hpp>
#include <tinyhhg_core/edgedofspace/EdgeDoFPetsc.hpp>

namespace hhg {
namespace petsc {

inline void createVectorFromFunction(const P2Function<PetscScalar> &function,
                                     const P2Function<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createVectorFromFunction(function.getVertexDoFFunction(), numerator.getVertexDoFFunction(), vec, level, flag);
  edgedof::createVectorFromFunction(function.getEdgeDoFFunction(), numerator.getEdgeDoFFunction(), vec, level, flag);
}

inline void createFunctionFromVector(const P2Function<PetscScalar> &function,
                                     const P2Function<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createFunctionFromVector(function.getVertexDoFFunction(), numerator.getVertexDoFFunction(), vec, level, flag);
  edgedof::createFunctionFromVector(function.getEdgeDoFFunction(), numerator.getEdgeDoFFunction(), vec, level, flag);
}

inline void applyDirichletBC(const P2Function<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {
  applyDirichletBC(numerator.getVertexDoFFunction(), mat, level);
  edgedof::applyDirichletBC(numerator.getEdgeDoFFunction(), mat, level);
}

template<class OperatorType>
inline void createMatrix(const OperatorType &opr,
                         const P2Function<PetscInt> &src,
                         const P2Function<PetscInt> &dst,
                         Mat &mat,
                         uint_t level,
                         DoFType flag) {


  createMatrix(opr.getVertexToVertexOpr(), src.getVertexDoFFunction(), dst.getVertexDoFFunction(), mat, level, flag);
  EdgeDoFToVertexDoF::createMatrix(opr.getEdgeToVertexOpr(), src.getEdgeDoFFunction(), dst.getVertexDoFFunction(), mat, level, flag);
  VertexDoFToEdgeDoF::createMatrix(opr.getVertexToEdgeOpr(), src.getVertexDoFFunction(), dst.getEdgeDoFFunction(), mat, level, flag);
  edgedof::createMatrix(opr.getEdgeToEdgeOpr(), src.getEdgeDoFFunction(), dst.getEdgeDoFFunction(), mat, level, flag);

}

}
}
