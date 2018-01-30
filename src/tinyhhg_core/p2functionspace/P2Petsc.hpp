#pragma once

#include <tinyhhg_core/p2functionspace/P2Function.hpp>

#include <tinyhhg_core/p1functionspace/P1Petsc.hpp>
#include <tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFPetsc.hpp>
#include <tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFPetsc.hpp>
#include <tinyhhg_core/edgedofspace/EdgeDoFPetsc.hpp>

namespace hhg {
namespace petsc {

inline void createVectorFromFunction(P2Function<PetscScalar> &function,
                                     P2Function<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createVectorFromFunction(*function.getVertexDoFFunction(), *numerator.getVertexDoFFunction(), vec, level, flag);
  EdgeDoF::createVectorFromFunction(*function.getEdgeDoFFunction(), *numerator.getEdgeDoFFunction(), vec, level, flag);
}

inline void createFunctionFromVector(P2Function<PetscScalar> &function,
                                     P2Function<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createFunctionFromVector(*function.getVertexDoFFunction(), *numerator.getVertexDoFFunction(), vec, level, flag);
  EdgeDoF::createFunctionFromVector(*function.getEdgeDoFFunction(), *numerator.getEdgeDoFFunction(), vec, level, flag);
}

inline void applyDirichletBC(P2Function<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {
  applyDirichletBC(*numerator.getVertexDoFFunction(), mat, level);
  EdgeDoF::applyDirichletBC(*numerator.getEdgeDoFFunction(), mat, level);
}

template<class OperatorType>
inline void createMatrix(OperatorType &opr,
                         P2Function<PetscInt> &src,
                         P2Function<PetscInt> &dst,
                         Mat &mat,
                         uint_t level,
                         DoFType flag) {


  createMatrix(opr.getVertexToVertexOpr(), *src.getVertexDoFFunction(), *dst.getVertexDoFFunction(), mat, level, flag);
  EdgeDoFToVertexDoF::createMatrix(opr.getEdgeToVertexOpr(), *src.getEdgeDoFFunction(), *dst.getVertexDoFFunction(), mat, level, flag);
  VertexDoFToEdgeDoF::createMatrix(opr.getVertexToEdgeOpr(), *src.getVertexDoFFunction(), *dst.getEdgeDoFFunction(), mat, level, flag);
  EdgeDoF::createMatrix(opr.getEdgeToEdgeOpr(), *src.getEdgeDoFFunction(), *dst.getEdgeDoFFunction(), mat, level, flag);

}

}
}
