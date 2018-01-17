#pragma once

#include <tinyhhg_core/p2functionspace/P2Function.hpp>

#include <tinyhhg_core/p1functionspace/P1Petsc.hpp>
#include <tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFPetsc.hpp>
#include <tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFPetsc.hpp>

namespace hhg {
namespace petsc {

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
//  createMatrix(opr.getEdgeToEdgeOpr(), *src.getEdgeDoFFunction(), *dst.getEdgeDoFFunction(), mat, level, flag);

}

}
}
