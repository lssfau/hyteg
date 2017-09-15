#pragma once

namespace hhg {
namespace petsc {

inline void createVectorFromFunction(BubbleFunction<PetscScalar> &function,
                                     BubbleFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      BubbleFace::createVectorFromFunction<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }
}

inline void createFunctionFromVector(BubbleFunction<PetscScalar> &function,
                                     BubbleFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      BubbleFace::createFunctionFromVector<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }
}

template<class OperatorType>
inline void createMatrix(OperatorType& opr, BubbleFunction< PetscInt > & src, BubbleFunction< PetscInt > & dst, Mat& mat, uint_t level, DoFType flag)
{
  for (auto& it : opr.getStorage()->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      BubbleFace::saveOperator<OperatorType>(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
    }
  }
}

}
}