#pragma once

namespace hhg {
namespace petsc {

template<class OperatorType>
inline void createMatrix(OperatorType& opr, P1Function< PetscInt > & src, BubbleFunction< PetscInt > & dst, Mat& mat, uint_t level, DoFType flag)
{
  for (auto& it : opr.getStorage()->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      P1ToBubbleFace::saveOperator(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
    }
  }
}

}
}