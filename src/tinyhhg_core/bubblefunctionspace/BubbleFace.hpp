#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "BubbleMemory.hpp"
#include "BubbleFaceIndex.hpp"

namespace hhg {
namespace BubbleFace {

template< typename ValueType, size_t Level >
inline void assign_tmpl(Face &face,
                        const std::vector<ValueType> &scalars,
                        const std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face>> &srcIds,
                        const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      ValueType tmp = scalars[0] * face.getData(srcIds[0])->getPointer( Level )[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)];

      for (size_t k = 1; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->getPointer( Level )[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)];
      }
      face.getData(dstId)->getPointer( Level )[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] = tmp;
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      ValueType tmp = scalars[0] * face.getData(srcIds[0])->getPointer( Level )[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)];

      for (size_t k = 1; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->getPointer( Level )[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)];
      }
      face.getData(dstId)->getPointer( Level )[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, assign_tmpl, assign)

template< typename ValueType, size_t Level >
inline void add_tmpl(Face &face,
                     const std::vector<ValueType> &scalars,
                     const std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face>> &srcIds,
                     const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      ValueType tmp = 0.0;

      for (size_t k = 0; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->getPointer( Level )[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)];
      }
      face.getData(dstId)->getPointer( Level )[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] += tmp;
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      ValueType tmp = 0.0;

      for (size_t k = 0; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->getPointer( Level )[indexFaceFromBlueFace<Level>(i, j,
                                                                                                                  stencilDirection::CELL_BLUE_C)];
      }
      face.getData(dstId)->getPointer( Level )[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] += tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, add_tmpl, add)

template< typename ValueType, size_t Level >
inline real_t dot_tmpl(Face &face,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &lhsId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &rhsId) {
  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto lhs_data = face.getData(lhsId)->getPointer( Level );
  auto rhs_data = face.getData(rhsId)->getPointer( Level );

  for (size_t i = 0; i < rowsize - 1; ++i) {
    for (size_t j = 0; j < inner_rowsize - 1; ++j) {
      sp += lhs_data[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)]
          *rhs_data[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)];
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i) {
    for (size_t j = 0; j < inner_rowsize - 2; ++j) {
      sp += lhs_data[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)]
          *rhs_data[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)];
    }
    --inner_rowsize;
  }

  return sp;
}

SPECIALIZE_WITH_VALUETYPE(real_t, dot_tmpl, dot)

template< typename ValueType, size_t Level >
inline void apply_tmpl(Face& face, const PrimitiveDataID<FaceBubbleStencilMemory, Face>& operatorId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &srcId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId, UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = face.getData(operatorId)->data[Level];

  auto& face_gray_stencil = opr_data[0];
  auto& face_blue_stencil = opr_data[1];

  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  ValueType tmp;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      tmp = face_gray_stencil[FaceCoordsCellGray::CELL_GRAY_C] * src[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)];

      if (update == Replace) {
        dst[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] = tmp;
      } else if (update == Add) {
        dst[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] += tmp;
      }
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      tmp = face_blue_stencil[FaceCoordsCellBlue::CELL_BLUE_C] * src[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)];

      if (update == Replace) {
        dst[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] = tmp;
      } else if (update == Add) {
        dst[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, apply_tmpl, apply)

template< typename ValueType, size_t Level >
inline void enumerate_tmpl(Face &face, const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId, uint_t& num) {
  using walberla::real_c;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      face.getData(dstId)->getPointer( Level )[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] = real_c(num++);
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      face.getData(dstId)->getPointer( Level )[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] = real_c(num++);
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, enumerate_tmpl, enumerate)

template< typename ValueType, size_t Level >
inline void printFunctionMemory(Face& face, const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId){
  using namespace std;
  ValueType* faceMemory = face.getData(dstId)->getPointer(Level);
  uint_t verticesPerDge = hhg::levelinfo::num_microvertices_per_edge(Level);
  cout << setfill('=') << setw(100) << "" << endl;
  cout << face << std::left << setprecision(1) << fixed << setfill(' ') << endl;
  std::cout << "Cell Blue: " << std::endl;
  for (size_t i = 0; i < verticesPerDge-2; ++i) {
    for (size_t j = 0; j < verticesPerDge-2 - i; ++j) {
      cout << setw(5) << faceMemory[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] << "|";
    }
    std::cout << std::endl;
  }
  cout << "Cell Gray: " << std::endl;
  for (size_t i = 0; i < verticesPerDge-1; ++i) {
    for (size_t j = 0; j < verticesPerDge-1 - i; ++j) {
      cout << setw(5) << faceMemory[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] << "|";
    }
    std::cout << std::endl;
  }
  cout << setw(100) << setfill(' ') << endl;

}

#ifdef HHG_BUILD_WITH_PETSC

template< typename ValueType, uint_t Level >
inline void createVectorFromFunctionTmpl(Face &face,
                                     const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &srcId,
                                     const PrimitiveDataID<FaceBubbleFunctionMemory< PetscInt >, Face> &numeratorId,
                                     Vec& vec) {
  PetscInt dofs = (PetscInt) levelinfo::num_microfaces_per_face(Level);

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );

  VecSetValues(vec, dofs, &numerator[0], &src[0], INSERT_VALUES);
}

SPECIALIZE_WITH_VALUETYPE(void, createVectorFromFunctionTmpl, createVectorFromFunction)

template< typename ValueType, uint_t Level >
inline void createFunctionFromVectorTmpl(Face &face,
                                         const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &srcId,
                                         const PrimitiveDataID<FaceBubbleFunctionMemory< PetscInt >, Face> &numeratorId,
                                         Vec& vec) {
  PetscInt dofs = (PetscInt) levelinfo::num_microfaces_per_face(Level);

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );

  VecGetValues(vec, dofs, &numerator[0], &src[0]);
}

SPECIALIZE_WITH_VALUETYPE(void, createFunctionFromVectorTmpl, createFunctionFromVector)

template< typename ValueType, size_t Level >
inline void saveOperator_tmpl(Face& face, const PrimitiveDataID<FaceBubbleStencilMemory, Face>& operatorId,
                              const PrimitiveDataID<FaceBubbleFunctionMemory< PetscInt >, Face> &srcId,
                              const PrimitiveDataID<FaceBubbleFunctionMemory< PetscInt >, Face> &dstId, Mat& mat)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = face.getData(operatorId)->data[Level];

  auto& face_gray_stencil = opr_data[0];
  auto& face_blue_stencil = opr_data[1];

  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      MatSetValues(mat, 1, &dst[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)], 1, &src[indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)], &face_gray_stencil[FaceCoordsCellGray::CELL_GRAY_C], INSERT_VALUES);
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      MatSetValues(mat,
                   1,
                   &dst[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)],
                   1,
                   &src[indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)],
                   &face_blue_stencil[FaceCoordsCellBlue::CELL_BLUE_C], INSERT_VALUES);
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, saveOperator_tmpl, saveOperator)

#endif


}// namespace P1BubbleFace
}// namespace hhg
