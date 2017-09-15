#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "P1ToBubbleMemory.hpp"

#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"
#include "tinyhhg_core/p1functionspace/P1FaceIndex.hpp"

namespace hhg {
namespace P1ToBubbleFace {
template<size_t Level>
inline void apply_tmpl(Face &face, const PrimitiveDataID<FaceP1ToBubbleStencilMemory, Face> &operatorId,
                       const PrimitiveDataID<FaceP1FunctionMemory< real_t >, Face> &srcId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &dstId, UpdateType update) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& face_stencil_stack = face.getData(operatorId)->data[Level];
  auto& src = face.getData(srcId)->data[Level];
  auto& dst = face.getData(dstId)->data[Level];

  auto& face_gray_stencil = face_stencil_stack[0];
  auto& face_blue_stencil = face_stencil_stack[1];

  real_t tmp;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      tmp = 0.0;

      for (auto neighbor : P1Face::FaceCoordsCellGray::neighbors)
      {
        tmp += face_gray_stencil[neighbor] * src[P1Face::FaceCoordsCellGray::index<Level>(i, j, neighbor)];
      }

      if (update == Replace) {
        dst[BubbleFace::FaceCoordsCellGray::index<Level>(i, j, BubbleFace::FaceCoordsCellGray::CELL_GRAY_C)] = tmp;
      } else if (update == Add) {
        dst[BubbleFace::FaceCoordsCellGray::index<Level>(i, j, BubbleFace::FaceCoordsCellGray::CELL_GRAY_C)] += tmp;
      }
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      tmp = 0.0;

      for (auto neighbor : P1Face::FaceCoordsCellBlue::neighbors)
      {
        tmp += face_blue_stencil[neighbor] * src[P1Face::FaceCoordsCellBlue::index<Level>(i, j, neighbor)];
      }

      if (update == Replace) {
        dst[BubbleFace::FaceCoordsCellBlue::index<Level>(i, j, BubbleFace::FaceCoordsCellBlue::CELL_BLUE_C)] = tmp;
      } else if (update == Add) {
        dst[BubbleFace::FaceCoordsCellBlue::index<Level>(i, j, BubbleFace::FaceCoordsCellBlue::CELL_BLUE_C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, apply_tmpl, apply)


#ifdef HHG_BUILD_WITH_PETSC
template<size_t Level>
inline void saveOperator_tmpl(Face &face, const PrimitiveDataID<FaceP1ToBubbleStencilMemory, Face> &operatorId,
                              const PrimitiveDataID<FaceP1FunctionMemory< PetscInt >, Face> &srcId,
                              const PrimitiveDataID<FaceBubbleFunctionMemory< PetscInt >, Face> &dstId, Mat& mat) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& face_stencil_stack = face.getData(operatorId)->data[Level];
  auto& src = face.getData(srcId)->data[Level];
  auto& dst = face.getData(dstId)->data[Level];

  auto& face_gray_stencil = face_stencil_stack[0];
  auto& face_blue_stencil = face_stencil_stack[1];

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      PetscInt dst_id = dst[BubbleFace::FaceCoordsCellGray::index<Level>(i, j, BubbleFace::FaceCoordsCellGray::CELL_GRAY_C)];

      for (auto neighbor : P1Face::FaceCoordsCellGray::neighbors)
      {
        MatSetValues(mat, 1, &dst_id, 1, &src[P1Face::FaceCoordsCellGray::index<Level>(i, j, neighbor)], &face_gray_stencil[neighbor], INSERT_VALUES);
      }
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      PetscInt dst_id = dst[BubbleFace::FaceCoordsCellBlue::index<Level>(i, j, BubbleFace::FaceCoordsCellBlue::CELL_BLUE_C)];

      for (auto neighbor : P1Face::FaceCoordsCellBlue::neighbors)
      {
        MatSetValues(mat, 1, &dst_id, 1, &src[P1Face::FaceCoordsCellBlue::index<Level>(i, j, neighbor)], &face_blue_stencil[neighbor], INSERT_VALUES);
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)
#endif

}
}
