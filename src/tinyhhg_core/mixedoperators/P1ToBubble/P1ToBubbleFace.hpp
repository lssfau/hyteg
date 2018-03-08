#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "P1ToBubbleMemory.hpp"

#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"

namespace hhg {
namespace P1ToBubbleFace {
template<size_t Level>
inline void apply_tmpl(Face &face, const PrimitiveDataID<FaceP1ToBubbleStencilMemory, Face> &operatorId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Face> &srcId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &dstId, UpdateType update) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& face_stencil_stack = face.getData(operatorId)->data[Level];
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  auto& face_gray_stencil = face_stencil_stack[0];
  auto& face_blue_stencil = face_stencil_stack[1];

  real_t tmp;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      tmp = 0.0;

      for ( const auto & neighbor : vertexdof::macroface::neighborsFromGrayFace )
      {
        tmp += face_gray_stencil[ vertexdof::stencilIndexFromGrayFace(neighbor)] * src[vertexdof::macroface::indexFromGrayFace( Level, i, j, neighbor )];
      }

      if (update == Replace) {
        dst[BubbleFace::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )] = tmp;
      } else if (update == Add) {
        dst[BubbleFace::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )] += tmp;
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

      for ( const auto neighbor : vertexdof::macroface::neighborsFromBlueFace )
      {
        tmp += face_blue_stencil[ vertexdof::stencilIndexFromBlueFace(neighbor)] * src[vertexdof::macroface::indexFromBlueFace( Level, i, j, neighbor )];
      }

      if (update == Replace) {
        dst[BubbleFace::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )] = tmp;
      } else if (update == Add) {
        dst[BubbleFace::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, apply_tmpl, apply)


#ifdef HHG_BUILD_WITH_PETSC
template<size_t Level>
inline void saveOperator_tmpl(Face &face, const PrimitiveDataID<FaceP1ToBubbleStencilMemory, Face> &operatorId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &srcId,
                              const PrimitiveDataID<FaceBubbleFunctionMemory< PetscInt >, Face> &dstId, Mat& mat) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& face_stencil_stack = face.getData(operatorId)->data[Level];
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  auto& face_gray_stencil = face_stencil_stack[0];
  auto& face_blue_stencil = face_stencil_stack[1];

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      PetscInt dst_id = dst[BubbleFace::indexFaceFromGrayFace( Level, i, j, stencilDirection ::CELL_GRAY_C)];

      for ( const auto & neighbor : vertexdof::macroface::neighborsFromGrayFace )
      {
        MatSetValues(mat, 1, &dst_id, 1, &src[vertexdof::macroface::indexFromGrayFace(Level, i, j, neighbor)], &face_gray_stencil[vertexdof::stencilIndexFromGrayFace(neighbor)], INSERT_VALUES);
      }
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      PetscInt dst_id = dst[BubbleFace::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C)];

      for ( const auto & neighbor : vertexdof::macroface::neighborsFromBlueFace )
      {
        MatSetValues(mat, 1, &dst_id, 1, &src[vertexdof::macroface::indexFromBlueFace( Level, i, j, neighbor)], &face_blue_stencil[vertexdof::stencilIndexFromBlueFace(neighbor)], INSERT_VALUES);
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)
#endif

}
}
