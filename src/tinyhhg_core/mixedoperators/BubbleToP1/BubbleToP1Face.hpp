#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "BubbleToP1Memory.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"

#include "tinyhhg_core/p1functionspace/P1FaceIndex.hpp"

namespace hhg {
namespace BubbleToP1Face {
template<size_t Level>
inline void apply_tmpl(Face &face, const PrimitiveDataID<FaceBubbleToP1StencilMemory, Face> &operatorId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &srcId,
                       const PrimitiveDataID<FaceP1FunctionMemory< real_t >, Face> &dstId, UpdateType update) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = face.getData(operatorId)->data[Level];
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  real_t tmp;

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = 0.0;

      for (auto neighbor : BubbleFace::neighbors) {
        tmp += opr_data[BubbleFace::indexFaceStencil(neighbor)]*src[BubbleFace::indexFaceFromVertex<Level>(i, j, neighbor)];
      }

      if (update==Replace) {
        dst[ vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C) ] = tmp;
      } else if (update==Add) {
        dst[ vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C) ] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, apply_tmpl, apply)

#ifdef HHG_BUILD_WITH_PETSC
template<size_t Level>
inline void saveOperator_tmpl(Face &face, const PrimitiveDataID<FaceBubbleToP1StencilMemory, Face> &operatorId,
                              const PrimitiveDataID<FaceBubbleFunctionMemory< PetscInt >, Face> &srcId,
                              const PrimitiveDataID<FaceP1FunctionMemory< PetscInt >, Face> &dstId, Mat& mat) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = face.getData(operatorId)->data[Level];
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {

      PetscInt dst_id = dst[P1Face::FaceCoordsVertex::index<Level>(i, j, stencilDirection::VERTEX_C)];

      for (auto neighbor : BubbleFace::neighbors) {
        MatSetValues(mat, 1, &dst_id, 1, &src[BubbleFace::indexFaceFromVertex<Level>(i, j, neighbor)],
                     &opr_data[BubbleFace::indexFaceStencil(neighbor)], INSERT_VALUES);
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)
#endif

}
}
