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
                       const PrimitiveDataID<FaceBubbleFunctionMemory, Face> &srcId,
                       const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId, UpdateType update) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = face.getData(operatorId)->data[Level];
  auto& src = face.getData(srcId)->data[Level];
  auto& dst = face.getData(dstId)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = 0.0;

      for (auto neighbor : BubbleFace::CoordsVertex::neighbors) {
        tmp += opr_data[neighbor]*src[BubbleFace::CoordsVertex::index<Level>(i, j, neighbor)];
      }

      if (update==Replace) {
        dst[P1Face::CoordsVertex::index<Level>(i, j, P1Face::CoordsVertex::VERTEX_C)] = tmp;
      } else if (update==Add) {
        dst[P1Face::CoordsVertex::index<Level>(i, j, P1Face::CoordsVertex::VERTEX_C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, apply_tmpl, apply)
}
}