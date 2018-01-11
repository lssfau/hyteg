#pragma once

#include <unordered_map>

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"

namespace hhg {
namespace P2Elements {

// Fenics P2 DoF ordering
// 2         1---5--0
// | \        \     |
// |  \        \    |
// 4   3        3   4
// |    \        \  |
// |     \        \ |
// 0--5---1         2

const uint_t ElementSize = 6;

typedef stencilDirection SD;
typedef std::array<SD, ElementSize> P2Element;

namespace P2Face {

// See absolute_indexing/edge.pdf for naming of the DoFs
const P2Element elementSW = {{SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S, SD::EDGE_HO_W, SD::EDGE_DI_SW, SD::EDGE_VE_S}};
const P2Element elementS = {{SD::VERTEX_C, SD::VERTEX_S, SD::VERTEX_SE, SD::EDGE_VE_S, SD::EDGE_HO_SE, SD::EDGE_DI_SE}};
const P2Element elementSE = {{SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E, SD::EDGE_DI_SE, SD::EDGE_VE_SE, SD::EDGE_HO_E}};
const P2Element elementNE = {{SD::VERTEX_C, SD::VERTEX_E, SD::VERTEX_N, SD::EDGE_HO_E, SD::EDGE_DI_NE, SD::EDGE_VE_N}};
const P2Element elementN = {{SD::VERTEX_C, SD::VERTEX_N, SD::VERTEX_NW, SD::EDGE_VE_N, SD::EDGE_HO_NW, SD::EDGE_DI_NW}};
const P2Element elementNW = {{SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_W, SD::EDGE_DI_NW, SD::EDGE_VE_NW, SD::EDGE_HO_W}};


static const std::array<P2Element, 3> P2GrayElements =
    {{
         elementS,
         elementNE,
         elementNW
     }};

static const std::array<P2Element, 3> P2BlueElements =
    {{
         elementSW,
         elementSE,
         elementN
     }};

namespace VertexToVertex {

typedef std::array<uint_t, 3> DoFMap;
typedef std::array<uint_t, 3> StencilMap;

static const std::array<StencilMap, 3> P2GrayStencilMaps =
    {{
         {{vertexdof::stencilIndexFromVertex(elementS[0]), vertexdof::stencilIndexFromVertex(elementS[1]), vertexdof::stencilIndexFromVertex(elementS[2])}},
         {{vertexdof::stencilIndexFromVertex(elementNE[0]), vertexdof::stencilIndexFromVertex(elementNE[1]), vertexdof::stencilIndexFromVertex(elementNE[2])}},
         {{vertexdof::stencilIndexFromVertex(elementNW[0]), vertexdof::stencilIndexFromVertex(elementNW[1]), vertexdof::stencilIndexFromVertex(elementNW[2])}}
     }};

static const std::array<StencilMap, 3> P2BlueStencilMaps =
    {{
         {{vertexdof::stencilIndexFromVertex(elementSW[0]), vertexdof::stencilIndexFromVertex(elementSW[1]), vertexdof::stencilIndexFromVertex(elementSW[2])}},
         {{vertexdof::stencilIndexFromVertex(elementSE[0]), vertexdof::stencilIndexFromVertex(elementSE[1]), vertexdof::stencilIndexFromVertex(elementSE[2])}},
         {{vertexdof::stencilIndexFromVertex(elementN[0]), vertexdof::stencilIndexFromVertex(elementN[1]), vertexdof::stencilIndexFromVertex(elementN[2])}}
     }};

static const std::array<DoFMap, 3> P2GrayDoFMaps =
    {{
         {{2, 0, 1}},
         {{0, 1, 2}},
         {{1, 2, 0}}
     }};


static const std::array<DoFMap, 3> P2BlueDoFMaps =
    {{
         {{0, 1, 2}},
         {{1, 2, 0}},
         {{2, 0, 1}}
     }};

template<typename StencilMemory>
inline void assembleStencil(const StencilMap &stencilMap, const DoFMap &dofMap, const Matrix6r &localMatrix,
                              StencilMemory &stencil) {
  for (uint_t j = 0; j < 3; ++j) {
    stencil[stencilMap[j]] += localMatrix(dofMap[0], dofMap[j]);
  }
}

} // VertexToVertex

} // P2Face

} // P2Elements

} // hhg