#pragma once

#include <unordered_map>

#include "tinyhhg_core/StencilDirections.hpp"

namespace P1Elements {

// Fenics P1 DoF ordering
// 2       1---0
// |\       \  |
// | \       \ |
// |  \       \|
// 0---1       2

const uint_t ElementSize = 3;

typedef stencilDirection SD;
typedef std::array<SD, ElementSize> P1Element;
typedef std::array<uint_t, ElementSize> DoFMap;
typedef std::array<uint_t, ElementSize> StencilMap;

namespace FaceVertexDoF {

const P1Element elementSW = {{SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S}};
const P1Element elementS = {{SD::VERTEX_C, SD::VERTEX_S, SD::VERTEX_SE}};
const P1Element elementSE = {{SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E}};
const P1Element elementNE = {{SD::VERTEX_C, SD::VERTEX_E, SD::VERTEX_N}};
const P1Element elementN = {{SD::VERTEX_C, SD::VERTEX_N, SD::VERTEX_NW}};
const P1Element elementNW = {{SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_W}};

constexpr inline uint_t stencilMap_(const stencilDirection dir) {
  switch (dir) {
    case SD::VERTEX_S:
      return 0;
    case SD::VERTEX_SE:
      return 1;
    case SD::VERTEX_W:
      return 2;
    case SD::VERTEX_C:
      return 3;
    case SD::VERTEX_E:
      return 4;
    case SD::VERTEX_NW:
      return 5;
    case SD::VERTEX_N:
      return 6;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

static const std::array<P1Element, 3> P1GrayElements =
    {{
         elementS,
         elementNE,
         elementNW
     }};

static const std::array<P1Element, 3> P1BlueElements =
    {{
         elementSW,
         elementSE,
         elementN
     }};

static const std::array<StencilMap, 3> P1GrayStencilMaps =
    {{
         {{stencilMap_(elementS[0]), stencilMap_(elementS[1]), stencilMap_(elementS[2])}},
         {{stencilMap_(elementNE[0]), stencilMap_(elementNE[1]), stencilMap_(elementNE[2])}},
         {{stencilMap_(elementNW[0]), stencilMap_(elementNW[1]), stencilMap_(elementNW[2])}}
     }};

static const std::array<StencilMap, 3> P1BlueStencilMaps =
    {{
         {{stencilMap_(elementSW[0]), stencilMap_(elementSW[1]), stencilMap_(elementSW[2])}},
         {{stencilMap_(elementSE[0]), stencilMap_(elementSE[1]), stencilMap_(elementSE[2])}},
         {{stencilMap_(elementN[0]), stencilMap_(elementN[1]), stencilMap_(elementN[2])}}
     }};

static const std::array<DoFMap, 3> P1GrayDoFMaps =
    {{
         {{2, 0, 1}},
         {{0, 1, 2}},
         {{1, 2, 0}}
     }};


static const std::array<DoFMap, 3> P1BlueDoFMaps =
    {{
         {{0, 1, 2}},
         {{1, 2, 0}},
         {{2, 0, 1}}
     }};
}

template<typename StencilMemory>
void assembleP1LocalStencil(const StencilMap& stencilMap, const DoFMap& dofMap, const Matrix3r& localMatrix, StencilMemory& stencil) {
  for (uint_t j = 0; j < 3; ++j) {
    stencil[stencilMap[j]] += localMatrix(dofMap[0], dofMap[j]);
  }
}

}