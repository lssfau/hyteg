#pragma once

#include <unordered_map>

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/types/matrix.hpp"

namespace hhg {
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

// ordered
const P1Element elementSWOrd = {{SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S}};
const P1Element elementSOrd = {{SD::VERTEX_S, SD::VERTEX_SE, SD::VERTEX_C}};
const P1Element elementSEOrd = {{SD::VERTEX_E, SD::VERTEX_C, SD::VERTEX_SE}};
const P1Element elementNEOrd = {{SD::VERTEX_C, SD::VERTEX_E, SD::VERTEX_N}};
const P1Element elementNOrd = {{SD::VERTEX_N, SD::VERTEX_NW, SD::VERTEX_C}};
const P1Element elementNWOrd = {{SD::VERTEX_W, SD::VERTEX_C, SD::VERTEX_NW}};

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
         {{vertexdof::stencilIndexFromVertex(elementS[0]), vertexdof::stencilIndexFromVertex(elementS[1]), vertexdof::stencilIndexFromVertex(elementS[2])}},
         {{vertexdof::stencilIndexFromVertex(elementNE[0]), vertexdof::stencilIndexFromVertex(elementNE[1]), vertexdof::stencilIndexFromVertex(elementNE[2])}},
         {{vertexdof::stencilIndexFromVertex(elementNW[0]), vertexdof::stencilIndexFromVertex(elementNW[1]), vertexdof::stencilIndexFromVertex(elementNW[2])}}
     }};

static const std::array<StencilMap, 3> P1BlueStencilMaps =
    {{
         {{vertexdof::stencilIndexFromVertex(elementSW[0]), vertexdof::stencilIndexFromVertex(elementSW[1]), vertexdof::stencilIndexFromVertex(elementSW[2])}},
         {{vertexdof::stencilIndexFromVertex(elementSE[0]), vertexdof::stencilIndexFromVertex(elementSE[1]), vertexdof::stencilIndexFromVertex(elementSE[2])}},
         {{vertexdof::stencilIndexFromVertex(elementN[0]), vertexdof::stencilIndexFromVertex(elementN[1]), vertexdof::stencilIndexFromVertex(elementN[2])}}
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

inline StencilMap convertStencilDirectionsToIndices( const P1Element & element )
{
  return {{ vertexdof::stencilIndexFromVertex( element[0] ), vertexdof::stencilIndexFromVertex( element[1] ), vertexdof::stencilIndexFromVertex( element[2] ) }};
}

template<typename StencilMemory>
inline void assembleP1LocalStencil(const StencilMap &stencilMap, const DoFMap &dofMap, const Matrix3r &localMatrix,
                            StencilMemory &stencil, double coeffWeight = 1.0) {
  for (uint_t j = 0; j < 3; ++j) {
    stencil[stencilMap[j]] += coeffWeight * localMatrix(dofMap[0], dofMap[j]);
  }
}

}
}
