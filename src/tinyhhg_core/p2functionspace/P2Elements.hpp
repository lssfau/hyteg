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
inline void assembleStencil(const Matrix6r &grayMatrix, const Matrix6r &blueMatrix, StencilMemory &stencil) {
  for (uint_t i = 0; i < P2Face::P2GrayElements.size(); ++i) {
    for (uint_t j = 0; j < 3; ++j) {
      stencil[P2GrayStencilMaps[i][j]] += grayMatrix(P2GrayDoFMaps[i][0], P2GrayDoFMaps[i][j]);
    }
  }

  for (uint_t i = 0; i < P2Face::P2BlueElements.size(); ++i) {
    for (uint_t j = 0; j < 3; ++j) {
      stencil[P2BlueStencilMaps[i][j]] += blueMatrix(P2BlueDoFMaps[i][0], P2BlueDoFMaps[i][j]);
    }
  }
}

} // VertexToVertex

namespace EdgeToVertex {

typedef std::array<uint_t, 4> DoFMap;
typedef std::array<uint_t, 3> StencilMap;

static const std::array<StencilMap, 3> P2GrayStencilMaps =
    {{
         {{indexing::edgedof::stencilIndexFromVertex(elementS[3]), indexing::edgedof::stencilIndexFromVertex(elementS[4]), indexing::edgedof::stencilIndexFromVertex(elementS[5])}},
         {{indexing::edgedof::stencilIndexFromVertex(elementNE[3]), indexing::edgedof::stencilIndexFromVertex(elementNE[4]), indexing::edgedof::stencilIndexFromVertex(elementNE[5])}},
         {{indexing::edgedof::stencilIndexFromVertex(elementNW[3]), indexing::edgedof::stencilIndexFromVertex(elementNW[4]), indexing::edgedof::stencilIndexFromVertex(elementNW[5])}}
     }};

static const std::array<StencilMap, 3> P2BlueStencilMaps =
    {{
         {{indexing::edgedof::stencilIndexFromVertex(elementSW[3]), indexing::edgedof::stencilIndexFromVertex(elementSW[4]), indexing::edgedof::stencilIndexFromVertex(elementSW[5])}},
         {{indexing::edgedof::stencilIndexFromVertex(elementSE[3]), indexing::edgedof::stencilIndexFromVertex(elementSE[4]), indexing::edgedof::stencilIndexFromVertex(elementSE[5])}},
         {{indexing::edgedof::stencilIndexFromVertex(elementN[3]), indexing::edgedof::stencilIndexFromVertex(elementN[4]), indexing::edgedof::stencilIndexFromVertex(elementN[5])}}
     }};

// First DoF is center vertex DoF
static const std::array<DoFMap, 3> P2GrayDoFMaps =
    {{
         {{2, 4, 5, 3}},
         {{0, 5, 3, 4}},
         {{1, 3, 4, 5}}
     }};


static const std::array<DoFMap, 3> P2BlueDoFMaps =
    {{
         {{0, 5, 3, 4}},
         {{1, 3, 4, 5}},
         {{2, 4, 5, 3}}
     }};

template<typename StencilMemory>
inline void assembleStencil(const Matrix6r &grayMatrix, const Matrix6r &blueMatrix, StencilMemory &stencil) {

  for (uint_t i = 0; i < P2Face::P2GrayElements.size(); ++i) {
    for (uint_t j = 0; j < 3; ++j) {
      stencil[P2GrayStencilMaps[i][j]] += grayMatrix(P2GrayDoFMaps[i][0], P2GrayDoFMaps[i][j+1]);
    }
  }

  for (uint_t i = 0; i < P2Face::P2BlueElements.size(); ++i) {
    for (uint_t j = 0; j < 3; ++j) {
      stencil[P2BlueStencilMaps[i][j]] += blueMatrix(P2BlueDoFMaps[i][0], P2BlueDoFMaps[i][j+1]);
    }
  }
}

} // EdgeToVertex

namespace VertexToEdge {

template<typename StencilMemory>
inline void assembleStencil(const Matrix6r &grayMatrix, const Matrix6r &blueMatrix, StencilMemory &stencil) {

  // Horizontal
  stencil[vertexdof::stencilIndexFromHorizontalEdge(SD::VERTEX_W)] = grayMatrix(5, 0) + blueMatrix(5, 1);
  stencil[vertexdof::stencilIndexFromHorizontalEdge(SD::VERTEX_E)] = grayMatrix(5, 1) + blueMatrix(5, 0);
  stencil[vertexdof::stencilIndexFromHorizontalEdge(SD::VERTEX_SE)] = blueMatrix(5, 2);
  stencil[vertexdof::stencilIndexFromHorizontalEdge(SD::VERTEX_NW)] = grayMatrix(5, 2);

  // Diagonal
  stencil[vertexdof::stencilIndexFromDiagonalEdge(SD::VERTEX_SE)] = grayMatrix(3, 1) + blueMatrix(3, 2);
  stencil[vertexdof::stencilIndexFromDiagonalEdge(SD::VERTEX_NE)] = blueMatrix(3, 0);
  stencil[vertexdof::stencilIndexFromDiagonalEdge(SD::VERTEX_NW)] = grayMatrix(3, 2) + blueMatrix(3, 1);
  stencil[vertexdof::stencilIndexFromDiagonalEdge(SD::VERTEX_SW)] = grayMatrix(3, 0);

  // Vertical
  stencil[vertexdof::stencilIndexFromVerticalEdge(SD::VERTEX_S)] = grayMatrix(4, 0) + blueMatrix(4, 2);
  stencil[vertexdof::stencilIndexFromVerticalEdge(SD::VERTEX_SE)] = grayMatrix(4, 1);
  stencil[vertexdof::stencilIndexFromVerticalEdge(SD::VERTEX_N)] = grayMatrix(4, 2) + blueMatrix(4, 0);
  stencil[vertexdof::stencilIndexFromVerticalEdge(SD::VERTEX_NW)] = blueMatrix(4, 1);

}

} // VertexToEdge

namespace EdgeToEdge {

template<typename StencilMemory>
inline void assembleStencil(const Matrix6r &grayMatrix, const Matrix6r &blueMatrix, StencilMemory &stencil) {

  // Horizontal
  stencil[indexing::edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_HO_C)] = grayMatrix(5, 5) + blueMatrix(5, 5);
  stencil[indexing::edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_DI_S)] = blueMatrix(5, 3);
  stencil[indexing::edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_VE_SE)] = blueMatrix(5, 4);
  stencil[indexing::edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_DI_N)] = grayMatrix(5, 3);
  stencil[indexing::edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_VE_NW)] = grayMatrix(5, 4);

  // Diagonal
  stencil[indexing::edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_DI_C)] = grayMatrix(3, 3) + blueMatrix(3, 3);
  stencil[indexing::edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_HO_S)] = grayMatrix(3, 5);
  stencil[indexing::edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_VE_E)] = blueMatrix(3, 4);
  stencil[indexing::edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_HO_N)] = blueMatrix(3, 5);
  stencil[indexing::edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_VE_W)] = grayMatrix(3, 4);

  // Vertical
  stencil[indexing::edgedof::stencilIndexFromVerticalEdge(SD::EDGE_VE_C)] = grayMatrix(4, 4) + blueMatrix(4, 4);
  stencil[indexing::edgedof::stencilIndexFromVerticalEdge(SD::EDGE_HO_SE)] = grayMatrix(4, 5);
  stencil[indexing::edgedof::stencilIndexFromVerticalEdge(SD::EDGE_DI_E)] = grayMatrix(4, 3);
  stencil[indexing::edgedof::stencilIndexFromVerticalEdge(SD::EDGE_HO_NW)] = blueMatrix(4, 5);
  stencil[indexing::edgedof::stencilIndexFromVerticalEdge(SD::EDGE_DI_W)] = blueMatrix(4, 3);

}

} // EdgeToEdge

} // P2Face

} // P2Elements

} // hhg