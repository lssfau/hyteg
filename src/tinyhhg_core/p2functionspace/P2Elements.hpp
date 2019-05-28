#pragma once

#include <unordered_map>

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

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

template<typename StencilMemory, uint_t M, uint_t N>
inline void assembleStencil(const Matrixr<M,N> &grayMatrix, const Matrixr<M,N> &blueMatrix, StencilMemory &stencil) {
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
         {{edgedof::stencilIndexFromVertex(elementS[3]), edgedof::stencilIndexFromVertex(elementS[4]), edgedof::stencilIndexFromVertex(elementS[5])}},
         {{edgedof::stencilIndexFromVertex(elementNE[3]), edgedof::stencilIndexFromVertex(elementNE[4]), edgedof::stencilIndexFromVertex(elementNE[5])}},
         {{edgedof::stencilIndexFromVertex(elementNW[3]), edgedof::stencilIndexFromVertex(elementNW[4]), edgedof::stencilIndexFromVertex(elementNW[5])}}
     }};

static const std::array<StencilMap, 3> P2BlueStencilMaps =
    {{
         {{edgedof::stencilIndexFromVertex(elementSW[3]), edgedof::stencilIndexFromVertex(elementSW[4]), edgedof::stencilIndexFromVertex(elementSW[5])}},
         {{edgedof::stencilIndexFromVertex(elementSE[3]), edgedof::stencilIndexFromVertex(elementSE[4]), edgedof::stencilIndexFromVertex(elementSE[5])}},
         {{edgedof::stencilIndexFromVertex(elementN[3]), edgedof::stencilIndexFromVertex(elementN[4]), edgedof::stencilIndexFromVertex(elementN[5])}}
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

template<typename StencilMemory, uint_t M, uint_t N>
inline void assembleStencil(const Matrixr<M,N> &grayMatrix, const Matrixr<M,N> &blueMatrix, StencilMemory &stencil) {

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

template<typename StencilMemory, uint_t M, uint_t N>
inline void assembleStencil(const Matrixr<M,N> &grayMatrix, const Matrixr<M,N> &blueMatrix, StencilMemory &stencil) {

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
  stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_HO_C)] = grayMatrix(5, 5) + blueMatrix(5, 5);
  stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_DI_S)] = blueMatrix(5, 3);
  stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_VE_SE)] = blueMatrix(5, 4);
  stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_DI_N)] = grayMatrix(5, 3);
  stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_VE_NW)] = grayMatrix(5, 4);

  // Diagonal
  stencil[edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_DI_C)] = grayMatrix(3, 3) + blueMatrix(3, 3);
  stencil[edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_HO_S)] = grayMatrix(3, 5);
  stencil[edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_VE_E)] = blueMatrix(3, 4);
  stencil[edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_HO_N)] = blueMatrix(3, 5);
  stencil[edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_VE_W)] = grayMatrix(3, 4);

  // Vertical
  stencil[edgedof::stencilIndexFromVerticalEdge(SD::EDGE_VE_C)] = grayMatrix(4, 4) + blueMatrix(4, 4);
  stencil[edgedof::stencilIndexFromVerticalEdge(SD::EDGE_HO_SE)] = grayMatrix(4, 5);
  stencil[edgedof::stencilIndexFromVerticalEdge(SD::EDGE_DI_E)] = grayMatrix(4, 3);
  stencil[edgedof::stencilIndexFromVerticalEdge(SD::EDGE_HO_NW)] = blueMatrix(4, 5);
  stencil[edgedof::stencilIndexFromVerticalEdge(SD::EDGE_DI_W)] = blueMatrix(4, 3);

}

} // EdgeToEdge

} // P2Face

namespace P2Edge {

inline void fillEdgeDoFMap(uint_t start_id, uint_t end_id, uint_t opposite_id, uint_t& edge_start_end_id, uint_t& edge_start_opposite_id, uint_t& edge_end_opposite_id) {
  if (start_id == 0 && end_id == 1 && opposite_id == 2) {
    edge_start_end_id = 5;
    edge_start_opposite_id = 4;
    edge_end_opposite_id = 3;
  }

  if (start_id == 0 && end_id == 2 && opposite_id == 1) {
    edge_start_end_id = 4;
    edge_start_opposite_id = 5;
    edge_end_opposite_id = 3;
  }

  if (start_id == 1 && end_id == 2 && opposite_id == 0) {
    edge_start_end_id = 3;
    edge_start_opposite_id = 5;
    edge_end_opposite_id = 4;
  }

  if (start_id == 1 && end_id == 0 && opposite_id == 2) {
    edge_start_end_id = 5;
    edge_start_opposite_id = 3;
    edge_end_opposite_id = 4;
  }

  if (start_id == 2 && end_id == 0 && opposite_id == 1) {
    edge_start_end_id = 4;
    edge_start_opposite_id = 3;
    edge_end_opposite_id = 5;
  }

  if (start_id == 2 && end_id == 1 && opposite_id == 0) {
    edge_start_end_id = 3;
    edge_start_opposite_id = 4;
    edge_end_opposite_id = 5;
  }
}

namespace VertexToVertex {

typedef std::array<uint_t, 3> DoFMap;
typedef std::array<uint_t, 3> StencilMap;

inline StencilMap convertStencilDirectionsToIndices(const P2Element& element)
{
  return {{ vertexdof::stencilIndexFromVertex( element[0] ), vertexdof::stencilIndexFromVertex( element[1] ), vertexdof::stencilIndexFromVertex( element[2] ) }};
}

template<typename StencilMemory, uint_t M, uint_t N>
inline void assembleStencil(const Edge& edge, const Face& face, const Matrixr<M,N> &grayMatrix, const Matrixr<M,N> &blueMatrix,
                            StencilMemory &stencil, bool south) {

  uint_t start_id = face.vertex_index(edge.neighborVertices()[0]);
  uint_t end_id = face.vertex_index(edge.neighborVertices()[1]);
  uint_t opposite_id = face.vertex_index(face.get_vertex_opposite_to_edge(edge.getID()));

  DoFMap dofMap;
  StencilMap stencilMap;

  if (south) {

    dofMap = DoFMap({{end_id, start_id, opposite_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementSW );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += grayMatrix(dofMap[0], dofMap[j]);
    }

    dofMap = DoFMap({{opposite_id, end_id, start_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementS );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += blueMatrix(dofMap[0], dofMap[j]);
    }

    dofMap = DoFMap({{start_id, opposite_id, end_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementSE );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += grayMatrix(dofMap[0], dofMap[j]);
    }

  } else {

    dofMap = DoFMap({{start_id, end_id, opposite_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementNE );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += grayMatrix(dofMap[0], dofMap[j]);
    }

    dofMap = DoFMap({{opposite_id, start_id, end_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementN );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += blueMatrix(dofMap[0], dofMap[j]);
    }

    dofMap = DoFMap({{end_id, opposite_id, start_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementNW );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += grayMatrix(dofMap[0], dofMap[j]);
    }

  }

}

} // VertexToVertex

namespace EdgeToVertex {

typedef std::array<uint_t, 4> DoFMap;
typedef std::array<uint_t, 3> StencilMap;

inline StencilMap convertStencilDirectionsToIndices(const P2Element& element)
{
  return {{ edgedof::stencilIndexFromVertex( element[3] ), edgedof::stencilIndexFromVertex( element[4] ), edgedof::stencilIndexFromVertex( element[5] ) }};
}

template<typename StencilMemory, uint_t M, uint_t N>
inline void assembleStencil(const Edge& edge, const Face& face, const Matrixr<M,N> &grayMatrix, const Matrixr<M,N> &blueMatrix,
                            StencilMemory &stencil, bool south) {

  uint_t start_id = face.vertex_index(edge.neighborVertices()[0]);
  uint_t end_id = face.vertex_index(edge.neighborVertices()[1]);
  uint_t opposite_id = face.vertex_index(face.get_vertex_opposite_to_edge(edge.getID()));

  uint_t edge_start_end_id;
  uint_t edge_start_opposite_id;
  uint_t edge_end_opposite_id;

  fillEdgeDoFMap(start_id, end_id, opposite_id, edge_start_end_id, edge_start_opposite_id, edge_end_opposite_id);

  DoFMap dofMap;
  StencilMap stencilMap;

  if (south) {

    dofMap = DoFMap({{end_id, edge_start_end_id, edge_start_opposite_id, edge_end_opposite_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementSW );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += grayMatrix(dofMap[0], dofMap[j+1]);
    }

    dofMap = DoFMap({{opposite_id, edge_end_opposite_id, edge_start_end_id, edge_start_opposite_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementS );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += blueMatrix(dofMap[0], dofMap[j+1]);
    }

    dofMap = DoFMap({{start_id, edge_start_opposite_id, edge_end_opposite_id, edge_start_end_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementSE );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += grayMatrix(dofMap[0], dofMap[j+1]);
    }

  } else {

    dofMap = DoFMap({{start_id, edge_start_end_id, edge_end_opposite_id, edge_start_opposite_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementNE );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += grayMatrix(dofMap[0], dofMap[j+1]);
    }

    dofMap = DoFMap({{opposite_id, edge_start_opposite_id, edge_start_end_id, edge_end_opposite_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementN );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += blueMatrix(dofMap[0], dofMap[j+1]);
    }

    dofMap = DoFMap({{end_id, edge_end_opposite_id, edge_start_opposite_id, edge_start_end_id}});
    stencilMap = convertStencilDirectionsToIndices( P2Face::elementNW );
    for (uint_t j = 0; j < 3; ++j) {
      stencil[stencilMap[j]] += grayMatrix(dofMap[0], dofMap[j+1]);
    }

  }

}

} // EdgeToVertex

namespace VertexToEdge {

template<typename StencilMemory, uint_t M, uint_t N>
inline void assembleStencil(const Edge& edge, const Face& face, const Matrixr<M,N> &grayMatrix, const Matrixr<M,N> &blueMatrix,
                            StencilMemory &stencil, bool south) {

  uint_t start_id = face.vertex_index(edge.neighborVertices()[0]);
  uint_t end_id = face.vertex_index(edge.neighborVertices()[1]);
  uint_t opposite_id = face.vertex_index(face.get_vertex_opposite_to_edge(edge.getID()));

  uint_t edge_start_end_id;
  uint_t edge_start_opposite_id;
  uint_t edge_end_opposite_id;

  fillEdgeDoFMap(start_id, end_id, opposite_id, edge_start_end_id, edge_start_opposite_id, edge_end_opposite_id);

  stencil[vertexdof::stencilIndexFromHorizontalEdge(SD::VERTEX_W)] += grayMatrix(edge_start_end_id, start_id);
  stencil[vertexdof::stencilIndexFromHorizontalEdge(SD::VERTEX_E)] += grayMatrix(edge_start_end_id, end_id);

  if (south) {
    stencil[vertexdof::stencilIndexFromHorizontalEdge(SD::VERTEX_SE)] += grayMatrix(edge_start_end_id, opposite_id);

  } else {
    stencil[vertexdof::stencilIndexFromHorizontalEdge(SD::VERTEX_NW)] += grayMatrix(edge_start_end_id, opposite_id);
  }

}

} // VertexToEdge

namespace EdgeToEdge {

template<typename StencilMemory, uint_t M, uint_t N>
inline void assembleStencil(const Edge& edge, const Face& face, const Matrixr<M,N> &grayMatrix, const Matrixr<M,N> &blueMatrix,
                            StencilMemory &stencil, bool south) {

  uint_t start_id = face.vertex_index(edge.neighborVertices()[0]);
  uint_t end_id = face.vertex_index(edge.neighborVertices()[1]);
  uint_t opposite_id = face.vertex_index(face.get_vertex_opposite_to_edge(edge.getID()));

  uint_t edge_start_end_id;
  uint_t edge_start_opposite_id;
  uint_t edge_end_opposite_id;

  fillEdgeDoFMap(start_id, end_id, opposite_id, edge_start_end_id, edge_start_opposite_id, edge_end_opposite_id);

  stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_HO_C)] += grayMatrix(edge_start_end_id, edge_start_end_id);

  if (south) {
    stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_DI_S)] = grayMatrix(edge_start_end_id, edge_start_opposite_id);
    stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_VE_SE)] = grayMatrix(edge_start_end_id, edge_end_opposite_id);
  } else {
    stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_DI_N)] = grayMatrix(edge_start_end_id, edge_end_opposite_id);
    stencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_VE_NW)] = grayMatrix(edge_start_end_id, edge_start_opposite_id);
  }

}

} // EdgeToEdge

} // P2Edge

namespace P2Vertex {

namespace VertexToVertex {

template<typename StencilMemory, uint_t M, uint_t N>
inline void assembleStencil(const Vertex& vertex, const Face& face, const Matrixr<M,N> &grayMatrix,
                            StencilMemory &stencil, const std::shared_ptr< PrimitiveStorage > & storage) {

  uint_t v_i = face.vertex_index(vertex.getID());

  std::vector<PrimitiveID> adj_edges = face.adjacent_edges(vertex.getID());

  std::array<uint_t, 3> stencilMap;
  stencilMap[0] = 0;

  std::array<uint_t, 3> dofMap;
  dofMap[0] = v_i;

  // iterate over adjacent edges
  for (uint_t i = 0; i < adj_edges.size(); ++i)
  {
    uint_t edge_idx = vertex.edge_index(adj_edges[i]) + 1;
    Edge* edge = storage->getEdge(adj_edges[i]);
    PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

    stencilMap[i+1] = edge_idx;
    dofMap[i+1] = face.vertex_index(vertex_j);
  }

  for (uint_t j = 0; j < 3; ++j) {
    stencil[stencilMap[j]] += grayMatrix(dofMap[0], dofMap[j]);
  }
}

} // VertexToVertex

namespace EdgeToVertex {

template<typename StencilMemory, uint_t M, uint_t N>
inline void assembleStencil(const Vertex& vertex, const Face& face, const Matrixr<M,N> &grayMatrix,
                            StencilMemory &stencil, const std::shared_ptr< PrimitiveStorage > & storage) {

  uint_t v_i = face.vertex_index(vertex.getID());

  std::vector<PrimitiveID> adj_edges = face.adjacent_edges(vertex.getID());

  std::array<uint_t, 3> stencilMap;

  std::array<uint_t, 3> vertexDofs;
  vertexDofs[0] = v_i;

  // iterate over adjacent edges
  for (uint_t i = 0; i < adj_edges.size(); ++i)
  {
    uint_t edge_stencil_idx = vertex.edge_index(adj_edges[i]);
    Edge* edge = storage->getEdge(adj_edges[i]);
    PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

    stencilMap[i] = edge_stencil_idx;
    vertexDofs[i+1] = face.vertex_index(vertex_j);
  }

  uint_t face_stencil_idx = vertex.getNumNeighborEdges() + vertex.face_index(face.getID());
  stencilMap[2] = face_stencil_idx;

  std::array<uint_t, 3> dofMap;
  P2Edge::fillEdgeDoFMap(vertexDofs[0], vertexDofs[1], vertexDofs[2], dofMap[0], dofMap[1], dofMap[2]);

  for (uint_t j = 0; j < 3; ++j) {
    stencil[stencilMap[j]] += grayMatrix(v_i, dofMap[j]);
  }
}

} // EdgeToVertex

} // P2Vertex

} // P2Elements

} // hhg
