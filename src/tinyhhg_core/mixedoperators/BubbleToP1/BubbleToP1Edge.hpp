#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "BubbleToP1Memory.hpp"
#include "tinyhhg_core/p1functionspace/P1EdgeIndex.hpp"

namespace hhg {
namespace BubbleToP1Edge {
template<size_t Level>
inline void apply_tmpl(Edge &edge, const PrimitiveDataID<EdgeBubbleToP1StencilMemory, Edge> &operatorId,
                       const PrimitiveDataID<EdgeBubbleFunctionMemory< real_t >, Edge> &srcId,
                       const PrimitiveDataID<EdgeP1FunctionMemory< real_t >, Edge> &dstId, UpdateType update) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &edge_vertex_stencil = edge.getData(operatorId)->data[Level];
  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );

  real_t tmp;

  for (size_t i = 1; i < rowsize - 1; ++i) {
    tmp = 0.0;

    for (auto neighbor : BubbleEdge::EdgeCoordsVertex::neighbors_south) {
      tmp += edge_vertex_stencil[neighbor]*src[BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (auto neighbor : BubbleEdge::EdgeCoordsVertex::neighbors_north) {
        tmp += edge_vertex_stencil[neighbor]*src[BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)];
      }
    }

    if (update == Replace) {
      dst[P1Edge::EdgeCoordsVertex::index<Level>(i, P1Edge::EdgeCoordsVertex::VERTEX_C)] = tmp;
    } else if (update == Add) {
      dst[P1Edge::EdgeCoordsVertex::index<Level>(i, P1Edge::EdgeCoordsVertex::VERTEX_C)] += tmp;
    }
  }
}

SPECIALIZE(void, apply_tmpl, apply)

template<size_t Level>
inline void saveOperator_tmpl(Edge &edge, const PrimitiveDataID<EdgeBubbleToP1StencilMemory, Edge> &operatorId,
                              const PrimitiveDataID<EdgeBubbleFunctionMemory< real_t >, Edge> &srcId,
                              const PrimitiveDataID<EdgeP1FunctionMemory< real_t >, Edge> &dstId, std::ostream& out) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &edge_vertex_stencil = edge.getData(operatorId)->data[Level];
  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );

  for (size_t i = 1; i < rowsize - 1; ++i) {

    uint_t dst_id = dst[P1Edge::EdgeCoordsVertex::index<Level>(i, P1Edge::EdgeCoordsVertex::VERTEX_C)];

    for (auto neighbor : BubbleEdge::EdgeCoordsVertex::neighbors_south) {
      out << fmt::format("{}\t{}\t{}\n", dst_id, src[BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)], edge_vertex_stencil[neighbor]);
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (auto neighbor : BubbleEdge::EdgeCoordsVertex::neighbors_north) {
        out << fmt::format("{}\t{}\t{}\n", dst_id, src[BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)], edge_vertex_stencil[neighbor]);
      }
    }
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)
}
}
