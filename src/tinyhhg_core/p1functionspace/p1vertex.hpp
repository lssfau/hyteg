#ifndef P1VERTEX_HPP
#define P1VERTEX_HPP

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/p1memory.hpp"

namespace hhg {

namespace P1Vertex {

inline void interpolate(Vertex &vertex,
                        const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &vertexMemoryId,
                        std::function<real_t(const hhg::Point3D &)> &expr,
                        size_t level) {
  VertexP1FunctionMemory *vertexMemory = vertex.getData(vertexMemoryId);
  vertexMemory->data[level][0] = expr(vertex.coords);
}

inline void assign(Vertex &vertex,
                   const std::vector<real_t> &scalars,
                   const std::vector<PrimitiveDataID<VertexP1FunctionMemory, Vertex>> &srcIds,
                   const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &dstId,
                   size_t level) {
  real_t tmp = scalars[0]*vertex.getData(srcIds[0])->data[level][0];

  for (size_t i = 1; i < srcIds.size(); ++i) {
    tmp += scalars[i]*vertex.getData(srcIds[i])->data[level][0];
  }

  vertex.getData(dstId)->data[level][0] = tmp;
}

inline void add(Vertex &vertex,
                const std::vector<real_t> &scalars,
                const std::vector<PrimitiveDataID<VertexP1FunctionMemory, Vertex>> &srcIds,
                const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &dstId,
                size_t level) {
  real_t tmp = 0.0;

  for (size_t i = 0; i < srcIds.size(); ++i) {
    tmp += scalars[i]*vertex.getData(srcIds[i])->data[level][0];
  }

  vertex.getData(dstId)->data[level][0] += tmp;
}

inline real_t dot(Vertex &vertex,
                  const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &lhsMemoryId,
                  const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &rhsMemoryId,
                  size_t level) {
  return vertex.getData(lhsMemoryId)->data[level][0]*vertex.getData(rhsMemoryId)->data[level][0];
}

inline void apply(Vertex &vertex,
                  const PrimitiveDataID<VertexP1StencilMemory, Vertex> &operatorId,
                  const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &srcId,
                  const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &dstId,
                  size_t level,
                  UpdateType update) {
  auto &opr_data = vertex.getData(operatorId)->data[level];
  auto &src = vertex.getData(srcId)->data[level];
  auto &dst = vertex.getData(dstId)->data[level];

  if (update==Replace) {
    dst[0] = opr_data[0]*src[0];
  } else if (update==Add) {
    dst[0] += opr_data[0]*src[0];
  }

  for (size_t i = 0; i < vertex.edges.size(); ++i) {
    dst[0] += opr_data[i + 1]*src[i + 1];
  }
}

inline void smooth_gs(Vertex &vertex, const PrimitiveDataID<VertexP1StencilMemory, Vertex> &operatorId,
                      const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &dstId,
                      const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &rhsId, size_t level) {
  auto &opr_data = vertex.getData(operatorId)->data[level];
  auto &dst = vertex.getData(dstId)->data[level];
  auto &rhs = vertex.getData(rhsId)->data[level];

  dst[0] = rhs[0];

  for (size_t i = 0; i < vertex.edges.size(); ++i) {
    dst[0] -= opr_data[i + 1]*dst[i + 1];
  }

  dst[0] /= opr_data[0];
}

inline void prolongate(Vertex &vertex, const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &memoryId, size_t level) {
  vertex.getData(memoryId)->data[level + 1][0] =
      vertex.getData(memoryId)->data[level][0];
}

inline void prolongateQuadratic(Vertex &vertex,
                                const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &memoryId,
                                size_t level) {
  prolongate(vertex, memoryId, level);
}

inline void restrict(Vertex &vertex, const PrimitiveDataID<VertexP1FunctionMemory, Vertex> &memoryId, size_t level) {
  auto &vertex_data_f = vertex.getData(memoryId)->data[level];
  auto &vertex_data_c = vertex.getData(memoryId)->data[level - 1];

  vertex_data_c[0] = vertex_data_f[0];

  size_t i = 1;
  for (Edge *edge : vertex.edges) {
    vertex_data_c[0] += 0.5*vertex_data_f[i];
    i += 1;
  }
}
}
}

#endif /* P1VERTEX_HPP */
