#ifndef P1EDGE_HPP
#define P1EDGE_HPP

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/p1memory.hpp"

namespace hhg {

namespace P1Edge {

inline void interpolate(Edge &edge,
                        const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &edgeMemoryId,
                        std::function<real_t(const hhg::Point3D &)> &expr,
                        size_t level) {
  EdgeP1FunctionMemory *edgeMemory = edge.getData(edgeMemoryId);

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  Point3D x = edge.coords[0];
  Point3D dx = edge.direction/(real_t) (rowsize - 1);
  x += dx;

  for (size_t i = 1; i < rowsize - 1; ++i) {
    edgeMemory->data[level][i] = expr(x);
    x += dx;
  }
}

inline void assign(Edge &edge,
                   const std::vector<real_t> &scalars,
                   const std::vector<PrimitiveDataID<EdgeP1FunctionMemory, Edge>> &srcIds,
                   const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId,
                   size_t level) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    real_t tmp = scalars[0]*edge.getData(srcIds[0])->data[level][i];

    for (size_t k = 1; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->data[level][i];
    }

    edge.getData(dstId)->data[level][i] = tmp;
  }
}

inline void add(Edge &edge,
                const std::vector<real_t> &scalars,
                const std::vector<PrimitiveDataID<EdgeP1FunctionMemory, Edge>> &srcIds,
                const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId,
                size_t level) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    real_t tmp = 0.0;

    for (size_t k = 0; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->data[level][i];
    }

    edge.getData(dstId)->data[level][i] += tmp;
  }
}

inline real_t dot(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &lhsMemoryId,
                  const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &rhsMemoryId, size_t level) {
  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    sp += edge.getData(lhsMemoryId)->data[level][i]*edge.getData(rhsMemoryId)->data[level][i];
  }

  return sp;
}

inline void apply(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                  const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &srcId,
                  const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId, size_t level, UpdateType update) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto &opr_data = edge.getData(operatorId)->data[level];
  auto &src = edge.getData(srcId)->data[level];
  auto &dst = edge.getData(dstId)->data[level];

  for (size_t i = 1; i < rowsize - 1; ++i) {
    if (update==Replace) {
      dst[i] = opr_data[2]*src[i - 1] + opr_data[3]*src[i] + opr_data[4]*src[i + 1];
    } else if (update==Add) {
      dst[i] += opr_data[2]*src[i - 1] + opr_data[3]*src[i] + opr_data[4]*src[i + 1];
    }
    dst[i] += opr_data[0]*src[rowsize + i - 1] + opr_data[1]*src[rowsize + i];

    if (edge.faces.size()==2) {
      dst[i] += opr_data[5]*src[rowsize + rowsize - 1 + i - 1] + opr_data[6]*src[rowsize + rowsize - 1 + i];
    }
  }
}

inline void smooth_gs(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                      const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId,
                      const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &rhsId, size_t level) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto &opr_data = edge.getData(operatorId)->data[level];
  auto &dst = edge.getData(dstId)->data[level];
  auto &rhs = edge.getData(rhsId)->data[level];

  for (size_t i = 1; i < rowsize - 1; ++i) {
    dst[i] = rhs[i] - opr_data[2]*dst[i - 1] - opr_data[4]*dst[i + 1];
    dst[i] -= opr_data[0]*dst[rowsize + i - 1] + opr_data[1]*dst[rowsize + i];

    if (edge.faces.size()==2) {
      dst[i] -= opr_data[5]*dst[rowsize + rowsize - 1 + i - 1] + opr_data[6]*dst[rowsize + rowsize - 1 + i];
    }

    dst[i] /= opr_data[3];
  }
}

inline void prolongate(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &memoryId, size_t level) {
  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level);
  size_t i_fine = 1;

  auto &edge_data_f = edge.getData(memoryId)->data[level + 1];
  auto &edge_data_c = edge.getData(memoryId)->data[level];

  for (size_t i_coarse = 0; i_coarse < rowsize_coarse - 1; ++i_coarse) {
    edge_data_f[i_fine] = 0.5*(edge_data_c[i_coarse] + edge_data_c[i_coarse + 1]);
    edge_data_f[i_fine + 1] = edge_data_c[i_coarse + 1];
    i_fine += 2;
  }
}

inline void prolongateQuadratic(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &memoryId, size_t level) {

  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level);
  size_t i_fine = 1;
  real_t invtemp = 1/8.;
  const real_t s1[3] = {3*invtemp, 6*invtemp, -invtemp};
  const real_t s2[3] = {-invtemp, 6*invtemp, 3*invtemp};

  auto &edge_data_f = edge.getData(memoryId)->data[level + 1];
  auto &edge_data_c = edge.getData(memoryId)->data[level];
  size_t i_coarse;
  for (i_coarse = 0; i_coarse < rowsize_coarse - 2; ++i_coarse) {
    edge_data_f[i_fine] =
        (s1[0]*edge_data_c[i_coarse] + s1[1]*edge_data_c[i_coarse + 1] + s1[2]*edge_data_c[i_coarse + 2]);
    edge_data_f[i_fine + 1] = edge_data_c[i_coarse + 1];
    i_fine += 2;
  }
  i_coarse--;
  edge_data_f[i_fine] =
      (s2[0]*edge_data_c[i_coarse] + s2[1]*edge_data_c[i_coarse + 1] + s2[2]*edge_data_c[i_coarse + 2]);

}

inline void restrict(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &memoryId, size_t level) {
  size_t rowsize_fine = levelinfo::num_microvertices_per_edge(level);
  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level - 1);

  auto &edge_data_f = edge.getData(memoryId)->data[level];
  auto &edge_data_c = edge.getData(memoryId)->data[level - 1];

  size_t i_fine = 2;
  size_t i_off = 1;

  for (size_t i_coarse = 1; i_coarse < rowsize_coarse - 1; ++i_coarse) {
    // mid edge
    edge_data_c[i_coarse] = 0.5*edge_data_f[i_fine - 1] + edge_data_f[i_fine] + 0.5*edge_data_f[i_fine + 1];

    for (size_t off_edge = 0; off_edge < edge.faces.size(); ++off_edge) {
      edge_data_c[i_coarse] += 0.5*edge_data_f[rowsize_fine + off_edge*(rowsize_fine - 1) + i_off]
          + 0.5*edge_data_f[rowsize_fine + off_edge*(rowsize_fine - 1) + i_off + 1];
    }

    i_fine += 2;
    i_off += 2;
  }
}

}
}

#endif /* P1EDGE_HPP */
