#ifndef P1FACE_HPP
#define P1FACE_HPP

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/p1functionspace/p1memory.hpp"

namespace hhg {
namespace P1Face {

enum Dir {
  S = 0,
  SE = 1,
  W = 2,
  C = 3,
  E = 4,
  NW = 5,
  N = 6
};

const Dir neighbors_with_center[] = {S, SE, W, C, E, NW, N};
const Dir neighbors[] = {S, SE, W, E, NW, N};

template<size_t Level>
inline size_t index(size_t row, size_t col, Dir dir) {
  size_t h = levelinfo::num_microvertices_per_edge(Level);
  size_t n = h*(h + 1)/2;
  size_t center = (n - (h - row)*(h - row + 1)/2) + col;
  switch (dir) {
    case C:return center;
    case N:return center + h - row;
    case E:return center + 1;
    case S:return center - h - 1 + row;
    case W:return center - 1;
    case SE:return center - h + row;
    case NW:return center + h - row - 1;
  }
  return 0;
}

template<size_t Level>
inline void interpolateTmpl(Face &face,
                            const PrimitiveDataID<FaceP1FunctionMemory, Face> &faceMemoryId,
                            std::function<real_t(const hhg::Point3D &)> &expr) {
  FaceP1FunctionMemory *faceMemory = face.getData(faceMemoryId);

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0;

  x0 = face.coords[0];

  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));

  size_t inner_rowsize = rowsize;

  for (size_t i = 1; i < rowsize - 2; ++i) {
    x = x0;
    x += i*d2 + d0;

    for (size_t j = 1; j < inner_rowsize - 2; ++j) {
      faceMemory->data[Level][index<Level>(i, j, C)] = expr(x);
      x += d0;
    }

    inner_rowsize -= 1;
  }
}

SPECIALIZE(void, interpolateTmpl, interpolate)

inline void assign(Face &face,
                   const std::vector<real_t> &scalars,
                   const std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>> &srcIds,
                   const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId,
                   size_t level) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i) {
    for (size_t j = 0; j < inner_rowsize - 3; ++j) {
      real_t tmp = scalars[0]*face.getData(srcIds[0])->data[level][mr];

      for (size_t k = 1; k < srcIds.size(); ++k) {
        tmp += scalars[k]*face.getData(srcIds[k])->data[level][mr];
      }
      face.getData(dstId)->data[level][mr] = tmp;

      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }
}

inline void add(Face &face,
                const std::vector<real_t> &scalars,
                const std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>> &srcIds,
                const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId,
                size_t level) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i) {
    for (size_t j = 0; j < inner_rowsize - 3; ++j) {
      real_t tmp = 0.0;

      for (size_t k = 0; k < srcIds.size(); ++k) {
        tmp += scalars[k]*face.getData(srcIds[k])->data[level][mr];
      }

      face.getData(dstId)->data[level][mr] += tmp;

      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }
}

inline real_t dot(Face &face,
                  const PrimitiveDataID<FaceP1FunctionMemory, Face> &lhsId,
                  const PrimitiveDataID<FaceP1FunctionMemory, Face> &rhsId,
                  size_t level) {
  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i) {
    for (size_t j = 0; j < inner_rowsize - 3; ++j) {
      sp += face.getData(lhsId)->data[level][mr]
          *face.getData(rhsId)->data[level][mr];
      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }

  return sp;
}

template<size_t Level>
inline void apply_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face> &operatorId,
                       const PrimitiveDataID<FaceP1FunctionMemory, Face> &srcId,
                       const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId, UpdateType update) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto &src = face.getData(srcId)->data[Level];
  auto &dst = face.getData(dstId)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = opr_data[C]*src[index<Level>(i, j, C)];

      for (auto neighbor : neighbors) {
        tmp += opr_data[neighbor]*src[index<Level>(i, j, neighbor)];
      }

      if (update==Replace) {
        dst[index<Level>(i, j, C)] = tmp;
      } else if (update==Add) {
        dst[index<Level>(i, j, C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, apply_tmpl, apply)

template<size_t Level>
inline void smooth_gs_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face> &operatorId,
                           const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId,
                           const PrimitiveDataID<FaceP1FunctionMemory, Face> &rhsId) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto &dst = face.getData(dstId)->data[Level];
  auto &rhs = face.getData(rhsId)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = rhs[index<Level>(i, j, C)];

      for (auto neighbor : neighbors) {
        tmp -= opr_data[neighbor]*dst[index<Level>(i, j, neighbor)];
      }

      dst[index<Level>(i, j, C)] = tmp/opr_data[C];
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, smooth_gs_tmpl, smooth_gs)

template<size_t Level>
inline void prolongate_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face> &memoryId) {
  size_t N_c = levelinfo::num_microvertices_per_edge(Level);
  size_t N_c_i = N_c;

  auto &v_f = face.getData(memoryId)->data[Level + 1];
  auto &v_c = face.getData(memoryId)->data[Level];

  size_t j;

  for (size_t i = 1; i < N_c - 1; ++i) {
    for (j = 1; j < N_c_i - 2; ++j) {
      v_f[index<Level + 1>(2*i, 2*j, C)] = v_c[index<Level>(i, j, C)];
      v_f[index<Level + 1>(2*i - 1, 2*j - 1, C)] =
          0.5*(v_c[index<Level>(i - 1, j, C)] + v_c[index<Level>(i, j - 1, C)]);
      v_f[index<Level + 1>(2*i - 1, 2*j, C)] = 0.5*(v_c[index<Level>(i, j, C)] + v_c[index<Level>(i - 1, j, C)]);
      v_f[index<Level + 1>(2*i, 2*j - 1, C)] = 0.5*(v_c[index<Level>(i, j, C)] + v_c[index<Level>(i, j - 1, C)]);
    }

    v_f[index<Level + 1>(2*i - 1, 2*j - 1, C)] = 0.5*(v_c[index<Level>(i - 1, j, C)] + v_c[index<Level>(i, j - 1, C)]);
    v_f[index<Level + 1>(2*i - 1, 2*j, C)] = 0.5*(v_c[index<Level>(i, j, C)] + v_c[index<Level>(i - 1, j, C)]);
    v_f[index<Level + 1>(2*i, 2*j - 1, C)] = 0.5*(v_c[index<Level>(i, j, C)] + v_c[index<Level>(i, j - 1, C)]);

    --N_c_i;
  }
}

SPECIALIZE(void, prolongate_tmpl, prolongate)

template<size_t Level>
inline void prolongateQuadratic_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face> &memoryId) {
  size_t N_c = levelinfo::num_microvertices_per_edge(Level);
  size_t N_c_i = N_c;
  auto &v_f = face.getData(memoryId)->data[Level + 1];
  auto &v_c = face.getData(memoryId)->data[Level];

  size_t i, j;
  real_t linearx, lineary, linearxy, offx, offy, offxy;
  i = 0;
  for (j = 2; j <= N_c - 1; j += 2) {
// upper triangle inner points
//calculate offsets
    linearx = 0.5*(v_c[index<Level>(i, j - 2, C)] + v_c[index<Level>(i, j, C)]);
    lineary = 0.5*(v_c[index<Level>(i + 2, j - 2, C)] + v_c[index<Level>(i, j - 2, C)]);
    linearxy = 0.5*(v_c[index<Level>(i + 2, j - 2, C)] + v_c[index<Level>(i, j, C)]);

    offx = v_c[index<Level>(i, j - 1, C)] - linearx;
    offy = v_c[index<Level>(i + 1, j - 2, C)] - lineary;
    offxy = v_c[index<Level>(i + 1, j - 1, C)] - linearxy;

// left bottom corner
    v_f[index<Level + 1>(2*i + 1, 2*j - 3, C)] = 0.5*(linearx + lineary) + 0.5*offx + 0.5*offy + 0.25*offxy;
// right bottom corner
    v_f[index<Level + 1>(2*i + 1, 2*j - 2, C)] = 0.5*(linearx + linearxy) + 0.5*offx + 0.25*offy + 0.5*offxy;
// top corner
    v_f[index<Level + 1>(2*i + 2, 2*j - 3, C)] = 0.5*(linearxy + lineary) + 0.25*offx + 0.5*offy + 0.5*offxy;
  }

  N_c_i -= 1;

  for (i = 2; i < N_c - 1; i += 2) {
    for (j = 2; j < N_c_i - 1; j += 2) {
// upper triangle inner points
//calculate offsets
      linearx = 0.5*(v_c[index<Level>(i, j - 2, C)] + v_c[index<Level>(i, j, C)]);
      lineary = 0.5*(v_c[index<Level>(i + 2, j - 2, C)] + v_c[index<Level>(i, j - 2, C)]);
      linearxy = 0.5*(v_c[index<Level>(i + 2, j - 2, C)] + v_c[index<Level>(i, j, C)]);

      offx = v_c[index<Level>(i, j - 1, C)] - linearx;
      offy = v_c[index<Level>(i + 1, j - 2, C)] - lineary;
      offxy = v_c[index<Level>(i + 1, j - 1, C)] - linearxy;
// left bottom corner
      v_f[index<Level + 1>(2*i + 1, 2*j - 3, C)] = 0.5*(linearx + lineary) + 0.5*offx + 0.5*offy + 0.25*offxy;
// right bottom corner
      v_f[index<Level + 1>(2*i + 1, 2*j - 2, C)] = 0.5*(linearx + linearxy) + 0.5*offx + 0.25*offy + 0.5*offxy;
// top corner
      v_f[index<Level + 1>(2*i + 2, 2*j - 3, C)] = 0.5*(linearxy + lineary) + 0.25*offx + 0.5*offy + 0.5*offxy;

// lower triangle all points
//calculate offsets
      lineary = 0.5*(v_c[index<Level>(i - 2, j, C)] + v_c[index<Level>(i, j, C)]);
      linearxy = 0.5*(v_c[index<Level>(i - 2, j, C)] + v_c[index<Level>(i, j - 2, C)]);

      offy = v_c[index<Level>(i - 1, j, C)] - lineary;
      offxy = v_c[index<Level>(i - 1, j - 1, C)] - linearxy;
// first inner points
// left bottom corner
      v_f[index<Level + 1>(2*i - 1, 2*j - 1, C)] = 0.5*(linearx + lineary) + 0.5*offx + 0.5*offy + 0.25*offxy;
// right bottom corner
      v_f[index<Level + 1>(2*i - 1, 2*j - 2, C)] = 0.5*(linearx + linearxy) + 0.5*offx + 0.25*offy + 0.5*offxy;
// top corner
      v_f[index<Level + 1>(2*i - 2, 2*j - 1, C)] = 0.5*(linearxy + lineary) + 0.25*offx + 0.5*offy + 0.5*offxy;

// boundary points
// x-direction
      v_f[index<Level + 1>(2*i, 2*j - 1, C)] = 0.5*(linearx + v_c[index<Level>(i, j, C)]) + 0.75*offx;
      v_f[index<Level + 1>(2*i, 2*j - 3, C)] = 0.5*(linearx + v_c[index<Level>(i, j - 2, C)]) + 0.75*offx;
//y-direction
      v_f[index<Level + 1>(2*i - 1, 2*j, C)] = 0.5*(v_c[index<Level>(i, j, C)] + lineary) + 0.75*offy;
      v_f[index<Level + 1>(2*i - 3, 2*j, C)] = 0.5*(v_c[index<Level>(i - 2, j, C)] + lineary) + 0.75*offy;
//xy-direction
      v_f[index<Level + 1>(2*i - 1, 2*j - 3, C)] = 0.5*(v_c[index<Level>(i, j - 2, C)] + linearxy) + 0.75*offxy;
      v_f[index<Level + 1>(2*i - 3, 2*j - 1, C)] = 0.5*(v_c[index<Level>(i - 2, j, C)] + linearxy) + 0.75*offxy;
// coarse points
      v_f[index<Level + 1>(2*i, 2*j, C)] = v_c[index<Level>(i, j, C)];
      v_f[index<Level + 1>(2*i, 2*j - 2, C)] = v_c[index<Level>(i, j - 1, C)];
      v_f[index<Level + 1>(2*i - 2, 2*j, C)] = v_c[index<Level>(i - 1, j, C)];
      v_f[index<Level + 1>(2*i - 2, 2*j - 2, C)] = v_c[index<Level>(i - 1, j - 1, C)];
    }
    N_c_i -= 2;

  }
}

SPECIALIZE(void, prolongateQuadratic_tmpl, prolongateQuadratic)

template<size_t Level>
inline void restrict_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face> &memoryId) {
  size_t N_c = levelinfo::num_microvertices_per_edge(Level - 1);
  size_t N_c_i = N_c;

  auto &v_f = face.getData(memoryId)->data[Level];
  auto &v_c = face.getData(memoryId)->data[Level - 1];

  real_t tmp;

  for (size_t i = 1; i < N_c - 2; ++i) {
    for (size_t j = 1; j < N_c_i - 2; ++j) {
      tmp = v_f[index<Level>(2*i, 2*j, C)];

      for (auto neighbor : neighbors) {
        tmp += 0.5*v_f[index<Level>(2*i, 2*j, neighbor)];
      }

      v_c[index<Level - 1>(i, j, C)] = tmp;
    }

    --N_c_i;
  }
}

SPECIALIZE(void, restrict_tmpl, restrict)

/// Checks if a given index is a the boundary of the face
/// \param index The index which should be checked
/// \param length Size of the triangle in the first dimension
bool is_boundary(size_t index, size_t length) {
  if (index < length) return true;
  while (index >= length) {
    index -= length;
    length--;
  }
  return (index==0 || index==(length - 1));
}


}// namespace P1Face
}// namespace hhg


#endif /* P1FACE_HPP */
