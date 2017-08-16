#pragma once

#include "core/debug/all.h"

#include "tinyhhg_core/primitives/face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "P1Memory.hpp"
#include "P1FaceIndex.hpp"

namespace hhg {
namespace P1Face {

using walberla::uint_t;
using walberla::real_c;

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

template<uint_t Level>
inline uint_t index(uint_t col, uint_t row, Dir dir)
{
  uint_t h = levelinfo::num_microvertices_per_edge(Level);
  WALBERLA_ASSERT_LESS_EQUAL(row,h);
  WALBERLA_ASSERT_LESS_EQUAL(col,h);
  uint_t n = h*(h + 1)/2;
  uint_t center = (n - (h - row)*(h - row + 1)/2) + col;
  switch (dir) {
    case C:return center;
    case N:return center + h - row;
    case E:return center + 1;
    case S:return center - h - 1 + row;
    case W:return center - 1;
    case SE:return center - h + row;
    case NW:return center + h - row - 1;
  }
  return std::numeric_limits<uint_t>::max();
}

template<uint_t Level>
inline void interpolateTmpl(Face &face,
                            const PrimitiveDataID<FaceP1FunctionMemory, Face>& faceMemoryId,
                            std::function<real_t(const hhg::Point3D &)> &expr) {
  FaceP1FunctionMemory *faceMemory = face.getData(faceMemoryId);

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0;

  x0 = face.coords[0];

  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));

  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    x = x0;
    x += real_c(i)*d2 + d0;

    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      faceMemory->data[Level][index<Level>(j, i, C)] = expr(x);
      x += d0;
    }

    inner_rowsize -= 1;
  }
}

SPECIALIZE(void, interpolateTmpl, interpolate)

inline void assign(Face &face,
                   const std::vector<real_t>& scalars,
                   const std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>> &srcIds,
                   const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId,
                   uint_t level) {
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  uint_t inner_rowsize = rowsize;

  uint_t mr = 1 + rowsize;

  for (uint_t i = 0; i < rowsize - 3; ++i) {
    for (uint_t j = 0; j < inner_rowsize - 3; ++j) {
      real_t tmp = scalars[0]*face.getData(srcIds[0])->data[level][mr];

      for (uint_t k = 1; k < srcIds.size(); ++k) {
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
                const std::vector<real_t>& scalars,
                const std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>> &srcIds,
                const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId,
                uint_t level) {
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  uint_t inner_rowsize = rowsize;

  uint_t mr = 1 + rowsize;

  for (uint_t i = 0; i < rowsize - 3; ++i) {
    for (uint_t j = 0; j < inner_rowsize - 3; ++j) {
      real_t tmp = 0.0;

      for (uint_t k = 0; k < srcIds.size(); ++k) {
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
                  const PrimitiveDataID<FaceP1FunctionMemory, Face>& lhsId,
                  const PrimitiveDataID<FaceP1FunctionMemory, Face>& rhsId,
                  uint_t level) {
  real_t sp = 0.0;
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  uint_t inner_rowsize = rowsize;

  uint_t mr = 1 + rowsize;

  for (uint_t i = 0; i < rowsize - 3; ++i) {
    for (uint_t j = 0; j < inner_rowsize - 3; ++j) {
      sp += face.getData(lhsId)->data[level][mr]
          *face.getData(rhsId)->data[level][mr];
      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }

  return sp;
}

template<uint_t Level>
inline void apply_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                       const PrimitiveDataID<FaceP1FunctionMemory, Face> &srcId,
                       const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId, UpdateType update) {
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto &src = face.getData(srcId)->data[Level];
  auto &dst = face.getData(dstId)->data[Level];

  real_t tmp;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
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

template<uint_t Level>
inline void smooth_gs_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                           const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId,
                           const PrimitiveDataID<FaceP1FunctionMemory, Face> &rhsId) {
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto &dst = face.getData(dstId)->data[Level];
  auto &rhs = face.getData(rhsId)->data[Level];

  real_t tmp;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
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

template<uint_t Level>
inline void prolongate_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face>& memoryId) {
  uint_t N_c = levelinfo::num_microvertices_per_edge(Level);
  uint_t N_c_i = N_c;

  auto &v_f = face.getData(memoryId)->data[Level + 1];
  auto &v_c = face.getData(memoryId)->data[Level];

  uint_t j;

  for (uint_t i = 1; i < N_c - 1; ++i) {
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

template<uint_t Level>
inline void prolongateQuadratic_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face>& memoryId) {
  uint_t N_c = levelinfo::num_microvertices_per_edge(Level);
  uint_t N_c_i = N_c;
  auto &v_f = face.getData(memoryId)->data[Level + 1];
  auto &v_c = face.getData(memoryId)->data[Level];

  uint_t i, j;
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

template<uint_t Level>
inline void restrict_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face> &memoryId) {
  uint_t N_c = levelinfo::num_microvertices_per_edge(Level - 1);
  uint_t N_c_i = N_c;

  auto &v_f = face.getData(memoryId)->data[Level];
  auto &v_c = face.getData(memoryId)->data[Level - 1];

  real_t tmp;

  for (uint_t i = 1; i < N_c - 2; ++i) {
    for (uint_t j = 1; j < N_c_i - 2; ++j) {
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
inline bool is_boundary(uint_t index, uint_t length) {
  if (index < length) return true;
  while (index >= length) {
    index -= length;
    length--;
  }
  return (index==0 || index==(length - 1));
}

//
//template<uint_t Level>
//inline void printFunctionMemory(Face& face, PrimitiveDataID<FaceP1FunctionMemory,Face> memory_id){
//  using namespace std;
//  auto& faceMemory = hhg::P1Bubble::getFaceFunctionMemory(face, 0)->data[Level];
//  uint_t verticesPerDge = hhg::levelinfo::num_microvertices_per_edge(Level);
//  cout << setfill('=') << setw(100) << "" << endl;
//  cout << face << std::left << setprecision(1) << fixed << setfill(' ') << endl;
//  cout << "Vertices: " << std::endl;
//  for (uint_t i = 0; i < verticesPerDge; ++i) {
//    for (uint_t j = 0; j < verticesPerDge - i; ++j) {
//      cout << setw(5) << faceMemory[CoordsVertex::index<Level>(i, j, CoordsVertex::VERTEX_C)] << "|";
//    }
//    std::cout << std::endl;
//  }
//  cout << setw(100) << setfill(' ') << endl;
//
//}

inline void enumerate(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId, size_t level, uint_t& num) {
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  uint_t inner_rowsize = rowsize;

  uint_t mr = 1 + rowsize;

  for (uint_t i = 0; i < rowsize - 3; ++i) {
    for (uint_t j = 0; j < inner_rowsize - 3; ++j) {

      face.getData(dstId)->data[level][mr] = walberla::real_c(num++);

      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }
}

template<uint_t Level>
inline void saveOperator_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                              const PrimitiveDataID<FaceP1FunctionMemory, Face> &srcId,
                              const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId, std::ostream& out) {
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto &src = face.getData(srcId)->data[Level];
  auto &dst = face.getData(dstId)->data[Level];


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, C)], src[index<Level>(i, j, C)], opr_data[C]);

      for (auto neighbor : neighbors) {
        out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, C)], src[index<Level>(i, j, neighbor)], opr_data[neighbor]);
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)

}// namespace P1Face
}// namespace hhg

