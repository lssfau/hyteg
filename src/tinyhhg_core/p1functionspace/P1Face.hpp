#pragma once

#include "core/debug/all.h"

#include "tinyhhg_core/primitives/face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "P1Memory.hpp"
#include "P1FaceIndex.hpp"
#include <petscmat.h>

namespace hhg {
namespace P1Face {

using walberla::uint_t;
using walberla::real_c;

template<uint_t Level>
inline void interpolateTmpl(Face &face,
                            const PrimitiveDataID<FaceP1FunctionMemory, Face>& faceMemoryId,
                            std::function<real_t(const hhg::Point3D &)> &expr) {
  using namespace CoordsVertex;

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
      faceMemory->data[Level][index<Level>(j, i, VERTEX_C)] = expr(x);
      x += d0;
    }

    inner_rowsize -= 1;
  }
}

SPECIALIZE(void, interpolateTmpl, interpolate)

template<uint_t Level>
inline void assignTmpl(Face &face,
                   const std::vector<real_t>& scalars,
                   const std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>> &srcIds,
                   const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId) {
  using namespace CoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      real_t tmp = scalars[0]*face.getData(srcIds[0])->data[Level][index<Level>(i, j, VERTEX_C)];

      for (uint_t k = 1; k < srcIds.size(); ++k) {
        tmp += scalars[k]*face.getData(srcIds[k])->data[Level][index<Level>(i, j, VERTEX_C)];
      }
      face.getData(dstId)->data[Level][index<Level>(i, j, VERTEX_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, assignTmpl, assign)

template<uint_t Level>
inline void addTmpl(Face &face,
                const std::vector<real_t>& scalars,
                const std::vector<PrimitiveDataID<FaceP1FunctionMemory, Face>> &srcIds,
                const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId) {
  using namespace CoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      real_t tmp = 0.0;

      for (uint_t k = 0; k < srcIds.size(); ++k) {
        tmp += scalars[k]*face.getData(srcIds[k])->data[Level][index<Level>(i, j, VERTEX_C)];
      }

      face.getData(dstId)->data[Level][index<Level>(i, j, VERTEX_C)] += tmp;
    }

    --inner_rowsize;
  }
}

SPECIALIZE(void, addTmpl, add)

template<uint_t Level>
inline real_t dotTmpl(Face &face,
                  const PrimitiveDataID<FaceP1FunctionMemory, Face>& lhsId,
                  const PrimitiveDataID<FaceP1FunctionMemory, Face>& rhsId) {
  using namespace CoordsVertex;

  real_t sp = 0.0;
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      sp += face.getData(lhsId)->data[Level][index<Level>(i, j, VERTEX_C)]
          *face.getData(rhsId)->data[Level][index<Level>(i, j, VERTEX_C)];
    }
    --inner_rowsize;
  }

  return sp;
}

SPECIALIZE(real_t, dotTmpl, dot)

template<uint_t Level>
inline void apply_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                       const PrimitiveDataID<FaceP1FunctionMemory, Face> &srcId,
                       const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId, UpdateType update) {
  using namespace CoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto &src = face.getData(srcId)->data[Level];
  auto &dst = face.getData(dstId)->data[Level];

  real_t tmp;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = opr_data[VERTEX_C]*src[index<Level>(i, j, VERTEX_C)];

      for (auto neighbor : neighbors) {
        tmp += opr_data[neighbor]*src[index<Level>(i, j, neighbor)];
      }

      if (update==Replace) {
        dst[index<Level>(i, j, VERTEX_C)] = tmp;
      } else if (update==Add) {
        dst[index<Level>(i, j, VERTEX_C)] += tmp;
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
  using namespace CoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto &dst = face.getData(dstId)->data[Level];
  auto &rhs = face.getData(rhsId)->data[Level];

  real_t tmp;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = rhs[index<Level>(i, j, VERTEX_C)];

      for (auto neighbor : neighbors) {
        tmp -= opr_data[neighbor]*dst[index<Level>(i, j, neighbor)];
      }

      dst[index<Level>(i, j, VERTEX_C)] = tmp/opr_data[VERTEX_C];
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, smooth_gs_tmpl, smooth_gs)

template<uint_t Level>
inline void smooth_jac_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                            const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId,
                            const PrimitiveDataID<FaceP1FunctionMemory, Face> &rhsId,
                            const PrimitiveDataID<FaceP1FunctionMemory, Face> &tmpId) {
  using namespace CoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto &dst = face.getData(dstId)->data[Level];
  auto &rhs = face.getData(rhsId)->data[Level];
  auto &tmpVar = face.getData(tmpId)->data[Level];

  real_t tmp;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = rhs[index<Level>(i, j, VERTEX_C)];

      for (auto neighbor : neighbors) {
        tmp -= opr_data[neighbor]*tmpVar[index<Level>(i, j, neighbor)];
      }

      dst[index<Level>(i, j, VERTEX_C)] = tmp/opr_data[VERTEX_C];
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, smooth_jac_tmpl, smooth_jac)

template<uint_t Level>
inline void prolongate_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face>& memoryId) {
  using namespace CoordsVertex;

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level);
  uint_t N_c_i = N_c;

  auto &v_f = face.getData(memoryId)->data[Level + 1];
  auto &v_c = face.getData(memoryId)->data[Level];

  uint_t j;

  for (uint_t i = 1; i < N_c - 1; ++i) {
    for (j = 1; j < N_c_i - 2; ++j) {
      v_f[index<Level + 1>(2*i, 2*j, VERTEX_C)] = v_c[index<Level>(i, j, VERTEX_C)];
      v_f[index<Level + 1>(2*i - 1, 2*j - 1, VERTEX_C)] =
          0.5*(v_c[index<Level>(i - 1, j, VERTEX_C)] + v_c[index<Level>(i, j - 1, VERTEX_C)]);
      v_f[index<Level + 1>(2*i - 1, 2*j, VERTEX_C)] = 0.5*(v_c[index<Level>(i, j, VERTEX_C)] + v_c[index<Level>(i - 1, j, VERTEX_C)]);
      v_f[index<Level + 1>(2*i, 2*j - 1, VERTEX_C)] = 0.5*(v_c[index<Level>(i, j, VERTEX_C)] + v_c[index<Level>(i, j - 1, VERTEX_C)]);
    }

    v_f[index<Level + 1>(2*i - 1, 2*j - 1, VERTEX_C)] = 0.5*(v_c[index<Level>(i - 1, j, VERTEX_C)] + v_c[index<Level>(i, j - 1, VERTEX_C)]);
    v_f[index<Level + 1>(2*i - 1, 2*j, VERTEX_C)] = 0.5*(v_c[index<Level>(i, j, VERTEX_C)] + v_c[index<Level>(i - 1, j, VERTEX_C)]);
    v_f[index<Level + 1>(2*i, 2*j - 1, VERTEX_C)] = 0.5*(v_c[index<Level>(i, j, VERTEX_C)] + v_c[index<Level>(i, j - 1, VERTEX_C)]);

    --N_c_i;
  }
}

SPECIALIZE(void, prolongate_tmpl, prolongate)

template<uint_t Level>
inline void prolongateQuadratic_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face>& memoryId) {
  using namespace CoordsVertex;

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
    linearx = 0.5*(v_c[index<Level>(i, j - 2, VERTEX_C)] + v_c[index<Level>(i, j, VERTEX_C)]);
    lineary = 0.5*(v_c[index<Level>(i + 2, j - 2, VERTEX_C)] + v_c[index<Level>(i, j - 2, VERTEX_C)]);
    linearxy = 0.5*(v_c[index<Level>(i + 2, j - 2, VERTEX_C)] + v_c[index<Level>(i, j, VERTEX_C)]);

    offx = v_c[index<Level>(i, j - 1, VERTEX_C)] - linearx;
    offy = v_c[index<Level>(i + 1, j - 2, VERTEX_C)] - lineary;
    offxy = v_c[index<Level>(i + 1, j - 1, VERTEX_C)] - linearxy;

// left bottom corner
    v_f[index<Level + 1>(2*i + 1, 2*j - 3, VERTEX_C)] = 0.5*(linearx + lineary) + 0.5*offx + 0.5*offy + 0.25*offxy;
// right bottom corner
    v_f[index<Level + 1>(2*i + 1, 2*j - 2, VERTEX_C)] = 0.5*(linearx + linearxy) + 0.5*offx + 0.25*offy + 0.5*offxy;
// top corner
    v_f[index<Level + 1>(2*i + 2, 2*j - 3, VERTEX_C)] = 0.5*(linearxy + lineary) + 0.25*offx + 0.5*offy + 0.5*offxy;
  }

  N_c_i -= 1;

  for (i = 2; i < N_c - 1; i += 2) {
    for (j = 2; j < N_c_i - 1; j += 2) {
// upper triangle inner points
//calculate offsets
      linearx = 0.5*(v_c[index<Level>(i, j - 2, VERTEX_C)] + v_c[index<Level>(i, j, VERTEX_C)]);
      lineary = 0.5*(v_c[index<Level>(i + 2, j - 2, VERTEX_C)] + v_c[index<Level>(i, j - 2, VERTEX_C)]);
      linearxy = 0.5*(v_c[index<Level>(i + 2, j - 2, VERTEX_C)] + v_c[index<Level>(i, j, VERTEX_C)]);

      offx = v_c[index<Level>(i, j - 1, VERTEX_C)] - linearx;
      offy = v_c[index<Level>(i + 1, j - 2, VERTEX_C)] - lineary;
      offxy = v_c[index<Level>(i + 1, j - 1, VERTEX_C)] - linearxy;
// left bottom corner
      v_f[index<Level + 1>(2*i + 1, 2*j - 3, VERTEX_C)] = 0.5*(linearx + lineary) + 0.5*offx + 0.5*offy + 0.25*offxy;
// right bottom corner
      v_f[index<Level + 1>(2*i + 1, 2*j - 2, VERTEX_C)] = 0.5*(linearx + linearxy) + 0.5*offx + 0.25*offy + 0.5*offxy;
// top corner
      v_f[index<Level + 1>(2*i + 2, 2*j - 3, VERTEX_C)] = 0.5*(linearxy + lineary) + 0.25*offx + 0.5*offy + 0.5*offxy;

// lower triangle all points
//calculate offsets
      lineary = 0.5*(v_c[index<Level>(i - 2, j, VERTEX_C)] + v_c[index<Level>(i, j, VERTEX_C)]);
      linearxy = 0.5*(v_c[index<Level>(i - 2, j, VERTEX_C)] + v_c[index<Level>(i, j - 2, VERTEX_C)]);

      offy = v_c[index<Level>(i - 1, j, VERTEX_C)] - lineary;
      offxy = v_c[index<Level>(i - 1, j - 1, VERTEX_C)] - linearxy;
// first inner points
// left bottom corner
      v_f[index<Level + 1>(2*i - 1, 2*j - 1, VERTEX_C)] = 0.5*(linearx + lineary) + 0.5*offx + 0.5*offy + 0.25*offxy;
// right bottom corner
      v_f[index<Level + 1>(2*i - 1, 2*j - 2, VERTEX_C)] = 0.5*(linearx + linearxy) + 0.5*offx + 0.25*offy + 0.5*offxy;
// top corner
      v_f[index<Level + 1>(2*i - 2, 2*j - 1, VERTEX_C)] = 0.5*(linearxy + lineary) + 0.25*offx + 0.5*offy + 0.5*offxy;

// boundary points
// x-direction
      v_f[index<Level + 1>(2*i, 2*j - 1, VERTEX_C)] = 0.5*(linearx + v_c[index<Level>(i, j, VERTEX_C)]) + 0.75*offx;
      v_f[index<Level + 1>(2*i, 2*j - 3, VERTEX_C)] = 0.5*(linearx + v_c[index<Level>(i, j - 2, VERTEX_C)]) + 0.75*offx;
//y-direction
      v_f[index<Level + 1>(2*i - 1, 2*j, VERTEX_C)] = 0.5*(v_c[index<Level>(i, j, VERTEX_C)] + lineary) + 0.75*offy;
      v_f[index<Level + 1>(2*i - 3, 2*j, VERTEX_C)] = 0.5*(v_c[index<Level>(i - 2, j, VERTEX_C)] + lineary) + 0.75*offy;
//xy-direction
      v_f[index<Level + 1>(2*i - 1, 2*j - 3, VERTEX_C)] = 0.5*(v_c[index<Level>(i, j - 2, VERTEX_C)] + linearxy) + 0.75*offxy;
      v_f[index<Level + 1>(2*i - 3, 2*j - 1, VERTEX_C)] = 0.5*(v_c[index<Level>(i - 2, j, VERTEX_C)] + linearxy) + 0.75*offxy;
// coarse points
      v_f[index<Level + 1>(2*i, 2*j, VERTEX_C)] = v_c[index<Level>(i, j, VERTEX_C)];
      v_f[index<Level + 1>(2*i, 2*j - 2, VERTEX_C)] = v_c[index<Level>(i, j - 1, VERTEX_C)];
      v_f[index<Level + 1>(2*i - 2, 2*j, VERTEX_C)] = v_c[index<Level>(i - 1, j, VERTEX_C)];
      v_f[index<Level + 1>(2*i - 2, 2*j - 2, VERTEX_C)] = v_c[index<Level>(i - 1, j - 1, VERTEX_C)];
    }
    N_c_i -= 2;

  }
}

SPECIALIZE(void, prolongateQuadratic_tmpl, prolongateQuadratic)

template<uint_t Level>
inline void restrict_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory, Face> &memoryId) {
  using namespace CoordsVertex;

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level - 1);
  uint_t N_c_i = N_c;

  auto &v_f = face.getData(memoryId)->data[Level];
  auto &v_c = face.getData(memoryId)->data[Level - 1];

  real_t tmp;

  for (uint_t i = 1; i < N_c - 2; ++i) {
    for (uint_t j = 1; j < N_c_i - 2; ++j) {
      tmp = v_f[index<Level>(2*i, 2*j, VERTEX_C)];

      for (auto neighbor : neighbors) {
        tmp += 0.5*v_f[index<Level>(2*i, 2*j, neighbor)];
      }

      v_c[index<Level - 1>(i, j, VERTEX_C)] = tmp;
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
                              const PrimitiveDataID<FaceP1FunctionMemory, Face> &dstId, Mat& mat) {
  using namespace CoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto &src = face.getData(srcId)->data[Level];
  auto &dst = face.getData(dstId)->data[Level];


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt srcInt = (PetscInt)src[index<Level>(i, j, VERTEX_C)];
      PetscInt dstInt = (PetscInt)dst[index<Level>(i, j, VERTEX_C)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, VERTEX_C)], opr_data[VERTEX_C]);
      MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[VERTEX_C] ,INSERT_VALUES);

      for (auto neighbor : neighbors) {
        srcInt = (PetscInt)src[index<Level>(i, j, neighbor)];
        //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, neighbor)], opr_data[neighbor]);
        MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[neighbor] ,INSERT_VALUES);
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)

}// namespace P1Face
}// namespace hhg

