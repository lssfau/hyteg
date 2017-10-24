#pragma once

#include "core/debug/all.h"

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "tinyhhg_core/p1functionspace/P1FaceIndex.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFaceIndex.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

namespace hhg {
namespace P1Face {

using walberla::uint_t;
using walberla::real_c;

template<typename ValueType, uint_t Level>
inline ValueType assembleLocal(uint_t i, uint_t j, const Matrix3r& localMatrix,
                               double* src,
                               double* coeff,
                               const std::array<FaceCoordsVertex::DirVertex,3>& vertices,
                               const std::array<uint_t,3>& idx)
{
  using namespace FaceCoordsVertex;

  ValueType meanCoeff = 1.0/3.0 * (coeff[index<Level>(i, j, vertices[0])]
                                 + coeff[index<Level>(i, j, vertices[1])]
                                 + coeff[index<Level>(i, j, vertices[2])]);

  ValueType tmp;
  tmp  = localMatrix(idx[0],idx[0]) * src[index<Level>(i, j, vertices[0])]
         + localMatrix(idx[0],idx[1]) * src[index<Level>(i, j, vertices[1])]
         + localMatrix(idx[0],idx[2]) * src[index<Level>(i, j, vertices[2])];
  return meanCoeff * tmp;
}

template<typename ValueType, uint_t Level>
inline ValueType assembleLocalDG(uint_t i, uint_t j, const Matrix3r& localMatrix,
                               double* src,
                               const std::array<FaceCoordsVertex::DirVertex,3>& vertices,
                               const std::array<uint_t,3>& idx)
{
  using namespace FaceCoordsVertex;

  ValueType tmp;
  tmp  = localMatrix(idx[0],idx[0]) * src[index<Level>(i, j, vertices[0])]
      + localMatrix(idx[0],idx[1]) * src[index<Level>(i, j, vertices[1])]
      + localMatrix(idx[0],idx[2]) * src[index<Level>(i, j, vertices[2])];
  return tmp;
}

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Face &face,
                            const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>& faceMemoryId,
                            std::function<ValueType(const hhg::Point3D &)> &expr) {
  using namespace FaceCoordsVertex;

  FaceP1FunctionMemory< ValueType > *faceMemory = face.getData(faceMemoryId);

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0;

  auto dstPtr = faceMemory->getPointer( Level );

  x0 = face.coords[0];

  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));

  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    x = x0;
    x += real_c(i)*d2 + d0;

    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      dstPtr[index<Level>(j, i, VERTEX_C)] = expr(x);
      x += d0;
    }

    inner_rowsize -= 1;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, interpolateTmpl, interpolate)

template< typename ValueType, uint_t Level >
inline void assignTmpl(Face &face,
                   const std::vector<ValueType>& scalars,
                   const std::vector<PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>> &srcIds,
                   const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &dstId) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* dst = face.getData(dstId)->getPointer( Level );
  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }
  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      ValueType tmp = scalars[0]*srcPtr[0][index<Level>(i, j, VERTEX_C)];

      for (uint_t k = 1; k < srcIds.size(); ++k) {
        tmp += scalars[k]*srcPtr[k][index<Level>(i, j, VERTEX_C)];
      }
      dst[index<Level>(i, j, VERTEX_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, assignTmpl, assign)

template< typename ValueType, uint_t Level >
inline void addTmpl(Face &face,
                const std::vector<ValueType>& scalars,
                const std::vector<PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>> &srcIds,
                const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &dstId) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* dstPtr = face.getData(dstId)->getPointer( Level );
  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      ValueType tmp = 0.0;

      for (uint_t k = 0; k < srcIds.size(); ++k) {
        tmp += scalars[k] * srcPtr[k][index<Level>(i, j, VERTEX_C)];
      }

      dstPtr[index<Level>(i, j, VERTEX_C)] += tmp;
    }

    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, addTmpl, add)

template< typename ValueType, uint_t Level >
inline real_t dotTmpl(Face &face,
                  const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>& lhsId,
                  const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>& rhsId) {
  using namespace FaceCoordsVertex;

  real_t sp = 0.0;
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* lhsPtr = face.getData(lhsId)->getPointer( Level );
  ValueType* rhsPtr = face.getData(rhsId)->getPointer( Level );

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      sp += lhsPtr[index<Level>(i, j, VERTEX_C)]
          * rhsPtr[index<Level>(i, j, VERTEX_C)];
    }
    --inner_rowsize;
  }

  return sp;
}

SPECIALIZE_WITH_VALUETYPE(real_t, dotTmpl, dot)

template< typename ValueType, uint_t Level >
inline void apply_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                       const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &srcId,
                       const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &dstId, UpdateType update) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  ValueType tmp;

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

SPECIALIZE_WITH_VALUETYPE(void, apply_tmpl, apply)

template< typename ValueType, uint_t Level >
inline void applyCoefficientTmpl(Face &face, const PrimitiveDataID<FaceP1LocalMatrixMemory, Face>& operatorId,
                                 const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &srcId,
                                 const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &dstId,
                                 const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &coeffId,
                                 UpdateType update) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto localMatrices = face.getData(operatorId);
  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);
  auto coeff = face.getData(coeffId)->getPointer(Level);

  ValueType tmp;

  std::array<DirVertex,3> triangleBlueSW = { VERTEX_C, VERTEX_W,  VERTEX_S  };
  std::array<DirVertex,3> triangleGrayS  = { VERTEX_C, VERTEX_S,  VERTEX_SE };
  std::array<DirVertex,3> triangleBlueSE = { VERTEX_C, VERTEX_SE, VERTEX_E  };
  std::array<DirVertex,3> triangleGrayNW = { VERTEX_C, VERTEX_W,  VERTEX_NW };
  std::array<DirVertex,3> triangleBlueN  = { VERTEX_C, VERTEX_NW, VERTEX_N  };
  std::array<DirVertex,3> triangleGrayNE = { VERTEX_C, VERTEX_N,  VERTEX_E  };

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {

      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[index<Level>(i, j, VERTEX_C)];
      }

      tmp += assembleLocal<ValueType, Level>(i, j, localMatrices->getGrayMatrix(Level), src, coeff, triangleGrayS, {2,0,1});
      tmp += assembleLocal<ValueType, Level>(i, j, localMatrices->getBlueMatrix(Level), src, coeff, triangleBlueSE, {1,2,0});
      tmp += assembleLocal<ValueType, Level>(i, j, localMatrices->getBlueMatrix(Level), src, coeff, triangleBlueSW, {0,1,2});
      tmp += assembleLocal<ValueType, Level>(i, j, localMatrices->getGrayMatrix(Level), src, coeff, triangleGrayNW, {1,0,2});
      tmp += assembleLocal<ValueType, Level>(i, j, localMatrices->getBlueMatrix(Level), src, coeff, triangleBlueN, {2,1,0});
      tmp += assembleLocal<ValueType, Level>(i, j, localMatrices->getGrayMatrix(Level), src, coeff, triangleGrayNE, {0,2,1});

      dst[index<Level>(i, j, VERTEX_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, applyCoefficientTmpl, applyCoefficient)

template< typename ValueType, uint_t Level >
inline void applyCoefficientDGTmpl(Face &face, const PrimitiveDataID<FaceP1LocalMatrixMemory, Face>& operatorId,
                                 const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &srcId,
                                 const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &dstId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &coeffId,
                                 UpdateType update) {
  using namespace FaceCoordsVertex;
  typedef stencilDirection sD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto localMatrices = face.getData(operatorId);
  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);
  auto coeff = face.getData(coeffId)->getPointer(Level);

  ValueType tmp;

  std::array<DirVertex,3> triangleBlueSW = { VERTEX_C, VERTEX_W,  VERTEX_S  };
  std::array<DirVertex,3> triangleGraySE = { VERTEX_C, VERTEX_S,  VERTEX_SE };
  std::array<DirVertex,3> triangleBlueSE = { VERTEX_C, VERTEX_SE, VERTEX_E  };
  std::array<DirVertex,3> triangleGrayNW = { VERTEX_C, VERTEX_W,  VERTEX_NW };
  std::array<DirVertex,3> triangleBlueNW = { VERTEX_C, VERTEX_NW, VERTEX_N  };
  std::array<DirVertex,3> triangleGrayNE = { VERTEX_C, VERTEX_N,  VERTEX_E  };



  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {

      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[index<Level>(i, j, VERTEX_C)];
      }

      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_GRAY_SE)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getGrayMatrix(Level), src, triangleGraySE, {2,0,1});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_BLUE_SE)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getBlueMatrix(Level), src, triangleBlueSE, {1,2,0});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_BLUE_SW)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getBlueMatrix(Level), src, triangleBlueSW, {0,1,2});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_GRAY_NW)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getGrayMatrix(Level), src, triangleGrayNW, {1,0,2});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_BLUE_NW)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getBlueMatrix(Level), src, triangleBlueNW, {2,1,0});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_GRAY_NE)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getGrayMatrix(Level), src, triangleGrayNE, {0,2,1});

      dst[index<Level>(i, j, VERTEX_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, applyCoefficientDGTmpl, applyCoefficientDG)

template< typename ValueType, uint_t Level >
inline void smooth_gs_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                           const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &dstId,
                           const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &rhsId) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp = rhs[index<Level>(i, j, VERTEX_C)];

      //for (auto neighbor : neighbors) {
      for(uint_t k = 0; k < neighbors.size(); ++k){
        tmp -= opr_data[neighbors[k]]*dst[index<Level>(i, j, neighbors[k])];
      }

      dst[index<Level>(i, j, VERTEX_C)] = tmp/opr_data[VERTEX_C];
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, smooth_gs_tmpl, smooth_gs)

template< typename ValueType, uint_t Level >
inline void smooth_sor_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                            const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &dstId,
                            const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &rhsId,
                            ValueType relax) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp = rhs[index<Level>(i, j, VERTEX_C)];

      //for (auto neighbor : neighbors) {
      for(uint_t k = 0; k < neighbors.size(); ++k){
        tmp -= opr_data[neighbors[k]]*dst[index<Level>(i, j, neighbors[k])];
      }

      dst[index<Level>(i, j, VERTEX_C)] = (1.0-relax) * dst[index<Level>(i, j, VERTEX_C)] + relax * tmp/opr_data[VERTEX_C];
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, smooth_sor_tmpl, smooth_sor)

template< typename ValueType, uint_t Level >
inline void smooth_jac_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                            const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &dstId,
                            const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &rhsId,
                            const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &tmpId) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );
  auto tmpVar = face.getData(tmpId)->getPointer( Level );

  ValueType tmp;

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

SPECIALIZE_WITH_VALUETYPE(void, smooth_jac_tmpl, smooth_jac)

template< typename ValueType, uint_t Level >
inline void prolongate_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>& memoryId) {
  using namespace FaceCoordsVertex;

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level);
  uint_t N_c_i = N_c;

  auto v_f = face.getData(memoryId)->getPointer( Level + 1 );
  auto v_c = face.getData(memoryId)->getPointer( Level );

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

SPECIALIZE_WITH_VALUETYPE(void, prolongate_tmpl, prolongate)

template< typename ValueType, uint_t Level >
inline void prolongateQuadratic_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>& memoryId) {
  using namespace FaceCoordsVertex;

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level);
  uint_t N_c_i = N_c;
  auto v_f = face.getData(memoryId)->getPointer( Level + 1 );
  auto v_c = face.getData(memoryId)->getPointer( Level );

  uint_t i, j;
  ValueType linearx, lineary, linearxy, offx, offy, offxy;
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

SPECIALIZE_WITH_VALUETYPE(void, prolongateQuadratic_tmpl, prolongateQuadratic)

template< typename ValueType, uint_t Level >
inline void restrict_tmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &memoryId) {
  using namespace FaceCoordsVertex;

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level - 1);
  uint_t N_c_i = N_c;

  auto v_f = face.getData(memoryId)->getPointer( Level );
  auto v_c = face.getData(memoryId)->getPointer( Level - 1 );

  ValueType tmp;

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

SPECIALIZE_WITH_VALUETYPE(void, restrict_tmpl, restrict)

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

template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Face &face, const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &dstId, uint_t& num) {
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  uint_t mr = 1 + rowsize;

  ValueType* dstPtr = face.getData(dstId)->getPointer( Level );

  for (uint_t i = 0; i < rowsize - 3; ++i) {
    for (uint_t j = 0; j < inner_rowsize - 3; ++j) {

      dstPtr[mr] = static_cast<ValueType>(num);
      num++;

      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate )

template< typename ValueType, uint_t Level >
inline void integrateDGTmpl(Face &face,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId) {
  using namespace FaceCoordsVertex;
  typedef stencilDirection sD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto rhs = face.getData(rhsId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);

  real_t faceArea = std::pow(4.0, -walberla::real_c(Level)) * face.area;
  real_t weightedFaceArea = faceArea / 3.0;

  ValueType tmp;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {

      tmp =   rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_GRAY_SE)]
            + rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_BLUE_SE)]
            + rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_BLUE_SW)]
            + rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_GRAY_NW)]
            + rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_BLUE_NW)]
            + rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, sD::CELL_GRAY_NE)];

      dst[index<Level>(i, j, VERTEX_C)] = weightedFaceArea  * tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, integrateDGTmpl, integrateDG )

#ifdef HHG_BUILD_WITH_PETSC
template< uint_t Level >
inline void saveOperator_tmpl(Face &face, const PrimitiveDataID<FaceP1StencilMemory, Face>& operatorId,
                              const PrimitiveDataID<FaceP1FunctionMemory< PetscInt >, Face> &srcId,
                              const PrimitiveDataID<FaceP1FunctionMemory< PetscInt >, Face> &dstId, Mat& mat) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto &opr_data = face.getData(operatorId)->data[Level];
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt srcInt = src[index<Level>(i, j, VERTEX_C)];
      PetscInt dstInt = dst[index<Level>(i, j, VERTEX_C)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, VERTEX_C)], opr_data[VERTEX_C]);
      MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[VERTEX_C] ,INSERT_VALUES);

      for (auto neighbor : neighbors) {
        srcInt = src[index<Level>(i, j, neighbor)];
        //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, neighbor)], opr_data[neighbor]);
        MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[neighbor] ,INSERT_VALUES);
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)


template< typename ValueType, uint_t Level >
inline void createVectorFromFunctionTmpl(Face &face,
                              const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &srcId,
                              const PrimitiveDataID<FaceP1FunctionMemory< PetscInt >, Face> &numeratorId,
                              Vec& vec) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt numeratorInt = numerator[index<Level>(i, j, VERTEX_C)];
      VecSetValues(vec,1,&numeratorInt,&src[index<Level>(i, j, VERTEX_C)],INSERT_VALUES);
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, createVectorFromFunctionTmpl, createVectorFromFunction)



template< typename ValueType, uint_t Level >
inline void createFunctionFromVectorTmpl(Face &face,
                                         const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face> &srcId,
                                         const PrimitiveDataID<FaceP1FunctionMemory< PetscInt >, Face> &numeratorId,
                                         Vec& vec) {
  using namespace FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt numeratorInt = numerator[index<Level>(i, j, VERTEX_C)];
      VecGetValues(vec,1,&numeratorInt,&src[index<Level>(i, j, VERTEX_C)]);
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, createFunctionFromVectorTmpl, createFunctionFromVector)
#endif

}// namespace P1Face
}// namespace hhg

