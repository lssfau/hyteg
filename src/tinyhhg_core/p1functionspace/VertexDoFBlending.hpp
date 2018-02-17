#pragma once

#include "core/debug/all.h"

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFaceIndex.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/polynomial/PolynomialEvaluator.hpp"

namespace hhg {
namespace vertexdof {
namespace blending {
namespace macroface {

template<typename ValueType, uint_t Level>
inline void assembleLocalStencil(uint_t i, uint_t j, const Matrix3r& localMatrix,
                                 double* opr_data,
                                 real_t coeff,
                                 const std::array<stencilDirection,3>& vertices,
                                 const std::array<uint_t,3>& idx)
{
  opr_data[vertexdof::stencilIndexFromVertex(vertices[0])] += coeff * localMatrix(idx[0],idx[0]);
  opr_data[vertexdof::stencilIndexFromVertex(vertices[1])] += coeff * localMatrix(idx[0],idx[1]);
  opr_data[vertexdof::stencilIndexFromVertex(vertices[2])] += coeff * localMatrix(idx[0],idx[2]);
}

template<typename ValueType, uint_t Level>
inline void applyBlendingTmpl(Face &face,
                                 const std::vector<PrimitiveDataID<FaceP1LocalMatrixMemory, Face>> &operatorIds,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Face> &srcId,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Face> &dstId,
                                 UpdateType update) {
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);

  Point3D x, x0;
  x0 = face.coords[0];
  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));

  Matrix2r mappingTensor;

  std::vector<FaceP1LocalMatrixMemory *> localMatricesVector;
  for (auto operatorId : operatorIds) {
    localMatricesVector.push_back(face.getData(operatorId));
  }

  ValueType tmp;

  std::array<SD, 3> triangleBlueSW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S};
  std::array<SD, 3> triangleGrayS = {SD::VERTEX_C, SD::VERTEX_S, SD::VERTEX_SE};
  std::array<SD, 3> triangleBlueSE = {SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E};
  std::array<SD, 3> triangleGrayNW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_NW};
  std::array<SD, 3> triangleBlueN = {SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N};
  std::array<SD, 3> triangleGrayNE = {SD::VERTEX_C, SD::VERTEX_N, SD::VERTEX_E};

  Point3D offsetSW = -1.0/3.0 * d0 - 1.0/3.0 * d2;
  Point3D offsetS = 1.0/3.0 * d0 - 2.0/3.0 * d2;
  Point3D offsetSE = 2.0/3.0 * d0 - 1.0/3.0 * d2;

  Point3D offsetNE = 1.0/3.0 * d0 + 1.0/3.0 * d2;
  Point3D offsetN = -1.0/3.0 * d0 + 2.0/3.0 * d2;
  Point3D offsetNW = -2.0/3.0 * d0 + 1.0/3.0 * d2;

  real_t* coeffs[] = { &mappingTensor(0,0), &mappingTensor(0,1), &mappingTensor(1,1) };
  std::vector<real_t> opr_data(vertexDoFMacroFaceStencilMemorySize(Level, 0));

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    x = x0;
    x += real_c(j)*d2 + d0;

    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      std::fill(opr_data.begin(), opr_data.end(), 0.0);

      face.blendingMap->evalTensorCoeff(x + offsetS, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleGrayS,
                                               {2, 0, 1});
      }

      face.blendingMap->evalTensorCoeff(x + offsetSE, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleBlueSE,
                                               {1, 2, 0});
      }

      face.blendingMap->evalTensorCoeff(x + offsetSW, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleBlueSW,
                                               {0, 1, 2});
      }

      face.blendingMap->evalTensorCoeff(x + offsetNW, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleGrayNW,
                                               {1, 0, 2});
      }

      face.blendingMap->evalTensorCoeff(x + offsetN, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleBlueN,
                                               {2, 1, 0});
      }

      face.blendingMap->evalTensorCoeff(x + offsetNE, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleGrayNE,
                                               {0, 2, 1});
      }

      auto stencil = PointND<real_t,7>(opr_data.data());

      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)];
      }

      tmp += opr_data[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)] = tmp;

      x += d0;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, applyBlendingTmpl, applyBlending)

}
}
}
}