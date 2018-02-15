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
namespace macroface {

using walberla::uint_t;
using walberla::real_c;

template<typename ValueType, uint_t Level>
inline ValueType assembleLocal(uint_t i, uint_t j, const Matrix3r& localMatrix,
                               double* src,
                               double* coeff,
                               const std::array<stencilDirection,3>& vertices,
                               const std::array<uint_t,3>& idx)
{
  ValueType meanCoeff = 1.0/3.0 * (coeff[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[0])]
                                 + coeff[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[1])]
                                 + coeff[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[2])]);

  ValueType tmp;
  tmp  = localMatrix(idx[0],idx[0]) * src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[0])]
         + localMatrix(idx[0],idx[1]) * src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[1])]
         + localMatrix(idx[0],idx[2]) * src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[2])];
  return meanCoeff * tmp;
}

template<typename ValueType, uint_t Level>
inline void assembleLocalStencil(uint_t i, uint_t j, const Matrix3r& localMatrix,
                                double* opr_data,
                                double* coeff,
                                const std::array<stencilDirection,3>& vertices,
                                const std::array<uint_t,3>& idx)
{
  ValueType meanCoeff = 1.0/3.0 * (coeff[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[0])]
                                   + coeff[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[1])]
                                   + coeff[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[2])]);

  opr_data[vertexdof::stencilIndexFromVertex(vertices[0])] += meanCoeff * localMatrix(idx[0],idx[0]);
  opr_data[vertexdof::stencilIndexFromVertex(vertices[1])] += meanCoeff * localMatrix(idx[0],idx[1]);
  opr_data[vertexdof::stencilIndexFromVertex(vertices[2])] += meanCoeff * localMatrix(idx[0],idx[2]);
}

template<typename ValueType, uint_t Level>
inline ValueType assembleLocalDG(uint_t i, uint_t j, const Matrix3r& localMatrix,
                               double* src,
                               const std::array<stencilDirection,3>& vertices,
                               const std::array<uint_t,3>& idx)
{
  ValueType tmp;
  tmp  = localMatrix(idx[0],idx[0]) * src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[0])]
      + localMatrix(idx[0],idx[1]) * src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[1])]
      + localMatrix(idx[0],idx[2]) * src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertices[2])];
  return tmp;
}

template<uint_t Level>
inline void fillLocalCoords(uint_t i, uint_t j, const P1Elements::P1Element& element, const std::array<real_t*, 2>& coords, real_t localCoords[6]) {

  localCoords[0] = coords[0][vertexdof::macroface::indexFromVertex<Level>(i, j, element[0])];
  localCoords[1] = coords[1][vertexdof::macroface::indexFromVertex<Level>(i, j, element[0])];
  localCoords[2] = coords[0][vertexdof::macroface::indexFromVertex<Level>(i, j, element[1])];
  localCoords[3] = coords[1][vertexdof::macroface::indexFromVertex<Level>(i, j, element[1])];
  localCoords[4] = coords[0][vertexdof::macroface::indexFromVertex<Level>(i, j, element[2])];
  localCoords[5] = coords[1][vertexdof::macroface::indexFromVertex<Level>(i, j, element[2])];
}

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Face &face,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face>& faceMemoryId,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &srcIds,
                            std::function<ValueType(const hhg::Point3D &, const std::vector<ValueType>&)> &expr) {

  FunctionMemory< ValueType > *faceMemory = face.getData(faceMemoryId);

  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }

  std::vector<ValueType> srcVector(srcIds.size());

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0, xBlend;

  auto dstPtr = faceMemory->getPointer( Level );

  x0 = face.coords[0];

  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));

  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    x = x0;
    x += real_c(i)*d2 + d0;

    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {

      for (uint_t k = 0; k < srcPtr.size(); ++k) {
        srcVector[k] = srcPtr[k][vertexdof::macroface::indexFromVertex<Level>(j, i, stencilDirection::VERTEX_C)];
      }

      face.blendingMap->evalF(x, xBlend);
      dstPtr[vertexdof::macroface::indexFromVertex<Level>(j, i, stencilDirection::VERTEX_C)] = expr(xBlend, srcVector);
      x += d0;
    }

    inner_rowsize -= 1;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, interpolateTmpl, interpolate)

template< typename ValueType, uint_t Level >
inline void assignTmpl(Face &face,
                   const std::vector<ValueType>& scalars,
                   const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &srcIds,
                   const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* dst = face.getData(dstId)->getPointer( Level );
  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }
  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      ValueType tmp = scalars[0]*srcPtr[0][vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

      for (uint_t k = 1; k < srcIds.size(); ++k) {
        tmp += scalars[k]*srcPtr[k][vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];
      }
      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, assignTmpl, assign)

template< typename ValueType, uint_t Level >
inline void addTmpl(Face &face,
                const std::vector<ValueType>& scalars,
                const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &srcIds,
                const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* dstPtr = face.getData(dstId)->getPointer( Level );
  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      ValueType tmp = 0.0;

      for (uint_t k = 0; k < srcIds.size(); ++k) {
        tmp += scalars[k] * srcPtr[k][vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];
      }

      dstPtr[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] += tmp;
    }

    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, addTmpl, add)

template< typename ValueType, uint_t Level >
inline real_t dotTmpl(Face &face,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Face>& lhsId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Face>& rhsId) {

  real_t sp = 0.0;
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* lhsPtr = face.getData(lhsId)->getPointer( Level );
  ValueType* rhsPtr = face.getData(rhsId)->getPointer( Level );

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      sp += lhsPtr[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)]
          * rhsPtr[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }

  return sp;
}

SPECIALIZE_WITH_VALUETYPE(real_t, dotTmpl, dot)

template< typename ValueType, uint_t Level >
inline void apply_tmpl(Face &face, const PrimitiveDataID< StencilMemory< ValueType >, Face>& operatorId,
                       const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                       const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId, UpdateType update) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* opr_data = face.getData(operatorId)->getPointer( Level );
  ValueType* src = face.getData(srcId)->getPointer( Level );
  ValueType* dst = face.getData(dstId)->getPointer( Level );

  ValueType tmp;

  if( update == Replace ) {
    for (uint_t j = 1; j < rowsize - 2; ++j) {
      for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
        tmp = opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

        //strangely the intel compiler cant handle this if it is a loop
        static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

        dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = tmp;
      }
      --inner_rowsize;
    }
  } else {
    for (uint_t j = 1; j < rowsize - 2; ++j) {
      for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
        tmp = opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

        //strangely the intel compiler cant handle this if it is a loop
        static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

        dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] += tmp;
      }
      --inner_rowsize;
    }
  }
}

SPECIALIZE_WITH_VALUETYPE(void, apply_tmpl, apply)

template< typename ValueType, uint_t Level >
inline void applyCoefficientTmpl(Face &face, const std::vector<PrimitiveDataID<FaceP1LocalMatrixMemory, Face>>& operatorIds,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                 const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &coeffIds,
                                 UpdateType update) {
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);

  std::vector<FaceP1LocalMatrixMemory*> localMatricesVector;
  for(auto operatorId : operatorIds) {
    localMatricesVector.push_back(face.getData(operatorId));
  }

  std::vector<real_t*> coeffs;
  for(auto coeffId : coeffIds) {
    coeffs.push_back(face.getData(coeffId)->getPointer( Level ));
  }

  ValueType tmp;

  std::array<SD,3> triangleBlueSW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_S  };
  std::array<SD,3> triangleGrayS  = { SD::VERTEX_C, SD::VERTEX_S,  SD::VERTEX_SE };
  std::array<SD,3> triangleBlueSE = { SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E  };
  std::array<SD,3> triangleGrayNW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_NW };
  std::array<SD,3> triangleBlueN  = { SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N  };
  std::array<SD,3> triangleGrayNE = { SD::VERTEX_C, SD::VERTEX_N,  SD::VERTEX_E  };

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {


      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)];
      }

      for (uint_t coeffIdx = 0; coeffIdx < coeffIds.size(); ++coeffIdx) {
        tmp += assembleLocal<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getGrayMatrix(Level), src, coeffs[coeffIdx], triangleGrayS,
                                               {2, 0, 1});
        tmp += assembleLocal<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getBlueMatrix(Level), src, coeffs[coeffIdx], triangleBlueSE,
                                               {1, 2, 0});
        tmp += assembleLocal<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getBlueMatrix(Level), src, coeffs[coeffIdx], triangleBlueSW,
                                               {0, 1, 2});
        tmp += assembleLocal<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getGrayMatrix(Level), src, coeffs[coeffIdx], triangleGrayNW,
                                               {1, 0, 2});
        tmp += assembleLocal<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getBlueMatrix(Level), src, coeffs[coeffIdx], triangleBlueN,
                                               {2, 1, 0});
        tmp += assembleLocal<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getGrayMatrix(Level), src, coeffs[coeffIdx], triangleGrayNE,
                                               {0, 2, 1});
      }

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, applyCoefficientTmpl, applyCoefficient)

template< typename ValueType, uint_t Level >
inline void applyCoefficientDGTmpl(Face &face, const PrimitiveDataID<FaceP1LocalMatrixMemory, Face>& operatorId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &coeffId,
                                 UpdateType update) {
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto localMatrices = face.getData(operatorId);
  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);
  auto coeff = face.getData(coeffId)->getPointer(Level);

  ValueType tmp;

  std::array<SD,3> triangleBlueSW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_S  };
  std::array<SD,3> triangleGraySE = { SD::VERTEX_C, SD::VERTEX_S,  SD::VERTEX_SE };
  std::array<SD,3> triangleBlueSE = { SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E  };
  std::array<SD,3> triangleGrayNW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_NW };
  std::array<SD,3> triangleBlueNW = { SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N  };
  std::array<SD,3> triangleGrayNE = { SD::VERTEX_C, SD::VERTEX_N,  SD::VERTEX_E  };



  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {


      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)];
      }

      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_GRAY_SE)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getGrayMatrix(Level), src, triangleGraySE, {2,0,1});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_BLUE_SE)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getBlueMatrix(Level), src, triangleBlueSE, {1,2,0});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_BLUE_SW)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getBlueMatrix(Level), src, triangleBlueSW, {0,1,2});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_GRAY_NW)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getGrayMatrix(Level), src, triangleGrayNW, {1,0,2});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_BLUE_NW)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getBlueMatrix(Level), src, triangleBlueNW, {2,1,0});
      tmp += coeff[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_GRAY_NE)] * assembleLocalDG<ValueType, Level>(i, j, localMatrices->getGrayMatrix(Level), src, triangleGrayNE, {0,2,1});

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, applyCoefficientDGTmpl, applyCoefficientDG)

template< typename ValueType, uint_t Level >
inline void applyElementwiseTmpl(Face &face, std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                 std::array<const PrimitiveDataID<FunctionMemory< ValueType >, Face>, 2> &coordIds,
                                 UpdateType update) {
  using namespace P1Elements;
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);
  std::array<ValueType*, 2> globalCoords{{face.getData(coordIds[0])->getPointer(Level), face.getData(coordIds[1])->getPointer(Level)}};

  ValueType tmp;
  real_t localCoords[6];
  Matrix3r localStiffness;
  std::vector<real_t> faceStencil(vertexDoFMacroFaceStencilMemorySize(Level, 0));

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      std::fill(faceStencil.begin(), faceStencil.end(), walberla::real_c(0.0));

      for (uint_t k = 0; k < FaceVertexDoF::P1GrayElements.size(); ++k) {

        // fill local coords
        fillLocalCoords<Level>(i, j, FaceVertexDoF::P1GrayElements[k], globalCoords, localCoords);

        // compute stencil
        computeElementMatrix(localStiffness, localCoords);
        assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[k], {{0,1,2}}, localStiffness, faceStencil);
      }

      for (uint_t k = 0; k < FaceVertexDoF::P1BlueElements.size(); ++k) {

        // fill local coords
        fillLocalCoords<Level>(i, j, FaceVertexDoF::P1BlueElements[k], globalCoords, localCoords);

        // fill coords
        computeElementMatrix(localStiffness, localCoords);
        assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[k], {{0,1,2}}, localStiffness, faceStencil);
      }

//      WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("FACE.id = {}:faceStencil = {}", face.getID().getID(), PointND<real_t, 7>(&faceStencil[0])));

      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)];
      }

      tmp += faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)];

      //strangely the intel compiler cant handle this if it is a loop
      static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );
      tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
      tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
      tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
      tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
      tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
      tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, applyElementwiseTmpl, applyElementwise)

template<typename ValueType, uint_t Level, uint_t PolyDegree>
inline void applyPolynomialTmpl(Face &face, const PrimitiveDataID<FaceP1PolynomialMemory, Face>& polynomialId,
                                const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId, UpdateType update) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto polynomials = face.getData(polynomialId);
  ValueType* src = face.getData(srcId)->getPointer( Level );
  ValueType* dst = face.getData(dstId)->getPointer( Level );

  ValueType tmp;
  std::vector<real_t> faceStencil(vertexDoFMacroFaceStencilMemorySize(Level, 0));
  Point2D x;
  real_t h = real_c(1.0) / real_c(rowsize-1);

  auto horiPoly = polynomials->getHoriPolynomial(PolyDegree);
  auto vertPoly = polynomials->getVertPolynomial(PolyDegree);
  auto diagPoly = polynomials->getDiagPolynomial(PolyDegree);

  Polynomial2DEvaluator evalHoriPoly(horiPoly);
  Polynomial2DEvaluator evalVertPolyS(vertPoly);
  Polynomial2DEvaluator evalVertPolyN(vertPoly);
  Polynomial2DEvaluator evalDiagPolySE(diagPoly);
  Polynomial2DEvaluator evalDiagPolyNW(diagPoly);

  if( update == Replace ) {
    for (uint_t j = 1; j < rowsize - 2; ++j) {
      x[1] = j * h;

      // Set new Y values
      evalHoriPoly.setY(x[1]);
      evalVertPolyS.setY(x[1] - 0.5 * h);
      evalVertPolyN.setY(x[1] + 0.5 * h);
      evalDiagPolySE.setY(x[1] - 0.5 * h);
      evalDiagPolyNW.setY(x[1] + 0.5 * h);

      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = evalHoriPoly.setStartX<PolyDegree>(-0.5 * h, h);
      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval<PolyDegree>();

      evalVertPolyS.setStartX<PolyDegree>(0.0, h);
      evalVertPolyN.setStartX<PolyDegree>(0.0, h);

      evalDiagPolySE.setStartX<PolyDegree>(0.5 * h, h);
      evalDiagPolyNW.setStartX<PolyDegree>(-0.5 * h, h);

      for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)];
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval<PolyDegree>();

        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)] = evalVertPolyS.incrementEval<PolyDegree>();
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)] = evalVertPolyN.incrementEval<PolyDegree>();

        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)] = evalDiagPolySE.incrementEval<PolyDegree>();
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)] = evalDiagPolyNW.incrementEval<PolyDegree>();

        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]
                        - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]
                        - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]
                        - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]
                        - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]
                        - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])];

        tmp = faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

        //strangely the intel compiler cant handle this if it is a loop
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

        dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = tmp;
      }
      --inner_rowsize;
    }
  } else {
    for (uint_t j = 1; j < rowsize - 2; ++j) {
      x[1] = j * h;

      // Set new Y values
      evalHoriPoly.setY(x[1]);
      evalVertPolyS.setY(x[1] - 0.5 * h);
      evalVertPolyN.setY(x[1] + 0.5 * h);
      evalDiagPolySE.setY(x[1] - 0.5 * h);
      evalDiagPolyNW.setY(x[1] + 0.5 * h);

      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = evalHoriPoly.setStartX<PolyDegree>(-0.5 * h, h);
      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval<PolyDegree>();

      evalVertPolyS.setStartX<PolyDegree>(0.0, h);
      evalVertPolyN.setStartX<PolyDegree>(0.0, h);

      evalDiagPolySE.setStartX<PolyDegree>(0.5 * h, h);
      evalDiagPolyNW.setStartX<PolyDegree>(-0.5 * h, h);

      for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)];
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval<PolyDegree>();

        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)] = evalVertPolyS.incrementEval<PolyDegree>();
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)] = evalVertPolyN.incrementEval<PolyDegree>();

        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)] = evalDiagPolySE.incrementEval<PolyDegree>();
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)] = evalDiagPolyNW.incrementEval<PolyDegree>();

        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]
                                                                                     - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]
                                                                                     - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]
                                                                                     - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]
                                                                                     - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]
                                                                                     - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])];

        tmp = faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

        //strangely the intel compiler cant handle this if it is a loop
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
        tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

        dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] += tmp;
      }
      --inner_rowsize;
    }
  }
}

SPECIALIZE_POLYNOMIAL(void, applyPolynomialTmpl, applyPolynomial)

template< typename ValueType, uint_t Level >
inline void smooth_gs_tmpl(Face &face, const PrimitiveDataID<StencilMemory< ValueType >, Face>& operatorId,
                           const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                           const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData(operatorId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp = rhs[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

      //for (auto neighbor : neighbors) {
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        tmp -= opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])]*dst[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[k])];
      }

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = tmp/opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, smooth_gs_tmpl, smooth_gs)

template< typename ValueType, uint_t Level >
inline void smooth_gs_coefficient_tmpl(Face &face,
                                       const std::vector<PrimitiveDataID<FaceP1LocalMatrixMemory, Face>>& operatorIds,
                                       const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                       const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId,
                                       const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &coeffIds) {

  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto dst = face.getData(dstId)->getPointer(Level);
  auto rhs = face.getData(rhsId)->getPointer( Level );

  std::vector<FaceP1LocalMatrixMemory*> localMatricesVector;
  for(auto operatorId : operatorIds) {
    localMatricesVector.push_back(face.getData(operatorId));
  }

  std::vector<real_t*> coeffs;
  for(auto coeffId : coeffIds) {
    coeffs.push_back(face.getData(coeffId)->getPointer( Level ));
  }

  std::array<SD,3> triangleBlueSW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_S  };
  std::array<SD,3> triangleGrayS  = { SD::VERTEX_C, SD::VERTEX_S,  SD::VERTEX_SE };
  std::array<SD,3> triangleBlueSE = { SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E  };
  std::array<SD,3> triangleGrayNW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_NW };
  std::array<SD,3> triangleBlueN  = { SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N  };
  std::array<SD,3> triangleGrayNE = { SD::VERTEX_C, SD::VERTEX_N,  SD::VERTEX_E  };

  ValueType tmp;
  std::vector<real_t> opr_data(vertexDoFMacroFaceStencilMemorySize(Level, 0));

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      std::fill(opr_data.begin(), opr_data.end(), 0.0);

      for (uint_t coeffIdx = 0; coeffIdx < coeffIds.size(); ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getGrayMatrix(Level), opr_data.data(), coeffs[coeffIdx],
                                               triangleGrayS, {2, 0, 1});
        assembleLocalStencil<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getBlueMatrix(Level), opr_data.data(), coeffs[coeffIdx],
                                               triangleBlueSE, {1, 2, 0});
        assembleLocalStencil<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getBlueMatrix(Level), opr_data.data(), coeffs[coeffIdx],
                                               triangleBlueSW, {0, 1, 2});
        assembleLocalStencil<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getGrayMatrix(Level), opr_data.data(), coeffs[coeffIdx],
                                               triangleGrayNW, {1, 0, 2});
        assembleLocalStencil<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getBlueMatrix(Level), opr_data.data(), coeffs[coeffIdx],
                                               triangleBlueN, {2, 1, 0});
        assembleLocalStencil<ValueType, Level>(i, j, localMatricesVector[coeffIdx]->getGrayMatrix(Level), opr_data.data(), coeffs[coeffIdx],
                                               triangleGrayNE, {0, 2, 1});
      }

      tmp = rhs[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

      //for (auto neighbor : neighbors) {
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        tmp -= opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])]*dst[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[k])];
      }

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = tmp/opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, smooth_gs_coefficient_tmpl, smooth_gs_coefficient)

template<typename ValueType, uint_t Level, uint_t PolyDegree>
inline void smooth_gs_polynomial_tmpl(Face &face, const PrimitiveDataID<FaceP1PolynomialMemory, Face>& polynomialId,
                                      const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                      const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto polynomials = face.getData(polynomialId);
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );

  std::vector<real_t> opr_data(vertexDoFMacroFaceStencilMemorySize(Level, 0));
  Point2D x;
  real_t h = real_c(1.0) / real_c(rowsize-1);

  auto horiPoly = polynomials->getHoriPolynomial(PolyDegree);
  auto vertPoly = polynomials->getVertPolynomial(PolyDegree);
  auto diagPoly = polynomials->getDiagPolynomial(PolyDegree);

  Polynomial2DEvaluator evalHoriPoly(horiPoly);
  Polynomial2DEvaluator evalVertPolyS(vertPoly);
  Polynomial2DEvaluator evalVertPolyN(vertPoly);
  Polynomial2DEvaluator evalDiagPolySE(diagPoly);
  Polynomial2DEvaluator evalDiagPolyNW(diagPoly);

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    x[1] = j * h;

    // Set new Y values
    evalHoriPoly.setY(x[1]);
    evalVertPolyS.setY(x[1] - 0.5 * h);
    evalVertPolyN.setY(x[1] + 0.5 * h);
    evalDiagPolySE.setY(x[1] - 0.5 * h);
    evalDiagPolyNW.setY(x[1] + 0.5 * h);

    opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = evalHoriPoly.setStartX<PolyDegree>(-0.5 * h, h);
    opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval<PolyDegree>();

    evalVertPolyS.setStartX<PolyDegree>(0.0, h);
    evalVertPolyN.setStartX<PolyDegree>(0.0, h);

    evalDiagPolySE.setStartX<PolyDegree>(0.5 * h, h);
    evalDiagPolyNW.setStartX<PolyDegree>(-0.5 * h, h);

    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)];
      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval<PolyDegree>();

      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)] = evalVertPolyS.incrementEval<PolyDegree>();
      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)] = evalVertPolyN.incrementEval<PolyDegree>();

      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)] = evalDiagPolySE.incrementEval<PolyDegree>();
      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)] = evalDiagPolyNW.incrementEval<PolyDegree>();

      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = - opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]
                                                                                   - opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]
                                                                                   - opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]
                                                                                   - opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]
                                                                                   - opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]
                                                                                   - opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])];

      tmp = rhs[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

      //for (auto neighbor : neighbors) {
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        tmp -= opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])]*dst[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[k])];
      }

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = tmp/opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }
}

SPECIALIZE_POLYNOMIAL(void, smooth_gs_polynomial_tmpl, smooth_gs_polynomial)

template< typename ValueType, uint_t Level >
inline void smooth_sor_tmpl(Face &face, const PrimitiveDataID<StencilMemory< ValueType >, Face>& operatorId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId,
                            ValueType relax) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData(operatorId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp = rhs[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

      //for (auto neighbor : neighbors) {
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        tmp -= opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])]*dst[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[k])];
      }

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = (1.0-relax) * dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] + relax * tmp/opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, smooth_sor_tmpl, smooth_sor)

template< typename ValueType, uint_t Level >
inline void smooth_jac_tmpl(Face &face, const PrimitiveDataID<StencilMemory< ValueType >, Face>& operatorId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &tmpId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData(operatorId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );
  auto tmpVar = face.getData(tmpId)->getPointer( Level );

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp = rhs[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];

      for (auto neighbor : vertexdof::macroface::neighborsWithoutCenter ) {
        tmp -= opr_data[vertexdof::stencilIndexFromVertex(neighbor)]*tmpVar[vertexdof::macroface::indexFromVertex<Level>(i, j, neighbor)];
      }

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = tmp/opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, smooth_jac_tmpl, smooth_jac)

template< typename ValueType, uint_t Level >
inline void prolongate_tmpl(Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face>& memoryId) {

  typedef stencilDirection SD;
  using namespace vertexdof::macroface;

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level);
  uint_t N_c_i = N_c;

  auto v_f = face.getData(memoryId)->getPointer( Level + 1 );
  auto v_c = face.getData(memoryId)->getPointer( Level );

  uint_t j;

  for (uint_t i = 1; i < N_c - 1; ++i) {
    for (j = 1; j < N_c_i - 2; ++j) {
      v_f[indexFromVertex<Level + 1>(2*i, 2*j, SD::VERTEX_C)] = v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)];
      v_f[indexFromVertex<Level + 1>(2*i - 1, 2*j - 1, SD::VERTEX_C)] =
          0.5*(v_c[indexFromVertex<Level>(i - 1, j, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j - 1, SD::VERTEX_C)]);
      v_f[indexFromVertex<Level + 1>(2*i - 1, 2*j, SD::VERTEX_C)] = 0.5*(v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i - 1, j, SD::VERTEX_C)]);
      v_f[indexFromVertex<Level + 1>(2*i, 2*j - 1, SD::VERTEX_C)] = 0.5*(v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j - 1, SD::VERTEX_C)]);
    }

    v_f[indexFromVertex<Level + 1>(2*i - 1, 2*j - 1, SD::VERTEX_C)] = 0.5*(v_c[indexFromVertex<Level>(i - 1, j, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j - 1, SD::VERTEX_C)]);
    v_f[indexFromVertex<Level + 1>(2*i - 1, 2*j, SD::VERTEX_C)] = 0.5*(v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i - 1, j, SD::VERTEX_C)]);
    v_f[indexFromVertex<Level + 1>(2*i, 2*j - 1, SD::VERTEX_C)] = 0.5*(v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j - 1, SD::VERTEX_C)]);

    --N_c_i;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, prolongate_tmpl, prolongate)

template< typename ValueType, uint_t Level >
inline void prolongateQuadratic_tmpl(Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face>& memoryId) {

  typedef stencilDirection SD;
  using namespace vertexdof::macroface;

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
    linearx = 0.5*(v_c[indexFromVertex<Level>(i, j - 2, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)]);
    lineary = 0.5*(v_c[indexFromVertex<Level>(i + 2, j - 2, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j - 2, SD::VERTEX_C)]);
    linearxy = 0.5*(v_c[indexFromVertex<Level>(i + 2, j - 2, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)]);

    offx = v_c[indexFromVertex<Level>(i, j - 1, SD::VERTEX_C)] - linearx;
    offy = v_c[indexFromVertex<Level>(i + 1, j - 2, SD::VERTEX_C)] - lineary;
    offxy = v_c[indexFromVertex<Level>(i + 1, j - 1, SD::VERTEX_C)] - linearxy;

// left bottom corner
    v_f[indexFromVertex<Level + 1>(2*i + 1, 2*j - 3, SD::VERTEX_C)] = 0.5*(linearx + lineary) + 0.5*offx + 0.5*offy + 0.25*offxy;
// right bottom corner
    v_f[indexFromVertex<Level + 1>(2*i + 1, 2*j - 2, SD::VERTEX_C)] = 0.5*(linearx + linearxy) + 0.5*offx + 0.25*offy + 0.5*offxy;
// top corner
    v_f[indexFromVertex<Level + 1>(2*i + 2, 2*j - 3, SD::VERTEX_C)] = 0.5*(linearxy + lineary) + 0.25*offx + 0.5*offy + 0.5*offxy;
  }

  N_c_i -= 1;

  for (j = 2; j < N_c - 1; j += 2) {
    for (i = 2; i < N_c_i - 1; i += 2) {

// upper triangle inner points
//calculate offsets
      linearx = 0.5*(v_c[indexFromVertex<Level>(i, j - 2, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)]);
      lineary = 0.5*(v_c[indexFromVertex<Level>(i + 2, j - 2, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j - 2, SD::VERTEX_C)]);
      linearxy = 0.5*(v_c[indexFromVertex<Level>(i + 2, j - 2, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)]);

      offx = v_c[indexFromVertex<Level>(i, j - 1, SD::VERTEX_C)] - linearx;
      offy = v_c[indexFromVertex<Level>(i + 1, j - 2, SD::VERTEX_C)] - lineary;
      offxy = v_c[indexFromVertex<Level>(i + 1, j - 1, SD::VERTEX_C)] - linearxy;
// left bottom corner
      v_f[indexFromVertex<Level + 1>(2*i + 1, 2*j - 3, SD::VERTEX_C)] = 0.5*(linearx + lineary) + 0.5*offx + 0.5*offy + 0.25*offxy;
// right bottom corner
      v_f[indexFromVertex<Level + 1>(2*i + 1, 2*j - 2, SD::VERTEX_C)] = 0.5*(linearx + linearxy) + 0.5*offx + 0.25*offy + 0.5*offxy;
// top corner
      v_f[indexFromVertex<Level + 1>(2*i + 2, 2*j - 3, SD::VERTEX_C)] = 0.5*(linearxy + lineary) + 0.25*offx + 0.5*offy + 0.5*offxy;

// lower triangle all points
//calculate offsets
      lineary = 0.5*(v_c[indexFromVertex<Level>(i - 2, j, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)]);
      linearxy = 0.5*(v_c[indexFromVertex<Level>(i - 2, j, SD::VERTEX_C)] + v_c[indexFromVertex<Level>(i, j - 2, SD::VERTEX_C)]);

      offy = v_c[indexFromVertex<Level>(i - 1, j, SD::VERTEX_C)] - lineary;
      offxy = v_c[indexFromVertex<Level>(i - 1, j - 1, SD::VERTEX_C)] - linearxy;
// first inner points
// left bottom corner
      v_f[indexFromVertex<Level + 1>(2*i - 1, 2*j - 1, SD::VERTEX_C)] = 0.5*(linearx + lineary) + 0.5*offx + 0.5*offy + 0.25*offxy;
// right bottom corner
      v_f[indexFromVertex<Level + 1>(2*i - 1, 2*j - 2, SD::VERTEX_C)] = 0.5*(linearx + linearxy) + 0.5*offx + 0.25*offy + 0.5*offxy;
// top corner
      v_f[indexFromVertex<Level + 1>(2*i - 2, 2*j - 1, SD::VERTEX_C)] = 0.5*(linearxy + lineary) + 0.25*offx + 0.5*offy + 0.5*offxy;

// boundary points
// x-direction
      v_f[indexFromVertex<Level + 1>(2*i, 2*j - 1, SD::VERTEX_C)] = 0.5*(linearx + v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)]) + 0.75*offx;
      v_f[indexFromVertex<Level + 1>(2*i, 2*j - 3, SD::VERTEX_C)] = 0.5*(linearx + v_c[indexFromVertex<Level>(i, j - 2, SD::VERTEX_C)]) + 0.75*offx;
//y-direction
      v_f[indexFromVertex<Level + 1>(2*i - 1, 2*j, SD::VERTEX_C)] = 0.5*(v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + lineary) + 0.75*offy;
      v_f[indexFromVertex<Level + 1>(2*i - 3, 2*j, SD::VERTEX_C)] = 0.5*(v_c[indexFromVertex<Level>(i - 2, j, SD::VERTEX_C)] + lineary) + 0.75*offy;
//xy-direction
      v_f[indexFromVertex<Level + 1>(2*i - 1, 2*j - 3, SD::VERTEX_C)] = 0.5*(v_c[indexFromVertex<Level>(i, j - 2, SD::VERTEX_C)] + linearxy) + 0.75*offxy;
      v_f[indexFromVertex<Level + 1>(2*i - 3, 2*j - 1, SD::VERTEX_C)] = 0.5*(v_c[indexFromVertex<Level>(i - 2, j, SD::VERTEX_C)] + linearxy) + 0.75*offxy;
// coarse points
      v_f[indexFromVertex<Level + 1>(2*i, 2*j, SD::VERTEX_C)] = v_c[indexFromVertex<Level>(i, j, SD::VERTEX_C)];
      v_f[indexFromVertex<Level + 1>(2*i, 2*j - 2, SD::VERTEX_C)] = v_c[indexFromVertex<Level>(i, j - 1, SD::VERTEX_C)];
      v_f[indexFromVertex<Level + 1>(2*i - 2, 2*j, SD::VERTEX_C)] = v_c[indexFromVertex<Level>(i - 1, j, SD::VERTEX_C)];
      v_f[indexFromVertex<Level + 1>(2*i - 2, 2*j - 2, SD::VERTEX_C)] = v_c[indexFromVertex<Level>(i - 1, j - 1, SD::VERTEX_C)];
    }
    N_c_i -= 2;

  }
}

SPECIALIZE_WITH_VALUETYPE(void, prolongateQuadratic_tmpl, prolongateQuadratic)

template< typename ValueType, uint_t Level >
inline void restrict_tmpl(Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face> &memoryId) {

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level - 1);
  uint_t N_c_i = N_c;

  auto v_f = face.getData(memoryId)->getPointer( Level );
  auto v_c = face.getData(memoryId)->getPointer( Level - 1 );

  ValueType tmp;

  for (uint_t j = 1; j < N_c - 2; ++j) {
    for (uint_t i = 1; i < N_c_i - 2; ++i) {

      tmp = v_f[vertexdof::macroface::indexFromVertex<Level>(2*i, 2*j, stencilDirection::VERTEX_C)];

      for ( const auto & neighbor : vertexdof::macroface::neighborsWithoutCenter ) {
        tmp += 0.5*v_f[vertexdof::macroface::indexFromVertex<Level>(2*i, 2*j, neighbor)];
      }

      v_c[vertexdof::macroface::indexFromVertex<Level - 1>(i, j, stencilDirection::VERTEX_C)] = tmp;
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
inline void enumerateTmpl(Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId, uint_t& num) {
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
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsP1Id,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId) {
  using namespace vertexdof::macroface;
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto rhs = face.getData(rhsId)->getPointer(Level);
  auto rhsP1 = face.getData(rhsP1Id)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);

  real_t faceArea = std::pow(4.0, -walberla::real_c(Level)) * face.area;
  real_t weightedFaceArea = faceArea / 3.0;

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp  = rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_BLUE_SW)] * (0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_W)]) + 0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_S)]));
      tmp += rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_GRAY_SE)] * (0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_S)]) + 0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_SE)]));
      tmp += rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_BLUE_SE)] * (0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_SE)]) + 0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_E)]));

      tmp += rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_GRAY_NW)] * (0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_W)]) + 0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_NW)]));
      tmp += rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_BLUE_NW)] * (0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_NW)]) + 0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_N)]));
      tmp += rhs[DGFace::indexDGFaceFromVertex<Level>(i, j, SD::CELL_GRAY_NE)] * (0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_N)]) + 0.5 * 0.5 * (rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_C)] + rhsP1[indexFromVertex<Level>(i, j, SD::VERTEX_E)]));

      dst[indexFromVertex<Level>(i, j, SD::VERTEX_C)] = weightedFaceArea * tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, integrateDGTmpl, integrateDG )

template< typename ValueType, uint_t Level >
inline real_t getMaxValueTmpl(Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer( Level );
  real_t localMax = std::numeric_limits<real_t>::min();

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
      localMax = std::max(localMax, src[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)]);
    }
    --inner_rowsize;
  }

  return localMax;
}

SPECIALIZE_WITH_VALUETYPE( real_t, getMaxValueTmpl, getMaxValue )

#ifdef HHG_BUILD_WITH_PETSC
template< uint_t Level >
inline void saveOperator_tmpl(Face &face, const PrimitiveDataID<StencilMemory< real_t >, Face>& operatorId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &srcId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &dstId, Mat& mat) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData(operatorId)->getPointer( Level );
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt srcInt = src[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];
      PetscInt dstInt = dst[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, VERTEX_C)], opr_data[VERTEX_C]);
      MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] ,INSERT_VALUES);

      for ( const auto & neighbor : vertexdof::macroface::neighborsWithoutCenter ) {
        srcInt = src[vertexdof::macroface::indexFromVertex<Level>(i, j, neighbor)];
        //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, neighbor)], opr_data[neighbor]);
        MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[vertexdof::stencilIndexFromVertex(neighbor)] ,INSERT_VALUES);
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)


template< typename ValueType, uint_t Level >
inline void createVectorFromFunctionTmpl(Face &face,
                              const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &numeratorId,
                              Vec& vec) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt numeratorInt = numerator[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];
      VecSetValues(vec,1,&numeratorInt,&src[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)],INSERT_VALUES);
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, createVectorFromFunctionTmpl, createVectorFromFunction)



template< typename ValueType, uint_t Level >
inline void createFunctionFromVectorTmpl(Face &face,
                                         const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                         const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &numeratorId,
                                         Vec& vec) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt numeratorInt = numerator[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)];
      VecGetValues(vec,1,&numeratorInt,&src[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)]);
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, createFunctionFromVectorTmpl, createFunctionFromVector)
#endif

template< typename ValueType, size_t Level >
inline void printFunctionMemory(Face& face, const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId){
  ValueType* faceMemory = face.getData(dstId)->getPointer( Level );
  using namespace std;
  cout << setfill('=') << setw(100) << "" << endl;
  cout << face << std::left << setprecision(1) << fixed << setfill(' ') << endl;
  for ( const auto & it : vertexdof::macroface::Iterator( Level, 0 ) ){
    if(it.col() == 0) std::cout << std::endl;
    cout << setw(5) << faceMemory[hhg::vertexdof::macroface::indexFromVertex< Level >(it.col(),it.row(), stencilDirection::VERTEX_C)] << "|";
  }
  cout << endl << setfill('=') << setw(100) << "" << endl << setfill(' ');

}

}// namespace macroface
}// namespace vertexdof
}// namespace hhg

