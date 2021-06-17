/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "core/debug/all.h"

namespace hyteg {
enum class OperatorType {
  MASS,
  EVEN,
  ODD
};
}

#include "hyteg/primitives/Face.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/Macros.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/p1functionspace/P1Elements.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/polynomial/PolynomialEvaluator.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {

using walberla::uint_t;
using walberla::real_c;
using indexing::Index;

template<typename ValueType, OperatorType OprType, uint_t PolyDegree>
inline void applyPolynomialTmpl(uint_t Level, Face &face, const PrimitiveDataID<FaceP1PolynomialMemory, Face> polynomialId,
                                const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId, UpdateType update) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto polynomials = face.getData(polynomialId);
  ValueType* src = face.getData(srcId)->getPointer( Level );
  ValueType* dst = face.getData(dstId)->getPointer( Level );

  ValueType tmp;
  std::vector<real_t> faceStencil(7);
  Point2D x;
  real_t h = real_c(1.0) / real_c(rowsize-1);

  auto centerPoly = polynomials->getPolynomialC(PolyDegree);
  auto horiPoly = polynomials->getPolynomialW(PolyDegree);
  auto vertPoly = polynomials->getPolynomialS(PolyDegree);
  auto diagPoly = polynomials->getPolynomialSE(PolyDegree);

  Polynomial2DEvaluator evalCenterPoly(centerPoly);
  Polynomial2DEvaluator evalHoriPoly(horiPoly);
  Polynomial2DEvaluator evalVertPolyS(vertPoly);
  Polynomial2DEvaluator evalVertPolyN(vertPoly);
  Polynomial2DEvaluator evalDiagPolySE(diagPoly);
  Polynomial2DEvaluator evalDiagPolyNW(diagPoly);

 for (uint_t j = 1; j < rowsize - 2; ++j) {
   x[1] = j * h;

   // Set new Y values
   if (OprType != OperatorType::EVEN) {
      evalCenterPoly.setY(x[1]);
   }
   evalHoriPoly.setY(x[1]);
   evalVertPolyS.setY(x[1] - 0.5 * h);
   evalVertPolyN.setY(x[1] + 0.5 * h);
   evalDiagPolySE.setY(x[1] - 0.5 * h);
   evalDiagPolyNW.setY(x[1] + 0.5 * h);

   faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = evalHoriPoly.setStartX(-0.5 * h, h);
   faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval();

   if (OprType != OperatorType::EVEN) {
      evalCenterPoly.setStartX(0.0, h);
   }
   evalVertPolyS.setStartX(0.0, h);
   evalVertPolyN.setStartX(0.0, h);

   evalDiagPolySE.setStartX(0.5 * h, h);
   evalDiagPolyNW.setStartX(-0.5 * h, h);

   for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
     faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)];
     faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval();

     faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)] = evalVertPolyS.incrementEval();
     faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)] = evalVertPolyN.incrementEval();

     faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)] = evalDiagPolySE.incrementEval();
     faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)] = evalDiagPolyNW.incrementEval();

//     faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = evalCenterPoly.incrementEval();

     if (OprType == OperatorType::MASS) {
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = evalCenterPoly.incrementEval();
     } else if (OprType == OperatorType::EVEN) {
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = -faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)]
                                                                                     -faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)]
                                                                                     -faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)]
                                                                                     -faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)]
                                                                                     -faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)]
                                                                                     -faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)];
     }

     // if (i == 1 && j == 1) {
     //    PointND<real_t, 7> test(faceStencil.data());
     //    WALBERLA_LOG_INFO("stencil = " << test);
     // }

     tmp = faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex(Level, i, j, stencilDirection::VERTEX_C)];

     //strangely the intel compiler cant handle this if it is a loop
     tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
     tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
     tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
     tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
     tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
     tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

     if (update == Replace) {
        dst[vertexdof::macroface::indexFromVertex(Level, i, j, stencilDirection::VERTEX_C)] = tmp;
     } else {
        dst[vertexdof::macroface::indexFromVertex(Level, i, j, stencilDirection::VERTEX_C)] += tmp;
     }
   }
   --inner_rowsize;
 }
}

SPECIALIZE_OPRTYPE_POLYNOMIAL(void, applyPolynomialTmpl, applyPolynomial)

template<typename ValueType, uint_t PolyDegree>
inline void applyPolynomialOddTmpl(uint_t Level, Face &face, const PrimitiveDataID<FaceP1PolynomialMemory, Face>& polynomialId,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                    const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId, UpdateType update) {

   uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
   uint_t inner_rowsize = rowsize;

   auto polynomials = face.getData(polynomialId);
   ValueType* src = face.getData(srcId)->getPointer( Level );
   ValueType* dst = face.getData(dstId)->getPointer( Level );

   ValueType tmp;
   std::vector<real_t> faceStencil(7);
   Point2D x;
   real_t h = real_c(1.0) / real_c(rowsize-1);

   auto polyS = polynomials->getPolynomialS(PolyDegree);
   auto polySE = polynomials->getPolynomialSE(PolyDegree);
   auto polyW = polynomials->getPolynomialW(PolyDegree);
//   auto polyC = polynomials->getPolynomialC(PolyDegree);
   auto polyE = polynomials->getPolynomialE(PolyDegree);
   auto polyNW = polynomials->getPolynomialNW(PolyDegree);
   auto polyN = polynomials->getPolynomialN(PolyDegree);

   Polynomial2DEvaluator evalPolyS(polyS);
   Polynomial2DEvaluator evalPolySE(polySE);
   Polynomial2DEvaluator evalPolyW(polyW);
//   Polynomial2DEvaluator evalPolyC(polyC);
   Polynomial2DEvaluator evalPolyE(polyE);
   Polynomial2DEvaluator evalPolyNW(polyNW);
   Polynomial2DEvaluator evalPolyN(polyN);

   for (uint_t j = 1; j < rowsize - 2; ++j) {
      x[1] = j * h;

      // Set new Y values
      evalPolyS.setY(x[1]);
      evalPolySE.setY(x[1]);
      evalPolyW.setY(x[1]);
//      evalPolyC.setY(x[1]);
      evalPolyE.setY(x[1]);
      evalPolyNW.setY(x[1]);
      evalPolyN.setY(x[1]);

      evalPolyS.setStartX(0.0, h);
      evalPolySE.setStartX(0.0, h);
      evalPolyW.setStartX(0.0, h);
//      evalPolyC.setStartX(0.0, h);
      evalPolyE.setStartX(0.0, h);
      evalPolyNW.setStartX(0.0, h);
      evalPolyN.setStartX(0.0, h);

      for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = evalPolyW.incrementEval();
         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalPolyE.incrementEval();

         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)] = evalPolyS.incrementEval();
         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)] = evalPolyN.incrementEval();

         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)] = evalPolySE.incrementEval();
         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)] = evalPolyNW.incrementEval();

//         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = evalPolyC.incrementEval();

//         if (i == 1 && j == 1) {
//            PointND<real_t, 7> test(faceStencil.data());
//            WALBERLA_LOG_INFO("stencil = " << test);
//         }

         tmp = faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex(Level, i, j, stencilDirection::VERTEX_C)];

         //strangely the intel compiler cant handle this if it is a loop
         tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
         tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
         tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
         tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
         tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
         tmp += faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

         if (update == Replace) {
            dst[vertexdof::macroface::indexFromVertex(Level, i, j, stencilDirection::VERTEX_C)] = tmp;
         } else {
            dst[vertexdof::macroface::indexFromVertex(Level, i, j, stencilDirection::VERTEX_C)] += tmp;
         }
      }
      --inner_rowsize;
   }
}

SPECIALIZE_POLYNOMIAL(void, applyPolynomialOddTmpl, applyPolynomialOdd)

template<typename ValueType, uint_t PolyDegree>
inline void smooth_gs_polynomial_even_tmpl(uint_t Level, Face &face, const PrimitiveDataID<FaceP1PolynomialMemory, Face>& polynomialId,
                                           const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                           const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto polynomials = face.getData(polynomialId);
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );

  std::vector<real_t> opr_data(7);
  Point2D x;
  real_t h = real_c(1.0) / real_c(rowsize-1);

//  auto centerPoly = polynomials->getPolynomialC(PolyDegree);
  auto horiPoly = polynomials->getPolynomialW(PolyDegree);
  auto vertPoly = polynomials->getPolynomialS(PolyDegree);
  auto diagPoly = polynomials->getPolynomialSE(PolyDegree);

//  Polynomial2DEvaluator evalCenterPoly(centerPoly);
  Polynomial2DEvaluator evalHoriPoly(horiPoly);
  Polynomial2DEvaluator evalVertPolyS(vertPoly);
  Polynomial2DEvaluator evalVertPolyN(vertPoly);
  Polynomial2DEvaluator evalDiagPolySE(diagPoly);
  Polynomial2DEvaluator evalDiagPolyNW(diagPoly);

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    x[1] = j * h;

    // Set new Y values
//    evalCenterPoly.setY(x[1]);
    evalHoriPoly.setY(x[1]);
    evalVertPolyS.setY(x[1] - 0.5 * h);
    evalVertPolyN.setY(x[1] + 0.5 * h);
    evalDiagPolySE.setY(x[1] - 0.5 * h);
    evalDiagPolyNW.setY(x[1] + 0.5 * h);

    opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = evalHoriPoly.setStartX(-0.5 * h, h);
    opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval();

//    evalCenterPoly.setStartX(0.0, h);
    evalVertPolyS.setStartX(0.0, h);
    evalVertPolyN.setStartX(0.0, h);

    evalDiagPolySE.setStartX(0.5 * h, h);
    evalDiagPolyNW.setStartX(-0.5 * h, h);

    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)];
      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval();

      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)] = evalVertPolyS.incrementEval();
      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)] = evalVertPolyN.incrementEval();

      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)] = evalDiagPolySE.incrementEval();
      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)] = evalDiagPolyNW.incrementEval();

      opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = -opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)]
                                                                                -opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)]
                                                                                -opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)]
                                                                                -opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)]
                                                                                -opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)]
                                                                                -opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)];

      tmp = rhs[vertexdof::macroface::indexFromVertex(Level, i, j, stencilDirection::VERTEX_C)];

      //for (auto neighbor : neighbors) {
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        tmp -= opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])]*dst[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[k])];
      }

      dst[vertexdof::macroface::indexFromVertex(Level, i, j, stencilDirection::VERTEX_C)] = tmp/opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }
}

SPECIALIZE_POLYNOMIAL(void, smooth_gs_polynomial_even_tmpl, smooth_gs_polynomial_even)

}// namespace macroface
}// namespace vertexdof
}// namespace hyteg

