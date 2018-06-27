#pragma once

#include "core/debug/all.h"

namespace hhg {
enum class OperatorType {
  MASS,
  EVEN,
  ODD
};
}

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/Macros.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/polynomial/PolynomialEvaluator.hpp"

namespace hhg {
namespace vertexdof {
namespace macroface {

using walberla::uint_t;
using walberla::real_c;
using indexing::Index;

template<typename ValueType, OperatorType OprType, uint_t PolyDegree>
inline void applyPolynomialTmpl(uint_t Level, Face &face, const PrimitiveDataID<FaceP1PolynomialMemory, Face>& polynomialId,
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

   faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = evalHoriPoly.setStartX<PolyDegree>(-0.5 * h, h);
   faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval<PolyDegree>();

   if (OprType != OperatorType::EVEN) {
      evalCenterPoly.setStartX<PolyDegree>(0.0, h);
   }
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

//     faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = evalCenterPoly.incrementEval<PolyDegree>();

     if (OprType == OperatorType::MASS) {
        faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = evalCenterPoly.incrementEval<PolyDegree>();
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

      evalPolyS.setStartX<PolyDegree>(0.0, h);
      evalPolySE.setStartX<PolyDegree>(0.0, h);
      evalPolyW.setStartX<PolyDegree>(0.0, h);
//      evalPolyC.setStartX<PolyDegree>(0.0, h);
      evalPolyE.setStartX<PolyDegree>(0.0, h);
      evalPolyNW.setStartX<PolyDegree>(0.0, h);
      evalPolyN.setStartX<PolyDegree>(0.0, h);

      for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = evalPolyW.incrementEval<PolyDegree>();
         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalPolyE.incrementEval<PolyDegree>();

         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)] = evalPolyS.incrementEval<PolyDegree>();
         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)] = evalPolyN.incrementEval<PolyDegree>();

         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)] = evalPolySE.incrementEval<PolyDegree>();
         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)] = evalPolyNW.incrementEval<PolyDegree>();

//         faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = evalPolyC.incrementEval<PolyDegree>();

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

    opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = evalHoriPoly.setStartX<PolyDegree>(-0.5 * h, h);
    opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = evalHoriPoly.incrementEval<PolyDegree>();

//    evalCenterPoly.setStartX<PolyDegree>(0.0, h);
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
}// namespace hhg

