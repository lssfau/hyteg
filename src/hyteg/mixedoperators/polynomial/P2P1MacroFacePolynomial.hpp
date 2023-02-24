/*
 * Copyright (c) 2017-2019 Benjamin Mann.
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

#include "hyteg/p1functionspace/polynomial/VertexDoFMacroFacePolynomial.hpp"
#include "hyteg/mixedoperators/polynomial/P2P1PolynomialMemory.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/mixedoperators/polynomial/P2P1PolynomialStencil.hpp"


namespace hyteg {
namespace P2toP1 {
namespace variablestencil {
namespace macroface {

using walberla::uint_t;
using walberla::real_c;
using indexing::Index;

template<uint_t PolyDegree>
inline void applyPolynomialTmpl(const PrimitiveDataID<P2toP1::FacePolynomialMemory, Face>& polynomialId,
                                uint_t level, Face& face,
                                const PrimitiveDataID<FunctionMemory<real_t>, Face>& srcVertexDoFID,
                                const PrimitiveDataID<FunctionMemory<real_t>, Face>& srcEdgeDoFID,
                                const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstVertexDoFID,
                                UpdateType update)
{
   const uint_t N = levelinfo::num_microedges_per_edge(level);

   const real_t* srcVertexDoF = face.getData(srcVertexDoFID)->getPointer(level);
   const real_t* srcEdgeDoF   = face.getData(srcEdgeDoFID)->getPointer(level);
   real_t* dstVertexDoF = face.getData(dstVertexDoFID)->getPointer(level);

   indexing::Index idx;
   real_t  h = 1.0 / walberla::real_c(N);

   P2toP1PolynomialStencil<PolyDegree> stencil(face.getData(polynomialId)->getPolynomial(PolyDegree));

   // loop over all DoFs on reference element
   for (idx.y() = 0; idx.y() < N; ++idx.y())
   {
      stencil.setY(idx.y() * h);
      stencil.setStartX(0.0, h);

      for (idx.x() = 0; idx.x() < N - idx.y(); ++idx.x())
      {
         P2::macroface::applyStencil_DoF(level, idx,
                                         stencil.vertexToVertex(), stencil.edgeToVertex(),
                                         nullptr, nullptr,
                                         srcVertexDoF, srcEdgeDoF, dstVertexDoF, nullptr, update);

         stencil.incrementX();
      }
   }
}

SPECIALIZE_R_POLYNOMIAL(void, applyPolynomialTmpl, applyPolynomial)

} // namespace macroface
} // namespace variablestencil
} // namespace P2toP1

namespace P1toP2 {
namespace variablestencil {
namespace macroface {

using walberla::uint_t;
using walberla::real_c;
using indexing::Index;

template<uint_t PolyDegree>
inline void applyPolynomialTmpl(const PrimitiveDataID<P1toP2::FacePolynomialMemory, Face>& polynomialId,
                                uint_t level, Face& face,
                                const PrimitiveDataID<FunctionMemory<real_t>, Face>& srcVertexDoFID,
                                const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstVertexDoFID,
                                const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstEdgeDoFID,
                                UpdateType update)
{
   const uint_t N = levelinfo::num_microedges_per_edge(level);

   const real_t* srcVertexDoF = face.getData(srcVertexDoFID)->getPointer(level);
   real_t* dstVertexDoF = face.getData(dstVertexDoFID)->getPointer(level);
   real_t* dstEdgeDoF   = face.getData(dstEdgeDoFID)->getPointer(level);

   indexing::Index idx;
   real_t  h = 1.0 / walberla::real_c(N);

   P1toP2PolynomialStencil<PolyDegree> stencil(face.getData(polynomialId)->getPolynomial(PolyDegree));

   // loop over all DoFs on reference element
   for (idx.y() = 0; idx.y() < N; ++idx.y())
   {
      stencil.setY(idx.y() * h);
      stencil.setStartX(0.0, h);

      for (idx.x() = 0; idx.x() < N - idx.y(); ++idx.x())
      {
         P2::macroface::applyStencil_DoF(level, idx,
                                         stencil.vertexToVertex(), nullptr,
                                         stencil.vertexToEdge(), nullptr,
                                         srcVertexDoF, nullptr, dstVertexDoF, dstEdgeDoF, update);

         stencil.incrementX();
      }
   }
}

SPECIALIZE_R_POLYNOMIAL(void, applyPolynomialTmpl, applyPolynomial)

} // namespace macroface
} // namespace variablestencil
} // namespace P1toP2
} // namespace hyteg
