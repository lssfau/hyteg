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

#include <hyteg/p1functionspace/VertexDoFMemory.hpp>
#include <hyteg/p2functionspace/variablestencil/P2VariableStencilCommon.hpp>

namespace hyteg {
namespace P2 {

class FacePolynomialMemory
{
 public:

   template <NumStencilentries2D N>
   using StencilPolynomial = std::array<std::shared_ptr<GeneralPolynomial2D>, N>;

   struct FacePolynomials
   {
      StencilPolynomial<NumStencilentries2D::VtV> VtV;
      StencilPolynomial<NumStencilentries2D::EtV> EtV;
      StencilPolynomial<NumStencilentries2D::VtE> VtE;
      StencilPolynomial<NumStencilentries2D::EtE> EtE;
   };

   FacePolynomialMemory() {}

   inline FacePolynomials& addDegree(uint_t degree)
   {
      if (polynomialDegreeExists(degree))
      {
         WALBERLA_LOG_WARNING("Degree already exists.");
      }

      FacePolynomials& poly = polynomials_[degree];

      initializePolynomial<NumStencilentries2D::VtV>(poly.VtV, degree);
      initializePolynomial<NumStencilentries2D::EtV>(poly.EtV, degree);
      initializePolynomial<NumStencilentries2D::VtE>(poly.VtE, degree);
      initializePolynomial<NumStencilentries2D::EtE>(poly.EtE, degree);

      return poly;
   }

   inline bool polynomialDegreeExists(uint_t degree)
   {
      return polynomials_.count(degree) > 0;
   }

   inline FacePolynomials& getPolynomial(uint_t degree)
   {
      return polynomials_[degree];
   }


 private:

   template <NumStencilentries2D N>
   inline void initializePolynomial(StencilPolynomial<N>& poly, uint_t degree)
   {
      for (uint_t i = 0; i < N; ++i)
      {
         poly[i] = std::make_shared<GeneralPolynomial2D>(degree);
      }
   }

   std::map<uint_t, FacePolynomials> polynomials_;
};

} // namespace P2
} // namespace hyteg