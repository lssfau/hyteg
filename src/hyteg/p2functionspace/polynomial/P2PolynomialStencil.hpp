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
#include "hyteg/p2functionspace/polynomial/P2PolynomialMemory.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/polynomial/PolynomialEvaluator.hpp"


namespace hyteg {
namespace P2 {
namespace variablestencil {

using walberla::uint_t;
using walberla::real_c;
using indexing::Index;

template <uint_t PolyDegree, NumStencilentries2D N>
class PolynomialStencil
{
 public:

   using StencilEvaluator = std::array<std::shared_ptr<Polynomial2DEvaluator>, N>;
   using StencilEntries = std::array<real_t, N>;

   PolynomialStencil(const P2::StencilPolynomial<N>& poly)
   {
      for (uint_t i = 0; i < N; ++i)
      {
         evaluator_[i] = std::make_shared<Polynomial2DEvaluator>(poly[i]);
      }
   }

   inline const real_t* operator&() const {return stencil_.data();}

   inline void setY(const real_t y)
   {
      for (uint_t i = 0; i < N; ++i)
      {
         evaluator_[i]->setY(y);
      }
   }

   // set x = x0 and use polynomial values for point (y,x) as stencilentries
   inline void setStartX(const real_t x0, const real_t h)
   {
      for (uint_t i = 0; i < N; ++i)
      {
         Polynomial2DEvaluator& eval = *evaluator_[i];
         stencil_[i] = eval.setStartX(x0, h);
      }
   }

   // set x += h and use polynomial values for point (y,x) as stencilentries
   inline void incrementX()
   {
      for (uint_t i = 0; i < N; ++i)
      {
         Polynomial2DEvaluator& eval = *evaluator_[i];
         stencil_[i] = eval.incrementEval();
      }
   }

 private:
   StencilEvaluator evaluator_;
   StencilEntries stencil_;
};

template <uint_t PolyDegree>
class P2PolynomialStencil
{
 public:
   P2PolynomialStencil(const FacePolynomialMemory::FacePolynomials& poly)
   : VtV_( poly.VtV )
   , VtE_( poly.VtE )
   , EtV_( poly.EtV )
   , EtE_( poly.EtE )
   {}

   inline void setY(const real_t y)
   {
      VtV_.setY(y);
      VtE_.setY(y);
      EtV_.setY(y);
      EtE_.setY(y);
   }

   // set x = x0 and use polynomial values for point (y,x) as stencilentries
   inline void setStartX(const real_t x0, const real_t h)
   {
      VtV_.setStartX(x0, h);
      VtE_.setStartX(x0, h);
      EtV_.setStartX(x0, h);
      EtE_.setStartX(x0, h);
   }

   // set x += h and use polynomial values for point (y,x) as stencilentries
   inline void incrementX()
   {
      VtV_.incrementX();
      VtE_.incrementX();
      EtV_.incrementX();
      EtE_.incrementX();
   }

   // getters for pointers to corresponding stencil data

   inline const real_t* vertexToVertex() const {return &VtV_;}

   inline const real_t* edgeToVertex() const {return &EtV_;}

   inline const real_t* vertexToEdge() const {return &VtE_;}

   inline const real_t* edgeToEdge() const {return &EtE_;}

 private:
   PolynomialStencil<PolyDegree, NumStencilentries2D::VtV> VtV_;
   PolynomialStencil<PolyDegree, NumStencilentries2D::VtE> VtE_;
   PolynomialStencil<PolyDegree, NumStencilentries2D::EtV> EtV_;
   PolynomialStencil<PolyDegree, NumStencilentries2D::EtE> EtE_;
};


} // namespace variablestencil
} // namespace P2
} // namespace hyteg
