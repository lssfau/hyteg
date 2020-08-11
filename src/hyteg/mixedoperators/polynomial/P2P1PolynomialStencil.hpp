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


#include <hyteg/p2functionspace/polynomial/P2PolynomialStencil.hpp>
#include <hyteg/mixedoperators/polynomial/P2P1PolynomialMemory.hpp>


namespace hyteg {
namespace P2toP1 {
namespace variablestencil {

using walberla::uint_t;
using walberla::real_c;
using indexing::Index;
using P2::variablestencil::PolynomialStencil;

template <uint_t PolyDegree>
class P2toP1PolynomialStencil
{
 public:
   P2toP1PolynomialStencil(const FacePolynomialMemory::FacePolynomials& poly)
      : VtV_(poly.VtV), EtV_(poly.EtV)
   {}

   inline void setY(const real_t y)
   {
      VtV_.setY(y);
      EtV_.setY(y);
   }

   // set x = x0 and use polynomial values for point (y,x) as stencilentries
   inline void setStartX(const real_t x0, const real_t h)
   {
      VtV_.setStartX(x0, h);
      EtV_.setStartX(x0, h);
   }

   // set x += h and use polynomial values for point (y,x) as stencilentries
   inline void incrementX()
   {
      VtV_.incrementX();
      EtV_.incrementX();
   }

   // getters for pointers to corresponding stencil data
   inline const real_t* vertexToVertex() const {return &VtV_;}
   inline const real_t* edgeToVertex() const {return &EtV_;}

 private:
   PolynomialStencil<PolyDegree, NumStencilentries2D::VtV> VtV_;
   PolynomialStencil<PolyDegree, NumStencilentries2D::EtV> EtV_;
};


} // namespace variablestencil
} // namespace P2toP1

namespace P1toP2 {
namespace variablestencil {

using walberla::uint_t;
using walberla::real_c;
using indexing::Index;
using P2::variablestencil::PolynomialStencil;

template <uint_t PolyDegree>
class P1toP2PolynomialStencil
{
 public:
   P1toP2PolynomialStencil(const FacePolynomialMemory::FacePolynomials& poly)
      : VtV_(poly.VtV), VtE_(poly.VtE)
   {}

   inline void setY(const real_t y)
   {
      VtV_.setY(y);
      VtE_.setY(y);
   }

   // set x = x0 and use polynomial values for point (y,x) as stencilentries
   inline void setStartX(const real_t x0, const real_t h)
   {
      VtV_.setStartX(x0, h);
      VtE_.setStartX(x0, h);
   }

   // set x += h and use polynomial values for point (y,x) as stencilentries
   inline void incrementX()
   {
      VtV_.incrementX();
      VtE_.incrementX();
   }

   // getters for pointers to corresponding stencil data
   inline const real_t* vertexToVertex() const {return &VtV_;}
   inline const real_t* vertexToEdge() const {return &VtE_;}

 private:
   PolynomialStencil<PolyDegree, NumStencilentries2D::VtV> VtV_;
   PolynomialStencil<PolyDegree, NumStencilentries2D::VtE> VtE_;
};


} // namespace variablestencil
} // namespace P1toP2
} // namespace hyteg
