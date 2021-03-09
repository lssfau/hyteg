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

#include <hyteg/polynomial/Polynomial.hpp>

namespace hyteg {
namespace P2 {

template <NumStencilentries2D N>
class StencilPolynomial
{
 public:

   inline void initialize(uint_t degree)
   {
      for (uint_t i = 0; i < N; ++i)
      {
         data_[i] = std::make_shared<GeneralPolynomial2D>(degree);
      }
   }

   inline GeneralPolynomial2D& operator[](uint_t i) {return *data_[i];}

   inline const GeneralPolynomial2D& operator[](uint_t i) const {return *data_[i];}

 private:
   std::array<std::shared_ptr<GeneralPolynomial2D>, N> data_;
};

} // namespace P2
} // namespace hyteg