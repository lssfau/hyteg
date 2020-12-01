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

#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/mixedoperators/polynomial/P2P1PolynomialMemory.hpp"
#include "hyteg/primitivedata/PrimitiveDataHandling.hpp"

namespace hyteg {

class FaceP2toP1PolynomialMemoryDataHandling : public OnlyInitializeDataHandling<P2toP1::FacePolynomialMemory, Face>
{
 public:
   FaceP2toP1PolynomialMemoryDataHandling() {}

   std::shared_ptr<P2toP1::FacePolynomialMemory> initialize(const Face* const) const
   {
      return std::make_shared<P2toP1::FacePolynomialMemory>();
   }
};

class FaceP1toP2PolynomialMemoryDataHandling : public OnlyInitializeDataHandling<P1toP2::FacePolynomialMemory, Face>
{
 public:
   FaceP1toP2PolynomialMemoryDataHandling() {}

   std::shared_ptr<P1toP2::FacePolynomialMemory> initialize(const Face* const) const
   {
      return std::make_shared<P1toP2::FacePolynomialMemory>();
   }
};

} // namespace hyteg
