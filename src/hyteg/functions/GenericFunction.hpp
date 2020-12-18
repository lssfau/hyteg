/*
 * Copyright (c) 2020 Marcus Mohr.
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

#include "core/DataTypes.h"

namespace hyteg {

using walberla::uint_t;

/// Non-templated base class for all function classes in HyTeG
class GenericFunction
{
 public:

   /// Returns the dimension of the field represented by the function
   ///
   /// A function in HyTeG
   ///
   ///  Value  | Explanation
   /// :-----: | :-----------------
   ///     1   | a scalar field represented by a Function
   ///   d>1   | a vector field with components of dimension d represented by a VectorFunction
   ///     0   | a composite function such as e.g. from a Taylor-Hood discretisation, represented by a BlockFunction
   ///
   virtual uint_t getDimension() const = 0;
};

} // namespace hyteg
