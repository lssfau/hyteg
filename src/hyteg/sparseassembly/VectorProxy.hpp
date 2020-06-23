/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// \brief This class serves as a proxy for the assembly of vectors, especially for external formats/libraries (e.g. sparse
/// solvers).
class VectorProxy
{
 public:
   /// \brief Sets the passed value in the vector.
   virtual void setValue( uint_t idx, real_t value ) = 0;

   /// \brief Returns the passed value of the vector.
   virtual real_t getValue( uint_t idx ) const = 0;
};

} // namespace hyteg