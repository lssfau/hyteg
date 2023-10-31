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

#include <vector>

#include "core/DataTypes.h"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// \brief This class serves as a proxy for the assembly of sparse matrices, especially for external formats/libraries.
class SparseMatrixProxy
{
 public:
   virtual std::shared_ptr< SparseMatrixProxy > createCopy() const = 0;

   /// \brief Adds the passed value on the existing value in the matrix, or sets it to the value if no value exists.
   virtual void addValue( uint_t row, uint_t col, real_t value ) = 0;

   /// \brief Adds a "block" of values to the sparse matrix at once.
   ///
   /// The values vector is expected to be of size rows.size() * cols.size().
   virtual void
       addValues( const std::vector< uint_t >& rows, const std::vector< uint_t >& cols, const std::vector< real_t >& values ) = 0;

   /// \brief Stores the product of the passed matrices (as ordered in the vector) in this matrix.
   ///        Can be used to assemble concatenated operators.
   virtual void createFromMatrixProduct( const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) = 0;

   /// \brief Stores the sum of the passed matrices (as ordered in the vector) in this matrix.
   ///        Can be used to assemble concatenated operators.
   virtual void createFromMatrixLinComb( const std::vector< real_t >&                               scalars,
                                         const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) = 0;
};

} // namespace hyteg