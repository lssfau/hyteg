/*
* Copyright (c) 2023 Nils Kohl.
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

#include "core/Abort.h"

#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

using walberla::int_c;

/// \brief This class can be used to conveniently create an Eigen sparse matrix object from a HyTeG operator.
///
/// The matrix is internally represented as a triplet of row/col indices and the value.
/// This allows efficient assembly of the sparse matrix.
/// Once finished, the current state of the proxy can be converted to an Eigen sparse matrix object.
class EigenSparseMatrixProxy : public SparseMatrixProxy
{
 public:
   /// \brief Creates an Eigen sparse matrix proxy object that can be used to build an Eigen sparse matrix from a HyTeG operator.
   EigenSparseMatrixProxy( uint_t rows, uint_t cols )
   : rows_( rows )
   , cols_( cols )
   {}

   virtual ~EigenSparseMatrixProxy() = default;

   std::shared_ptr< SparseMatrixProxy > createCopy() const override { WALBERLA_ABORT( "Not implemented." ); }
   std::shared_ptr< SparseMatrixProxy > createEmptyCopy() const override { WALBERLA_ABORT( "Not implemented." ); }
   std::shared_ptr< SparseMatrixProxy > createMatrix( uint_t          localRows,
                                                      uint_t          localCols,
                                                      uint_t          globalRows,
                                                      uint_t          globalCols,
                                                      const MPI_Comm& MpiCommunicator ) const override
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   /// \brief Resets the proxy such that it can be reused.
   void clear() { tripletList_.clear(); }

   /// \brief Adds the passed value on the existing value in the matrix, or sets it to the value if no value exists.
   void addValue( uint_t row, uint_t col, real_t value ) override
   {
      tripletList_.push_back( Triplet( int_c( row ), int_c( col ), value ) );
   }

   /// \brief Adds a "block" of values to the sparse matrix at once.
   ///
   /// The values vector is expected to be of size rows.size() * cols.size().
   void addValues( const std::vector< uint_t >& rows,
                   const std::vector< uint_t >& cols,
                   const std::vector< real_t >& values ) override
   {
      WALBERLA_ASSERT_EQUAL( values.size(), rows.size() * cols.size() );

      for ( uint_t i = 0; i < rows.size(); i++ )
      {
         for ( uint_t j = 0; j < cols.size(); j++ )
         {
            addValue( rows[i], cols[j], values[i * cols.size() + j] );
         }
      }
   }

   /// \brief Stores the product of the passed matrices (as ordered in the vector) in this matrix.
   ///        Can be used to assemble concatenated operators.
   void createFromMatrixProduct( const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) override
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   /// \brief Stores the sum of the passed matrices (as ordered in the vector) in this matrix.
   ///        Can be used to assemble concatenated operators.
   void createFromMatrixLinComb( const std::vector< real_t >&                               scalars,
                                 const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) override
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   /// \brief Builds a compressed Eigen::SparseMatrix from the current state of the proxy.
   Eigen::SparseMatrix< real_t, Eigen::RowMajor > getSparseMatrix() const
   {
      Eigen::SparseMatrix< real_t, Eigen::RowMajor > mat( static_cast< Eigen::Index >( rows_ ),
                                                          static_cast< Eigen::Index >( cols_ ) );
      mat.setFromTriplets( tripletList_.begin(), tripletList_.end() );
      return mat;
   }

 private:
   uint_t rows_;
   uint_t cols_;

   using Triplet = Eigen::Triplet< real_t >;
   std::vector< Triplet > tripletList_;
};

} // namespace hyteg
