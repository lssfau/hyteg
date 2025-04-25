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

#include "core/Abort.h"

#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/trilinos/TpetraWrapper.hpp"

namespace hyteg {

using Teuchos::RCP;
using walberla::uint_t;

class TrilinosSparseMatrixProxy : public SparseMatrixProxy
{
 public:
   explicit TrilinosSparseMatrixProxy( const RCP< Tpetra::CrsMatrix<> >& mat )
   : mat_( mat )
   {}

   std::shared_ptr< SparseMatrixProxy > createCopy() const override
   {
      WALBERLA_ABORT( "Trilinos sparse matrix copy not implemented." );
   }
   std::shared_ptr< SparseMatrixProxy > createEmptyCopy() const override
   {
      WALBERLA_ABORT( "Trilinos sparse matrix copy not implemented." );
   }
   std::shared_ptr< SparseMatrixProxy > createMatrix( uint_t          localRows,
                                                      uint_t          localCols,
                                                      uint_t          globalRows,
                                                      uint_t          globalCols,
                                                      const MPI_Comm& MpiCommunicator ) const override
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   void addValue( uint_t row, uint_t col, real_t value ) override
   {
      mat_->insertGlobalValues( row,
                                Teuchos::tuple< Tpetra::Vector<>::global_ordinal_type >( col ),
                                Teuchos::tuple< Tpetra::Vector<>::scalar_type >( value ) );
   }

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

   void createFromMatrixProduct( const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) override
   {
      WALBERLA_ABORT( "Trilinos sparse matrix construction from matrix product not implemented." );
   };

   void createFromMatrixLinComb( const std::vector< real_t >&                               scalars,
                                 const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) override
   {
      WALBERLA_ABORT( "Not implemented yet!" );
   }

 private:
   RCP< Tpetra::CrsMatrix<> > mat_;
};

} // namespace hyteg