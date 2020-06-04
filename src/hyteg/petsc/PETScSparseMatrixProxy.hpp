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

#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

using walberla::uint_t;

class PETScSparseMatrixProxy : public SparseMatrixProxy
{
 public:
   PETScSparseMatrixProxy( Mat mat )
   : mat_( mat )
   {}

   void addValue( uint_t row, uint_t col, real_t value )
   {
      MatSetValue( mat_, static_cast< PetscInt >( row ), static_cast< PetscInt >( col ), value, ADD_VALUES );
   }

   void addValues( const std::vector< uint_t >& rows, const std::vector< uint_t >& cols, const std::vector< real_t >& values )
   {
      WALBERLA_ASSERT_EQUAL( rows.size(), cols.size() );
      WALBERLA_ASSERT_EQUAL( values.size(), rows.size() * cols.size() );
      std::vector< PetscInt > petscRows( rows.size() );
      std::vector< PetscInt > petscCols( cols.size() );
      for ( uint_t i = 0; i < rows.size(); i++ )
      {
         petscRows[i] = static_cast< PetscInt >( rows[i] );
         petscCols[i] = static_cast< PetscInt >( cols[i] );
      }
      MatSetValues( mat_,
                    static_cast< PetscInt >( petscRows.size() ),
                    petscRows.data(),
                    static_cast< PetscInt >( petscCols.size() ),
                    petscCols.data(),
                    values.data(),
                    ADD_VALUES );
   }

 private:
   Mat mat_;
};

} // namespace hyteg