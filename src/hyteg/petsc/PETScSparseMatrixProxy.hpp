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
   {
      MatSetOption( mat_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );
   }

   virtual std::shared_ptr< SparseMatrixProxy > createCopy() const
   {
      Mat matCopy;

      MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );

      MatDuplicate( mat_, MAT_COPY_VALUES, &matCopy );
      return std::make_shared< PETScSparseMatrixProxy >( matCopy );
   }

#ifndef PETSC_USE_COMPLEX
   void addValue( uint_t row, uint_t col, real_t value ) override
   {
      PetscReal petscVal;

      // check whether we need to convert between PetscReal and real_t
      if constexpr ( std::is_same< PetscReal, real_t >::value )
      {
         petscVal = value;
      }
      else
      {
         static_cast< PetscReal >( value );
      }

      MatSetValue( mat_, static_cast< PetscInt >( row ), static_cast< PetscInt >( col ), petscVal, ADD_VALUES );
   }

   void addValues( const std::vector< uint_t >& rows,
                   const std::vector< uint_t >& cols,
                   const std::vector< real_t >& values ) override
   {
      WALBERLA_ASSERT_EQUAL( values.size(), rows.size() * cols.size() );
      std::vector< PetscInt > petscRows( rows.size() );
      std::vector< PetscInt > petscCols( cols.size() );
      for ( uint_t i = 0; i < rows.size(); i++ )
      {
         petscRows[i] = static_cast< PetscInt >( rows[i] );
      }
      for ( uint_t i = 0; i < cols.size(); i++ )
      {
         petscCols[i] = static_cast< PetscInt >( cols[i] );
      }

      // check whether we need to convert between PetscReal and real_t
      if constexpr ( std::is_same< PetscReal, real_t >::value )
      {
         MatSetValues( mat_,
                       static_cast< PetscInt >( petscRows.size() ),
                       petscRows.data(),
                       static_cast< PetscInt >( petscCols.size() ),
                       petscCols.data(),
                       values.data(),
                       ADD_VALUES );
      }
      else
      {
         std::vector< PetscReal > petscVals( values.size() );
         for ( uint_t k = 0; k < values.size(); k++ )
         {
            petscVals[k] = static_cast< PetscReal >( values[k] );
         }
         MatSetValues( mat_,
                       static_cast< PetscInt >( petscRows.size() ),
                       petscRows.data(),
                       static_cast< PetscInt >( petscCols.size() ),
                       petscCols.data(),
                       petscVals.data(),
                       ADD_VALUES );
      }
   }

#else

   void addValue( uint_t row, uint_t col, real_t value )
   {
      PetscComplex petscVal = static_cast< PetscReal >( value );
      MatSetValue( mat_, static_cast< PetscInt >( row ), static_cast< PetscInt >( col ), petscVal, ADD_VALUES );
   }

   void addValues( const std::vector< uint_t >& rows, const std::vector< uint_t >& cols, const std::vector< real_t >& values )
   {
      WALBERLA_ASSERT_EQUAL( values.size(), rows.size() * cols.size() );

      std::vector< PetscInt >     petscRows( rows.size() );
      std::vector< PetscInt >     petscCols( cols.size() );
      std::vector< PetscComplex > petscVals( rows.size() * cols.size() );

      // convert between datatypes
      for ( uint_t i = 0; i < rows.size(); i++ )
      {
         petscRows[i] = static_cast< PetscInt >( rows[i] );
      }
      for ( uint_t i = 0; i < cols.size(); i++ )
      {
         petscCols[i] = static_cast< PetscInt >( cols[i] );
      }
      for ( uint_t i = 0; i < rows.size() * cols.size(); i++ )
      {
         petscVals[i] = static_cast< PetscReal >( values[i] );
      }

      MatSetValues( mat_,
                    static_cast< idx_t >( petscRows.size() ),
                    petscRows.data(),
                    static_cast< idx_t >( petscCols.size() ),
                    petscCols.data(),
                    petscVals.data(),
                    ADD_VALUES );
   }

#endif

   void createFromMatrixProduct( const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) override
   {
      Mat tmp;

      PetscErrorCode err;

      err = MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );
      err = MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );

      err = MatDuplicate( mat_, MAT_DO_NOT_COPY_VALUES, &tmp );
      WALBERLA_CHECK( !err );

      WALBERLA_CHECK_GREATER( matrices.size(), 0 );

      auto petscProxy = std::dynamic_pointer_cast< PETScSparseMatrixProxy >( matrices.at( 0 ) );
      WALBERLA_CHECK_NOT_NULLPTR( petscProxy );

      err = MatAssemblyBegin( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );
      err = MatAssemblyEnd( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );

      err = MatCopy( petscProxy->mat_, mat_, DIFFERENT_NONZERO_PATTERN );
      WALBERLA_CHECK( !err );

      for ( uint_t i = 1; i < matrices.size(); i++ )
      {
         petscProxy = std::dynamic_pointer_cast< PETScSparseMatrixProxy >( matrices.at( i ) );

         err = MatAssemblyBegin( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
         err = MatAssemblyEnd( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );

         err = MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
         err = MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );

         err = MatMatMult( mat_, petscProxy->mat_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmp );
         WALBERLA_CHECK( !err );
         err = MatCopy( tmp, mat_, DIFFERENT_NONZERO_PATTERN );
         WALBERLA_CHECK( !err );

         err = MatAssemblyBegin( tmp, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
         err = MatAssemblyEnd( tmp, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );

         err = MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
         err = MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
      }
   }

 private:
   Mat mat_;
};

} // namespace hyteg
