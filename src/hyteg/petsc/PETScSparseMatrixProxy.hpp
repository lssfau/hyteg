/*
 * Copyright (c) 2017-2025 Nils Kohl, Daniel Bauer, Andreas Burkhart, Marcus Mohr.
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

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"

#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

using walberla::uint_t;

class PETScSparseMatrixProxy : public SparseMatrixProxy
{
 public:
   PETScSparseMatrixProxy( Mat mat, bool autoDestroy = false )
   : mat_( mat )
   , autoDestroy_( autoDestroy )
   {
      PETScManager::ensureIsInitialized();
      MatSetOption( mat_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );
      MatMPIAIJSetPreallocation( mat, 500, NULL, 500, NULL );
   }

   ~PETScSparseMatrixProxy()
   {
      if ( autoDestroy_ )
      {
         MatDestroy( &mat_ );
      }
   }

   virtual std::shared_ptr< SparseMatrixProxy > createCopy() const override
   {
      Mat matCopy;

      MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );

      MatDuplicate( mat_, MAT_COPY_VALUES, &matCopy );
      return std::make_shared< PETScSparseMatrixProxy >( matCopy, true );
   }

   virtual std::shared_ptr< SparseMatrixProxy > createEmptyCopy() const override
   {
      Mat matCopy;

      MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );

      MatDuplicate( mat_, MAT_DO_NOT_COPY_VALUES, &matCopy );
      return std::make_shared< PETScSparseMatrixProxy >( matCopy, true );
   }

   virtual std::shared_ptr< SparseMatrixProxy > createMatrix( uint_t          localRows,
                                                              uint_t          localCols,
                                                              uint_t          globalRows,
                                                              uint_t          globalCols,
                                                              const MPI_Comm& MpiCommunicator ) const override
   {
      Mat matCreate;

      MatCreate( MpiCommunicator, &matCreate );
      MatSetType( matCreate, MATMPIAIJ );
      MatSetSizes( matCreate,
                   static_cast< PetscInt >( localRows ),
                   static_cast< PetscInt >( localCols ),
                   static_cast< PetscInt >( globalRows ),
                   static_cast< PetscInt >( globalCols ) );

      // Roughly overestimate number of non-zero entries for faster assembly of matrix
      MatMPIAIJSetPreallocation( matCreate, 500, NULL, 500, NULL );

      return std::make_shared< PETScSparseMatrixProxy >( matCreate, true );
   }

#ifndef PETSC_USE_COMPLEX
   void addValue( uint_t row, uint_t col, real_t value ) override
   {
      PetscReal petscVal = static_cast< PetscReal >( value );
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
                       // This cast is only here to prevent compiler errors
                       reinterpret_cast< const PetscReal* >( values.data() ),
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
      auto res = std::make_shared< Mat >();

      PetscErrorCode err;

      err = MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );
      err = MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );

      WALBERLA_CHECK_GREATER( matrices.size(), 0 );

      auto petscProxy = std::dynamic_pointer_cast< PETScSparseMatrixProxy >( matrices.at( 0 ) );
      WALBERLA_CHECK_NOT_NULLPTR( petscProxy );

      err = MatAssemblyBegin( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );
      err = MatAssemblyEnd( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );

      err = MatDuplicate( petscProxy->mat_, MAT_COPY_VALUES, res.get() );
      WALBERLA_CHECK( !err );

      for ( uint_t i = 1; i < matrices.size(); i++ )
      {
         Mat tmp;

         petscProxy = std::dynamic_pointer_cast< PETScSparseMatrixProxy >( matrices.at( i ) );

         err = MatAssemblyBegin( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
         err = MatAssemblyEnd( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );

         err = MatMatMult( *res, petscProxy->mat_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmp );
         WALBERLA_CHECK( !err );

         err = MatDestroy( res.get() );
         WALBERLA_CHECK( !err );
         res = std::make_shared< Mat >();

         err = MatDuplicate( tmp, MAT_COPY_VALUES, res.get() );
         WALBERLA_CHECK( !err );

         err = MatAssemblyBegin( tmp, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
         err = MatAssemblyEnd( tmp, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );

         err = MatAssemblyBegin( *res, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
         err = MatAssemblyEnd( *res, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );

         err = MatDestroy( &tmp );
         WALBERLA_CHECK( !err );
      }

      err = MatCopy( *res, mat_, DIFFERENT_NONZERO_PATTERN );
      WALBERLA_CHECK( !err );

      err = MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );
      err = MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );

      err = MatDestroy( res.get() );
      WALBERLA_CHECK( !err );
   }

   void createFromMatrixLinComb( const std::vector< real_t >&                               scalars,
                                 const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) override
   {
      auto tmp = std::make_shared< Mat >();

      PetscErrorCode err;

      err = MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );
      err = MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );

      WALBERLA_CHECK_GREATER( matrices.size(), 0 );

      auto petscProxy = std::dynamic_pointer_cast< PETScSparseMatrixProxy >( matrices.at( 0 ) );
      WALBERLA_CHECK_NOT_NULLPTR( petscProxy );

      err = MatAssemblyBegin( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );
      err = MatAssemblyEnd( petscProxy->mat_, MAT_FINAL_ASSEMBLY );
      WALBERLA_CHECK( !err );

      MatScale( petscProxy->mat_, scalars.at( 0 ) );
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

         Mat mats[2];

         err = MatDuplicate( mat_, MAT_COPY_VALUES, &mats[0] );
         WALBERLA_CHECK( !err );

         err = MatDuplicate( petscProxy->mat_, MAT_COPY_VALUES, &mats[1] );
         WALBERLA_CHECK( !err );

         err = MatScale( mats[1], scalars.at( i ) );
         WALBERLA_CHECK( !err );

         err = MatDestroy( tmp.get() );
         WALBERLA_CHECK( !err );
         tmp = std::make_shared< Mat >();

         err = MatCreateComposite( PETSC_COMM_WORLD, 2, mats, tmp.get() );
         WALBERLA_CHECK( !err );
         err = MatCompositeMerge( *tmp );
         WALBERLA_CHECK( !err );
         err = MatCopy( *tmp, mat_, DIFFERENT_NONZERO_PATTERN );
         WALBERLA_CHECK( !err );

         err = MatAssemblyBegin( *tmp, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
         err = MatAssemblyEnd( *tmp, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );

         err = MatAssemblyBegin( mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );
         err = MatAssemblyEnd( mat_, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !err );

         err = MatDestroy( &mats[0] );
         WALBERLA_CHECK( !err );

         err = MatDestroy( &mats[1] );
         WALBERLA_CHECK( !err );
      }

      err = MatDestroy( tmp.get() );
      WALBERLA_CHECK( !err );
   }

 private:
   Mat  mat_;
   bool autoDestroy_;
};

} // namespace hyteg
