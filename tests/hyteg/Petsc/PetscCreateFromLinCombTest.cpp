/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
#include <core/timing/Timer.h>

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScSparseMatrixProxy.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
#error "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON"
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   PETScManager petscManager( &argc, &argv );

   // setup 3 matrices
   Mat M1, M2, M3;
   MatCreate( PETSC_COMM_WORLD, &M1 );
   MatCreate( PETSC_COMM_WORLD, &M2 );
   MatCreate( PETSC_COMM_WORLD, &M3 );
   MatSetSizes( M1, 2, 2, 2, 2 );
   MatSetSizes( M2, 2, 2, 2, 2 );
   MatSetSizes( M3, 2, 2, 2, 2 );
   MatSetUp( M1 );
   MatSetUp( M2 );
   MatSetUp( M3 );
   MatZeroEntries( M1 );
   MatZeroEntries( M2 );
   MatZeroEntries( M3 );

   // and corresponding hyteg proxy objects
   PETScSparseMatrixProxy pM1( M1 );
   PETScSparseMatrixProxy pM2( M2 );
   PETScSparseMatrixProxy pM3( M3 );

   // fill operand matrices
   pM1.addValue( 0, 0, 1 );
   pM1.addValue( 0, 1, 2 );
   pM1.addValue( 1, 0, 3 );
   pM1.addValue( 1, 1, 4 );

   pM2.addValue( 0, 0, 2 );
   pM2.addValue( 0, 1, 4 );
   pM2.addValue( 1, 0, 6 );
   pM2.addValue( 1, 1, 8 );

   // petsc assembly stuff
   MatAssemblyBegin( M1, MAT_FINAL_ASSEMBLY );
   MatAssemblyEnd( M1, MAT_FINAL_ASSEMBLY );
   MatAssemblyBegin( M2, MAT_FINAL_ASSEMBLY );
   MatAssemblyEnd( M2, MAT_FINAL_ASSEMBLY );

   // M3 = 2*M1 - M2
   pM3.createFromMatrixLinComb(
       { 2, -1 }, { std::make_shared< PETScSparseMatrixProxy >( pM1 ), std::make_shared< PETScSparseMatrixProxy >( pM2 ) } );

   // check result
   PetscReal val = 0;
   for ( int i = 0; i < 2; ++i )
   {
      for ( int j = 0; j < 2; ++j )
      {
         MatGetValue( M3, i, j, &val );
         WALBERLA_CHECK_FLOAT_EQUAL( val, 0.0 );
      }
   }
   return EXIT_SUCCESS;
}
