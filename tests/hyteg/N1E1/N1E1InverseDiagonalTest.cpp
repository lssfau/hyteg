/*
* Copyright (c) 2022-2023 Daniel Bauer.
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

#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_curl_curl_affine_q0.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_mass_affine_qe.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/petsc/PETScVectorProxy.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using namespace hyteg;
using walberla::real_t;

void test( uint_t level, MeshInfo meshInfo )
{
   using namespace n1e1;

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   N1E1VectorFunction< idx_t > numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   forms::n1e1_curl_curl_affine_q0          curlCurlForm;
   forms::n1e1_mass_affine_qe               massForm;
   N1E1ElementwiseLinearCombinationOperator A( storage, level, level, { { 1.0, 1.0 }, { &curlCurlForm, &massForm } } );

   PETScSparseMatrix< N1E1ElementwiseLinearCombinationOperator > matrix;
   matrix.createMatrixFromOperator( A, level, numerator );

   PETScVector< real_t, N1E1VectorFunction > vector( *A.getInverseDiagonalValues(), numerator, level );
   PETScVectorProxy                          inverseDiagonal( vector.get() );

   const uint_t n = numberOfGlobalDoFs( numerator, level );
   for ( uint_t k = 0; k < n; ++k )
   {
      PetscScalar diagVal;
      MatGetValue( matrix.get(), numeric_cast< PetscInt >( k ), numeric_cast< PetscInt >( k ), &diagVal );

      real_t invDiagVal = 0.0;
#if ( defined( PETSC_USE_COMPLEX ) && !defined( PETSC_SKIP_COMPLEX ) )
      invDiagVal = 1.0 / diagVal.real();
#else
      invDiagVal = 1.0 / diagVal;
#endif

      WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( inverseDiagonal.getValue( k ), invDiagVal, 1.0e-12, "k = " << k )
   }
}

int main( int argc, char** argv )
{
   using std::sin;

   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   test( 3, MeshInfo::fromGmshFile( prependHyTeGMeshDir( "3D/tet_1el.msh" ) ) );
   test( 1, MeshInfo::meshSymmetricCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), 1, 1, 1 ) );

   return EXIT_SUCCESS;
}
