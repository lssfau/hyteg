/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "hyteg/gridtransferoperators/N1E1toP1Lifting.hpp"

#include <Eigen/Sparse>

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1toN1E1Gradient.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

void test( const uint_t lvl, const MeshInfo meshInfo, const bool print = false )
{
   using namespace n1e1;

   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t nN1E1 = numberOfGlobalDoFs< N1E1VectorFunctionTag >( *storage, lvl );
   const uint_t nP1   = numberOfGlobalDoFs< P1FunctionTag >( *storage, lvl );

   N1E1VectorFunction< real_t > fN1E1( "fN1E1", storage, lvl, lvl );
   N1E1VectorFunction< idx_t >  numeratorN1E1( "numeratorN1E1", storage, lvl, lvl );
   numeratorN1E1.enumerate( lvl );

   P1Function< real_t > fP1( "fP1", storage, lvl, lvl );
   P1Function< idx_t >  numeratorP1( "numeratorP1", storage, lvl, lvl );
   numeratorP1.enumerate( lvl );

   PETScVector< real_t, N1E1VectorFunction > uVectorN1E1( fN1E1, numeratorN1E1, lvl );
   PETScVector< real_t, P1Function >         uVectorP1( fP1, numeratorP1, lvl );
   PETScVectorProxy                          unitVectorN1E1( uVectorN1E1.get() );
   PETScVectorProxy                          unitVectorP1( uVectorP1.get() );

   Eigen::SparseMatrix< real_t > gradientMatrix( numeric_cast< Eigen::Index >( nN1E1 ), numeric_cast< Eigen::Index >( nP1 ) );
   Eigen::SparseMatrix< real_t > liftingMatrix( numeric_cast< Eigen::Index >( nP1 ), numeric_cast< Eigen::Index >( nN1E1 ) );

   // list of non-zeros for efficient insertion in sparse matrices
   std::vector< Eigen::Triplet< double > > tripletList;

   // assemble gradient matrix
   for ( uint_t j = 0; j < nP1; ++j )
   {
      // make `unitVectorP1` and `fP1` the j-th unit vector
      if ( j != 0 )
      {
         unitVectorP1.setValue( j - 1, 0.0 );
      }
      unitVectorP1.setValue( j, 1.0 );
      uVectorP1.createFunctionFromVector( fP1, numeratorP1, lvl );

      // apply operator to unit vector, result is j-th column of matrix
      P1toN1E1Gradient( fP1, fN1E1, lvl, All );
      uVectorN1E1.createVectorFromFunction( fN1E1, numeratorN1E1, lvl );

      // copy j-th column to matrix
      for ( uint_t i = 0; i < nN1E1; ++i )
      {
         const real_t val = unitVectorN1E1.getValue( i );
         if ( walberla::debug::check_functions_detail::check_float_unequal( val, 0.0 ) )
         {
            tripletList.emplace_back( i, j, val );
         }
      }
   }
   gradientMatrix.setFromTriplets( tripletList.begin(), tripletList.end() );

   tripletList.clear();
   fN1E1.setToZero( lvl );
   uVectorN1E1.createVectorFromFunction( fN1E1, numeratorN1E1, lvl );

   // assemble lifting matrix
   for ( uint_t j = 0; j < nN1E1; ++j )
   {
      // make `unitVectorN1E1` and `fN1E1` the j-th unit vector
      if ( j != 0 )
      {
         unitVectorN1E1.setValue( j - 1, 0.0 );
      }
      unitVectorN1E1.setValue( j, 1.0 );
      uVectorN1E1.createFunctionFromVector( fN1E1, numeratorN1E1, lvl );

      // apply operator to unit vector, result is j-th column of matrix
      N1E1toP1Lifting( fN1E1, fP1, lvl, All );
      uVectorP1.createVectorFromFunction( fP1, numeratorP1, lvl );

      // copy j-th column to matrix
      for ( uint_t i = 0; i < nP1; ++i )
      {
         const real_t val = unitVectorP1.getValue( i );
         if ( walberla::debug::check_functions_detail::check_float_unequal( val, 0.0 ) )
         {
            tripletList.emplace_back( i, j, val );
         }
      }
   }
   liftingMatrix.setFromTriplets( tripletList.begin(), tripletList.end() );

   // check lifting == transpose(gradient)
   Eigen::SparseMatrix< real_t > gradientTransp = gradientMatrix.transpose();

   if ( print )
   {
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( gradientTransp );
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( liftingMatrix );
   }

   WALBERLA_CHECK_EQUAL( liftingMatrix.outerSize(), gradientTransp.outerSize() );
   WALBERLA_CHECK_EQUAL( liftingMatrix.innerSize(), gradientTransp.innerSize() );

   for ( int i = 0; i < liftingMatrix.outerSize(); ++i )
   {
      Eigen::SparseMatrix< real_t >::InnerIterator itL( liftingMatrix, i );
      Eigen::SparseMatrix< real_t >::InnerIterator itGT( gradientTransp, i );

      for ( ; itL && itGT; ++itL, ++itGT )
      {
         WALBERLA_CHECK_EQUAL( itL.row(), itGT.row() );
         WALBERLA_CHECK_EQUAL( itL.col(), itGT.col() );
         WALBERLA_CHECK_FLOAT_EQUAL(
             itL.value(), itGT.value(), "Mismatch at lifting (row, col): (" << itL.row() << ", " << itL.col() << ")" );
      }
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );

   test( 1, MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), true );
   test( 4, MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ) );
   test( 3, MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ) );
   test( 3, MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ) );
   test( 3, MeshInfo::fromGmshFile( "../../data/meshes/3D/regular_octahedron_8el.msh" ) );

   return EXIT_SUCCESS;
}
