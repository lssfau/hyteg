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

// clang-format off
#include "hyteg/eigen/EigenWrapper.hpp" // for hyteg/eigen/EigenMatrixPlugin.hpp
#include <Eigen/Sparse>
// clang-format on

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "hyteg/gridtransferoperators/N1E1toN1E1Prolongation.hpp"
#include "hyteg/gridtransferoperators/N1E1toN1E1Restriction.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

void test( const uint_t fineLevel, const MeshInfo meshInfo, const bool print = false )
{
   using namespace n1e1;

   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t coarseLevel = fineLevel - 1;
   const uint_t nCoarse     = numberOfGlobalDoFs< N1E1VectorFunctionTag >( *storage, coarseLevel );
   const uint_t nFine       = numberOfGlobalDoFs< N1E1VectorFunctionTag >( *storage, fineLevel );

   N1E1VectorFunction< real_t > f( "f", storage, coarseLevel, fineLevel );
   N1E1VectorFunction< idx_t >  numerator( "numerator", storage, coarseLevel, fineLevel );
   numerator.enumerate( coarseLevel );
   numerator.enumerate( fineLevel );

   PETScVector< real_t, N1E1VectorFunction > uVectorC( f, numerator, coarseLevel );
   PETScVector< real_t, N1E1VectorFunction > uVectorF( f, numerator, fineLevel );
   PETScVectorProxy                          unitVectorCoarse( uVectorC.get() );
   PETScVectorProxy                          unitVectorFine( uVectorF.get() );

   N1E1toN1E1Prolongation        prolongation;
   N1E1toN1E1Restriction         restriction;
   Eigen::SparseMatrix< real_t > prolongationMatrix( numeric_cast< Eigen::Index >( nFine ),
                                                     numeric_cast< Eigen::Index >( nCoarse ) );
   Eigen::SparseMatrix< real_t > restrictionMatrix( numeric_cast< Eigen::Index >( nCoarse ),
                                                    numeric_cast< Eigen::Index >( nFine ) );

   // list of non-zeros for efficient insertion in sparse matrices
   std::vector< Eigen::Triplet< double > > tripletList;

   // fill vectors with junk, in particular ghost layers!
   f.interpolate( 3.14, coarseLevel );
   f.interpolate( 3.14, fineLevel );
   f.communicate< Edge, Face >( coarseLevel );
   f.communicate< Face, Cell >( coarseLevel );
   f.communicate< Edge, Face >( fineLevel );
   f.communicate< Face, Cell >( fineLevel );

   // assemble prolongation matrix
   for ( uint_t j = 0; j < nCoarse; ++j )
   {
      // make `unitVectorCoarse` and `f` the j-th unit vector
      if ( j != 0 )
      {
         unitVectorCoarse.setValue( j - 1, 0.0 );
      }
      unitVectorCoarse.setValue( j, 1.0 );
      uVectorC.createFunctionFromVector( f, numerator, coarseLevel );

      // apply operator to unit vector, result is j-th column of matrix
      prolongation.prolongate( f, coarseLevel, All );
      uVectorF.createVectorFromFunction( f, numerator, fineLevel );

      // copy j-th column to matrix
      for ( uint_t i = 0; i < nFine; ++i )
      {
         const real_t val = unitVectorFine.getValue( i );
         if ( walberla::debug::check_functions_detail::check_float_unequal( val, 0.0 ) )
         {
            tripletList.emplace_back( i, j, val );
         }
      }
   }
   prolongationMatrix.setFromTriplets( tripletList.begin(), tripletList.end() );

   tripletList.clear();
   f.setToZero( fineLevel );
   uVectorF.createVectorFromFunction( f, numerator, fineLevel );

   // fill vectors with junk, in particular ghost layers!
   f.interpolate( 3.14, coarseLevel );
   f.interpolate( 3.14, fineLevel );
   f.communicate< Edge, Face >( coarseLevel );
   f.communicate< Face, Cell >( coarseLevel );
   f.communicate< Edge, Face >( fineLevel );
   f.communicate< Face, Cell >( fineLevel );

   // assemble restriction matrix
   for ( uint_t j = 0; j < nFine; ++j )
   {
      // make `unitVectorFine` and `f` the j-th unit vector
      if ( j != 0 )
      {
         unitVectorFine.setValue( j - 1, 0.0 );
      }
      unitVectorFine.setValue( j, 1.0 );
      uVectorF.createFunctionFromVector( f, numerator, fineLevel );

      // apply operator to unit vector, result is j-th column of matrix
      restriction.restrict( f, fineLevel, All );
      uVectorC.createVectorFromFunction( f, numerator, coarseLevel );

      // copy j-th column to matrix
      for ( uint_t i = 0; i < nCoarse; ++i )
      {
         const real_t val = unitVectorCoarse.getValue( i );
         if ( walberla::debug::check_functions_detail::check_float_unequal( val, 0.0 ) )
         {
            tripletList.emplace_back( i, j, val );
         }
      }
   }
   restrictionMatrix.setFromTriplets( tripletList.begin(), tripletList.end() );

   // check restriction == transpose(prolongation)
   Eigen::SparseMatrix< real_t > prolongationTransp = prolongationMatrix.transpose();

   if ( print )
   {
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( prolongationTransp );
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( restrictionMatrix );
   }

   WALBERLA_CHECK_EQUAL( restrictionMatrix.outerSize(), prolongationTransp.outerSize() );
   WALBERLA_CHECK_EQUAL( restrictionMatrix.innerSize(), prolongationTransp.innerSize() );

   for ( int i = 0; i < restrictionMatrix.outerSize(); ++i )
   {
      Eigen::SparseMatrix< real_t >::InnerIterator itR( restrictionMatrix, i );
      Eigen::SparseMatrix< real_t >::InnerIterator itPT( prolongationTransp, i );

      for ( ; itR && itPT; ++itR, ++itPT )
      {
         WALBERLA_CHECK_EQUAL( itR.row(), itPT.row() );
         WALBERLA_CHECK_EQUAL( itR.col(), itPT.col() );
         WALBERLA_CHECK_FLOAT_EQUAL(
             itR.value(), itPT.value(), "Mismatch at restriction (row, col): (" << itR.row() << ", " << itR.col() << ")" );
      }
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );

   test( 1, MeshInfo::fromGmshFile( "../../meshes/3D/tet_1el.msh" ), true );
   test( 4, MeshInfo::fromGmshFile( "../../meshes/3D/tet_1el.msh" ) );
   test( 3, MeshInfo::fromGmshFile( "../../meshes/3D/pyramid_2el.msh" ) );
   test( 3, MeshInfo::fromGmshFile( "../../meshes/3D/pyramid_4el.msh" ) );
   test( 3, MeshInfo::fromGmshFile( "../../meshes/3D/regular_octahedron_8el.msh" ) );
   test( 3, MeshInfo::meshSymmetricCuboid( Point3D(  0, 0, 0  ), Point3D(  1, 1, 1  ), 1, 1, 1 ) );

   return EXIT_SUCCESS;
}
