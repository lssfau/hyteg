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
#include "hyteg/functions/FunctionIterator.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"

#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

using walberla::uint_t;

namespace hyteg {

static void testFunctionIterator( const std::string& meshFileName, const uint_t& level )
{
   auto storage = PrimitiveStorage::createFromGmshFile( meshFileName );

   P1Function< int >      nVertex( "nVertex", storage, level, level );
   EdgeDoFFunction< int > nEdge( "nEdge", storage, level, level );

   nVertex.enumerate( level );
   nEdge.enumerate( level );

   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );

   std::set< int > testEnumerateSet;
   uint_t          numGlobalDoFs;
   uint_t          numLocalDoFs;

   WALBERLA_LOG_INFO_ON_ROOT( "--- Testing VertexDoFSpace ---" );
   WALBERLA_MPI_BARRIER();

   numGlobalDoFs = numberOfGlobalDoFs< P1FunctionTag >( *storage, level );
   numLocalDoFs  = numberOfLocalDoFs< P1FunctionTag >( *storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( "global dofs: " << numGlobalDoFs );
   WALBERLA_LOG_INFO( "local dofs: " << numLocalDoFs );

   for ( const auto& dof : FunctionIterator< P1Function< int > >( nVertex, level ) )
   {
      WALBERLA_CHECK( dof.isVertexDoF() )
      WALBERLA_LOG_INFO( dof );
      WALBERLA_CHECK_EQUAL( testEnumerateSet.count( dof.value() ), 0 );
      testEnumerateSet.insert( dof.value() );
   }
   WALBERLA_CHECK_EQUAL( testEnumerateSet.size(), numLocalDoFs );

   testEnumerateSet.clear();

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_INFO_ON_ROOT( "--- Testing EdgeDoFSpace ---" );
   WALBERLA_MPI_BARRIER();

   numGlobalDoFs = numberOfGlobalDoFs< EdgeDoFFunctionTag >( *storage, level );
   numLocalDoFs  = numberOfLocalDoFs< EdgeDoFFunctionTag >( *storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( "global dofs: " << numGlobalDoFs );
   WALBERLA_LOG_INFO( "local dofs: " << numLocalDoFs );

   for ( const auto& dof : FunctionIterator< EdgeDoFFunction< int > >( nEdge, level ) )
   {
      WALBERLA_CHECK( dof.isEdgeDoF() )
      WALBERLA_LOG_INFO( dof );
      WALBERLA_CHECK_EQUAL( testEnumerateSet.count( dof.value() ), 0 );
      testEnumerateSet.insert( dof.value() );
   }
   WALBERLA_CHECK_EQUAL( testEnumerateSet.size(), numLocalDoFs );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "2D/annulus_coarse.msh" ), 0 );
   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "2D/annulus_coarse.msh" ), 1 );
   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "2D/annulus_coarse.msh" ), 2 );
   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "2D/annulus_coarse.msh" ), 3 );

   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 0 );
   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 1 );

   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "3D/cube_24el.msh" ), 0 );
   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "3D/cube_24el.msh" ), 1 );
   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "3D/cube_24el.msh" ), 2 );
   hyteg::testFunctionIterator( hyteg::prependHyTeGMeshDir( "3D/cube_24el.msh" ), 3 );

   return EXIT_SUCCESS;
}
