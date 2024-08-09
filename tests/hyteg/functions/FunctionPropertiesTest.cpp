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
#include "hyteg/functions/FunctionProperties.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

using walberla::uint_t;

namespace hyteg {

template < typename FunctionTag_T >
void testFunctionProperties( const std::string& meshFileName,
                             const uint_t&      level,
                             const uint_t&      expectedGlobalDoFs,
                             const uint_t&      expectedInnerDoFs )
{
   auto storage = PrimitiveStorage::createFromGmshFile( meshFileName );

   auto storageInfo = storage->getGlobalInfo();

   const uint_t numGlobalDoFs = numberOfGlobalDoFs< FunctionTag_T >( *storage, level );
   const uint_t numLocalDoFs  = numberOfLocalDoFs< FunctionTag_T >( *storage, level );
   const uint_t numInnerDoFs  = numberOfGlobalInnerDoFs< FunctionTag_T >( *storage, level );

   bool verbose = false;
   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "numGlobalDoFs = " << numGlobalDoFs );
      WALBERLA_LOG_INFO_ON_ROOT( "numGlobalInnerDoFs = " << numInnerDoFs );
   }

   WALBERLA_CHECK_EQUAL( expectedGlobalDoFs, numGlobalDoFs, meshFileName )
   WALBERLA_CHECK_EQUAL( expectedInnerDoFs, numInnerDoFs, meshFileName )

   if ( walberla::mpi::MPIManager::instance()->numProcesses() == 1 )
   {
      //local and global DoFs are the same in non parallel case:
      WALBERLA_CHECK_EQUAL( numGlobalDoFs, numLocalDoFs, meshFileName )
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   // P1Function
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 2, 15, 3 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 3, 45, 21 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >(
       hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 20, 549757386753, 549754241025 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( hyteg::prependHyTeGMeshDir( "2D/tri_2el.msh" ), 2, 25, 9 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( hyteg::prependHyTeGMeshDir( "2D/tri_2el.msh" ), 3, 81, 49 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 2, 35, 1 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 3, 165, 35 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >(
       hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 17, 375317149057025, 375282789318655 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 7, 723905, 674751 );

   // EdgeDofFunction
   hyteg::testFunctionProperties< hyteg::EdgeDoFFunctionTag >(
       hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 2, 12 + 6 * 3, 6 * 3 );
   hyteg::testFunctionProperties< hyteg::EdgeDoFFunctionTag >(
       hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 3, 24 + 28 * 3, 28 * 3 );

   // P2Function
   hyteg::testFunctionProperties< hyteg::P2FunctionTag >(
       hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 19, 549757386753, 549754241025 );
   hyteg::testFunctionProperties< hyteg::P2FunctionTag >( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 2, 165, 35 );
   hyteg::testFunctionProperties< hyteg::P2FunctionTag >( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 6, 723905, 674751 );

   // P0Function
   hyteg::testFunctionProperties< hyteg::P0FunctionTag >( hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 2, 16, 16 );
   hyteg::testFunctionProperties< hyteg::P0FunctionTag >( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 2, 64, 64 );

   return EXIT_SUCCESS;
}
