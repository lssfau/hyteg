/*
 * Copyright (c) 2024 Nils Kohl, Marcus Mohr.
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"

#include "hyteg/dataexport/Preprocessing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

namespace hyteg {

template < typename FunctionType >
void testDuplicatePointCountingForOutput( std::string meshFile, uint_t level, typename FunctionType::valueType maxNeighbors )
{
   auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

   FunctionType u( "u", storage, level, level );

   interpolateNumberOfAdjacentMacroVolumes( u, level );

   auto min = u.getMinDoFValue( level );
   auto max = u.getMaxDoFValue( level );

   const bool vtk = false;

   if ( vtk )
   {
      writeDomainPartitioningVTK( *storage, ".", "domain" );

      VTKOutput vtkMesh( ".", "data", storage );
      vtkMesh.add( u );
      vtkMesh.write( level );
   }

   WALBERLA_CHECK_GREATER( min, 0 );
   WALBERLA_CHECK_LESS_EQUAL( max, maxNeighbors );
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::testDuplicatePointCountingForOutput< hyteg::P1Function< int32_t > >( hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 3, 1 );
   hyteg::testDuplicatePointCountingForOutput< hyteg::P2Function< int32_t > >( hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 3, 1 );

   hyteg::testDuplicatePointCountingForOutput< hyteg::P1Function< int32_t > >( hyteg::prependHyTeGMeshDir( "2D/bfs_126el.msh" ), 3, 8 );
   hyteg::testDuplicatePointCountingForOutput< hyteg::P2Function< int32_t > >( hyteg::prependHyTeGMeshDir( "2D/bfs_126el.msh" ), 3, 8 );

   hyteg::testDuplicatePointCountingForOutput< hyteg::P1Function< int32_t > >( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 3, 1 );
   hyteg::testDuplicatePointCountingForOutput< hyteg::P2Function< int32_t > >( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 3, 1 );

   hyteg::testDuplicatePointCountingForOutput< hyteg::P1Function< int32_t > >(
       hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 3, 2 );
   hyteg::testDuplicatePointCountingForOutput< hyteg::P2Function< int32_t > >(
       hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 3, 2 );

   hyteg::testDuplicatePointCountingForOutput< hyteg::P1Function< int32_t > >( hyteg::prependHyTeGMeshDir( "3D/cube_6el.msh" ), 3, 6 );
   hyteg::testDuplicatePointCountingForOutput< hyteg::P2Function< int32_t > >(
       hyteg::prependHyTeGMeshDir( "3D/cube_6el.msh" ), 3, 6 );

   return 0;
}
