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

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/eigen/typeAliases.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

void test3D()
{
   // order of vertices is important
   MeshInfo meshInfo =
       MeshInfo::singleTetrahedron( { Point3D{ { 0, 0, 0 } }, { { 1, 0, 0 } }, { { 0, 1, 0 } }, { { 0, 0, 1 } } } );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t minLevel = 2;
   const size_t maxLevel = 4;

   const Eigen::Vector3r c = { 1, 3, 6 };

   n1e1::N1E1VectorFunction< real_t > f( "f", storage, minLevel, maxLevel );

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Level: " << level )

      f.interpolate( c, level );
      f.communicate< Edge, Face >( level );

      auto dofs = f.getDoFs();

      for ( const auto& edgeIt : storage->getEdges() )
      {
         const Edge& edge     = *edgeIt.second;
         const auto  edgeData = edge.getData( dofs->getEdgeDataID() )->getPointer( level );
         const auto  nEdges   = real_c( levelinfo::num_microedges_per_edge( level ) );

         const auto correct = c.dot( edge.getDirection().vector_ ) / nEdges;

         for ( const auto& it : edgedof::macroedge::Iterator( level ) )
         {
            const uint_t idx    = edgedof::macroedge::index( level, it.col() );
            const auto   actual = edgeData[idx];

            WALBERLA_CHECK_FLOAT_EQUAL( actual,
                                        correct,
                                        "edge coordinates: (" << edge.getCoordinates()[0] << ", " << edge.getCoordinates()[1]
                                                              << "), edge index: " << idx );
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Edges passed" )

      for ( const auto& faceIt : storage->getFaces() )
      {
         const Face& face     = *faceIt.second;
         const auto  faceData = face.getData( dofs->getFaceDataID() )->getPointer( level );
         const auto  nEdges   = real_c( levelinfo::num_microedges_per_edge( level ) );

         const auto correctX  = c.dot( face.getCoordinates()[1].vector_ - face.getCoordinates()[0].vector_ ) / nEdges;
         const auto correctY  = c.dot( face.getCoordinates()[2].vector_ - face.getCoordinates()[0].vector_ ) / nEdges;
         const auto correctXY = c.dot( face.getCoordinates()[2].vector_ - face.getCoordinates()[1].vector_ ) / nEdges;

         for ( const auto& it : edgedof::macroface::Iterator( level ) )
         {
            const uint_t idxX  = edgedof::macroface::horizontalIndex( level, it.col(), it.row() );
            const uint_t idxY  = edgedof::macroface::verticalIndex( level, it.col(), it.row() );
            const uint_t idxXY = edgedof::macroface::diagonalIndex( level, it.col(), it.row() );

            const auto actualX  = faceData[idxX];
            const auto actualY  = faceData[idxY];
            const auto actualXY = faceData[idxXY];

            WALBERLA_CHECK_FLOAT_EQUAL( actualX, correctX, "face index: " << idxX );
            WALBERLA_CHECK_FLOAT_EQUAL( actualY, correctY, "face index: " << idxY );
            WALBERLA_CHECK_FLOAT_EQUAL( actualXY, correctXY, "face index: " << idxXY );
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Faces passed" )

      {
         // there is only one cell, equal to the reference tet
         const Cell& cell     = *storage->getCell( storage->getCellIDs()[0] );
         const auto  cellData = cell.getData( dofs->getCellDataID() )->getPointer( level );
         const auto  nEdges   = real_c( levelinfo::num_microedges_per_edge( level ) );

         const real_t correctX   = 1.0 / nEdges;
         const real_t correctY   = 3.0 / nEdges;
         const real_t correctZ   = 6.0 / nEdges;
         const real_t correctXY  = 2.0 / nEdges;
         const real_t correctXZ  = 5.0 / nEdges;
         const real_t correctYZ  = 3.0 / nEdges;
         const real_t correctXYZ = 4.0 / nEdges;

         for ( const auto& it : edgedof::macrocell::Iterator( level ) )
         {
            const uint_t idxX  = edgedof::macrocell::xIndex( level, it.x(), it.y(), it.z() );
            const uint_t idxY  = edgedof::macrocell::yIndex( level, it.x(), it.y(), it.z() );
            const uint_t idxZ  = edgedof::macrocell::zIndex( level, it.x(), it.y(), it.z() );
            const uint_t idxXY = edgedof::macrocell::xyIndex( level, it.x(), it.y(), it.z() );
            const uint_t idxXZ = edgedof::macrocell::xzIndex( level, it.x(), it.y(), it.z() );
            const uint_t idxYZ = edgedof::macrocell::yzIndex( level, it.x(), it.y(), it.z() );

            const auto actualX  = cellData[idxX];
            const auto actualY  = cellData[idxY];
            const auto actualZ  = cellData[idxZ];
            const auto actualXY = cellData[idxXY];
            const auto actualXZ = cellData[idxXZ];
            const auto actualYZ = cellData[idxYZ];

            WALBERLA_CHECK_FLOAT_EQUAL( actualX, correctX, "index: " << it );
            WALBERLA_CHECK_FLOAT_EQUAL( actualY, correctY, "index: " << it );
            WALBERLA_CHECK_FLOAT_EQUAL( actualZ, correctZ, "index: " << it );
            WALBERLA_CHECK_FLOAT_EQUAL( actualXY, correctXY, "index: " << it );
            WALBERLA_CHECK_FLOAT_EQUAL( actualXZ, correctXZ, "index: " << it );
            WALBERLA_CHECK_FLOAT_EQUAL( actualYZ, correctYZ, "index: " << it );
         }

         for ( const auto& it : edgedof::macrocell::IteratorXYZ( level ) )
         {
            const uint_t idxXYZ    = edgedof::macrocell::xyzIndex( level, it.x(), it.y(), it.z() );
            const auto   actualXYZ = cellData[idxXYZ];
            WALBERLA_CHECK_FLOAT_EQUAL( actualXYZ, correctXYZ, "cell index: " << idxXYZ );
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Cell passed" )
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   test3D();

   return EXIT_SUCCESS;
}
