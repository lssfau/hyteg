/*
 * Copyright (c) 2023 Daniel Bauer.
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

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/n1e1functionspace/N1E1MacroCell.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using namespace hyteg;

void testPyramid()
{
   const MeshInfo              pyramid = MeshInfo::fromGmshFile( "../../meshes/3D/pyramid_2el.msh" );
   const SetupPrimitiveStorage setupStorage( pyramid, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const PrimitiveStorage      storage( setupStorage );

   std::shared_ptr< Cell > cell0, cell1;
   for ( const auto& cell : storage.getCells() )
   {
      if ( cell.second->getCoordinates()[1] == Point3D{ 0.0, 0, 0 } )
      {
         cell0 = cell.second;
      }
      else if ( cell.second->getCoordinates()[2] == Point3D{ 1.0, 1, 0 } )
      {
         cell1 = cell.second;
      }
      else
      {
         WALBERLA_ABORT( "Unexpected mesh" );
      }
   }

   // lvl 0
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 0, *cell0, { 0, 0, 0 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, 1, 1, 1, -1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 0, *cell1, { 0, 0, 0 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ -1, 1, 1, 1, 1, 1 } ) )

   // lvl 1
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell0, { 0, 0, 0 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, 1, 1, 1, -1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell0, { 1, 0, 0 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, 1, 1, 1, -1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell0, { 0, 1, 0 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, 1, 1, 1, -1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell0, { 0, 0, 1 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, 1, 1, 1, -1 } ) )

   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell1, { 0, 0, 0 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ -1, 1, 1, 1, 1, 1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell1, { 1, 0, 0 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ -1, 1, 1, 1, 1, 1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell1, { 0, 1, 0 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ -1, 1, 1, 1, 1, 1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell1, { 0, 0, 1 }, celldof::CellType::WHITE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ -1, 1, 1, 1, 1, 1 } ) )

   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell0, { 0, 0, 0 }, celldof::CellType::BLUE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, -1, 1, 1, 1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell0, { 0, 0, 0 }, celldof::CellType::BLUE_DOWN ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, -1, 1, 1, 1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell0, { 0, 0, 0 }, celldof::CellType::GREEN_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ -1, 1, 1, 1, 1, 1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell0, { 0, 0, 0 }, celldof::CellType::GREEN_DOWN ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, 1, 1, 1, -1 } ) )

   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell1, { 0, 0, 0 }, celldof::CellType::BLUE_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ -1, 1, 1, 1, 1, 1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell1, { 0, 0, 0 }, celldof::CellType::BLUE_DOWN ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, 1, 1, 1, -1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell1, { 0, 0, 0 }, celldof::CellType::GREEN_UP ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, -1, 1, 1, 1 } ) )
   WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( 1, *cell1, { 0, 0, 0 }, celldof::CellType::GREEN_DOWN ).diagonal(),
                         ( Eigen::Matrix< real_t, 6, 1 >{ 1, 1, -1, 1, 1, 1 } ) )
}

void testCube( const uint_t minLevel, const uint_t maxLevel )
{
   const MeshInfo              pyramid = MeshInfo::meshSymmetricCuboid( { 0, 0, 0 }, { 1, 1, 1 }, 1, 1, 1 );
   const SetupPrimitiveStorage setupStorage( pyramid, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   n1e1::N1E1VectorFunction< real_t > edgeDirs( "edge directions", storage, minLevel, maxLevel );
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      edgeDirs.getDoFs()->interpolate( 1, level );
      edgeDirs.communicate< Face, Cell >( level );
      edgeDirs.communicate< Edge, Cell >( level );

      for ( auto& macroIter : storage->getCells() )
      {
         Cell&               cell         = *macroIter.second;
         const real_t* const edgeDirsData = cell.getData( edgeDirs.getDoFs()->getCellDataID() )->getPointer( level );

         for ( const auto& cellType : celldof::allCellTypes )
         {
            for ( const auto& microCell : celldof::macrocell::Iterator( level, cellType, 0 ) )
            {
               std::array< uint_t, 6 > edgeDoFIndices;
               n1e1::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cellType, level, edgeDoFIndices );

               Eigen::Matrix< real_t, 6, 1 > reference;
               for ( uint_t k = 0; k < 6; ++k )
               {
                  reference[numeric_cast< Eigen::Index >( k )] = edgeDirsData[edgeDoFIndices[k]];
               }

               WALBERLA_CHECK_EQUAL( n1e1::macrocell::basisTransformation( level, cell, microCell, cellType ).diagonal(),
                                     reference,
                                     "level = " << level << std::endl
                                                << "cell = " << cell.getCoordinates()[0] << ", " << cell.getCoordinates()[1]
                                                << ", " << cell.getCoordinates()[2] << ", " << cell.getCoordinates()[3] << ", "
                                                << std::endl
                                                << "microCell = " << microCell << std::endl
                                                << "cellType = " << celldof::CellTypeToStr.at( cellType ) )
            }
         }
      }
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   testPyramid();
   testCube( 0, 4 );

   return EXIT_SUCCESS;
}
