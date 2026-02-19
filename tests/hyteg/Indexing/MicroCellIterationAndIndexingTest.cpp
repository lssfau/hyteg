/*
 * Copyright (c) 2017-2026 Marcus Mohr.
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

#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/mpi/all.h"

#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

void runMicroCellTest()
{
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "TEST iterating over all micro-cells inside a macro-cell with double loop" );
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------------------------------------------" );

   bool verbose = false;

   uint_t level = 4;
   uint_t count = 0;

   for ( const auto& cType : celldof::allCellTypes )
   {
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "\nIterating over all micro-cells of type: " << celldof::CellTypeToStr.at( cType ) << "\n" );
         WALBERLA_LOG_INFO_ON_ROOT( "cells/row = " << celldof::macrocell::numCellsPerRowByType( level, cType ) );
      }

      std::array< uint_t, 4 > vertexDoFIndices;
      std::array< uint_t, 6 > edgeDoFIndices;

      for ( const auto& it : celldof::macrocell::Iterator( level, cType, 0 ) )
      {
         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "xIdx = " << it.x() << ", yIdx = " << it.y() << ", zIdx = " << it.z() );
         }
         auto microVertexIndices = celldof::macrocell::getMicroVerticesFromMicroCell( it, cType );

         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Vertex list of micro-cell:" );
         }
         for ( uint_t k = 0; k < 4; ++k )
         {
            indexing::Index vert = microVertexIndices[k];
            if ( verbose )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "(" << vert.x() << "," << vert.y() << "," << vert.z() << ")" );
            }
         }

         vertexdof::getVertexDoFDataIndicesFromMicroCell( it, cType, level, vertexDoFIndices );

         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Data indices of vertexDoFs of micro-cell: "
                                       << vertexDoFIndices[0] << " " << vertexDoFIndices[1] << " " << vertexDoFIndices[2] << " "
                                       << vertexDoFIndices[3] );
         }

         edgedof::getEdgeDoFDataIndicesFromMicroCell( it, cType, level, edgeDoFIndices );

         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Data indices of edgeDoFs of micro-cell:   "
                                       << edgeDoFIndices[0] << " " << edgeDoFIndices[1] << " " << edgeDoFIndices[2] << " "
                                       << edgeDoFIndices[3] << " " << edgeDoFIndices[4] << " " << edgeDoFIndices[5] );
         }
         count++;
      }
   }

   // Did we get the correct number of micro-cells?
   WALBERLA_LOG_INFO_ON_ROOT( "Found " << count << " micro-cells altogether" );
   WALBERLA_CHECK_EQUAL( count, ( 1ul << 3 * level ) );

   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "TEST computing data indices of micro-cell's vertex and edge DoFs" );
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------------------------------------------" );

   // Select the top WHITE_UP cell and check its data dofs
   uint_t          width         = ( 1u << level ) + 1;
   uint_t          numVertexDoFs = ( ( width + 2 ) * ( width + 1 ) * width ) / 6;
   indexing::Index cIdx( 0, 0, idx_t( width - 2 ) );

   // check its vertex dof data indices
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( cIdx, celldof::CellType::WHITE_UP, level, vertexDoFIndices );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of VertexDoFs = " << numVertexDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "Data indices of top WHITE_UP cell (vertexdofs): "
                              << vertexDoFIndices[0] << " " << vertexDoFIndices[1] << " " << vertexDoFIndices[2] << " "
                              << vertexDoFIndices[3] );
   WALBERLA_CHECK_EQUAL( vertexDoFIndices[0], numVertexDoFs - 4 );
   WALBERLA_CHECK_EQUAL( vertexDoFIndices[1], numVertexDoFs - 3 );
   WALBERLA_CHECK_EQUAL( vertexDoFIndices[2], numVertexDoFs - 2 );
   WALBERLA_CHECK_EQUAL( vertexDoFIndices[3], numVertexDoFs - 1 );

   // check its edge dof data indices
   std::array< uint_t, 6 > edgeDoFIndices;
   uint_t                  numEdgesPerType = ( ( width - 1 ) * ( width + 1 ) * width ) / 6;
   uint_t                  numEdgeDoFs     = 6 * numEdgesPerType;
   WALBERLA_LOG_INFO_ON_ROOT( "Number of edges per type = " << numEdgesPerType );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of edgeDoFs = " << numEdgeDoFs );
   edgedof::getEdgeDoFDataIndicesFromMicroCell( cIdx, celldof::CellType::WHITE_UP, level, edgeDoFIndices );
   WALBERLA_LOG_INFO_ON_ROOT( "Data indices of top WHITE_UP cell (edgedofs): "
                              << edgeDoFIndices[0] << " " << edgeDoFIndices[1] << " " << edgeDoFIndices[2] << " "
                              << edgeDoFIndices[3] << " " << edgeDoFIndices[4] << " " << edgeDoFIndices[5] );

   WALBERLA_CHECK_EQUAL( edgeDoFIndices[0], numEdgeDoFs - 5 * numEdgesPerType - 1 );
   WALBERLA_CHECK_EQUAL( edgeDoFIndices[1], numEdgeDoFs - 4 * numEdgesPerType - 1 );
   WALBERLA_CHECK_EQUAL( edgeDoFIndices[2], numEdgeDoFs - 3 * numEdgesPerType - 1 );
   WALBERLA_CHECK_EQUAL( edgeDoFIndices[3], numEdgeDoFs - 2 * numEdgesPerType - 1 );
   WALBERLA_CHECK_EQUAL( edgeDoFIndices[4], numEdgeDoFs - 1 * numEdgesPerType - 1 );
   WALBERLA_CHECK_EQUAL( edgeDoFIndices[5], numEdgeDoFs - 0 * numEdgesPerType - 1 );
}

void runMicroEdgeOrientationTest()
{
   const real_t a      = real_c( 1.0 );
   const real_t deltaX = real_c( 2.0 );
   const real_t deltaY = real_c( 0.8 );

   MeshInfo meshInfo = MeshInfo::singleTriangle( Point2D{ a, a }, Point2D{ a + deltaX, a }, Point2D{ a + deltaX, a + deltaY } );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   using indexing::Index;
   Index microFaceIndex{ 2, 4, 0 };

   std::array< Index, 3 > controlGray;
   std::array< Index, 3 > controlBlue;
   std::array< Index, 3 > listGray;
   std::array< Index, 3 > listBlue;

   for ( uint_t checkCase = 0; checkCase <= 1; ++checkCase )
   {
      if ( checkCase == 0 )
      {
         controlGray = { Index( 2, 4, 0 ), Index( 3, 4, 0 ), Index( 2, 5, 0 ) };
         controlBlue = { Index( 3, 4, 0 ), Index( 2, 5, 0 ), Index( 3, 5, 0 ) };
      }
      else
      {
         controlGray = { Index( 2, 4, 0 ), Index( 3, 4, 0 ), Index( 2, 5, 0 ) };
         controlBlue = { Index( 3, 4, 0 ), Index( 3, 5, 0 ), Index( 2, 5, 0 ) };
      }

      if ( checkCase == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "* Testing with consistent edge orientation:" );
         listGray = facedof::macroface::getMicroVerticesFromMicroFace< true >( microFaceIndex, facedof::FaceType::GRAY );
         listBlue = facedof::macroface::getMicroVerticesFromMicroFace< true >( microFaceIndex, facedof::FaceType::BLUE );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "* Testing with anti-clockwise edge orientation:" );
         listGray = facedof::macroface::getMicroVerticesFromMicroFace< false >( microFaceIndex, facedof::FaceType::GRAY );
         listBlue = facedof::macroface::getMicroVerticesFromMicroFace< false >( microFaceIndex, facedof::FaceType::BLUE );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "  -> testing gray type" );
      for ( uint_t k = 0; k < 3; ++k )
      {
         WALBERLA_CHECK_EQUAL( controlGray[k], listGray[k] );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "  -> testing blue type" );
      for ( uint_t k = 0; k < 3; ++k )
      {
         WALBERLA_CHECK_EQUAL( controlBlue[k], listBlue[k] );
      }
   }
}

int main( int argc, char* argv[] )
{
   walberla::mpi::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "=========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Running MICRO CELL TEST" );
   WALBERLA_LOG_INFO_ON_ROOT( "=========================" );
   runMicroCellTest();

   WALBERLA_LOG_INFO_ON_ROOT( "=========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Running MICRO EDGE TEST" );
   WALBERLA_LOG_INFO_ON_ROOT( "=========================" );
   runMicroEdgeOrientationTest();
}
