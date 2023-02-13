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
#include "core/Environment.h"

#include "hyteg/primitives/all.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

using namespace hyteg;

int main( int argc, char* argv[] )
{
   //this test is written for level 3
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   /// CHECK VERTEX this is the same as for bubble///
   std::vector< size_t > refOneOne  = { 1, 9, 8, 37, 43, 36 };
   std::vector< size_t > refTwoFive = { 28, 32, 31, 60, 62, 59 };
   std::vector< size_t > result;
   for ( auto n : facedof::macroface::neighbors )
   {
      size_t idx = facedof::macroface::indexFaceFromVertex( 3, 1, 1, n );
      result.push_back( idx );
      //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
   }
   for ( size_t i = 0; i < refOneOne.size(); ++i )
   {
      WALBERLA_CHECK_EQUAL( refOneOne[i], result[i], "i: " << i );
   }
   result.clear();
   for ( auto n : facedof::macroface::neighbors )
   {
      size_t idx = facedof::macroface::indexFaceFromVertex( 3, 2, 5, n );
      result.push_back( idx );
      //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
   }
   for ( size_t i = 0; i < refTwoFive.size(); ++i )
   {
      WALBERLA_CHECK_EQUAL( refTwoFive[i], result[i], "i: " << i );
   }
   /// END CHECK VERTEX ///
   result.clear();

   WALBERLA_CHECK_EQUAL( 50, facedof::macroface::indexFaceFromGrayFace( 3, 1, 3, stencilDirection::CELL_BLUE_S ) );
   WALBERLA_CHECK_EQUAL( 55, facedof::macroface::indexFaceFromGrayFace( 3, 2, 3, stencilDirection::CELL_BLUE_W ) );
   WALBERLA_CHECK_EQUAL( 48, facedof::macroface::indexFaceFromGrayFace( 3, 5, 1, stencilDirection::CELL_BLUE_E ) );

   WALBERLA_CHECK_EQUAL( 27, facedof::macroface::indexFaceFromBlueFace( 3, 1, 3, stencilDirection::CELL_GRAY_N ) );
   WALBERLA_CHECK_EQUAL( 23, facedof::macroface::indexFaceFromBlueFace( 3, 2, 3, stencilDirection::CELL_GRAY_W ) );
   WALBERLA_CHECK_EQUAL( 14, facedof::macroface::indexFaceFromBlueFace( 3, 5, 1, stencilDirection::CELL_GRAY_E ) );
   WALBERLA_CHECK_EQUAL( 13, facedof::macroface::indexFaceFromBlueFace( 3, 5, 0, stencilDirection::CELL_GRAY_N ) );
   WALBERLA_CHECK_EQUAL( 14, facedof::macroface::indexFaceFromBlueFace( 3, 6, 0, stencilDirection::CELL_GRAY_N ) );
   WALBERLA_CHECK_EQUAL( 7, facedof::macroface::indexFaceFromBlueFace( 3, 6, 0, stencilDirection::CELL_GRAY_E ) );
   WALBERLA_CHECK_EQUAL( 6, facedof::macroface::indexFaceFromBlueFace( 3, 6, 0, stencilDirection::CELL_GRAY_W ) );

   /// EDGE INDEXING ///
   std::vector< size_t > refOne  = { 2, 17, 15, 0, 1, 16 };
   std::vector< size_t > refFive = { 10, 25, 23, 8, 9, 24 };
   for ( auto n : facedof::macroedge::neighbors )
   {
      size_t idx = facedof::macroedge::indexFaceFromVertex( 3, 1, n );
      result.push_back( idx );
      //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
   }
   for ( size_t i = 0; i < refOne.size(); ++i )
   {
      WALBERLA_CHECK_EQUAL( refOne[i], result[i], "i: " << i );
   }
   result.clear();
   for ( auto n : facedof::macroedge::neighbors )
   {
      size_t idx = facedof::macroedge::indexFaceFromVertex( 3, 5, n );
      result.push_back( idx );
      //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
   }
   for ( size_t i = 0; i < refFive.size(); ++i )
   {
      WALBERLA_CHECK_EQUAL( refFive[i], result[i], "i: " << i );
   }
   result.clear();

   WALBERLA_CHECK_EQUAL( 0, facedof::macroedge::indexFaceFromVertex( 3, 1, stencilDirection::CELL_GRAY_SW ) );
   WALBERLA_CHECK_EQUAL( 1, facedof::macroedge::indexFaceFromVertex( 3, 1, stencilDirection::CELL_BLUE_SE ) );
   WALBERLA_CHECK_EQUAL( 2, facedof::macroedge::indexFaceFromVertex( 3, 1, stencilDirection::CELL_GRAY_SE ) );
   WALBERLA_CHECK_EQUAL( 15, facedof::macroedge::indexFaceFromVertex( 3, 1, stencilDirection::CELL_GRAY_NW ) );
   WALBERLA_CHECK_EQUAL( 16, facedof::macroedge::indexFaceFromVertex( 3, 1, stencilDirection::CELL_BLUE_NW ) );
   WALBERLA_CHECK_EQUAL( 17, facedof::macroedge::indexFaceFromVertex( 3, 1, stencilDirection::CELL_GRAY_NE ) );

   WALBERLA_CHECK_EQUAL( 2, facedof::macroedge::indexFaceFromVertex( 3, 2, stencilDirection::CELL_GRAY_SW ) );
   WALBERLA_CHECK_EQUAL( 3, facedof::macroedge::indexFaceFromVertex( 3, 2, stencilDirection::CELL_BLUE_SE ) );
   WALBERLA_CHECK_EQUAL( 4, facedof::macroedge::indexFaceFromVertex( 3, 2, stencilDirection::CELL_GRAY_SE ) );
   WALBERLA_CHECK_EQUAL( 17, facedof::macroedge::indexFaceFromVertex( 3, 2, stencilDirection::CELL_GRAY_NW ) );
   WALBERLA_CHECK_EQUAL( 18, facedof::macroedge::indexFaceFromVertex( 3, 2, stencilDirection::CELL_BLUE_NW ) );
   WALBERLA_CHECK_EQUAL( 19, facedof::macroedge::indexFaceFromVertex( 3, 2, stencilDirection::CELL_GRAY_NE ) );

   WALBERLA_CHECK_EQUAL( 4, facedof::macroedge::indexFaceFromVertex( 3, 3, stencilDirection::CELL_GRAY_SW ) );
   WALBERLA_CHECK_EQUAL( 5, facedof::macroedge::indexFaceFromVertex( 3, 3, stencilDirection::CELL_BLUE_SE ) );
   WALBERLA_CHECK_EQUAL( 6, facedof::macroedge::indexFaceFromVertex( 3, 3, stencilDirection::CELL_GRAY_SE ) );
   WALBERLA_CHECK_EQUAL( 19, facedof::macroedge::indexFaceFromVertex( 3, 3, stencilDirection::CELL_GRAY_NW ) );
   WALBERLA_CHECK_EQUAL( 20, facedof::macroedge::indexFaceFromVertex( 3, 3, stencilDirection::CELL_BLUE_NW ) );
   WALBERLA_CHECK_EQUAL( 21, facedof::macroedge::indexFaceFromVertex( 3, 3, stencilDirection::CELL_GRAY_NE ) );

   WALBERLA_CHECK_EQUAL( 6, facedof::macroedge::indexFaceFromVertex( 3, 4, stencilDirection::CELL_GRAY_SW ) );
   WALBERLA_CHECK_EQUAL( 7, facedof::macroedge::indexFaceFromVertex( 3, 4, stencilDirection::CELL_BLUE_SE ) );
   WALBERLA_CHECK_EQUAL( 8, facedof::macroedge::indexFaceFromVertex( 3, 4, stencilDirection::CELL_GRAY_SE ) );
   WALBERLA_CHECK_EQUAL( 21, facedof::macroedge::indexFaceFromVertex( 3, 4, stencilDirection::CELL_GRAY_NW ) );
   WALBERLA_CHECK_EQUAL( 22, facedof::macroedge::indexFaceFromVertex( 3, 4, stencilDirection::CELL_BLUE_NW ) );
   WALBERLA_CHECK_EQUAL( 23, facedof::macroedge::indexFaceFromVertex( 3, 4, stencilDirection::CELL_GRAY_NE ) );

   WALBERLA_CHECK_EQUAL( 8, facedof::macroedge::indexFaceFromVertex( 3, 5, stencilDirection::CELL_GRAY_SW ) );
   WALBERLA_CHECK_EQUAL( 9, facedof::macroedge::indexFaceFromVertex( 3, 5, stencilDirection::CELL_BLUE_SE ) );
   WALBERLA_CHECK_EQUAL( 10, facedof::macroedge::indexFaceFromVertex( 3, 5, stencilDirection::CELL_GRAY_SE ) );
   WALBERLA_CHECK_EQUAL( 23, facedof::macroedge::indexFaceFromVertex( 3, 5, stencilDirection::CELL_GRAY_NW ) );
   WALBERLA_CHECK_EQUAL( 24, facedof::macroedge::indexFaceFromVertex( 3, 5, stencilDirection::CELL_BLUE_NW ) );
   WALBERLA_CHECK_EQUAL( 25, facedof::macroedge::indexFaceFromVertex( 3, 5, stencilDirection::CELL_GRAY_NE ) );

   WALBERLA_CHECK_EQUAL( 10, facedof::macroedge::indexFaceFromVertex( 3, 6, stencilDirection::CELL_GRAY_SW ) );
   WALBERLA_CHECK_EQUAL( 11, facedof::macroedge::indexFaceFromVertex( 3, 6, stencilDirection::CELL_BLUE_SE ) );
   WALBERLA_CHECK_EQUAL( 12, facedof::macroedge::indexFaceFromVertex( 3, 6, stencilDirection::CELL_GRAY_SE ) );
   WALBERLA_CHECK_EQUAL( 25, facedof::macroedge::indexFaceFromVertex( 3, 6, stencilDirection::CELL_GRAY_NW ) );
   WALBERLA_CHECK_EQUAL( 26, facedof::macroedge::indexFaceFromVertex( 3, 6, stencilDirection::CELL_BLUE_NW ) );
   WALBERLA_CHECK_EQUAL( 27, facedof::macroedge::indexFaceFromVertex( 3, 6, stencilDirection::CELL_GRAY_NE ) );

   WALBERLA_CHECK_EQUAL( 12, facedof::macroedge::indexFaceFromVertex( 3, 7, stencilDirection::CELL_GRAY_SW ) );
   WALBERLA_CHECK_EQUAL( 13, facedof::macroedge::indexFaceFromVertex( 3, 7, stencilDirection::CELL_BLUE_SE ) );
   WALBERLA_CHECK_EQUAL( 14, facedof::macroedge::indexFaceFromVertex( 3, 7, stencilDirection::CELL_GRAY_SE ) );
   WALBERLA_CHECK_EQUAL( 27, facedof::macroedge::indexFaceFromVertex( 3, 7, stencilDirection::CELL_GRAY_NW ) );
   WALBERLA_CHECK_EQUAL( 28, facedof::macroedge::indexFaceFromVertex( 3, 7, stencilDirection::CELL_BLUE_NW ) );
   WALBERLA_CHECK_EQUAL( 29, facedof::macroedge::indexFaceFromVertex( 3, 7, stencilDirection::CELL_GRAY_NE ) );
}
