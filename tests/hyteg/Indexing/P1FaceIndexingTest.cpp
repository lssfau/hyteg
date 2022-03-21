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

#include <hyteg/p1functionspace/VertexDoFMacroFace.hpp>
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"

#include "core/mpi/all.h"

namespace hyteg {

using indexing::Index;

void checkIndices(idx_t col, idx_t row, std::vector<uint_t> ref, uint_t type){
  std::vector<size_t> result;
  switch(type){
    //vertex
    case 0:
      for( const auto & n : hyteg::vertexdof::macroface::neighborsWithCenter )
      {
        result.push_back( hyteg::vertexdof::macroface::indexFromVertex( 3, col, row, n ));
      }
      break;
    case 1:
      for( const auto & n : hyteg::vertexdof::macroface::neighborsFromGrayFace )
      {
        result.push_back( hyteg::vertexdof::macroface::indexFromGrayFace( 3, col, row, n ));
      }
      break;
    case 2:
      for( const auto & n : hyteg::vertexdof::macroface::neighborsFromBlueFace )
      {
        result.push_back( hyteg::vertexdof::macroface::indexFromBlueFace( 3, col, row, n ));
      }
      break;
    default:
      WALBERLA_ABORT("wrong type")
  }

  for(size_t i = 0; i < ref.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(ref[i],result[i],"col: " << col << " row: " << row <<" i: " << i);
  }
}

void testNeighborhood2( const Index & index, const std::map< stencilDirection, uint_t > & dirToLinearIndex )
{
  const uint_t level = 2;
  for ( const auto & it : dirToLinearIndex )
  {
    const auto dir         = it.first;
    const auto linearIndex = it.second;
    WALBERLA_CHECK_EQUAL( vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), dir ), linearIndex );
  }
}

void testNeighborhood3( const Index & index, const std::map< stencilDirection, uint_t > & dirToLinearIndex )
{
  const uint_t level = 3;
  for ( const auto & it : dirToLinearIndex )
  {
    const auto dir         = it.first;
    const auto linearIndex = it.second;
    WALBERLA_CHECK_EQUAL( vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), dir ), linearIndex );
  }
}

}

//const Dir neighbors_with_center[] = {S, SE, W, C, E, NW, N};
int main(int argc, char* argv[])
{
  typedef hyteg::stencilDirection sd;

  //this test is written for level 3
  walberla::mpi::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  /// CHECK VERTEX INDEXING ///
  std::vector<size_t> refOneOne = {10,1,2,11,18,17,9};
  hyteg::checkIndices(1,1,refOneOne,0);
  std::vector<size_t> refFiveTwo = {22,14,15,23,29,28,21};
  hyteg::checkIndices(5,2,refFiveTwo,0);

  /// CHECK CELL GRAY INDEXING ///
  hyteg::checkIndices(2,2,{19,20,26},1);
  hyteg::checkIndices(2,4,{32,33,37},1);

  /// CHECK CELL BLUE INDEXING ///
  hyteg::checkIndices(3,3,{28,33,34},2);
  hyteg::checkIndices(4,1,{14,21,22},2);

  std::map< sd, uint_t > level2FirstInner;

  level2FirstInner[ sd::VERTEX_C ]   =  6;
  level2FirstInner[ sd::VERTEX_E ]   =  7;
  level2FirstInner[ sd::VERTEX_W ]   =  5;
  level2FirstInner[ sd::VERTEX_N ]   = 10;
  level2FirstInner[ sd::VERTEX_S ]   =  1;

  level2FirstInner[ sd::VERTEX_NW ]  =  9;
  level2FirstInner[ sd::VERTEX_SE ]  =  2;

  level2FirstInner[ sd::VERTEX_TC ]  = 20;
  level2FirstInner[ sd::VERTEX_TW ]  = 19;
  level2FirstInner[ sd::VERTEX_TS ]  = 16;
  level2FirstInner[ sd::VERTEX_TSE ] = 17;

  level2FirstInner[ sd::VERTEX_BC ]  = 30;
  level2FirstInner[ sd::VERTEX_BE ]  = 31;
  level2FirstInner[ sd::VERTEX_BN ]  = 33;
  level2FirstInner[ sd::VERTEX_BNW ] = 32;

  hyteg::testNeighborhood2( hyteg::indexing::Index( 1, 1, 1 ), level2FirstInner );

  std::map< sd, uint_t > level2SecondInner;

  level2SecondInner[ sd::VERTEX_C ]   =  7;
  level2SecondInner[ sd::VERTEX_E ]   =  8;
  level2SecondInner[ sd::VERTEX_W ]   =  6;
  level2SecondInner[ sd::VERTEX_N ]   = 11;
  level2SecondInner[ sd::VERTEX_S ]   =  2;

  level2SecondInner[ sd::VERTEX_NW ]  = 10;
  level2SecondInner[ sd::VERTEX_SE ]  =  3;

  level2SecondInner[ sd::VERTEX_TC ]  = 21;
  level2SecondInner[ sd::VERTEX_TW ]  = 20;
  level2SecondInner[ sd::VERTEX_TS ]  = 17;
  level2SecondInner[ sd::VERTEX_TSE ] = 18;

  level2SecondInner[ sd::VERTEX_BC ]  = 31;
  level2SecondInner[ sd::VERTEX_BE ]  = 32;
  level2SecondInner[ sd::VERTEX_BN ]  = 34;
  level2SecondInner[ sd::VERTEX_BNW ] = 33;

  hyteg::testNeighborhood2( hyteg::indexing::Index( 2, 1, 1 ), level2SecondInner );

  std::map< sd, uint_t > level3FirstInner;

  level3FirstInner[ sd::VERTEX_C ]   = 10;
  level3FirstInner[ sd::VERTEX_E ]   = 11;
  level3FirstInner[ sd::VERTEX_W ]   =  9;
  level3FirstInner[ sd::VERTEX_N ]   = 18;
  level3FirstInner[ sd::VERTEX_S ]   =  1;

  level3FirstInner[ sd::VERTEX_NW ]  = 17;
  level3FirstInner[ sd::VERTEX_SE ]  =  2;

  level3FirstInner[ sd::VERTEX_TC ]  = 54;
  level3FirstInner[ sd::VERTEX_TW ]  = 53;
  level3FirstInner[ sd::VERTEX_TS ]  = 46;
  level3FirstInner[ sd::VERTEX_TSE ] = 47;

  level3FirstInner[ sd::VERTEX_BC ]  = 90;
  level3FirstInner[ sd::VERTEX_BE ]  = 91;
  level3FirstInner[ sd::VERTEX_BN ]  = 97;
  level3FirstInner[ sd::VERTEX_BNW ] = 96;

  hyteg::testNeighborhood3( hyteg::indexing::Index( 1, 1, 1 ), level3FirstInner );

}
