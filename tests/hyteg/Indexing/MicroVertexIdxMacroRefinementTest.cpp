/*
* Copyright (c) 2017-2022 Nils Kohl.
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
#include "hyteg/indexing/Common.hpp"
#include "hyteg/volumedofspace/VolumeDoFIndexing.hpp"

namespace hyteg {

using celldof::CellType;
using indexing::Index;

/// This function tests the search of corresponding micro-vertex indices between coarse and refined macro-tets.
void test()
{
   const uint_t levelFine = 2;

   const auto widthFine   = levelinfo::num_microvertices_per_edge( levelFine );
   const auto widthCoarse = levelinfo::num_microvertices_per_edge( levelFine + 1 );

   std::vector< Index > vertices = {
       Index( 0, 0, 0 ),
       Index( widthFine - 1, 0, 0 ),
       Index( widthFine - 1, 0, 0 ),
       Index( widthCoarse - 1, 0, 0 ),
       Index( widthFine - 1, 0, widthFine - 1 ),
       Index( widthFine - 1, 0, widthFine - 2 ),
       Index( 0, widthFine - 1, 0 ),
   };

   std::vector< CellType > cellTypes = { CellType::WHITE_UP,
                                         CellType::WHITE_UP,
                                         CellType::WHITE_UP,
                                         CellType::WHITE_UP,
                                         CellType::BLUE_UP,
                                         CellType::BLUE_UP,
                                         CellType::GREEN_UP };
   std::vector< uint_t >   cellIdx   = { 0, 0, 1, 1, 0, 0, 0, 0 };

   std::vector< Index > result = { Index( 0, 0, 0 ),
                                   Index( widthFine - 1, 0, 0 ),
                                   Index( 0, 0, 0 ),
                                   Index( widthFine - 1, 0, 0 ),
                                   Index( 0, 0, widthFine - 1 ),
                                   Index( 0, 0, widthFine - 2 ),
                                   Index( widthFine - 1, 0, 0 ) };

   for ( uint_t i = 0; i < vertices.size(); i++ )
   {
      auto idxFine =
          volumedofspace::indexing::getMicroVertexIdxOnRefinedMacro( vertices[i], levelFine, cellTypes[i], cellIdx[i] );

      WALBERLA_CHECK_GREATER_EQUAL( idxFine[0], 0 );
      WALBERLA_CHECK_GREATER_EQUAL( idxFine[1], 0 );
      WALBERLA_CHECK_GREATER_EQUAL( idxFine[2], 0 );
      WALBERLA_CHECK_LESS( idxFine[0] + idxFine[1] + idxFine[2], widthFine );

      WALBERLA_CHECK_EQUAL( idxFine, result[i] );

      auto idxCoarse = volumedofspace::indexing::getMicroVertexIdxOnCoarserMacro( idxFine, levelFine, cellTypes[i], cellIdx[i] );

      WALBERLA_CHECK_EQUAL( idxCoarse, vertices[i] );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::mpi::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::test();
}
