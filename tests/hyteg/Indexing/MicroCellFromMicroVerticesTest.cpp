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
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

namespace hyteg {

using celldof::CellType;
using indexing::Index;

/// This function tests the conversion from an array of four micro-vertex indices to the correct micro-cell (including type).
void test()
{
   std::vector< std::array< Index, 4 > > vertices = { {
       // WHITE_UP
       { Index( 0, 0, 0 ), Index( 1, 0, 0 ), Index( 0, 1, 0 ), Index( 0, 0, 1 ) },
       { Index( 1, 1, 1 ), Index( 2, 1, 1 ), Index( 1, 2, 1 ), Index( 1, 1, 2 ) },
       // WHITE_DOWN
       { Index( 0, 2, 1 ), Index( 1, 2, 1 ), Index( 1, 1, 1 ), Index( 1, 2, 0 ) },
       // BLUE_UP
       { Index( 0, 2, 0 ), Index( 1, 2, 0 ), Index( 1, 1, 0 ), Index( 1, 1, 1 ) },
       // BLUE_DOWN
       { Index( 0, 0, 3 ), Index( 1, 0, 3 ), Index( 0, 1, 3 ), Index( 0, 1, 2 ) },
       // GREEN_UP
       { Index( 3, 0, 1 ), Index( 2, 1, 0 ), Index( 2, 0, 1 ), Index( 3, 0, 0 ) },
       // GREEN_DOWN
       { Index( 0, 2, 1 ), Index( 0, 2, 2 ), Index( 1, 2, 1 ), Index( 1, 1, 2 ) },
   } };

   std::vector< Index > cellIdx = { { Index( 0, 0, 0 ),
                                      Index( 1, 1, 1 ),
                                      Index( 0, 1, 0 ),
                                      Index( 0, 1, 0 ),
                                      Index( 0, 0, 2 ),
                                      Index( 2, 0, 0 ),
                                      Index( 0, 1, 1 ) } };

   std::vector< CellType > cellType = { { CellType::WHITE_UP,
                                          CellType::WHITE_UP,
                                          CellType::WHITE_DOWN,
                                          CellType::BLUE_UP,
                                          CellType::BLUE_DOWN,
                                          CellType::GREEN_UP,
                                          CellType::GREEN_DOWN } };

   for ( uint_t i = 0; i < vertices.size(); i++ )
   {
      auto p = celldof::macrocell::getMicroCellFromMicroVertices( vertices[i] );

      auto idx = p.first;
      auto t   = p.second;

      WALBERLA_CHECK_EQUAL( idx, cellIdx[i] );
      WALBERLA_CHECK_EQUAL( t, cellType[i] );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::mpi::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::test();
}
