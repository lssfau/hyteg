/*
* Copyright (c) 2025 Benjamin Mann.
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
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/globalIndices.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

namespace hyteg {

/// This function tests the conversion between cube and tet indexing
void test()
{
   const uint_t lvl = 4;
   for ( const auto& micro : celldof::macrocell::Iterator( lvl, celldof::CellType::WHITE_DOWN, 0 ) )
   {
      std::array< uint_t, 8 > cubeIndices{};
      std::array< uint_t, 4 > tetIndices{};

      p1::getGlobalCubeIndices3D( lvl, micro, cubeIndices );

      for ( auto& ct : celldof::allCellTypes )
      {
         const auto& cts = celldof::CellTypeToStr.at( ct );
         vertexdof::getVertexDoFDataIndicesFromMicroCell( micro, ct, lvl, tetIndices );
         for ( uint_t i = 0; i < 4; ++i )
         {
            // check if tetIndices is a subset of cubeIndices
            bool found = std::find( cubeIndices.begin(), cubeIndices.end(), tetIndices[i] ) != cubeIndices.end();
            WALBERLA_CHECK( found, walberla::format( "%s not contained in the micro cube", cts.c_str() ) );
            // check if the mapping between tet and cube indices is correct
            auto j = p1::cubeIndicesFromCellIndices[ct][i];
            WALBERLA_CHECK_EQUAL( tetIndices[i],
                                  cubeIndices[j],
                                  walberla::format( "Mapping between %s-indices and cube-indices incorrect", cts.c_str() ) );
         }
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::mpi::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::test();
}
