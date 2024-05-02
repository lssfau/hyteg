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

#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"

namespace hyteg {

static void testVertexDoFFunctionMemorySize()
{
  using namespace edgedof;

  auto storage = PrimitiveStorage::createFromGmshFile( "../../meshes/3D/cube_24el.msh" );

  for ( const auto & it : storage->getEdges() )
  {
    const Edge & edge = *it.second;

    if ( edge.getNumNeighborFaces() == 1 )
    {
      WALBERLA_CHECK_EQUAL( vertexDoFMacroEdgeFunctionMemorySize( 2, edge ),  9 );
      WALBERLA_CHECK_EQUAL( vertexDoFMacroEdgeFunctionMemorySize( 3, edge ), 17 );
    }
    else if ( edge.getNumNeighborFaces() == 2 && edge.getNumNeighborCells() == 1 )
    {
      WALBERLA_CHECK_EQUAL( vertexDoFMacroEdgeFunctionMemorySize( 2, edge ), 16 );
      WALBERLA_CHECK_EQUAL( vertexDoFMacroEdgeFunctionMemorySize( 3, edge ), 32 );
    }
    else if ( edge.getNumNeighborFaces() == 3 && edge.getNumNeighborCells() == 3 )
    {
      WALBERLA_CHECK_EQUAL( vertexDoFMacroEdgeFunctionMemorySize( 2, edge ), 26 );
      WALBERLA_CHECK_EQUAL( vertexDoFMacroEdgeFunctionMemorySize( 3, edge ), 54 );
    }
    else if ( edge.getNumNeighborFaces() == 4 && edge.getNumNeighborCells() == 4 )
    {
      WALBERLA_CHECK_EQUAL( vertexDoFMacroEdgeFunctionMemorySize( 2, edge ), 33 );
      WALBERLA_CHECK_EQUAL( vertexDoFMacroEdgeFunctionMemorySize( 3, edge ), 69 );
    }
  }

  for ( const auto & it : storage->getFaces() )
  {
    const Face & face = *it.second;

    WALBERLA_CHECK_GREATER( face.getNumNeighborCells(), 0 );

    if ( face.getNumNeighborCells() == 1 )
    {
      WALBERLA_CHECK_EQUAL( vertexDoFMacroFaceFunctionMemorySize( 2, face ), 25 );
      WALBERLA_CHECK_EQUAL( vertexDoFMacroFaceFunctionMemorySize( 3, face ), 81 );
    }
    else if ( face.getNumNeighborCells() == 2 )
    {
      WALBERLA_CHECK_EQUAL( vertexDoFMacroFaceFunctionMemorySize( 2, face ),  35 );
      WALBERLA_CHECK_EQUAL( vertexDoFMacroFaceFunctionMemorySize( 3, face ), 117 );
    }
  }

  for ( const auto & it : storage->getCells() )
  {
    const Cell & cell = *it.second;
    WALBERLA_CHECK_EQUAL( vertexDoFMacroCellFunctionMemorySize( 2, cell ),  35 );
    WALBERLA_CHECK_EQUAL( vertexDoFMacroCellFunctionMemorySize( 3, cell ), 165 );
  }


}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testVertexDoFFunctionMemorySize();

   return EXIT_SUCCESS;
}
