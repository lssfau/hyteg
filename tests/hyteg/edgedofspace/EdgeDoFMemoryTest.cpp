/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr.
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
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"

namespace hyteg {

static void testEdgeDoFFunctionMemorySize()
{
   using namespace edgedof;

   auto storage = PrimitiveStorage::createFromGmshFile( prependHyTeGMeshDir( "2D/quad_8el.msh" ) );

   for ( const auto& it : storage->getEdges() )
   {
      const Edge& edge = *it.second;

      if ( edge.getNumNeighborFaces() == 1 )
      {
         WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 2, edge ), 15 );
         WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 3, edge ), 31 );
      }
      else if ( edge.getNumNeighborFaces() == 2 )
      {
         WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 2, edge ), 26 );
         WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 3, edge ), 54 );
      }
   }

   for ( const auto& it : storage->getFaces() )
   {
      const Face& face = *it.second;

      WALBERLA_CHECK_EQUAL( edgeDoFMacroFaceFunctionMemorySize( 2, face ), 30 );
      WALBERLA_CHECK_EQUAL( edgeDoFMacroFaceFunctionMemorySize( 3, face ), 108 );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testEdgeDoFFunctionMemorySize();

   return EXIT_SUCCESS;
}
