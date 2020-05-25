/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include <core/Environment.h>

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"

namespace hyteg {

void testInverse()
{
   const uint_t level = 4;
   auto meshInfo = MeshInfo::meshSphericalShell( 5, 2, 0.5, 1.0 );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   writeDomainPartitioningVTK( storage, "../../output", "IcosahedralShellMapTestDomain" );

   for ( const auto & it : storage->getCells() )
   {
      auto cell = it.second;
      IcosahedralShellMap geometryMap( *cell, *setupStorage );
      for ( const auto & idx : vertexdof::macrocell::Iterator( level ) )
      {
         const auto position = vertexdof::macrocell::coordinateFromIndex( level, *cell, idx );
         Point3D mappedPosition;
         Point3D invMappedPosition;
         geometryMap.evalF( position, mappedPosition );
         geometryMap.evalFinv( mappedPosition, invMappedPosition );
         auto error = (position - invMappedPosition).norm();
         WALBERLA_CHECK_LESS( error, 1e-14 );
      }
   }
}

}

int main( int argc, char ** argv )
{
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::testInverse();
   return 0;
}