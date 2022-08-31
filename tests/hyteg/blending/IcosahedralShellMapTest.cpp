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

#include "hyteg/geometry/IcosahedralShellMap.hpp"

#include <core/Environment.h>

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

namespace hyteg {

void testInverse( real_t radInnerShell, real_t radOuterShell )
{
   const uint_t level    = 4;
   auto         meshInfo = MeshInfo::meshSphericalShell( 5, 2, radInnerShell, radOuterShell );
   auto         setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   writeDomainPartitioningVTK( storage, "../../output", "IcosahedralShellMapTestDomain" );

   for ( const auto& it : storage->getCells() )
   {
      auto                cell = it.second;
      IcosahedralShellMap geometryMap( *cell, *setupStorage );
      for ( const auto& idx : vertexdof::macrocell::Iterator( level ) )
      {
         const auto position = vertexdof::macrocell::coordinateFromIndex( level, *cell, idx );
         Point3D    mappedPosition;
         Point3D    invMappedPosition;
         geometryMap.evalF( position, mappedPosition );
         geometryMap.evalFinv( mappedPosition, invMappedPosition );
         auto error = ( position - invMappedPosition ).norm();
         WALBERLA_CHECK_LESS( error, 5e-13 );
      }
   }
}

void testBlending( real_t radInnerShell, real_t radOuterShell )
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
   meshInfo          = MeshInfo::meshSphericalShell( 2, 3, radInnerShell, radOuterShell );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   IcosahedralShellMap::setMap( *setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );
}

} // namespace hyteg

using walberla::real_c;

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::testInverse( real_c( 0.5 ), real_c( 1.0 ) );

   // Check numerical stability for large radii, see issue #183
   hyteg::testBlending( real_c( 2e6 ), real_c( 4e6 ) );

   return 0;
}
