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

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::int_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void ARHangingNodesStorageTest()
{
   MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_184el.msh" );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t numRefinements = 3;

   writeDomainPartitioningVTK( *storage, "../../output/", "ARHangingNodesStorageTest_Domain_Refinement_0" );

   auto refinementLevelAtPoint = []( const Point3D& x ) -> uint_t {
      // refine near (0, 0)
      auto distance = x.norm();

      auto targetLevel = 3 - int_c( distance / 0.2 );
      if ( targetLevel < 0 )
         targetLevel = 0;
      return uint_c( targetLevel );
   };

   std::map< uint_t, uint_t > numFaces;

   for ( uint_t i = 1; i <= numRefinements; i++ )
   {
      // Collect all volume primitives that shall be refined.
      std::vector< PrimitiveID > refine, coarsen, refineResult, coarsenResult;

      for ( const auto& faceID : storage->getFaceIDs() )
      {
         auto face = storage->getFace( faceID );

         auto currentLevel = face->getID().numAncestors();

         uint_t targetLevel = currentLevel;

         for ( const auto& c : face->getCoordinates() )
         {
            targetLevel = std::max( targetLevel, refinementLevelAtPoint( c ) );
         }

         if ( targetLevel > currentLevel )
         {
            refine.push_back( face->getID() );
         }
      }

      // Refine all primitives
      storage->refinementAndCoarseningHanging( refine, coarsen, refineResult, coarsenResult );

      numFaces[i] = storage->getNumberOfGlobalFaces();
      WALBERLA_CHECK_GREATER( numFaces[i], numFaces[i - 1] );

      writeDomainPartitioningVTK(
          *storage, "../../output/", "ARHangingNodesStorageTest_Domain_Refinement_" + std::to_string( i ) );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::ARHangingNodesStorageTest();
}