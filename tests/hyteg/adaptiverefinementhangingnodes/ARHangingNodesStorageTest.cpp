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

#include "hyteg/Algorithms.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::int_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

/// Just testing the most basic features.
void smokeTest()
{
   MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_184el.msh" );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   const uint_t numRefinements = 3;

   writeDomainPartitioningVTK( *storage, "../../output/", "ARHangingNodesStorageTest_SmokeTest_Domain_Refinement_0" );

   auto refinementLevelAtPoint = []( const Point3D& x ) -> uint_t {
      // refine near (0, 0)
      auto distance = x.norm();

      auto targetLevel = 3 - int_c( distance / 0.2 );
      if ( targetLevel < 0 )
         targetLevel = 0;
      return uint_c( targetLevel );
   };

   std::map< uint_t, uint_t > numFaces;

   numFaces[0] = storage->getNumberOfGlobalFaces();

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
          *storage, "../../output/", "ARHangingNodesStorageTest_SmokeTest_Domain_Refinement_" + std::to_string( i ) );
   }
}

/// Checks indirect (volume-)neighborhood after refinement.
void indirectNeighborhoodTest()
{
   WALBERLA_CHECK_EQUAL(
       walberla::mpi::MPIManager::instance()->numProcesses(), 1, "Test currently only written for serial runs." );

   MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   const uint_t numRefinements = 5;

   WALBERLA_CHECK_LESS( numRefinements, storage->getNumberOfGlobalFaces(), "Cannot run test on this mesh." );

   writeDomainPartitioningVTK(
       *storage, "../../output/", "ARHangingNodesStorageTest_IndNeighTest_Domain_Refinement_" + std::to_string( 0 ) );

   for ( uint_t i = 1; i <= numRefinements; i++ )
   {
      // Refine a single macro-face on each currently available refinement level.
      const auto maxRefinementLevel = storage->getCurrentGlobalMaxRefinement();
      const auto faceIDs            = storage->getFaceIDs();

      std::vector< PrimitiveID > refine, coarsen;

      // Running over all available refinement levels.
      for ( uint_t refLevel = 0; refLevel < maxRefinementLevel; refLevel++ )
      {
         for ( const auto faceID : faceIDs )
         {
            if ( storage->getRefinementLevel( faceID ) == refLevel )
            {
               refine.push_back( faceID );
               break;
            }
         }
      }

      // Refine all selected primitives.
      storage->refinementAndCoarseningHanging( refine, coarsen );

      writeDomainPartitioningVTK(
          *storage, "../../output/", "ARHangingNodesStorageTest_IndNeighTest_Domain_Refinement_" + std::to_string( i ) );

      // Check neighborhood.
      for ( const auto& faceID : storage->getFaceIDs() )
      {
         const auto face = storage->getFace( faceID );

         std::vector< PrimitiveID > faceNeighbors;
         for ( const auto& [tmp0, nFaceID] : face->getIndirectNeighborFaceIDsOverEdges() )
         {
            faceNeighbors.push_back( nFaceID );
            const auto                 nface = storage->getFace( nFaceID );
            std::vector< PrimitiveID > nFaceNeighbors;
            for ( const auto& [tmp1, nFaceNFaceID] : nface->getIndirectNeighborFaceIDsOverEdges() )
            {
               nFaceNeighbors.push_back( nFaceNFaceID );
            }

            WALBERLA_LOG_INFO_ON_ROOT( "Face (pid = " << faceID << ") neighbors: " );
            for ( const auto& f : faceNeighbors )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "  " << f );
            }

            WALBERLA_LOG_INFO_ON_ROOT( "nFace (pid = " << nFaceID << ") neighbors: " );
            for ( const auto& f : nFaceNeighbors )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "  " << f );
            }

            WALBERLA_LOG_INFO_ON_ROOT( "" );

            WALBERLA_CHECK( algorithms::contains( nFaceNeighbors, faceID ), "Face is not a neighbor of its neighbors." );
         }
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   using namespace hyteg;

   smokeTest();
   indirectNeighborhoodTest();
}