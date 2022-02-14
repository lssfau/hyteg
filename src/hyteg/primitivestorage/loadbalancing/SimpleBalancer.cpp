/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include <queue>

#include "core/DataTypes.h"
#include "core/load_balancing/ParMetisWrapper.h"
#include "core/mpi/Broadcast.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/MPIWrapper.h"

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {
namespace loadbalancing {

using walberla::int64_c;
using walberla::int64_t;
using walberla::int_c;
using walberla::real_c;
using walberla::mpi::MPIRank;
using namespace walberla::mpistubs;

void parmetis( SetupPrimitiveStorage& setupStorage, uint_t subCommunicatorSize )
{
   const std::map< uint_t, int64_t > graphEdgeWeightsPerNumCommonVertices = {
       { 0, 0 },
       { 1, 1 },
       { 2, 10 },
       { 3, 100 },
   };

   MPI_Comm     communicator = walberla::mpi::MPIManager::instance()->comm();
   const uint_t rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );
   const uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   MPI_Comm  subCommunicator;
   const int inSubCommunicator = rank < subCommunicatorSize ? 1 : 0;
   MPI_Comm_split( communicator, inSubCommunicator, int_c( rank ), &subCommunicator );

   const bool hasGlobalCells = setupStorage.getNumberOfCells() > 0;

   roundRobin( setupStorage, subCommunicatorSize );

   if ( hasGlobalCells )
   {
      for ( const auto& it : setupStorage.getFaces() )
      {
         setupStorage.setTargetRank( it.first, numProcesses );
      }
      for ( const auto& it : setupStorage.getEdges() )
      {
         setupStorage.setTargetRank( it.first, numProcesses );
      }
      for ( const auto& it : setupStorage.getVertices() )
      {
         setupStorage.setTargetRank( it.first, numProcesses );
      }
   }

   // We calculate the mapping only on a subset of processes, and then broadcast them to all processes later on.
   std::vector< PrimitiveID::IDType > globalPrimitiveIDs;
   std::vector< uint_t >              globalRanks;

   if ( inSubCommunicator )
   {
      // Parameters needed by parmetis

      std::vector< int64_t > vtxdist;              // range of vertices that are local to each processor (identical on all ranks)
      std::vector< int64_t > xadj;                 // index and vtxdist indicate the vertex ID, entries the edge indices in adjncy
      std::vector< int64_t > adjncy;               // IDs of the neighboring vertices
      std::vector< int64_t > vwgt;                 // vertex weights
      std::vector< int64_t > adjwgt;               // edge weights
      int64_t                wgtflag;              // indicates types of weights
      int64_t                numflag;              // numbering scheme
      std::vector< int64_t > ndims;                // space dimensions
      std::vector< double >  xyz;                  // vertex coordinates
      int64_t                ncon;                 // number of weights per vertex
      int64_t                nparts;               // desired number of subdomains
      std::vector< double >  tpwgts;               // specifies vertex weight distribution
      std::vector< double >  ubvec;                // imbalance tolerance
      std::vector< int64_t > options;              // parmetis options
      std::vector< int64_t > edgecut;              // [out] edge cut
      std::vector< int64_t > part;                 // [out] resulting partition
      MPI_Comm*              parmetisCommunicator; // MPI communicator used by parmetis

      uint_t numLocalVolumePrimitives = 0;

      if ( hasGlobalCells )
      {
         for ( auto it : setupStorage.getCells() )
         {
            if ( setupStorage.getTargetRank( it.first ) == rank )
            {
               numLocalVolumePrimitives++;
            }
         }
      }
      else
      {
         WALBERLA_ABORT( "Not implemented for faces." )
      }

      //////////////////////
      // Building vtxdist //
      //////////////////////

      // Collecting number of primitives (== parmetis vertices) on each processor to build vtxdist
      std::vector< uint_t > numberOfLocalPrimitives;
      numberOfLocalPrimitives.push_back( numLocalVolumePrimitives );
      std::vector< uint_t > numberOfLocalPrimitivesOnProcesses =
          walberla::mpi::allGatherv( numberOfLocalPrimitives, subCommunicator );

      WALBERLA_ASSERT_EQUAL( numberOfLocalPrimitivesOnProcesses.size(), numProcesses );
      WALBERLA_ASSERT_EQUAL( numberOfLocalPrimitivesOnProcesses[rank], numLocalVolumePrimitives );

      int64_t sum = 0;
      for ( const auto& numberOfPrimitivesOnProcess : numberOfLocalPrimitivesOnProcesses )
      {
         vtxdist.push_back( sum );
         sum += int64_c( numberOfPrimitivesOnProcess );
      }
      vtxdist.push_back( sum );

      WALBERLA_ASSERT_EQUAL( vtxdist.size(), numProcesses + 1 );

      ///////////////////////////////////////////////////////////
      // Map PrimitiveIDs to global unique parmetis vertex IDs //
      ///////////////////////////////////////////////////////////

      // Creating a mapping from PrimitiveIDs of local primitives to a consecutive chunk of parmetis indices.
      // The chunks correspond to [ vtxdist[rank], vtxdist[rank+1] ).

      std::map< PrimitiveID::IDType, int64_t >
          localPrimitiveIDToGlobalParmetisIDMap; // contains all local PrimitiveIDs as keys and maps them to global parmetis IDs
      std::map< int64_t, PrimitiveID > globalParmetisIDToLocalPrimitiveIDMap; // reverse of the above map

      std::vector< PrimitiveID > localPrimitiveIDs;
      if ( hasGlobalCells )
      {
         for ( auto it : setupStorage.getCells() )
         {
            if ( setupStorage.getTargetRank( it.first ) == rank )
            {
               localPrimitiveIDs.push_back( it.first );
            }
         }
      }
      else
      {
         WALBERLA_ABORT( "Not implemented for faces." )
      }

      int64_t parmetisIDCounter = vtxdist[rank];
      for ( const auto& id : localPrimitiveIDs )
      {
         localPrimitiveIDToGlobalParmetisIDMap[id.getID()] = parmetisIDCounter;
         parmetisIDCounter++;
      }

      // Reverse the mapping (for convenience)
      for ( const auto& it : localPrimitiveIDToGlobalParmetisIDMap )
      {
         WALBERLA_ASSERT_EQUAL( globalParmetisIDToLocalPrimitiveIDMap.count( it.second ), 0 );
         globalParmetisIDToLocalPrimitiveIDMap[it.second] = PrimitiveID( it.first );
      }

      // To build the parmetis graph, we now need the mappings (PrimitiveID to parmetisID) from all neighboring processes

      std::set< MPIRank > neighboringRanks;
      for ( MPIRank r = 0; r < (MPIRank) subCommunicatorSize; r++ )
      {
         neighboringRanks.insert( r );
      }

      // Mapping neighboring process ranks to their ID mapping
      std::map< uint_t, std::map< PrimitiveID::IDType, int64_t > > neighboringPrimitiveIDToGlobalParmetisIDMaps;

      walberla::mpi::BufferSystem bufferSystem( communicator );
      bufferSystem.setReceiverInfo( neighboringRanks, true );

      for ( const MPIRank neighborRank : neighboringRanks )
      {
         bufferSystem.sendBuffer( neighborRank ) << localPrimitiveIDToGlobalParmetisIDMap;
      }
      bufferSystem.sendAll();
      for ( auto recv = bufferSystem.begin(); recv != bufferSystem.end(); ++recv )
      {
         recv.buffer() >> neighboringPrimitiveIDToGlobalParmetisIDMaps[uint_c( recv.rank() )];
      }

#ifndef NDEBUG
      for ( const MPIRank neighborRank : neighboringRanks )
      {
         WALBERLA_ASSERT_EQUAL( neighboringPrimitiveIDToGlobalParmetisIDMaps[uint_c( neighborRank )].size(),
                                numberOfLocalPrimitivesOnProcesses[uint_c( neighborRank )] );
      }
#endif

      //////////////////////////////
      // Building xadj and adjncy //
      //////////////////////////////

      numflag = 0;

      // Now that we got the assignment from PrimitiveIDs to parmetis IDs of the local and all neighboring processes, we can build the parmetis graph

      WALBERLA_ASSERT(
          std::is_sorted( globalParmetisIDToLocalPrimitiveIDMap.begin(), globalParmetisIDToLocalPrimitiveIDMap.end() ) );

      for ( const auto& it : globalParmetisIDToLocalPrimitiveIDMap )
      {
         const PrimitiveID          primitiveID = it.second;
         std::vector< PrimitiveID > neighborIDs;

         if ( hasGlobalCells )
         {
            const auto cell = setupStorage.getCell( primitiveID );
            neighborIDs.clear();
            for ( auto itt : cell->getIndirectNeighborCellIDs() )
            {
               neighborIDs.push_back( itt.second );
            }
            std::set< PrimitiveID > neighborIDsUnique( neighborIDs.begin(), neighborIDs.end() );
            neighborIDs = std::vector< PrimitiveID >( neighborIDsUnique.begin(), neighborIDsUnique.end() );
         }
         else
         {
            WALBERLA_ABORT( "Not implemented for faces." )
         }

         xadj.push_back( int64_c( adjncy.size() ) );

         for ( const auto& neighborID : neighborIDs )
         {
            uint_t  neighborRank;
            int64_t neighborParmetisID;

            neighborRank       = setupStorage.getTargetRank( neighborID );
            neighborParmetisID = neighboringPrimitiveIDToGlobalParmetisIDMaps[neighborRank][neighborID.getID()];

            // How are the neighboring volume primitives connected?
            // -> let's calculate number of common macro-vertices.

            std::set< PrimitiveID > adjVertices;

            const auto primitive         = setupStorage.getPrimitive( primitiveID );
            const auto neighborPrimitive = setupStorage.getPrimitive( neighborID );

            adjVertices.insert( primitive->neighborVertices().begin(), primitive->neighborVertices().end() );
            adjVertices.insert( neighborPrimitive->neighborVertices().begin(), neighborPrimitive->neighborVertices().end() );

            const auto commonVertices = 2 * primitive->getNumNeighborVertices() - adjVertices.size();

            WALBERLA_CHECK_LESS( commonVertices, 4 );

            adjncy.push_back( neighborParmetisID );

            // Depending on the connection (over macro-face, -edge, or -vertex) we add weights to the graph-edges.

            adjwgt.push_back( graphEdgeWeightsPerNumCommonVertices.at( commonVertices ) );
         }
      }
      xadj.push_back( int64_c( adjncy.size() ) );

      WALBERLA_ASSERT_EQUAL( xadj.size(), numLocalVolumePrimitives + 1 );

      //////////////////////////
      // Number of subdomains //
      //////////////////////////

      nparts = int64_c( numProcesses );

      /////////////////////////////
      // Vertex and edge weights //
      /////////////////////////////

      vwgt.resize( numLocalVolumePrimitives );
      std::fill( vwgt.begin(), vwgt.end(), 1 );

      wgtflag = int64_c( 1 );
      ncon    = int64_c( 1 );

      tpwgts.resize( uint_c( ncon * nparts ) );
      std::fill( tpwgts.begin(), tpwgts.end(), 1.0 / static_cast< double >( nparts ) );

      ubvec.resize( uint_c( ncon ) );
      std::fill( ubvec.begin(), ubvec.end(), 1.05 );

      //////////////////////
      // Parmetis options //
      //////////////////////

      options.push_back( 0 );
      options.push_back( 0 );
      options.push_back( 0 );

      ///////////////////
      // Output arrays //
      ///////////////////

      edgecut.resize( 1 );
      part.resize( numLocalVolumePrimitives );

      //////////////////////
      // MPI communicator //
      //////////////////////

      parmetisCommunicator = &subCommunicator;

      //////////////////////
      // Calling parmetis //
      //////////////////////

      int parmetisError = walberla::core::ParMETIS_V3_PartKway( vtxdist.data(),
                                                                xadj.data(),
                                                                adjncy.data(),
                                                                vwgt.data(),
                                                                adjwgt.data(),
                                                                &wgtflag,
                                                                &numflag,
                                                                &ncon,
                                                                &nparts,
                                                                tpwgts.data(),
                                                                ubvec.data(),
                                                                options.data(),
                                                                edgecut.data(),
                                                                part.data(),
                                                                parmetisCommunicator );

      WALBERLA_CHECK_EQUAL( parmetisError, walberla::core::METIS_OK );

      /////////////////////////
      // Primitive migration //
      /////////////////////////

      for ( uint_t partIdx = 0; partIdx < part.size(); partIdx++ )
      {
         const int64_t parmetisID  = vtxdist[rank] + int64_c( partIdx );
         const auto    primitiveID = globalParmetisIDToLocalPrimitiveIDMap[int64_c( parmetisID )];
         const auto    targetRank  = part[partIdx];
         globalPrimitiveIDs.push_back( primitiveID.getID() );
         globalRanks.push_back( uint_c( targetRank ) );
      }

      globalPrimitiveIDs = walberla::mpi::allGatherv( globalPrimitiveIDs, subCommunicator );
      globalRanks        = walberla::mpi::allGatherv( globalRanks, subCommunicator );

      WALBERLA_CHECK_EQUAL( globalPrimitiveIDs.size(), globalRanks.size() );
   }

   WALBERLA_MPI_BARRIER();

   // The volume primitive distribution is now finished on the process subset.
   // We distribute it now to all processes!

   walberla::mpi::broadcastObject( globalPrimitiveIDs, 0, communicator );
   walberla::mpi::broadcastObject( globalRanks, 0, communicator );

   WALBERLA_CHECK_EQUAL( globalPrimitiveIDs.size(), globalRanks.size() );

   for ( uint_t i = 0; i < globalPrimitiveIDs.size(); i++ )
   {
      setupStorage.setTargetRank( globalPrimitiveIDs[i], globalRanks[i] );
   }

   // Finally, we set the ranks of the lower dimensional primitives
   // depending on the ranks of the neighboring volume-primitives.

   if ( hasGlobalCells )
   {
      for ( auto it : setupStorage.getFaces() )
      {
         std::map< uint_t, uint_t > nbRanks;
         for ( const auto& nbCell : it.second->neighborCells() )
         {
            nbRanks[setupStorage.getTargetRank( nbCell )] += 1;
         }

         auto mostCommonNBRank = std::max_element( nbRanks.begin(),
                                                   nbRanks.end(),
                                                   []( std::pair< uint_t, uint_t > a, std::pair< uint_t, uint_t > b ) {
                                                      return a.second < b.second;
                                                   } )
                                     ->first;

         setupStorage.setTargetRank( it.first, mostCommonNBRank );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         std::map< uint_t, uint_t > nbRanks;
         for ( const auto& nbCell : it.second->neighborCells() )
         {
            nbRanks[setupStorage.getTargetRank( nbCell )] += 1;
         }

         auto mostCommonNBRank = std::max_element( nbRanks.begin(),
                                                   nbRanks.end(),
                                                   []( std::pair< uint_t, uint_t > a, std::pair< uint_t, uint_t > b ) {
                                                      return a.second < b.second;
                                                   } )
                                     ->first;

         setupStorage.setTargetRank( it.first, mostCommonNBRank );
      }

      for ( auto it : setupStorage.getVertices() )
      {
         std::map< uint_t, uint_t > nbRanks;
         for ( const auto& nbCell : it.second->neighborCells() )
         {
            nbRanks[setupStorage.getTargetRank( nbCell )] += 1;
         }

         auto mostCommonNBRank = std::max_element( nbRanks.begin(),
                                                   nbRanks.end(),
                                                   []( std::pair< uint_t, uint_t > a, std::pair< uint_t, uint_t > b ) {
                                                      return a.second < b.second;
                                                   } )
                                     ->first;

         setupStorage.setTargetRank( it.first, mostCommonNBRank );
      }
   }
   else
   {
      WALBERLA_ABORT( "Not implemented for faces." )
   }

   WALBERLA_LOG_DEVEL_ON_ROOT( setupStorage )
}

void parmetis( SetupPrimitiveStorage& setupStorage )
{
   parmetis( setupStorage, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
}

void allPrimitivesOnOneRank( SetupPrimitiveStorage& storage, const uint_t& targetRank )
{
   SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
   storage.getSetupPrimitives( setupPrimitives );
   for ( auto it : setupPrimitives )
   {
      storage.setTargetRank( it.first, uint_c( targetRank ) );
   }
}

void allPrimitivesOnRoot( SetupPrimitiveStorage& storage )
{
   allPrimitivesOnOneRank( storage, uint_c( 0 ) );
}

void roundRobin( SetupPrimitiveStorage& storage, uint_t numRanks )
{
   uint_t currentRank = 0;

   SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
   storage.getSetupPrimitives( setupPrimitives );

   for ( auto it : setupPrimitives )
   {
      storage.setTargetRank( it.first, uint_c( currentRank % numRanks ) );
      currentRank++;
   }
}

void roundRobin( SetupPrimitiveStorage& storage )
{
   roundRobin( storage, storage.getNumberOfProcesses() );
}

void roundRobinVolume( SetupPrimitiveStorage& storage, uint_t numRanks )
{
   SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
   storage.getSetupPrimitives( setupPrimitives );

   // Set every primitive to an invalid rank to allow counting primitives per rank later on.
   for ( const auto& it : setupPrimitives )
   {
      storage.setTargetRank( it.first, numRanks );
   }

   // Heuristic: let's assign primitives with consecutive IDs to the same process.
   // Often, meshes are generated in a way that primitives with consecutive IDs are located next to each other.
   if ( storage.getNumberOfCells() == 0 )
   {
      auto it = storage.getFaces().begin();
      for ( uint_t currentRank = 0; currentRank < numRanks; currentRank++ )
      {
         uint_t numVolumesOnThisRank = 0;
         while ( numVolumesOnThisRank < storage.getNumberOfFaces() / numRanks )
         {
            storage.setTargetRank( it->first, currentRank );
            numVolumesOnThisRank++;
            it++;
         }

         if ( currentRank < storage.getNumberOfFaces() % numRanks )
         {
            storage.setTargetRank( it->first, currentRank );
            numVolumesOnThisRank++;
            it++;
         }
      }
      WALBERLA_CHECK( it == storage.getFaces().end() );
   }
   else
   {
      auto it = storage.getCells().begin();
      for ( uint_t currentRank = 0; currentRank < numRanks; currentRank++ )
      {
         uint_t numVolumesOnThisRank = 0;
         while ( numVolumesOnThisRank < storage.getNumberOfCells() / numRanks )
         {
            storage.setTargetRank( it->first, currentRank );
            numVolumesOnThisRank++;
            it++;
         }

         if ( currentRank < storage.getNumberOfCells() % numRanks )
         {
            storage.setTargetRank( it->first, currentRank );
            numVolumesOnThisRank++;
            it++;
         }
      }
      WALBERLA_CHECK( it == storage.getCells().end() );
   }

   // Cache number of primitives per rank to optimize the performance of this function, by a lot.
   std::vector< uint_t > numPrimitivesPerRank( numRanks, 0 );

   // We assign lower-dimensional primitives to a rank of their higher-dimensional neighbors.
   // To equalize the weights a little, we choose the neighbor with least amount of primitives of the current type.
   if ( storage.getNumberOfCells() > 0 )
   {
      for ( const auto& it : storage.getFaces() )
      {
         const auto                 facePID = it.first;
         const auto                 face    = it.second;
         std::map< uint_t, uint_t > rankFaceCount;
         for ( const auto& neighborCell : face->neighborCells() )
         {
            const auto neighborRank           = storage.getTargetRank( neighborCell );
            const auto numFacesOnNeighborRank = numPrimitivesPerRank[neighborRank];
            rankFaceCount[neighborRank]       = numFacesOnNeighborRank;
         }

         auto leastFullNBRank = std::min_element( rankFaceCount.begin(),
                                                  rankFaceCount.end(),
                                                  []( std::pair< uint_t, uint_t > a, std::pair< uint_t, uint_t > b ) {
                                                     return a.second < b.second;
                                                  } )
                                    ->first;
         storage.setTargetRank( facePID, leastFullNBRank );
         numPrimitivesPerRank[leastFullNBRank]++;
      }
   }

   std::fill( numPrimitivesPerRank.begin(), numPrimitivesPerRank.end(), 0 );

   for ( const auto& it : storage.getEdges() )
   {
      const auto                 edgePID = it.first;
      const auto                 edge    = it.second;
      std::map< uint_t, uint_t > rankEdgeCount;
      for ( const auto& neighborFace : edge->neighborFaces() )
      {
         const auto neighborRank           = storage.getTargetRank( neighborFace );
         const auto numEdgesOnNeighborRank = numPrimitivesPerRank[neighborRank];
         rankEdgeCount[neighborRank]       = numEdgesOnNeighborRank;
      }

      auto leastFullNBRank =
          std::min_element( rankEdgeCount.begin(),
                            rankEdgeCount.end(),
                            []( std::pair< uint_t, uint_t > a, std::pair< uint_t, uint_t > b ) { return a.second < b.second; } )
              ->first;
      storage.setTargetRank( edgePID, leastFullNBRank );
      numPrimitivesPerRank[leastFullNBRank]++;
   }

   std::fill( numPrimitivesPerRank.begin(), numPrimitivesPerRank.end(), 0 );

   for ( const auto& it : storage.getVertices() )
   {
      const auto                 vertexPID = it.first;
      const auto                 vertex    = it.second;
      std::map< uint_t, uint_t > rankVertexCount;
      for ( const auto& neighborEdge : vertex->neighborEdges() )
      {
         const auto neighborRank              = storage.getTargetRank( neighborEdge );
         const auto numVerticesOnNeighborRank = numPrimitivesPerRank[neighborRank];
         rankVertexCount[neighborRank]        = numVerticesOnNeighborRank;
      }

      auto leastFullNBRank =
          std::min_element( rankVertexCount.begin(),
                            rankVertexCount.end(),
                            []( std::pair< uint_t, uint_t > a, std::pair< uint_t, uint_t > b ) { return a.second < b.second; } )
              ->first;
      storage.setTargetRank( vertexPID, leastFullNBRank );
      numPrimitivesPerRank[leastFullNBRank]++;
   }
}

void roundRobinVolume( SetupPrimitiveStorage& storage )
{
   roundRobinVolume( storage, storage.getNumberOfProcesses() );
}

void greedy( SetupPrimitiveStorage& storage )
{
   const uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   const uint_t maxVertices = storage.getNumberOfVertices() / numProcesses;
   const uint_t maxEdges    = storage.getNumberOfEdges() / numProcesses;
   const uint_t maxFaces    = storage.getNumberOfFaces() / numProcesses;
   const uint_t maxCells    = storage.getNumberOfCells() / numProcesses;

   // Set all target ranks to zero
   SetupPrimitiveStorage::PrimitiveMap primitives;
   storage.getSetupPrimitives( primitives );

   for ( const auto& it : primitives )
   {
      storage.setTargetRank( it.first, 0 );
   }

   // Main algorithm
   for ( uint_t rank = 1; rank < numProcesses; rank++ )
   {
      uint_t currentNumVertices = 0;
      uint_t currentNumEdges    = 0;
      uint_t currentNumFaces    = 0;
      uint_t currentNumCells    = 0;

      std::queue< PrimitiveID >             nextPrimitives;
      std::map< PrimitiveID::IDType, bool > wasAddedToQueue;

      while ( true )
      {
         if ( currentNumVertices >= maxVertices && currentNumEdges >= maxEdges && currentNumFaces >= maxFaces &&
              currentNumCells >= maxCells )
         {
            break;
         }

         // Push a random, unassigned primitive to the queue if the queue is empty and we are not finished
         if ( nextPrimitives.size() == 0 )
         {
            for ( const auto& it : primitives )
            {
               if ( storage.getTargetRank( it.first ) == 0 && !wasAddedToQueue[it.first] )
               {
                  if ( ( storage.vertexExists( it.first ) && currentNumVertices < maxVertices ) ||
                       ( storage.edgeExists( it.first ) && currentNumEdges < maxEdges ) ||
                       ( storage.faceExists( it.first ) && currentNumFaces < maxFaces ) ||
                       ( storage.cellExists( it.first ) && currentNumCells < maxCells ) )
                  {
                     nextPrimitives.push( it.first );
                     wasAddedToQueue[it.first] = true;
                     break;
                  }
               }
            }
         }

         WALBERLA_ASSERT_GREATER( nextPrimitives.size(), 0 );

         // Pop a primitive from the queue
         const Primitive* currentPrimitive = storage.getPrimitive( nextPrimitives.front() );
         nextPrimitives.pop();

         // Set the target rank to the current process if the process does not already carry enough primitives of that type.
         // Then set the primitive to visited. Otherwise, the next primitive is pulled from the queue.
         if ( storage.getTargetRank( currentPrimitive->getID().getID() ) == 0 )
         {
            // Check primitive type and assign rank
            if ( storage.vertexExists( currentPrimitive->getID() ) && currentNumVertices < maxVertices )
            {
               storage.setTargetRank( currentPrimitive->getID().getID(), rank );
               currentNumVertices++;
            }
            else if ( storage.edgeExists( currentPrimitive->getID() ) && currentNumEdges < maxEdges )
            {
               storage.setTargetRank( currentPrimitive->getID().getID(), rank );
               currentNumEdges++;
            }
            else if ( storage.faceExists( currentPrimitive->getID() ) && currentNumFaces < maxFaces )
            {
               storage.setTargetRank( currentPrimitive->getID().getID(), rank );
               currentNumFaces++;
            }
            else if ( storage.cellExists( currentPrimitive->getID() ) && currentNumCells < maxCells )
            {
               storage.setTargetRank( currentPrimitive->getID().getID(), rank );
               currentNumCells++;
            }
         }

         // Put neighboring primitives into queue
         std::vector< PrimitiveID > neighbors;
         currentPrimitive->getNeighborPrimitives( neighbors );

         for ( const auto& neighborID : neighbors )
         {
            if ( storage.getTargetRank( neighborID ) == 0 && !wasAddedToQueue[neighborID.getID()] )
            {
               nextPrimitives.push( neighborID );
               wasAddedToQueue[neighborID.getID()] = true;
            }
         }
      }
   }
}

} // namespace loadbalancing
} // namespace hyteg
