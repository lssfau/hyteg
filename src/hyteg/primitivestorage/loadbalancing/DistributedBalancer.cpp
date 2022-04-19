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

#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"

#include <algorithm>
#include <queue>

#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"
#include "core/load_balancing/ParMetisWrapper.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/MPIWrapper.h"

#include "hyteg/communication/PackageBufferSystem.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

namespace hyteg {
namespace loadbalancing {
namespace distributed {

using walberla::double_c;
using walberla::int64_c;
using walberla::int64_t;
using walberla::int_c;
using walberla::mpi::MPIRank;
using namespace walberla::mpistubs;

MigrationMap_T parmetis( PrimitiveStorage& storage )
{
   WALBERLA_CHECK_GREATER(
       storage.getNumberOfLocalPrimitives(), 0, "ParMeTis not supported (yet) for distributions with empty processes." );

   MPI_Comm communicator = walberla::mpi::MPIManager::instance()->comm();
   uint_t   rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );
   uint_t   numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

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

   //////////////////////
   // Building vtxdist //
   //////////////////////

   // Collecting number of primitives (== parmetis vertices) on each processor to build vtxdist
   std::vector< uint_t > numberOfLocalPrimitives;
   numberOfLocalPrimitives.push_back( storage.getNumberOfLocalPrimitives() );
   std::vector< uint_t > numberOfLocalPrimitivesOnProcesses = walberla::mpi::allGatherv( numberOfLocalPrimitives, communicator );

   WALBERLA_ASSERT_EQUAL( numberOfLocalPrimitivesOnProcesses.size(), numProcesses );
   WALBERLA_ASSERT_EQUAL( numberOfLocalPrimitivesOnProcesses[rank], storage.getNumberOfLocalPrimitives() );

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

   std::map< PrimitiveID, int64_t >
                                    localPrimitiveIDToGlobalParmetisIDMap; // contains all local PrimitiveIDs as keys and maps them to global parmetis IDs
   std::map< int64_t, PrimitiveID > globalParmetisIDToLocalPrimitiveIDMap; // reverse of the above map

   std::vector< PrimitiveID > localPrimitiveIDs;
   storage.getPrimitiveIDs( localPrimitiveIDs );

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
   storage.getNeighboringRanks( neighboringRanks );

   // Mapping neighboring process ranks to their ID mapping
   std::map< uint_t, std::map< PrimitiveID, int64_t > > neighboringPrimitiveIDToGlobalParmetisIDMaps;

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
      const PrimitiveID primitiveID = it.second;
      const Primitive*  primitive   = storage.getPrimitive( primitiveID );

      std::vector< PrimitiveID > neighborIDs;
      primitive->getNeighborPrimitives( neighborIDs );

      xadj.push_back( int64_c( adjncy.size() ) );

      for ( const auto& neighborID : neighborIDs )
      {
         uint_t  neighborRank;
         int64_t neighborParmetisID;
         ;
         if ( storage.primitiveExistsInNeighborhood( neighborID ) )
         {
            neighborRank       = storage.getNeighborPrimitiveRank( neighborID );
            neighborParmetisID = neighboringPrimitiveIDToGlobalParmetisIDMaps[neighborRank][neighborID.getID()];
         }
         else
         {
            WALBERLA_ASSERT( storage.primitiveExistsLocally( neighborID ) );
            neighborRank       = rank;
            neighborParmetisID = localPrimitiveIDToGlobalParmetisIDMap[neighborID.getID()];
         }

         adjncy.push_back( neighborParmetisID );
      }
   }
   xadj.push_back( int64_c( adjncy.size() ) );

   WALBERLA_ASSERT_EQUAL( xadj.size(), storage.getNumberOfLocalPrimitives() + 1 );

   //////////////////////////
   // Number of subdomains //
   //////////////////////////

   nparts = int64_c( numProcesses );

   /////////////////////////////
   // Vertex and edge weights //
   /////////////////////////////

   // TODO

   vwgt.resize( storage.getNumberOfLocalPrimitives() );
   std::fill( vwgt.begin(), vwgt.end(), 1 );

   wgtflag = int64_c( 2 );
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
   part.resize( storage.getNumberOfLocalPrimitives() );

   //////////////////////
   // MPI communicator //
   //////////////////////

   parmetisCommunicator = &communicator;

   //////////////////////
   // Calling parmetis //
   //////////////////////

   int parmetisError = walberla::core::ParMETIS_V3_PartKway( vtxdist.data(),
                                                             xadj.data(),
                                                             adjncy.data(),
                                                             vwgt.data(),
                                                             /* adjwgt */ NULL,
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

   MigrationMap_T migrationMap;
   for ( uint_t partIdx = 0; partIdx < part.size(); partIdx++ )
   {
      const int64_t     parmetisID  = vtxdist[rank] + int64_c( partIdx );
      const PrimitiveID primitiveID = globalParmetisIDToLocalPrimitiveIDMap[parmetisID];

      migrationMap[primitiveID.getID()] = uint_c( part[partIdx] );
   }

   const auto numReceivingPrimitives = getNumReceivingPrimitives( migrationMap );

   MigrationInfo migrationInfo( migrationMap, numReceivingPrimitives );
   storage.migratePrimitives( migrationInfo );

   return migrationMap;
}

MigrationInfo roundRobin( PrimitiveStorage& storage )
{
   return roundRobin( storage, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
}

MigrationInfo roundRobin( PrimitiveStorage& storage, uint_t numberOfTargetProcesses )
{
   return roundRobin( storage, 0, numberOfTargetProcesses - 1 );
}

MigrationInfo roundRobin( PrimitiveStorage& storage, uint_t minRank, uint_t maxRank )
{
   WALBERLA_CHECK_LESS_EQUAL( minRank, maxRank );
   WALBERLA_CHECK( maxRank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ),
                   "Cannot distribute to more than available processes." );

   const uint_t numReceivingProcesses = maxRank - minRank + 1;
   const uint_t rank                  = uint_c( walberla::mpi::MPIManager::instance()->rank() );

   MigrationMap_T migrationMap;
   uint_t         numReceivingPrimitives = 0;

   for ( const auto& primitiveID : storage.getPrimitiveIDs() )
   {
      const uint_t targetRank           = minRank + uint_c( primitiveID.getID() % numReceivingProcesses );
      migrationMap[primitiveID.getID()] = targetRank;
   }

   const uint_t numPrimitives = storage.getNumberOfGlobalPrimitives();

   if ( rank >= minRank && rank <= maxRank )
   {
      numReceivingPrimitives += numPrimitives / numReceivingProcesses;
      if ( numPrimitives % numReceivingProcesses > rank - minRank )
      {
         numReceivingPrimitives++;
      }
   }

   MigrationInfo migrationInfo( migrationMap, uint_c( numReceivingPrimitives ) );

   storage.migratePrimitives( migrationInfo );
   return migrationInfo;
}

MigrationInfo roundRobinInterval( PrimitiveStorage& storage, uint_t interval )
{
   return roundRobinInterval(
       storage, interval, ( uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) / interval ) * interval );
}

MigrationInfo roundRobinInterval( PrimitiveStorage& storage, uint_t interval, uint_t numProcesses )
{
   WALBERLA_CHECK_GREATER( interval, 0, "Process interval must be larger than zero." );
   WALBERLA_CHECK_LESS_EQUAL( interval * numProcesses,
                              uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ),
                              "Cannot fit " << numProcesses << " processes with interval " << interval << " into total of "
                                            << uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) << " processes." );

   const uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

   MigrationMap_T migrationMap;
   uint_t         numReceivingPrimitives = 0;

   const uint_t modulo = interval * numProcesses;

   for ( const auto& primitiveID : storage.getPrimitiveIDs() )
   {
      const uint_t targetRank           = ( uint_c( primitiveID.getID() ) * interval ) % modulo;
      migrationMap[primitiveID.getID()] = targetRank;
   }

   const uint_t numPrimitives = storage.getNumberOfGlobalPrimitives();

   if ( rank % interval == 0 && rank / interval < numProcesses )
   {
      numReceivingPrimitives += numPrimitives / numProcesses;
      if ( numPrimitives % numProcesses > rank / interval )
      {
         numReceivingPrimitives++;
      }
   }

   MigrationInfo migrationInfo( migrationMap, uint_c( numReceivingPrimitives ) );

   storage.migratePrimitives( migrationInfo );
   return migrationInfo;
}

MigrationInfo reverseDistribution( const MigrationInfo& originalMigrationInfo, PrimitiveStorage& storageToRedistribute )
{
   MigrationInfo migrationInfo = reverseDistributionDry( originalMigrationInfo );
   storageToRedistribute.migratePrimitives( migrationInfo );
   return migrationInfo;
}

MigrationInfo reverseDistributionDry( const MigrationInfo& originalMigrationInfo )
{
   MigrationMap_T migrationMap;
   const uint_t   numReceivingPrimitivesInverseMapping = originalMigrationInfo.getMap().size();

   // send messages along the original migration info to build the inverse mapping
   PackageBufferSystem bs( walberla::mpi::MPIManager::instance()->comm(), originalMigrationInfo.getNumReceivingPrimitives() );

   for ( const auto& it : originalMigrationInfo.getMap() )
   {
      bs.getPackageSendBuffer( it.second ) << it.first;
   }

   bs.sendAll();

   while ( !bs.allPackagesReceived() )
   {
      auto                package = bs.getNextPackage();
      auto&               buffer  = package.buffer();
      PrimitiveID pID;
      buffer >> pID;
      migrationMap[pID] = uint_c( package.senderRank() );
   }

   MigrationInfo migrationInfo( migrationMap, numReceivingPrimitivesInverseMapping );

   return migrationInfo;
}

void diffusiveSmooth( PrimitiveStorage& storage, uint_t outerIterations, uint_t smoothingIterations )
{
#define VISUALIZE_OUTER_ITERATIONS 1
#if VISUALIZE_OUTER_ITERATIONS
   writeDomainPartitioningVTK( storage, "/tmp/", walberla::format( "test_diffusion_%04d", 0 ) );
#endif

   const auto ownRank = walberla::mpi::MPIManager::instance()->rank();
   const auto comm    = walberla::mpi::MPIManager::instance()->comm();

   for ( uint_t outerIteration = 0; outerIteration < outerIterations; outerIteration++ )
   {
      uint_t numLocalVolumes = storage.getNumberOfLocalFaces();

      // Find all clusters.

      std::map< uint_t, std::set< PrimitiveID > > clusters;
      uint_t                                      currentClusterID   = 0;
      uint_t                                      assignedPrimitives = 0;

      while ( assignedPrimitives < numLocalVolumes )
      {
         // Choose any unassigned volume primitive.
         PrimitiveID unassignedVolume;
         for ( const auto& pid : storage.getFaceIDs() )
         {
            if ( clusters[currentClusterID].count( pid ) == 0 )
            {
               unassignedVolume = pid;
               break;
            }
         }

         std::set< PrimitiveID > nextPrimitives;
         nextPrimitives.insert( unassignedVolume );

         while ( !nextPrimitives.empty() )
         {
            auto currentPID = *nextPrimitives.begin();
            nextPrimitives.erase( currentPID );

            WALBERLA_CHECK( storage.faceExistsLocally( currentPID ) );

            WALBERLA_CHECK_EQUAL( clusters[currentClusterID].count( currentPID ), 0 );
            clusters[currentClusterID].insert( currentPID );
            assignedPrimitives++;

            const auto face = storage.getFace( currentPID );
            for ( const auto& epid : face->neighborEdges() )
            {
               const auto edge = storage.getEdge( epid );
               if ( edge->opposite_face_exists( currentPID ) )
               {
                  const auto npid = edge->get_opposite_face( currentPID );
                  if ( storage.faceExistsLocally( npid ) && clusters[currentClusterID].count( npid ) == 0 )
                  {
                     nextPrimitives.insert( npid );
                  }
               }
            }
         }

         currentClusterID++;
      }

      // Find largest cluster

      uint_t largestClusterID = 0;
      for ( const auto& it : clusters )
      {
         if ( it.second.size() > clusters[largestClusterID].size() )
         {
            largestClusterID = it.first;
         }
      }

      const auto cluster = clusters[largestClusterID];

      WALBERLA_CHECK_GREATER( cluster.size(), 0, "Cannot apply diffusive LB if a process is empty." )

      // Init load

      const double INITIAL_LOAD_TOTAL     = 1.0;
      const double INITIAL_LOAD_PRIMITIVE = INITIAL_LOAD_TOTAL / double_c( cluster.size() );

      std::map< PrimitiveID, std::map< MPIRank, double > > loadSrc;
      std::map< PrimitiveID, std::map< MPIRank, double > > loadDst;

      for ( const auto& pid : storage.getFaceIDs() )
      {
         loadSrc[pid][ownRank] = 0;
         loadDst[pid][ownRank] = 0;
      }

      for ( const auto& pid : cluster )
      {
         loadSrc[pid][ownRank] = INITIAL_LOAD_PRIMITIVE;
      }

      for ( uint_t smoothingIteration = 0; smoothingIteration < smoothingIterations; smoothingIteration++ )
      {
         // Communicate loads

         walberla::mpi::BufferSystem bs( comm );
         std::set< MPIRank >         neighborRanks;

         for ( const auto& pid : storage.getFaceIDs() )
         {
            const auto face = storage.getFace( pid );
            for ( const auto& epid : face->neighborEdges() )
            {
               if ( storage.getEdge( epid )->opposite_face_exists( pid ) )
               {
                  const auto npid  = storage.getEdge( epid )->get_opposite_face( pid );
                  const auto nrank = (MPIRank) storage.getPrimitiveRank( npid );
                  for ( const auto& it : loadSrc[pid] )
                  {
                     auto loadrank = it.first;
                     auto load     = it.second;
                     bs.sendBuffer( nrank ) << pid;
                     bs.sendBuffer( nrank ) << loadrank;
                     bs.sendBuffer( nrank ) << load;
                     neighborRanks.insert( nrank );
                  }
               }
            }
         }

         bs.setReceiverInfo( neighborRanks, true );

         bs.sendAll();

         for ( auto msg = bs.begin(); msg != bs.end(); ++msg )
         {
            while ( !msg.buffer().isEmpty() )
            {
               MPIRank     loadrank;
               PrimitiveID npid;
               double      npload;
               msg.buffer() >> npid;
               msg.buffer() >> loadrank;
               msg.buffer() >> npload;
               loadSrc[npid][loadrank] = npload;
            }
         }

         // Smooth

         for ( const auto& dstpid : storage.getFaceIDs() )
         {
            loadDst[dstpid].clear();

            for ( const auto& it : loadSrc[dstpid] )
            {
               const auto rank = it.first;
               const auto load = it.second;
               loadDst[dstpid][rank] += 0.25 * load;
            }

            const auto face = storage.getFace( dstpid );
            for ( const auto& epid : face->neighborEdges() )
            {
               auto srcpid = dstpid;
               if ( storage.getEdge( epid )->opposite_face_exists( dstpid ) )
               {
                  srcpid = storage.getEdge( epid )->get_opposite_face( dstpid );
               }

               for ( const auto& it : loadSrc[srcpid] )
               {
                  const auto rank = it.first;
                  const auto load = it.second;
                  loadDst[dstpid][rank] += 0.25 * load;
               }
            }
         }

         loadSrc = loadDst;
      }

      // Migrate all primitives to rank that has largest load.
      // If that rank is no neighbor, leave the primitive as is.

      MigrationMap_T             migrationMap;
      std::map< uint_t, uint_t > numPrimitivesToBeSent;
      uint_t                     numReceivingPrimitives = 0;

      for ( const auto& pid : storage.getFaceIDs() )
      {
         double  maxLoad         = 0;
         MPIRank rankWithMaxLoad = walberla::mpi::INVALID_RANK;

         for ( const auto& it : loadDst[pid] )
         {
            const auto rank = it.first;
            const auto load = it.second;

            if ( load > maxLoad )
            {
               maxLoad         = load;
               rankWithMaxLoad = rank;
            }
         }

         if ( storage.getNeighboringRanks().count( uint_c( rankWithMaxLoad ) ) == 1 )
         {
            migrationMap[pid.getID()] = uint_c( rankWithMaxLoad );
            numPrimitivesToBeSent[uint_c( rankWithMaxLoad )]++;
         }
      }

      for ( const auto& it : storage.getPrimitiveIDs() )
      {
         if ( migrationMap.count( it.getID() ) == 0 )
         {
            migrationMap[it.getID()] = uint_c( ownRank );
            numPrimitivesToBeSent[uint_c( ownRank )]++;
         }
      }

      walberla::mpi::BufferSystem bs( comm );
      std::set< MPIRank > nranks;
      for ( const auto & it : storage.getNeighboringRanks() )
      {
         nranks.insert( static_cast< MPIRank >( it ) );
      }
      nranks.insert( ownRank );
      bs.setReceiverInfo( nranks, true );

      for ( auto rank : nranks )
      {
         bs.sendBuffer( rank ) << numPrimitivesToBeSent[uint_c( rank )];
      }

      bs.sendAll();

      for ( auto msg = bs.begin(); msg != bs.end(); ++msg )
      {
         uint_t nprim;
         msg.buffer() >> nprim;
         numReceivingPrimitives += nprim;
      }

      MigrationInfo migrationInfo( migrationMap, numReceivingPrimitives );
      storage.migratePrimitives( migrationInfo );

#if VISUALIZE_OUTER_ITERATIONS
      writeDomainPartitioningVTK( storage, "/tmp/", walberla::format( "test_diffusion_%04d", outerIteration + 1 ) );
#endif
   }
}

} // namespace distributed
} // namespace loadbalancing
} // namespace hyteg
