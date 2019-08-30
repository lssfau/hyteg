
#include "tinyhhg_core/primitivestorage/loadbalancing/DistributedBalancer.hpp"

#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"
#include "core/load_balancing/ParMetisWrapper.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/MPIWrapper.h"

#include <algorithm>

namespace hyteg {
namespace loadbalancing {
namespace distributed {

using walberla::int64_t;
using walberla::int64_c;
using walberla::mpi::MPIRank;
using namespace walberla::mpistubs;

void parmetis( PrimitiveStorage & storage )
{
  WALBERLA_CHECK_GREATER( storage.getNumberOfLocalPrimitives(), 0, "ParMeTis not supported (yet) for distributions with empty processes." );

  MPI_Comm communicator = walberla::mpi::MPIManager::instance()->comm();
  uint_t   rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );
  uint_t   numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  // Parameters needed by parmetis

  std::vector< int64_t > vtxdist; // range of vertices that are local to each processor (identical on all ranks)
  std::vector< int64_t > xadj;    // index and vtxdist indicate the vertex ID, entries the edge indices in adjncy
  std::vector< int64_t > adjncy;  // IDs of the neighboring vertices
  std::vector< int64_t > vwgt;    // vertex weights
  std::vector< int64_t > adjwgt;  // edge weights
  int64_t                wgtflag; // indicates types of weights
  int64_t                numflag; // numbering scheme
  std::vector< int64_t > ndims;   // space dimensions
  std::vector< double  > xyz;     // vertex coordinates
  int64_t                ncon;    // number of weights per vertex
  int64_t                nparts;  // desired number of subdomains
  std::vector< double  > tpwgts;  // specifies vertex weight distribution
  std::vector< double  > ubvec;   // imbalance tolerance
  std::vector< int64_t > options; // parmetis options
  std::vector< int64_t > edgecut; // [out] edge cut
  std::vector< int64_t > part;    // [out] resulting partition
  MPI_Comm * parmetisCommunicator; // MPI communicator used by parmetis

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
  for ( const auto & numberOfPrimitivesOnProcess : numberOfLocalPrimitivesOnProcesses )
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

  std::map< PrimitiveID::IDType, int64_t     > localPrimitiveIDToGlobalParmetisIDMap; // contains all local PrimitiveIDs as keys and maps them to global parmetis IDs
  std::map< int64_t            , PrimitiveID > globalParmetisIDToLocalPrimitiveIDMap; // reverse of the above map

  std::vector< PrimitiveID >                   localPrimitiveIDs;
  storage.getPrimitiveIDs( localPrimitiveIDs );

  int64_t parmetisIDCounter = vtxdist[ rank ];
  for ( const auto & id : localPrimitiveIDs )
  {
    localPrimitiveIDToGlobalParmetisIDMap[ id.getID() ] = parmetisIDCounter;
    parmetisIDCounter++;
  }

  // Reverse the mapping (for convenience)
  for ( const auto it : localPrimitiveIDToGlobalParmetisIDMap )
  {
    WALBERLA_ASSERT_EQUAL( globalParmetisIDToLocalPrimitiveIDMap.count( it.second ), 0 );
    globalParmetisIDToLocalPrimitiveIDMap[ it.second ] = PrimitiveID( it.first );
  }

  // To build the parmetis graph, we now need the mappings (PrimitiveID to parmetisID) from all neighboring processes

  std::set< MPIRank > neighboringRanks;
  storage.getNeighboringRanks( neighboringRanks );

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
    recv.buffer() >> neighboringPrimitiveIDToGlobalParmetisIDMaps[ uint_c( recv.rank() ) ];
  }

#ifndef NDEBUG
  for ( const MPIRank neighborRank : neighboringRanks )
  {
    WALBERLA_ASSERT_EQUAL( neighboringPrimitiveIDToGlobalParmetisIDMaps[ uint_c( neighborRank ) ].size(), numberOfLocalPrimitivesOnProcesses[ uint_c( neighborRank ) ] );
  }
#endif

  //////////////////////////////
  // Building xadj and adjncy //
  //////////////////////////////

  numflag = 0;

  // Now that we got the assignment from PrimitiveIDs to parmetis IDs of the local and all neighboring processes, we can build the parmetis graph

  WALBERLA_ASSERT( std::is_sorted( globalParmetisIDToLocalPrimitiveIDMap.begin(), globalParmetisIDToLocalPrimitiveIDMap.end() ) );

  for ( const auto & it : globalParmetisIDToLocalPrimitiveIDMap )
  {
    const PrimitiveID primitiveID = it.second;
    const Primitive * primitive   = storage.getPrimitive( primitiveID );

    std::vector< PrimitiveID > neighborIDs;
    primitive->getNeighborPrimitives( neighborIDs );

    xadj.push_back( int64_c( adjncy.size() ) );

    for ( const auto & neighborID : neighborIDs )
    {
      uint_t  neighborRank;
      int64_t neighborParmetisID;;
      if ( storage.primitiveExistsInNeighborhood( neighborID ) )
      {
        neighborRank       = storage.getNeighborPrimitiveRank( neighborID );
        neighborParmetisID = neighboringPrimitiveIDToGlobalParmetisIDMaps[ neighborRank ][ neighborID.getID() ];
      }
      else
      {
        WALBERLA_ASSERT( storage.primitiveExistsLocally( neighborID ) );
        neighborRank = rank;
        neighborParmetisID = localPrimitiveIDToGlobalParmetisIDMap[ neighborID.getID() ];
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

  int parmetisError =
  walberla::core::ParMETIS_V3_PartKway( vtxdist.data(), xadj.data(), adjncy.data(), vwgt.data(), /* adjwgt */ NULL, &wgtflag,
                                        &numflag, &ncon, &nparts, tpwgts.data(), ubvec.data(), options.data(),
                                        edgecut.data(), part.data(), parmetisCommunicator);

  WALBERLA_CHECK_EQUAL( parmetisError, walberla::core::METIS_OK );

  /////////////////////////
  // Primitive migration //
  /////////////////////////

  std::map< PrimitiveID::IDType, uint_t > migrationMap;
  for ( uint_t partIdx = 0; partIdx < part.size(); partIdx++ )
  {
    const int64_t     parmetisID  = vtxdist[ rank ] + int64_c( partIdx );
    const PrimitiveID primitiveID = globalParmetisIDToLocalPrimitiveIDMap[ parmetisID ];

    migrationMap[ primitiveID.getID() ] = uint_c( part[ partIdx ] );
  }

  storage.migratePrimitives( migrationMap );



}

}
}
}
