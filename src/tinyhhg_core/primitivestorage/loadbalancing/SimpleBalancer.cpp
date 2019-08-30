
#include "core/DataTypes.h"

#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include <queue>

namespace hyteg {
namespace loadbalancing {

void allPrimitivesOnOneRank( SetupPrimitiveStorage & storage, const uint_t & targetRank )
{
  SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
  storage.getSetupPrimitives( setupPrimitives );
  for ( auto it : setupPrimitives )
  {
    storage.setTargetRank( it.first, uint_c( targetRank ) );
  }
}


void allPrimitivesOnRoot( SetupPrimitiveStorage & storage )
{
  allPrimitivesOnOneRank( storage, uint_c( 0 ) );
}


void roundRobin( SetupPrimitiveStorage & storage )
{
  uint_t currentRank = 0;

  SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
  storage.getSetupPrimitives( setupPrimitives );

  for ( auto it : setupPrimitives )
  {
    storage.setTargetRank( it.first, uint_c( currentRank % storage.getNumberOfProcesses() ) );
    currentRank++;
  }
}


void greedy( SetupPrimitiveStorage & storage )
{
  const uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  const uint_t maxVertices   = storage.getNumberOfVertices  () / numProcesses;
  const uint_t maxEdges      = storage.getNumberOfEdges     () / numProcesses;
  const uint_t maxFaces      = storage.getNumberOfFaces     () / numProcesses;
  const uint_t maxCells      = storage.getNumberOfCells     () / numProcesses;

  // Set all target ranks to zero
  SetupPrimitiveStorage::PrimitiveMap primitives;
  storage.getSetupPrimitives( primitives );

  for ( const auto & it : primitives )
  {
    storage.setTargetRank( it.first, 0 );
  }


  // Main algorithm
  for ( uint_t rank = 1; rank < numProcesses; rank++ )
  {
    uint_t currentNumVertices   = 0;
    uint_t currentNumEdges      = 0;
    uint_t currentNumFaces      = 0;
    uint_t currentNumCells      = 0;

    std::queue< PrimitiveID > nextPrimitives;
    std::map< PrimitiveID::IDType, bool > wasAddedToQueue;

    while ( true )
    {
      if (    currentNumVertices >= maxVertices
           && currentNumEdges    >= maxEdges
           && currentNumFaces    >= maxFaces
           && currentNumCells    >= maxCells )
      {
        break;
      }

      // Push a random, unassigned primitive to the queue if the queue is empty and we are not finished
      if ( nextPrimitives.size() == 0 )
      {
        for ( const auto & it : primitives )
        {
          if ( storage.getTargetRank( it.first ) == 0 && !wasAddedToQueue[ it.first ] )
          {
            if (   ( storage.vertexExists( it.first ) && currentNumVertices < maxVertices )
                || ( storage.edgeExists  ( it.first ) && currentNumEdges < maxEdges )
                || ( storage.faceExists  ( it.first ) && currentNumFaces < maxFaces )
                || ( storage.cellExists  ( it.first ) && currentNumCells < maxCells )
                )
            {
              nextPrimitives.push( it.first );
              wasAddedToQueue[ it.first ] = true;
              break;
            }
          }
        }
      }

      WALBERLA_ASSERT_GREATER( nextPrimitives.size(), 0 );

      // Pop a primitive from the queue
      const Primitive * currentPrimitive = storage.getPrimitive( nextPrimitives.front() );
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

      for ( const auto & neighborID : neighbors )
      {
        if ( storage.getTargetRank( neighborID ) == 0 && !wasAddedToQueue[ neighborID.getID() ] )
        {
          nextPrimitives.push( neighborID );
          wasAddedToQueue[ neighborID.getID() ] = true;
        }
      }
    }
  }
}


} // namespace loadbalancing
} // namespace hyteg
