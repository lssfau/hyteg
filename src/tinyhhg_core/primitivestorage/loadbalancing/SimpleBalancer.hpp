
#pragma once

#include "core/DataTypes.h"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hhg {
namespace loadbalancing {

/// \brief Load balancing function for \ref SetupPrimitiveStorage that locates all primitives on a certain rank
/// \author Nils Kohl (nils.kohl@fau.de)
void allPrimitivesOnOneRank( SetupPrimitiveStorage & storage, const uint_t & targetRank )
{
  SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
  storage.getSetupPrimitives( setupPrimitives );
  for ( auto it : setupPrimitives )
  {
    storage.setTargetRank( it.first, uint_c( targetRank ) );
  }
}


/// \brief Load balancing function for \ref SetupPrimitiveStorage that locates all primitives on root
/// \author Nils Kohl (nils.kohl@fau.de)
void allPrimitivesOnRoot( SetupPrimitiveStorage & storage )
{
  allPrimitivesOnOneRank( storage, uint_c( 0 ) );
}


/// \brief Load balancing function for \ref SetupPrimitiveStorage that distributed all primitives in a round robin fashion
/// \author Nils Kohl (nils.kohl@fau.de)
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


/// \brief Load balancing function for \ref SetupPrimitiveStorage that distributes all primitives in a greedy fashion.
///
/// This balancer ignores the type of the primitive (vertex, edge, ...) and therefore does not result in a distribution
/// that balances the number of the primitives per process if they are distinguished by type.
/// It's also pretty slow.
void greedyIgnoringPrimitiveType( SetupPrimitiveStorage & storage )
{
  const uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
  const uint_t maxPrimitives = storage.getNumberOfPrimitives() / numProcesses;

  // Set all target ranks to zero
  SetupPrimitiveStorage::PrimitiveMap primitives;
  storage.getSetupPrimitives( primitives );

  PrimitiveID       startPrimitiveID;
  const Primitive * startPrimitive;

  for ( const auto & it : primitives )
  {
    startPrimitiveID = it.first;
    startPrimitive   = it.second.get();
    storage.setTargetRank( it.first, 0 );
  }

  // Main algorithm
  const Primitive * currentPrimitive = startPrimitive;
  std::vector< PrimitiveID > visited;

  for ( uint_t rank = 1; rank < numProcesses; rank++ )
  {
    uint_t currentNumPrimitives = 0;

    while ( true )
    {
      visited.push_back( currentPrimitive->getID() );

      if ( storage.getTargetRank( currentPrimitive->getID().getID() ) == 0 )
      {
        storage.setTargetRank( currentPrimitive->getID().getID(), rank );
        currentNumPrimitives++;

        if ( currentNumPrimitives >= maxPrimitives )
        {
          break;
        }
      }

      // Choose next primitive
      const Primitive * nextPrimitive;
      for ( const auto & it : primitives )
      {
        PrimitiveID id = it.first;
        if ( std::find( visited.begin(), visited.end(), id ) == visited.end() )
        {
          nextPrimitive = storage.getPrimitive( id );
          break;
        }
      }

      std::vector< PrimitiveID > neighbors;
      currentPrimitive->getNeighborPrimitives( neighbors );

      for ( const auto & neighborID : neighbors )
      {
        if ( std::find( visited.begin(), visited.end(), neighborID ) == visited.end() )
        {
          nextPrimitive = storage.getPrimitive( neighborID );
          break;
        }
      }

      currentPrimitive = nextPrimitive;
    }
  }
}


} // namespace loadbalancing
} // namespace hhg
