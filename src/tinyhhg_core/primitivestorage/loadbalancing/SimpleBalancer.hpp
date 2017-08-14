
#pragma once

#include "core/DataTypes.h"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hhg {

/// \brief Loadbalanding callback for \ref SetupPrimitiveStorage that locates all blocks on a certain rank
/// \author Nils Kohl (nils.kohl@fau.de)
class AllBlocksOnOneRank
{
public:

  AllBlocksOnOneRank( uint_t targetRank ) : targetRank_( targetRank ) {}

  uint_t operator()( SetupPrimitiveStorage & storage, const memory_t & perProcessMemoryLimit ) const
  {
    SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
    storage.getSetupPrimitives( setupPrimitives );
    for ( auto it : setupPrimitives )
    {
      storage.setTargetRank( it.first, uint_c( targetRank_ ) );
    }
    return uint_c( 1 );
  }

private:

  uint_t targetRank_;

};


/// \brief Load balanding callback for \ref SetupPrimitiveStorage that locates all blocks on root
/// \author Nils Kohl (nils.kohl@fau.de)
class AllBlocksOnRoot : public AllBlocksOnOneRank
{
public:
  AllBlocksOnRoot() : AllBlocksOnOneRank( uint_c( 0 ) ) {}
};


/// \brief Load balancing callback for \ref SetupPrimitiveStorage that distributed all primitives in a round robin fashion
/// \author Nils Kohl (nils.kohl@fau.de)
class RoundRobin
{
public:
  uint_t operator()( SetupPrimitiveStorage & storage, const memory_t & perProcessMemoryLimit ) const
  {
    uint_t currentRank = 0;

    SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
    storage.getSetupPrimitives( setupPrimitives );

    for ( auto it : setupPrimitives )
    {
      storage.setTargetRank( it.first, uint_c( currentRank % storage.getNumberOfProcesses() ) );
      currentRank++;
    }

    return std::min( currentRank + 1, storage.getNumberOfProcesses() );
  }
};


/// \brief Load balancing callback for \ref SetupPrimitiveStorage that distributes all primitives in a greedy fashion.
///
/// This balancer ignores the type of the primitive (vertex, edge, ...) and therefore does not result in a distribution
/// that balances the number of the primitives per process if they are distinguished by type.
/// It's also pretty slow.
class GreedyIgnoringPrimitiveType
{
public:

  uint_t operator()( SetupPrimitiveStorage & storage, const memory_t & perProcessMemoryLimit ) const
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

    return 0;
  }

};


} // namespace hhg
