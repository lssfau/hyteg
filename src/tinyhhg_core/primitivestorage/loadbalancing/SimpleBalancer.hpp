
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
    SetupPrimitiveStorage::SetupPrimitiveMap setupPrimitives;
    storage.getSetupPrimitives( setupPrimitives );
    for ( auto it : setupPrimitives )
    {
      it.second->setTargetRank( uint_c( targetRank_ ) );
    }
    return uint_c( 1 );
  }

private:

  uint_t targetRank_;

};


/// \brief Loadbalanding callback for \ref SetupPrimitiveStorage that locates all blocks on root
/// \author Nils Kohl (nils.kohl@fau.de)
class AllBlocksOnRoot : public AllBlocksOnOneRank
{
public:
  AllBlocksOnRoot() : AllBlocksOnOneRank( uint_c( 0 ) ) {}
};


/// \brief Loadbalanding callback for \ref SetupPrimitiveStorage that distributed all primitives in a round robin fashion
/// \autho Nils Kohl (nils.kohl@fau.de)
class RoundRobin
{
public:
  uint_t operator()( SetupPrimitiveStorage & storage, const memory_t & perProcessMemoryLimit ) const
  {
    uint_t currentRank = 0;

    SetupPrimitiveStorage::SetupPrimitiveMap setupPrimitives;
    storage.getSetupPrimitives( setupPrimitives );

    for ( auto it : setupPrimitives )
    {
      it.second->setTargetRank( uint_c( currentRank % storage.getNumberOfProcesses() ) );
      currentRank++;
    }

    return std::min( currentRank + 1, storage.getNumberOfProcesses() );
  }
};


} // namespace hhg
