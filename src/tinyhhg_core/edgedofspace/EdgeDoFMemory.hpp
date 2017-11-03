#pragma once

#include "core/DataTypes.h"
#include "tinyhhg_core/levelinfo.hpp"

namespace hhg {

using walberla::uint_t;

inline uint_t EdgeDoFMacroVertexFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  WALBERLA_UNUSED( level );
  return 2 * numDependencies;
}

inline uint_t EdgeDoFMacroEdgeFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  return levelinfo::num_microedges_per_edge( level ) + numDependencies * ( 3 * ( levelinfo::num_microedges_per_edge( level ) ) - 1 );
}

inline uint_t EdgeDoFMacroFaceFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  WALBERLA_UNUSED( numDependencies );
  return 3 * ( ( ( levelinfo::num_microedges_per_edge( level ) + 1 ) * levelinfo::num_microedges_per_edge( level ) ) / 2 );
}

}
