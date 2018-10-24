
#pragma once

#include "core/DataTypes.h"

#include <map>
#include <set>

namespace hhg {
namespace indexing {

using walberla::uint_t;

const std::map< uint_t, std::set< uint_t > > cellLocalEdgeIDsToCellLocalNeighborFaceIDs = {
  { 0, std::set< uint_t >( { 0, 1 } ) },
  { 1, std::set< uint_t >( { 0, 2 } ) },
  { 2, std::set< uint_t >( { 0, 3 } ) },
  { 3, std::set< uint_t >( { 1, 2 } ) },
  { 4, std::set< uint_t >( { 1, 3 } ) },
  { 5, std::set< uint_t >( { 2, 3 } ) }
};

const std::map< uint_t, std::set< uint_t > > cellLocalFaceIDsToSpanningVertexIDs = {
  { 0, std::set< uint_t >( { 0, 1, 2 } ) },
  { 1, std::set< uint_t >( { 0, 1, 3 } ) },
  { 2, std::set< uint_t >( { 0, 2, 3 } ) },
  { 3, std::set< uint_t >( { 1, 2, 3 } ) }
};


}
}