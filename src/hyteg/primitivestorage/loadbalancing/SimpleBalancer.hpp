
#pragma once

#include "core/DataTypes.h"

namespace hyteg {
namespace loadbalancing {

using walberla::uint_t;

/// \brief Load balancing function for \ref SetupPrimitiveStorage that locates all primitives on a certain rank
void allPrimitivesOnOneRank( SetupPrimitiveStorage & storage, const uint_t & targetRank );


/// \brief Load balancing function for \ref SetupPrimitiveStorage that locates all primitives on root
void allPrimitivesOnRoot( SetupPrimitiveStorage & storage );


/// \brief Load balancing function for \ref SetupPrimitiveStorage that distributed all primitives in a round robin fashion
void roundRobin( SetupPrimitiveStorage & storage );


/// \brief Load balancing function for \ref SetupPrimitiveStorage that distributes all primitives in a greedy fashion.
/// It is expected to result in a much better edge-cut ratio than the roundRobin algorithm. But still not optimal.
void greedy( SetupPrimitiveStorage & storage );


} // namespace loadbalancing
} // namespace hyteg
