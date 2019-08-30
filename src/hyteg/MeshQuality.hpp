#pragma once

#include <memory>

#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {
namespace MeshQuality {

real_t getMinimalEdgeLength( const std::shared_ptr< hyteg::PrimitiveStorage >& storage, uint_t level )
{
   real_t localMin = std::numeric_limits< real_t >::max();

   for( auto& it : storage->getEdges() )
   {
      Edge& edge = *it.second;
      localMin   = std::min( localMin, edge.getLength() );
   }

   real_t globalMin = walberla::mpi::allReduce( localMin, walberla::mpi::MIN );

   return std::pow( 2.0, -walberla::real_c( level ) ) * globalMin;
}

} // namespace MeshQuality
} // namespace hyteg