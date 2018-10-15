#pragma once

#include <memory>

#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

namespace hhg {
namespace MeshQuality {

real_t getMinimalEdgeLength( const std::shared_ptr< hhg::PrimitiveStorage >& storage, uint_t level )
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
} // namespace hhg