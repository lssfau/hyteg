/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <memory>

#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {
namespace MeshQuality {

inline real_t
    getMinimalEdgeLength( const std::shared_ptr< hyteg::PrimitiveStorage >& storage, uint_t level, bool onRootOnly = false )
{
   real_t localMin = std::numeric_limits< real_t >::max();

   for ( auto& it : storage->getEdges() )
   {
      Edge& edge = *it.second;
      localMin   = std::min( localMin, edge.getLength() );
   }

   real_t globalMin;
   if ( onRootOnly )
   {
      globalMin = walberla::mpi::reduce( localMin, walberla::mpi::MIN );
   }
   else
   {
      globalMin = walberla::mpi::allReduce( localMin, walberla::mpi::MIN );
   }

   return std::pow( real_c( 2.0 ), -walberla::real_c( level ) ) * globalMin;
}

inline real_t
    getMaximalEdgeLength( const std::shared_ptr< hyteg::PrimitiveStorage >& storage, uint_t level, bool onRootOnly = false )
{
   real_t localMax = 0;

   for ( auto& it : storage->getEdges() )
   {
      Edge& edge = *it.second;
      localMax   = std::max( localMax, edge.getLength() );
   }

   real_t globalMax;
   if ( onRootOnly )
   {
      globalMax = walberla::mpi::reduce( localMax, walberla::mpi::MAX );
   }
   else
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return std::pow( real_c( 2.0 ), -walberla::real_c( level ) ) * globalMax;
}

} // namespace MeshQuality
} // namespace hyteg