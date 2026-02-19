/*
 * Copyright (c) 2026 Marcus Mohr.
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

#include "P3Function.hpp"

#include <algorithm>

namespace hyteg {

template < typename ValueType >
P3Function< ValueType >::P3Function( const std::string&                         name,
                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                     uint_t                                     minLevel,
                                     uint_t                                     maxLevel )
: P3Function( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
{}

template < typename ValueType >
P3Function< ValueType >::P3Function( const std::string&                         name,
                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                     uint_t                                     minLevel,
                                     uint_t                                     maxLevel,
                                     BoundaryCondition                          boundaryCondition )
: Function< P3Function< ValueType > >( name, storage, minLevel, maxLevel )
, vertexDoFFunction_(
      vertexdof::VertexDoFFunction< ValueType >( name + "_VertexDoF", storage, minLevel, maxLevel, boundaryCondition ) )
, edgeDoFFunctionBlue_( EdgeDoFFunction< ValueType >( name + "_EdgeDoFBlue", storage, minLevel, maxLevel, boundaryCondition ) )
, edgeDoFFunctionRed_( EdgeDoFFunction< ValueType >( name + "_EdgeDoFRed", storage, minLevel, maxLevel, boundaryCondition ) )
, faceDoFFunction_( volumedofspace::VolumeDoFFunction< ValueType >( name + "_FaceDoF", storage, minLevel, maxLevel, 1, volumedofspace::indexing::VolumeDoFMemoryLayout::SoA ) )
{
   if ( storage->hasGlobalCells() )
   {
      WALBERLA_ABORT( "This proof-of-concept implementation only works for 2D! We are missing a FaceDoFFunction for 3D!" );
   }
#if 0
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      /// one has to use the communicators of the vertexDoF and edgeDoF function to communicate
      /// TODO: find better solution
      communicators_[level] = nullptr;
   }
#endif
}

// ========================
//  explicit instantiation
// ========================
template class P3Function< double >;
template class P3Function< float >;
template class P3Function< int32_t >;
template class P3Function< int64_t >;

} //namespace hyteg
