/*
 * Copyright (c) 2023 Marcus Mohr.
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

#include "hyteg/mixedoperators/Embeddings.hpp"

namespace hyteg {

template < typename ValueType >
void embedP1IntoP2( const P1Function< ValueType >& src,
                    const P2Function< ValueType >& dst,
                    const uint_t&                  level,
                    const DoFType&                 flag )
{
   dst.prolongateP1ToP2( src, level, flag );
}

template < typename ValueType >
void embedP1IntoP2( const P1VectorFunction< ValueType >& src,
                    const P2VectorFunction< ValueType >& dst,
                    const uint_t&                        level,
                    const DoFType&                       flag )
{
   for ( uint_t k = 0; k < src.getDimension(); k++ )
   {
      dst[k].prolongateP1ToP2( src[k], level, flag );
   }
}

// -------------------------
//  Explicit instantiations
// -------------------------
template void embedP1IntoP2< real_t >( const P1Function< real_t >& src,
                                       const P2Function< real_t >& dst,
                                       const uint_t&               level,
                                       const DoFType&              flag );

template void embedP1IntoP2( const P1VectorFunction< real_t >& src,
                             const P2VectorFunction< real_t >& dst,
                             const uint_t&                     level,
                             const DoFType&                    flag );

} // namespace hyteg
