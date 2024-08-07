/*
* Copyright (c) 2024 Michael Zikeli.
*
* This file is part of HyTeG
* (see https://i10git.cs.fau.de/hyteg/hyteg).
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by * the Free Software Foundation, either version 3 of the License, or
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

#include "core/DataTypes.h"
#include "core/math/Random.h"

#include "hyteg/types/PointND.hpp"

namespace hyteg::mixedPrecisionMT {

using hyteg::Point3D;

namespace setupMasks {
template < typename T >
static constexpr uint_t n_mask( const T n )
{
   return ( uint_t( 1 ) << static_cast< uint_t >( n ) );
}
template < typename... T >
static constexpr uint_t n_maskSet( const T... args )
{
   return ( 0 | ... | n_mask( args ) );
}
enum stoppingCriterion
{
   Iterations        = n_mask( 1 ),
   ResidualRate      = n_mask( 2 ),
   ResidualThreshold = n_mask( 3 ),
   ErrorRate         = n_mask( 4 ),
   ErrorThreshold    = n_mask( 5 ),
   All               = n_maskSet( 1, 2, 3, 4, 5 )
};
enum normTypes
{
   Relative         = n_mask( 6 ),
   Absolute         = n_mask( 7 ),
   Vector           = n_mask( 8 ),
   NormalizedVector = n_mask( 9 ),
   Weak             = n_mask( 10 )
};
enum errorFlags
{
   ErrIsNaN = n_mask( 11 ),
   ResIsNaN = n_mask( 12 ),
   RunErr   = n_mask( 13 )
};
} // namespace setupMasks

template < typename ValueType >
struct AccuracyTrade
{
   static ValueType initZero( const Point3D& ) { return numeric_cast< ValueType >( 0.0 ); }

   static ValueType initOne( const Point3D& ) { return numeric_cast< ValueType >( 1.0 ); }

   static ValueType initRandom( const Point3D& )
   {
      return numeric_cast< ValueType >( walberla::math::realRandom< real_t >( 0.0, 1.0 ) );
   }

   static std::string name()
   {
      WALBERLA_ASSERT( false,
                       "" << typeid( ValueType ).name() << " is no known ValueType,"
                          << " therefore no user readable name can be specified." );
      return "unknown";
   }
};

template <>
std::string AccuracyTrade< walberla::float64 >::name()
{
   return "fp64";
}
template <>
std::string AccuracyTrade< walberla::float32 >::name()
{
   return "fp32";
}
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template <>
std::string AccuracyTrade< walberla::float16 >::name()
{
   return "fp16";
}
#endif // WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

} // namespace hyteg::mixedPrecisionMT
