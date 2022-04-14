#pragma once

#include "core/DataTypes.h"

#include "hyteg/types/types.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// \brief Returns the maximum velocity magnitude: max( sqrt( ux^2 + uy^2 + uz^2 ) )
template < typename vfunc_t >
real_t velocityMaxMagnitude( const vfunc_t&                                                velocity,
                             const typename FunctionTrait< vfunc_t >::VectorComponentType& tmp,
                             const typename FunctionTrait< vfunc_t >::VectorComponentType& magnitudeSquared,
                             const uint_t&                                                 level,
                             const DoFType&                                                flag )
{
   magnitudeSquared.interpolate( real_c( 0 ), level, flag );

   for ( uint_t k = 0; k < velocity.getDimension(); ++k )
   {
      tmp.assign( {1.0}, {velocity[k]}, level, flag );
      tmp.multElementwise( {tmp, tmp}, level, flag );
      magnitudeSquared.assign( {1.0, 1.0}, {magnitudeSquared, tmp}, level, flag );
   }

   return std::sqrt( magnitudeSquared.getMaxMagnitude( level, flag ) );
}

} // namespace hyteg
