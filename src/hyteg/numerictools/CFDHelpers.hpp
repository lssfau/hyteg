#pragma once

#include "core/DataTypes.h"

#include "hyteg/types/flags.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// \brief Returns the maximum velocity magnitude: max( sqrt( ux^2 + uy^2 + uz^2 ) )
template < typename FunctionType >
real_t velocityMaxMagnitude( const FunctionType& velocityX,
                             const FunctionType& velocityY,
                             const FunctionType& velocityZ,
                             const FunctionType& tmp1,
                             const FunctionType& tmp2,
                             const uint_t & level,
                             const DoFType & flag )
{
   // Could be done with interpolateExtended() and without tmp variables, however, this version should be much faster.

   tmp2.interpolate( 0, level, flag );

   tmp1.assign( {1.0}, {velocityX}, level, flag );
   tmp1.multElementwise( {tmp1, tmp1}, level, flag );
   tmp2.assign( {1.0, 1.0}, {tmp2, tmp1}, level, flag );

   tmp1.assign( {1.0}, {velocityY}, level, flag );
   tmp1.multElementwise( {tmp1, tmp1}, level, flag );
   tmp2.assign( {1.0, 1.0}, {tmp2, tmp1}, level, flag );

   tmp1.assign( {1.0}, {velocityZ}, level, flag );
   tmp1.multElementwise( {tmp1, tmp1}, level, flag );
   tmp2.assign( {1.0, 1.0}, {tmp2, tmp1}, level, flag );

   return std::sqrt( tmp2.getMaxMagnitude( level, flag ) );
}

} // namespace hyteg