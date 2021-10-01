#pragma once

#include "core/DataTypes.h"

#include "hyteg/types/types.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// \brief Returns the maximum velocity magnitude: max( sqrt( ux^2 + uy^2 + uz^2 ) )
template < typename FunctionType >
real_t velocityMaxMagnitude( const FunctionType& velocityX,
                             const FunctionType& velocityY,
                             const FunctionType& velocityZ,
                             const FunctionType& tmp,
                             const FunctionType& magnitudeSquared,
                             const uint_t & level,
                             const DoFType & flag )
{
   magnitudeSquared.interpolate( 0, level, flag );

   tmp.assign( {1.0}, {velocityX}, level, flag );
   tmp.multElementwise( {tmp, tmp}, level, flag );
   magnitudeSquared.assign( {1.0, 1.0}, {magnitudeSquared, tmp}, level, flag );

   tmp.assign( {1.0}, {velocityY}, level, flag );
   tmp.multElementwise( {tmp, tmp}, level, flag );
   magnitudeSquared.assign( {1.0, 1.0}, {magnitudeSquared, tmp}, level, flag );

   tmp.assign( {1.0}, {velocityZ}, level, flag );
   tmp.multElementwise( {tmp, tmp}, level, flag );
   magnitudeSquared.assign( {1.0, 1.0}, {magnitudeSquared, tmp}, level, flag );

   return std::sqrt( magnitudeSquared.getMaxMagnitude( level, flag ) );
}

template < typename FunctionType >
real_t velocityMaxMagnitude( const FunctionType& velocityX,
                             const FunctionType& velocityY,
                             const FunctionType& tmp,
                             const FunctionType& magnitudeSquared,
                             const uint_t & level,
                             const DoFType & flag )
{
   magnitudeSquared.interpolate( 0, level, flag );

   tmp.assign( {1.0}, {velocityX}, level, flag );
   tmp.multElementwise( {tmp, tmp}, level, flag );
   magnitudeSquared.assign( {1.0, 1.0}, {magnitudeSquared, tmp}, level, flag );

   tmp.assign( {1.0}, {velocityY}, level, flag );
   tmp.multElementwise( {tmp, tmp}, level, flag );
   magnitudeSquared.assign( {1.0, 1.0}, {magnitudeSquared, tmp}, level, flag );

   return std::sqrt( magnitudeSquared.getMaxMagnitude( level, flag ) );
}

} // namespace hyteg