/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/primitivedata/PrimitiveDataID.hpp"
#include "hyteg/primitivedata/PrimitiveDataHandling.hpp"

#include "core/DataTypes.h"
#include <core/DataTypes.h>
#include <core/logging/Logging.h>
#include <core/mpi/SendBuffer.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/BufferSizeTrait.h>
#include <core/mpi/BufferDataTypeExtensions.h>

#include <map>

namespace hyteg {

using walberla::uint_t;

template< typename DataType >
class LevelWiseMemory
{
public:

    LevelWiseMemory( const uint_t & minLevel, const uint_t & maxLevel )
    {
      for ( uint_t i = minLevel; i <= maxLevel; i++ )
      {
        data_[i] = DataType();
      }
    }

    bool hasLevel( const uint_t & level ) const
    {
      return data_.count( level ) > 0;
    }

    const DataType & getData( const uint_t & level ) const
    {
      WALBERLA_ASSERT( hasLevel( level ) );
      return data_.at( level );
    }

    DataType & getData( const uint_t & level )
    {
      WALBERLA_ASSERT( hasLevel( level ) );
      return data_[level];
    }

    inline void serialize( walberla::mpi::SendBuffer & sendBuffer ) const
    {
      sendBuffer << data_;
    }

    inline void deserialize( walberla::mpi::RecvBuffer & recvBuffer )
    {
      recvBuffer >> data_;
    }

private:

    std::map < uint_t, DataType > data_;
};


template< typename DataType, typename PrimitiveType >
class LevelWiseMemoryDataHandling : public PrimitiveDataHandling< DataType, PrimitiveType >
{
public:

    LevelWiseMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel )
      : minLevel_( minLevel ),
        maxLevel_( maxLevel )
    {}

    virtual ~LevelWiseMemoryDataHandling() {}

    std::shared_ptr< DataType > initialize( const PrimitiveType * const primitive ) const
    {
      return std::make_shared< DataType >( minLevel_, maxLevel_ );
    }

    virtual void serialize( const PrimitiveType * const primitive, const PrimitiveDataID< DataType, PrimitiveType > & id, walberla::mpi::SendBuffer & buffer ) const
    {
      DataType * data = primitive->getData( id );
      buffer << *data;
    }

    virtual void deserialize( const PrimitiveType * const primitive, const PrimitiveDataID< DataType, PrimitiveType > & id, walberla::mpi::RecvBuffer & buffer ) const
    {
      DataType * data = primitive->getData( id );
      buffer >> *data;
    }


private:

    const uint_t minLevel_;
    const uint_t maxLevel_;

};

}

namespace walberla {
namespace mpi {

template< typename T,  // Element type of SendBuffer
typename G,  // Growth policy of SendBuffer
typename ValueType>
inline mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const hyteg::LevelWiseMemory< ValueType > & levelWiseMemory )
{
  levelWiseMemory.serialize( buffer );
  return buffer;
}

template< typename T, // Element type  of RecvBuffer
typename ValueType >
inline mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buffer,
                                                hyteg::LevelWiseMemory< ValueType > & levelWiseMemory )
{
  levelWiseMemory.deserialize( buffer );
  return buffer;
}

template< typename ValueType >
struct BufferSizeTrait< hyteg::LevelWiseMemory< ValueType > > { static const bool constantSize = false; };

} // namespace mpi
} // namespace walberla

