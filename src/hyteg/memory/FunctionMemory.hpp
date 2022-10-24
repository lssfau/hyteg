/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include <core/DataTypes.h>
#include <core/debug/CheckFunctions.h>
#include <core/logging/Logging.h>
#include <core/mpi/BufferSizeTrait.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/Reduce.h>
#include <core/mpi/SendBuffer.h>
#include <map>
#include <memory>

#include "hyteg/misc/zeros.hpp"
#include "hyteg/primitivedata/PrimitiveDataHandling.hpp"
#include "hyteg/primitives/Primitive.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using walberla::mpi::RecvBuffer;
using walberla::mpi::SendBuffer;

template < typename ValueType >
class FunctionMemory
{
   static_assert( std::is_arithmetic< ValueType >::value, "Wrong ValueType template" );

 public:
   /// Constructs memory for a function
   explicit FunctionMemory( const ValueType fillValue = ValueType() )
   : fillValue_( fillValue )
   {}

   FunctionMemory( const std::function< uint_t( uint_t level, const Primitive& primitive ) >& sizeFunction,
                   const Primitive&                                                           primitive,
                   const uint_t&                                                              minLevel,
                   const uint_t&                                                              maxLevel,
                   const ValueType                                                            fillValue = ValueType() )
   : fillValue_( fillValue )
   {
      WALBERLA_ASSERT_LESS_EQUAL(
          minLevel, maxLevel, "minLevel should be equal or less than maxLevel during FunctionMemory allocation." );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         addData( level, sizeFunction( level, primitive ), fillValue );
      }
   }

   ~FunctionMemory()
   {
      for ( const auto& it : data_ )
      {
         totalAllocatedMemoryInBytes_ -= it.second->size() * sizeof( ValueType );
      }
   }

   /// Returns true if data is allocated at the specified level, false otherwise.
   inline bool hasLevel( const uint_t& level ) const { return data_.count( level ) > 0; }

   inline uint_t getSize( const uint_t& level ) const
   {
      WALBERLA_CHECK( hasLevel( level ), "Requested level not allocated" );
      return data_.at( level )->size();
   }

   /// Allocates an array of size size for a certain level
   inline void addData( const uint_t& level, const uint_t& size, const ValueType& fillValue )
   {
      WALBERLA_ASSERT( !hasLevel( level ),
                       "Attempting to overwrite already existing level (level == " << level << ") in function memory!" );
      data_[level] = std::unique_ptr< std::vector< ValueType > >( new std::vector< ValueType >( size, fillValue ) );
      totalAllocatedMemoryInBytes_ += size * sizeof( ValueType );
   }

   /// Deletes data of a certain level
   inline void deleteData( const uint_t& level )
   {
      if ( !hasLevel( level ) )
         return;
      totalAllocatedMemoryInBytes_ -= data_[level]->size() * sizeof( ValueType );
      data_.erase( level );
   }

   /// Returns a pointer to the first entry of the allocated array
   inline ValueType* getPointer( const uint_t& level ) const
   {
      WALBERLA_CHECK( hasLevel( level ), "Requested level " << level << " not allocated" );
      return data_.at( level )->data();
   }

   /// Copies the data of one leve from the other FunctionMemory.
   inline void copyFrom( const FunctionMemory& other, const uint_t& level ) { *data_[level] = *other.data_.at( level ); }

   inline void swap( const FunctionMemory< ValueType >& other, const uint_t& level ) const
   {
      WALBERLA_ASSERT( hasLevel( level ), "Requested level not allocated." );
      WALBERLA_ASSERT( other.hasLevel( level ), "Requested level not allocated." );
      WALBERLA_ASSERT_EQUAL( getSize( level ), other.getSize( level ), "Cannot swap FunctionMemory of different sizes." );
      data_.at( level )->swap( *( other.data_.at( level ) ) );
   }

   inline void setToZero( const uint_t& level ) const
   {
      WALBERLA_ASSERT( hasLevel( level ), "Requested level not allocated." );
      ValueType* ptr = data_.at( level )->data();
      for ( uint_t k = 0; k < data_.at( level )->size(); ++k )
      {
         ptr[k] = generateZero< ValueType >();
      }
   }

   inline static unsigned long long getLocalAllocatedMemoryInBytes() { return totalAllocatedMemoryInBytes_; }
   inline static unsigned long long getMinLocalAllocatedMemoryInBytes()
   {
      return walberla::mpi::allReduce(
          totalAllocatedMemoryInBytes_, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() );
   }
   inline static unsigned long long getMaxLocalAllocatedMemoryInBytes()
   {
      return walberla::mpi::allReduce(
          totalAllocatedMemoryInBytes_, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );
   }
   inline static unsigned long long getGlobalAllocatedMemoryInBytes()
   {
      return walberla::mpi::allReduce(
          totalAllocatedMemoryInBytes_, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   }

   /// Serializes the allocated data to a send buffer
   inline void serialize( SendBuffer& sendBuffer ) const
   {
      const uint_t numLevels = data_.size();
      sendBuffer << numLevels;

      for ( const auto& it : data_ )
      {
         const uint_t level     = it.first;
         const uint_t levelSize = it.second->size();

         sendBuffer << level;
         sendBuffer << levelSize;
         sendBuffer << getVector( level );
      }
   }

   /// Deserializes data from a recv buffer (clears all already allocated data and replaces it with the recv buffer's content)
   inline void deserialize( RecvBuffer& recvBuffer )
   {
      data_.clear();

      uint_t numLevels;

      recvBuffer >> numLevels;

      for ( uint_t levelCnt = 0; levelCnt < numLevels; levelCnt++ )
      {
         uint_t level;
         uint_t levelSize;

         recvBuffer >> level;
         recvBuffer >> levelSize;

         addData( level, levelSize, fillValue_ );

         std::vector< ValueType >& dataVector = getVector( level );
         recvBuffer >> dataVector;
      }
   }
   inline const std::vector< ValueType >& getVector( const uint_t& level ) const { return *( data_.at( level ) ); }
   inline std::vector< ValueType >&       getVector( const uint_t& level ) { return *( data_[level] ); }

 private:
   /// Maps a level to the respective allocated data
   std::map< uint_t, std::unique_ptr< std::vector< ValueType > > > data_;

   const ValueType fillValue_;

   static unsigned long long totalAllocatedMemoryInBytes_;
};

template < typename ValueType >
unsigned long long FunctionMemory< ValueType >::totalAllocatedMemoryInBytes_ = 0;

} // namespace hyteg

namespace walberla {
namespace mpi {

template < typename T, // Element type of SendBuffer
           typename G, // Growth policy of SendBuffer
           typename ValueType >
inline mpi::GenericSendBuffer< T, G >& operator<<( mpi::GenericSendBuffer< T, G >&           buffer,
                                                   const hyteg::FunctionMemory< ValueType >& functionMemory )
{
   functionMemory.serialize( buffer );
   return buffer;
}

template < typename T, // Element type  of RecvBuffer
           typename ValueType >
inline mpi::GenericRecvBuffer< T >& operator>>( mpi::GenericRecvBuffer< T >&        buffer,
                                                hyteg::FunctionMemory< ValueType >& functionMemory )
{
   functionMemory.deserialize( buffer );
   return buffer;
}

template < typename ValueType >
struct BufferSizeTrait< hyteg::FunctionMemory< ValueType > >
{
   static const bool constantSize = false;
};

} // namespace mpi
} // namespace walberla

namespace hyteg {

template < typename DataType, typename PrimitiveType >
class FunctionMemoryDataHandling : public PrimitiveDataHandling< DataType, PrimitiveType >
{
 public:
   virtual ~FunctionMemoryDataHandling() {}

   virtual void serialize( const PrimitiveType* const                        primitive,
                           const PrimitiveDataID< DataType, PrimitiveType >& id,
                           SendBuffer&                                       buffer ) const
   {
      DataType* data = primitive->getData( id );
      buffer << *data;
   }

   virtual void deserialize( const PrimitiveType* const                        primitive,
                             const PrimitiveDataID< DataType, PrimitiveType >& id,
                             RecvBuffer&                                       buffer ) const
   {
      DataType* data = primitive->getData( id );
      buffer >> *data;
   }
};

template < typename DataType, typename PrimitiveType >
class MemoryDataHandling : public PrimitiveDataHandling< DataType, PrimitiveType >
{
 public:
   explicit MemoryDataHandling()
   : initializeDefault_( true )
   {}

   MemoryDataHandling( const uint_t&                                                              minLevel,
                       const uint_t&                                                              maxLevel,
                       const std::function< uint_t( uint_t level, const Primitive& primitive ) >& sizeFunction )
   : initializeDefault_( false )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , sizeFunction_( sizeFunction )
   {}

   std::shared_ptr< DataType > initialize( const PrimitiveType* const primitive ) const
   {
      if ( initializeDefault_ )
         return std::make_shared< DataType >();
      else
         return std::make_shared< DataType >( sizeFunction_, *primitive, minLevel_, maxLevel_, 0 );
   }

   virtual void serialize( const PrimitiveType* const                        primitive,
                           const PrimitiveDataID< DataType, PrimitiveType >& id,
                           SendBuffer&                                       buffer ) const
   {
      DataType* data = primitive->getData( id );
      buffer << *data;
   }

   virtual void deserialize( const PrimitiveType* const                        primitive,
                             const PrimitiveDataID< DataType, PrimitiveType >& id,
                             RecvBuffer&                                       buffer ) const
   {
      DataType* data = primitive->getData( id );
      buffer >> *data;
   }

 private:
   bool                                                                initializeDefault_;
   uint_t                                                              minLevel_;
   uint_t                                                              maxLevel_;
   std::function< uint_t( uint_t level, const Primitive& primitive ) > sizeFunction_;
};

} // namespace hyteg
