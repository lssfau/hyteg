/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

#include <blockforest/Utility.h>
#include <core/debug/Debug.h>
#include <core/math/Uint.h>
#include <core/mpi/SendBuffer.h>
#include <core/mpi/RecvBuffer.h>
#include <domain_decomposition/IBlockID.h>

namespace hyteg {

// to be removed when moving to walberla namespace
using walberla::uint8_c;
using walberla::uint64_c;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::uintMSBPosition;
using walberla::numeric_cast;
using walberla::blockforest::uintToBitString;
using walberla::math::UINT_BYTES;

//**********************************************************************************************************************
/*!
*   \brief Implementation of the PrimitiveID (unique identifiers for primitives)
*
*   It serves as a (global) unique identifier of a primitive.
*/
//**********************************************************************************************************************
class PrimitiveID
{
public:

  typedef uint64_t IDType;
  typedef std::vector< PrimitiveID >::const_iterator const_iterator;

  inline PrimitiveID() : id_( uint64_c( 0 ) ) {}
  inline PrimitiveID( const PrimitiveID& id ) : id_( id.id_ ) {}
  inline PrimitiveID( const uint64_t id ) : id_( id ) {}

  bool operator< ( const PrimitiveID& rhs ) const
    { WALBERLA_ASSERT_EQUAL( dynamic_cast< const PrimitiveID* >( &rhs ), &rhs ); return id_ <  static_cast< const PrimitiveID* >( &rhs )->id_; }
  bool operator> ( const PrimitiveID& rhs ) const
    { WALBERLA_ASSERT_EQUAL( dynamic_cast< const PrimitiveID* >( &rhs ), &rhs ); return id_ >  static_cast< const PrimitiveID* >( &rhs )->id_; }
  bool operator==( const PrimitiveID& rhs ) const
    { WALBERLA_ASSERT_EQUAL( dynamic_cast< const PrimitiveID* >( &rhs ), &rhs ); return id_ == static_cast< const PrimitiveID* >( &rhs )->id_; }
  bool operator!=( const PrimitiveID& rhs ) const
    { WALBERLA_ASSERT_EQUAL( dynamic_cast< const PrimitiveID* >( &rhs ), &rhs ); return id_ != static_cast< const PrimitiveID* >( &rhs )->id_; }

  uint_t getUsedBits()  const { return uintMSBPosition( id_ ); }
  uint_t getUsedBytes() const { const uint_t bits( getUsedBits() ); return ( bits >> uint_c(3) ) + ( ( bits & uint_c(7) ) ? uint_c(1) : uint_c(0) ); }

  inline IDType getID() const;

  inline std::ostream& toStream( std::ostream& os ) const;

  template< typename Buffer_T > void   toBuffer( Buffer_T& buffer ) const;
  template< typename Buffer_T > void fromBuffer( Buffer_T& buffer );

private:

  uint64_t id_;

};


inline PrimitiveID::IDType PrimitiveID::getID() const
{
  return numeric_cast< PrimitiveID::IDType > ( id_ );
}


inline std::ostream& PrimitiveID::toStream( std::ostream& os ) const
{
  os << id_;
  return os;
}

inline std::ostream & operator<<( std::ostream & os, const PrimitiveID & id )
{
  id.toStream( os );
  return os;
}

template< typename Buffer_T >
void PrimitiveID::toBuffer( Buffer_T& buffer ) const {

  const uint8_t bytes = uint8_t( getUsedBytes() );
  buffer << bytes;
  for( uint8_t i = 0; i != bytes; ++i )
    buffer << uint8_c( ( id_ >> ( uint_c(i) * uint_c(8) ) ) & uint_c(255) );
}


template< typename Buffer_T >
void PrimitiveID::fromBuffer( Buffer_T& buffer ) {

  uint8_t bytes(0);
  buffer >> bytes;

  WALBERLA_ASSERT_LESS_EQUAL( bytes, UINT_BYTES );

  id_ = uint_c(0);
  for( uint8_t i = 0; i != bytes; ++i ) {
     uint8_t byte(0);
     buffer >> byte;
     id_ |= uint_c( byte ) << ( uint_c(i) * uint_c(8) );
  }
}


} // namespace hyteg


//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T,  // Element type of SendBuffer
          typename G > // Growth policy of SendBuffer
inline mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const hyteg::PrimitiveID & id )
{
   buffer.addDebugMarker( "bi" );
   id.toBuffer( buffer );
   return buffer;
}

template< typename T > // Element type  of RecvBuffer
inline mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buffer, hyteg::PrimitiveID & id )
{
   buffer.readDebugMarker( "bi" );
   id.fromBuffer( buffer );
   return buffer;
}

template<>
struct BufferSizeTrait< hyteg::PrimitiveID > { static const bool constantSize = false; };

} // namespace mpi
} // namespace walberla
