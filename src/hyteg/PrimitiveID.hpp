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

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/math/Uint.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "blockforest/Utility.h"
#include "domain_decomposition/IBlockID.h"

namespace hyteg {

using walberla::numeric_cast;
using walberla::uint64_c;
using walberla::uint8_c;
using walberla::uint_c;
using walberla::uint_t;
using walberla::blockforest::uintToBitString;
using walberla::math::UINT_BYTES;
using walberla::math::uintMSBPosition;

/// \brief Implementation of the PrimitiveID (unique identifiers for primitives).
///
/// It serves as a (global) unique identifier of a primitive.
///
/// To be able to create new, unique PrimitiveIDs during refinement without requiring global communication in a distributed
/// setting, the refinement is encoded in the bitstring. This is similar to what is done in walberla's BlockID and described e.g.
/// in
///
///     Schornbaum, F. (2018). Block-structured adaptive mesh refinement for simulations on extreme-scale supercomputers
///     (Doctoral dissertation, Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU)).
///
/// Example:
///
///   0000..000 1 00110011 000 010 000
///   |         | |        |
///   a         b c        d
///
/// a) leading zeros
/// b) marker bit (always one)
/// c) the ID of the coarsest primitive (here, 8 bits are reserved)
/// d) one additional, fixed-length bitstring for each refinement level (here, 3 bits are reserved)
///
/// The size of the bitstring can be either determined at compile time (e.g. 64 bit unsigned) or the PrimitiveID's representation
/// is implemented via a variable length bitstring. The latter case allows for infinite refinement.
///
/// The length of the ID for the coarsest Primitive is fixed at compile time.
///
/// The length of the short bitstrings for each refinement level is fixed at compile time. If set to 3, refinement of
/// Primitives is limited to 8 child primitives (which is sufficient for tets).
///
/// The marker bit is required to signal the start of the ID. Without the marker bit, refinement may produce duplicate IDs.
///
/// For example: if we encode the ID in a 64-bit unsigned integer, subtract one marker bit, reserve 30 bits for the coarse level
/// ID, we have 33 bits left for refinement. If each level requires 3 bits, this allows for 10 additional refinement levels.
///
class PrimitiveID
{
 private:
   /// The internal data type used to store the PrimitiveID.
   typedef uint64_t IDType;

   /// The number of bits that are reserved for each refinement step.
   static const uint_t BITS_REFINEMENT = 3;

   /// The number of bits that are available for primitives on the coarsest level.
   static const uint_t BITS_COARSE_LEVEL_ID = 30;

 public:
   typedef std::vector< PrimitiveID >::const_iterator const_iterator;

   /// Creates a unique coarse level PrimitiveID safely.
   inline static PrimitiveID create( const uint64_t& coarseID )
   {
      WALBERLA_CHECK_LESS( walberla::math::uintMSBPosition( coarseID ),
                           BITS_COARSE_LEVEL_ID + 1,
                           "Could not construct PrimitiveID with coarse ID " + std::to_string( coarseID ) + "." );

      PrimitiveID pid;
      pid.id_ = ( 1 << BITS_COARSE_LEVEL_ID ) | coarseID;
      return pid;
   }

   /// Creates an invalid PrimitiveID.
   inline PrimitiveID()
   : id_( IDType( 0 ) )
   {}

   /// Creates 2**BITS_REFINEMENT new PrimitiveIDs that can be assigned to new Primitives that result from mesh refinement.
   inline std::vector< PrimitiveID > createChildren() const
   {
      std::vector< PrimitiveID > children;
      const auto                 numChildren = walberla::math::uintPow2( BITS_REFINEMENT );
      for ( uint_t i = 0; i < numChildren; i++ )
      {
         auto childID = ( id_ << BITS_REFINEMENT ) | i;
         // children.push_back( childID );
      }
      return children;
   }

   inline uint_t numParents() const
   {
      return ( walberla::math::uintMSBPosition( id_ ) - BITS_COARSE_LEVEL_ID - 1 ) / BITS_REFINEMENT;
   }

   inline bool hasParents() const { return numParents() > 0; }

   bool operator<( const PrimitiveID& rhs ) const
   {
      WALBERLA_ASSERT_EQUAL( dynamic_cast< const PrimitiveID* >( &rhs ), &rhs );
      return id_ < static_cast< const PrimitiveID* >( &rhs )->id_;
   }

   bool operator>( const PrimitiveID& rhs ) const
   {
      WALBERLA_ASSERT_EQUAL( dynamic_cast< const PrimitiveID* >( &rhs ), &rhs );
      return id_ > static_cast< const PrimitiveID* >( &rhs )->id_;
   }

   bool operator==( const PrimitiveID& rhs ) const
   {
      WALBERLA_ASSERT_EQUAL( dynamic_cast< const PrimitiveID* >( &rhs ), &rhs );
      return id_ == static_cast< const PrimitiveID* >( &rhs )->id_;
   }

   bool operator!=( const PrimitiveID& rhs ) const
   {
      WALBERLA_ASSERT_EQUAL( dynamic_cast< const PrimitiveID* >( &rhs ), &rhs );
      return id_ != static_cast< const PrimitiveID* >( &rhs )->id_;
   }

   uint_t getUsedBits() const { return uintMSBPosition( id_ ); }

   uint_t getUsedBytes() const
   {
      const uint_t bits( getUsedBits() );
      return ( bits >> uint_c( 3 ) ) + ( ( bits & uint_c( 7 ) ) ? uint_c( 1 ) : uint_c( 0 ) );
   }

   inline std::ostream& toStream( std::ostream& os ) const;

   template < typename Buffer_T >
   void toBuffer( Buffer_T& buffer ) const;

   template < typename Buffer_T >
   void fromBuffer( Buffer_T& buffer );

 private:
   IDType id_;
};

inline std::ostream& PrimitiveID::toStream( std::ostream& os ) const
{
   os << id_;
   return os;
}

inline std::ostream& operator<<( std::ostream& os, const PrimitiveID& id )
{
   id.toStream( os );
   return os;
}

template < typename Buffer_T >
void PrimitiveID::toBuffer( Buffer_T& buffer ) const
{
   const uint8_t bytes = uint8_t( getUsedBytes() );
   buffer << bytes;
   for ( uint8_t i = 0; i != bytes; ++i )
      buffer << uint8_c( ( id_ >> ( uint_c( i ) * uint_c( 8 ) ) ) & uint_c( 255 ) );
}

template < typename Buffer_T >
void PrimitiveID::fromBuffer( Buffer_T& buffer )
{
   uint8_t bytes( 0 );
   buffer >> bytes;

   WALBERLA_ASSERT_LESS_EQUAL( bytes, UINT_BYTES );

   id_ = uint_c( 0 );
   for ( uint8_t i = 0; i != bytes; ++i )
   {
      uint8_t byte( 0 );
      buffer >> byte;
      id_ |= uint_c( byte ) << ( uint_c( i ) * uint_c( 8 ) );
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

template < typename T,  // Element type of SendBuffer
           typename G > // Growth policy of SendBuffer
inline mpi::GenericSendBuffer< T, G >& operator<<( mpi::GenericSendBuffer< T, G >& buffer, const hyteg::PrimitiveID& id )
{
   buffer.addDebugMarker( "bi" );
   id.toBuffer( buffer );
   return buffer;
}

template < typename T > // Element type  of RecvBuffer
inline mpi::GenericRecvBuffer< T >& operator>>( mpi::GenericRecvBuffer< T >& buffer, hyteg::PrimitiveID& id )
{
   buffer.readDebugMarker( "bi" );
   id.fromBuffer( buffer );
   return buffer;
}

template <>
struct BufferSizeTrait< hyteg::PrimitiveID >
{
   static const bool constantSize = false;
};

} // namespace mpi
} // namespace walberla
