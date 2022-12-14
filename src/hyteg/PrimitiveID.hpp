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

#include <bitset>

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
/// The length of the short bitstrings for each refinement level is fixed at compile time. If set to 5, refinement of
/// Primitives is limited to 32 child primitives (17 is required for tets).
///
/// The marker bit is required to signal the start of the ID. Without the marker bit, refinement may produce duplicate IDs.
///
/// For example: if we encode the ID in a 128-bit unsigned integer, subtract one marker bit, reserve 32 bits for the coarse level
/// ID, we have 95 bits left for refinement. If each level requires 5 bits, this allows for 19 additional refinement levels.
///
class PrimitiveID
{
 private:
   /// Total number of bits used to encode the id !MUST BE A MULTIPLE OF 64!
   static constexpr uint_t BITS_TOTAL = 128;

   /// The number of bits that are reserved for each refinement step.
   static constexpr uint_t BITS_REFINEMENT = 5;

   /// The number of bits that are available for primitives on the coarsest level.
   static constexpr uint_t BITS_COARSE_LEVEL_ID = 32;

   /// The internal data type used to store the PrimitiveID.
   using IDType = std::bitset< BITS_TOTAL >;

   /// Representation of IDType as an array of integers
   template < typename UINT = uint64_t >
   inline std::array< UINT, BITS_TOTAL / ( sizeof( UINT ) * 8 ) >& asIntArray()
   {
      return *( reinterpret_cast< std::array< UINT, BITS_TOTAL / ( sizeof( UINT ) * 8 ) >* >( &id_ ) );
   }
   /// Representation of IDType as an array of integers
   template < typename UINT = uint64_t >
   inline const std::array< UINT, BITS_TOTAL / ( sizeof( UINT ) * 8 ) >& asIntArray() const
   {
      return *( reinterpret_cast< std::array< UINT, BITS_TOTAL / ( sizeof( UINT ) * 8 ) > const* >( &id_ ) );
   }

   /// compute the position of the most significant bit
   inline uint_t msbPosition() const
   {
      for ( uint_t i = BITS_TOTAL; i > 0; --i )
      {
         if ( id_[i - 1] )
            return i;
      }
      return 0;
   }

 public:
   using const_iterator = std::vector< PrimitiveID >::const_iterator;

   /// Creates a unique coarse level PrimitiveID safely.
   inline static PrimitiveID create( const uint64_t& coarseID )
   {
      WALBERLA_CHECK_LESS( walberla::math::uintMSBPosition( coarseID ),
                           BITS_COARSE_LEVEL_ID + 1,
                           "Could not construct PrimitiveID with coarse ID " + std::to_string( coarseID ) + "." );

      PrimitiveID pid;
      pid.id_ = ( uint64_t( 1 ) << BITS_COARSE_LEVEL_ID ) | coarseID;
      return pid;
   }

   /// Creates an invalid PrimitiveID.
   inline PrimitiveID()
   : id_( IDType( 0 ) )
   {}

   /// Creates 2**BITS_REFINEMENT new PrimitiveIDs that can be assigned to new Primitives that result from mesh refinement.
   inline auto createChildren() const
   {
      WALBERLA_CHECK( id_ != IDType( 0 ), "Can't refineme invalid PrimitiveID" );
      WALBERLA_CHECK_LESS( numAncestors(), maxRefinement(), "Could not construct childIDs" );

      std::array< PrimitiveID, ( 1u << BITS_REFINEMENT ) > children;
      for ( uint_t i = 0; i < children.size(); ++i )
      {
         children[i].id_ = ( id_ << BITS_REFINEMENT ) | IDType( i );
      }
      return children;
   }

   // maximum number of times the function createChildren can be applied recursively on a coarse ID
   static inline uint_t maxRefinement() { return ( BITS_TOTAL - BITS_COARSE_LEVEL_ID - 1 ) / BITS_REFINEMENT; }

   inline uint_t numAncestors() const { return ( msbPosition() - BITS_COARSE_LEVEL_ID - 1 ) / BITS_REFINEMENT; }

   inline bool hasAncestors() const { return numAncestors() > 0; }

   inline PrimitiveID getParent() const
   {
      PrimitiveID pid;
      pid.id_ = id_ >> BITS_REFINEMENT;
      return pid;
   }

   bool operator<( const PrimitiveID& rhs ) const { return asIntArray() < rhs.asIntArray(); }

   bool operator>( const PrimitiveID& rhs ) const { return asIntArray() > rhs.asIntArray(); }

   bool operator==( const PrimitiveID& rhs ) const { return id_ == rhs.id_; }

   bool operator!=( const PrimitiveID& rhs ) const { return id_ != rhs.id_; }

   uint_t getUsedBits() const { return msbPosition(); }

   uint_t getUsedBytes() const
   {
      const uint_t bits = getUsedBits();
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
   for ( auto& byte : asIntArray< uint8_t >() )
      buffer << byte;
}

template < typename Buffer_T >
void PrimitiveID::fromBuffer( Buffer_T& buffer )
{
   for ( auto& byte : asIntArray< uint8_t >() )
      buffer >> byte;
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
   buffer.addDebugMarker( "pi" );
   id.toBuffer( buffer );
   return buffer;
}

template < typename T > // Element type  of RecvBuffer
inline mpi::GenericRecvBuffer< T >& operator>>( mpi::GenericRecvBuffer< T >& buffer, hyteg::PrimitiveID& id )
{
   buffer.readDebugMarker( "pi" );
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
