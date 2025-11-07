/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Benjamin Mann, Marcus Mohr
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

#include <array>
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

struct PrimitiveIDFormatter;

/// \brief Implementation of the PrimitiveID (unique identifiers for primitives).
///
/// A %PrimitiveID serves as a globally, i.e. across all MPI ranks, unique identifier of a primitive.
///
/// To be able to create new and unique PrimitiveIDs in the case of refinement without the need for
/// global communication in a distributed setting, the %PrimitiveID is encoded in the form of a
/// bitstream, which incorporates the refinement information. This is similar to what is done in
/// the BlockID of waLBerla and described in detail in
///
/// <b>Schornbaum, F. (2018):</b>
/// <a href="https://www10.cs.fau.de/publications/dissertations/Diss_2018-Schornbaum.pdf">Block-structured
/// adaptive mesh refinement for simulations on extreme-scale supercomputers</a>,<br/>
/// <b>doctoral dissertation, Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU)</b>
///
/// The bit pattern of a %PrimitiveID is composed of four parts which from (conceptual) left to right
/// are given by:
///
///  - a block of variable length containing zero bits
///  - a single control or most significant bit (MSB) marking the start of the actual ID by being always set to one
///  - a block of fixed length containing PrimitiveID::BITS_COARSE_LEVEL_ID bits giving either the primitive's ID,
///    if no refinement is used, or the ID of the primitive's grand-grand-parent on the coarsest level,
///    if refinement is used
///  - a block of variable length composed of one set of PrimitiveID::BITS_REFINEMENT bits per refinement level
///
/// The following shows three examples for a choice of BITS_TOTAL = 64, BITS_COARSE_LEVEL_ID = 32 and BITS_REFINEMENT = 5
/*! \htmlonly
  <center>
  <table border="0">
  <tr>
  <td>structure of the pattern for base mesh w/o refinement</td>
  </tr>
  <tr>
  <td><img src="PrimitiveID-noRefinement.png" width="869" height="148"/></td>
  </tr>
  <tr>
  <td>structure of the pattern after one refinement step</td>
  </tr>
  <tr>
  <td><img src="PrimitiveID-oneRefinementStep.png" width="869" height="148"/></td>
  </tr>
  <tr>
  <td>structure of the pattern after the maximally possible six refinement steps</td>
  </tr>
  <tr>
  <td><img src="PrimitiveID-sixRefinementSteps.png" width="869" height="148"/></td>
  </tr>
  </table>
  </center></br>
  \endhtmlonly
*/
/// Using a variable length bitstring will allow for an infinite number of refinement steps. The current implementation,
/// however, uses a fixed sized std::bitset. The length of the latter is fixed at compile time to PrimitiveID::BITS_TOTAL.
/// For internal technical reasons this number must be a multiple of 64.
///
/// The length of the ID for the coarsest primitive is also fixed at compile time by PrimitiveID::BITS_COARSE_LEVEL_ID,
/// while PrimitiveID::BITS_REFINEMENT gives the length of the block required for an individual refinement level. This
/// value limits the number of child primitives a parent primitive can be split into to 2<sup>BITS_REFINEMENT</sup>. Thus,
/// if set to 5, we can construct 32 child primitives. In the case of tetrahedra the number should be (at least) 17.
///
/// Example: Assume we choose the total number of bits to be 128 and reserve 32 bits for the coarse level ID. Then, after
/// subtracting the single control bit, we are left with 95 bits for refinement. If each new level requires 5 bits, this
/// will allow to perform a maximum of 19 refinement step.
///
/// The control bit is required to signal the start of actual ID part of the bit pattern. Without the control bit,
/// refinement may produce duplicate IDs.
///
class PrimitiveID
{
 private:
   friend struct PrimitiveIDFormatter;

   /// Total number of bits used to encode the id !MUST BE A MULTIPLE OF 64!
   static constexpr uint_t BITS_TOTAL = 128;

   /// The number of bits that are reserved for each refinement step.
   static constexpr uint_t BITS_REFINEMENT = 5;

   /// The number of bits that are available for primitives on the coarsest level.
   static constexpr uint_t BITS_COARSE_LEVEL_ID = 32;

   /// The internal data type used to store the PrimitiveID.
   using IDType = std::bitset< BITS_TOTAL >;

   ///@{
   ///
   /// \brief Representation of IDType as an array of integers
   template < typename UINT = uint64_t >
   inline std::array< UINT, BITS_TOTAL / ( sizeof( UINT ) * 8 ) >& asIntArray()
   {
      return *( reinterpret_cast< std::array< UINT, BITS_TOTAL / ( sizeof( UINT ) * 8 ) >* >( &id_ ) );
   }

   template < typename UINT = uint64_t >
   inline const std::array< UINT, BITS_TOTAL / ( sizeof( UINT ) * 8 ) > asIntArray() const
   {
      std::array< UINT, BITS_TOTAL / ( sizeof( UINT ) * 8 ) > result;
      std::memcpy( &result, &id_, sizeof( result ) );
      return result;
   }
   ///@}

   /// Determine the position of the most significant (a.k.a. control) bit
   ///
   /// The control bit marks the start of the bits in the bit pattern that are
   /// truly significant for the ID. The position is counted starting from the
   /// (conceptual) right and one-based.
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

   /// Creates 2<sup>BITS_REFINEMENT</sup> new PrimitiveIDs that can be assigned to new Primitives that result from mesh refinement
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

   // Maximal number of times the function createChildren can be applied recursively on a coarse ID
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

   /// Returns the number of used, i.e. significant, bits in the pattern representing the %PrimitiveID
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

/// Enum for changing output conversion of PrimitiveID objects
enum class PrimitiveIDFormat
{
   FULL,     ///< output complete bit-pattern of ID
   TRIM,     ///< output complete bit-pattern w/o leading turned-off bits and MSB
   COARSE_ID ///< output only COARSE_ID as integer value
};

/// Class for handling output conversion of PrimitiveID objects
struct PrimitiveIDFormatter
{
   static void              set( PrimitiveIDFormat format ) { currentFormat_ = format; }
   static PrimitiveIDFormat get() { return currentFormat_; }
   static std::ostream&     toStream( std::ostream& os, const PrimitiveID& pid )
   {
      switch ( currentFormat_ )
      {
      case PrimitiveIDFormat::FULL: {
         os << pid.id_;
         break;
      }
      case PrimitiveIDFormat::TRIM: {
         std::string pattern = pid.id_.to_string();
         os << pattern.substr( PrimitiveID::BITS_TOTAL - pid.msbPosition() + 1 );
         break;
      }
      case PrimitiveIDFormat::COARSE_ID: {
         PrimitiveID::IDType aux         = pid.id_;
         uint_t              msb         = pid.msbPosition(); // one-based
         uint_t              numChildren = ( msb - PrimitiveID::BITS_COARSE_LEVEL_ID - 1 ) / PrimitiveID::BITS_REFINEMENT;

         aux.reset( msb - 1 );
         aux >> ( numChildren * PrimitiveID::BITS_REFINEMENT );

         os << aux.to_ullong();
         break;
      }
      }

      return os;
   }

 private:
   static inline PrimitiveIDFormat currentFormat_ = PrimitiveIDFormat::FULL;
};

inline std::ostream& PrimitiveID::toStream( std::ostream& os ) const
{
   return PrimitiveIDFormatter::toStream( os, *this );
}

inline std::ostream& operator<<( std::ostream& os, const PrimitiveIDFormat& format )
{
   PrimitiveIDFormatter::set( format );
   return os;
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
