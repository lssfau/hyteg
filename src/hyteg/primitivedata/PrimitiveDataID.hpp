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

#include "core/DataTypes.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

namespace hyteg {

using walberla::uint_t;
using walberla::mpi::RecvBuffer;
using walberla::mpi::SendBuffer;

namespace internal {

/// Simple class that counts references.
class ReferenceCounter
{
 public:
   ReferenceCounter()
   : refs_( 0 )
   {}

   void increaseRefs() { refs_++; }

   void decreaseRefs() { refs_--; }

   uint_t refs() const { return refs_; }

 private:
   uint_t refs_;
};
} // namespace internal

/// \brief Identifier for data attached to \ref Primitive instances
/// \author Nils Kohl (nils.kohl@fau.de)
///
/// Identifies a data item attached to a \ref Primitive.
/// The data can be obtained by passing the respective \ref PrimitiveDataID to
/// the corresponding getter methods.
///
/// Since this class is templated with the data type of the passed data,
/// retrieving the data from primitives is more type-safe and checked via static casts
/// during compile-time.
///
/// PrimitiveDataIDs are generated when data is added to primitives via a storage instance
/// like \ref PrimitiveStorage.
///
/// PrimitiveDataID point to a reference counter. Copying a PrimitiveDataID increases this counter.
/// This allows to safely delete corresponding data when the reference counter drops to zero.
///
template < typename DataType, typename PrimitiveType >
class PrimitiveDataID
{
 public:
   friend class PrimitiveStorage;

   PrimitiveDataID()
   : id_( std::numeric_limits< uint_t >::max() )
   , refCounter_( new internal::ReferenceCounter() )
   {
      refCounter_->increaseRefs();
   }

   explicit PrimitiveDataID( const walberla::uint_t id )
   : id_( id )
   , refCounter_( new internal::ReferenceCounter() )
   {
      refCounter_->increaseRefs();
   }

   /// Copy-constructor
   PrimitiveDataID( const PrimitiveDataID& id )
   : id_( id.id_ )
   , refCounter_( id.refCounter_ )
   {
      refCounter_->increaseRefs();
   }

   /// Copy-assignment
   PrimitiveDataID& operator=( const PrimitiveDataID& id )
   {
      id_         = id.id_;
      refCounter_ = id.refCounter_;
      refCounter_->increaseRefs();
      return *this;
   }

   ~PrimitiveDataID() { refCounter_->decreaseRefs(); }

   bool operator==( const PrimitiveDataID& id ) const { return id_ == id.id_; }
   bool operator!=( const PrimitiveDataID& id ) const { return id_ != id.id_; }
   bool operator<( const PrimitiveDataID& id ) const { return id_ < id.id_; }

   /// Cast-operator
   operator uint_t() const { return id_; }

   /// Returns the number of copies of this ID that are still alive (including this).
   /// If this method returns 1, and this ID is going to be deleted, no other PrimitiveDataID with this ID exists.
   /// The corresponding data can then be deleted, too.
   uint_t numRefs() const { return refCounter_->refs(); }

 private:
   uint_t                                        id_;
   std::shared_ptr< internal::ReferenceCounter > refCounter_;

}; // class DataID

} // namespace hyteg
