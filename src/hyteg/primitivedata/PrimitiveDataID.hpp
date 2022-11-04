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
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"


namespace hyteg {

using walberla::mpi::SendBuffer;
using walberla::mpi::RecvBuffer;
using walberla::uint_t;

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
template< typename DataType, typename PrimitiveType >
class PrimitiveDataID
{
public:

  friend class PrimitiveStorage;

  PrimitiveDataID()                            : id_( std::numeric_limits< uint_t >::max() ) {}
  explicit PrimitiveDataID( const walberla::uint_t id ) : id_( id ) {}

  /// Copy-constructor
  PrimitiveDataID( const PrimitiveDataID& id ) : id_( id.id_ ) {}

  void pack( SendBuffer & buffer ) const { buffer << id_; }
  void unpack( RecvBuffer & buffer )     { buffer >> id_; }

  /// Copy-assignment
  PrimitiveDataID& operator=( const PrimitiveDataID& id ) { id_ = id.id_; return *this; }

  bool operator==( const PrimitiveDataID& id ) const { return id_ == id.id_; }
  bool operator!=( const PrimitiveDataID& id ) const { return id_ != id.id_; }
  bool operator< ( const PrimitiveDataID& id ) const { return id_ <  id.id_; }

  /// Cast-operator
  operator uint_t() const { return id_; }

private:

  uint_t id_;

}; // class DataID

}

