
#pragma once

#include "core/DataTypes.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"

namespace walberla {
namespace hhg {

template< typename DataType >
class PrimitiveDataID
{
public:

	   PrimitiveDataID()                   : id_( 0 ) {}
  explicit PrimitiveDataID( const uint_t id )  : id_( id ) {}
 	   PrimitiveDataID( const PrimitiveDataID& id ) : id_( id.id_ ) {}

  void pack( mpi::SendBuffer & buffer ) const { buffer << id_; }
  void unpack( mpi::RecvBuffer & buffer )     { buffer >> id_; }

  PrimitiveDataID& operator=( const PrimitiveDataID& id ) { id_ = id.id_; return *this; }

  bool operator==( const PrimitiveDataID& id ) const { return id_ == id.id_; }
  bool operator!=( const PrimitiveDataID& id ) const { return id_ != id.id_; }
  bool operator< ( const PrimitiveDataID& id ) const { return id_ <  id.id_; }

  operator uint_t() const { return id_; }

private:

  uint_t id_;

}; // class DataID


}
}

