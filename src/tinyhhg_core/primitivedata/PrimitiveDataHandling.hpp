
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/primitives/Primitive.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"

#include <memory>

namespace hhg {

using walberla::mpi::SendBuffer;
using walberla::mpi::RecvBuffer;

class Primitive;

template< typename DataType >
class PrimitiveDataHandling
{
public:

  typedef DataType value_type;

  virtual ~PrimitiveDataHandling() {}

  /// must be thread-safe !
  virtual DataType * initialize( Primitive * const primitive ) const = 0;

  /// must be thread-safe !
  virtual void serialize( Primitive * const primitive, const PrimitiveDataID< DataType > & id, SendBuffer & buffer ) const = 0;

  /// must be thread-safe !
  virtual void deserialize( Primitive * const primitive, const PrimitiveDataID< DataType > & id, RecvBuffer & buffer ) const = 0;
};


template< typename DataType >
class NoSerializePrimitiveDataHandling : public PrimitiveDataHandling< DataType >
{
public:

  ~NoSerializePrimitiveDataHandling() {}

  void serialize( Primitive * const primitive, const PrimitiveDataID< DataType > & id, SendBuffer & buffer ) const {};
  void deserialize( Primitive * const primitive, const PrimitiveDataID< DataType > & id, RecvBuffer & buffer ) const {};

};




} // namespace hhg

