
#pragma once

#include "core/debug/Debug.h"
#include "hyteg/primitivedata/PrimitiveDataID.hpp"

#include <memory>

namespace hyteg {

using walberla::mpi::SendBuffer;
using walberla::mpi::RecvBuffer;

// Forward declaration
class Primitive;


/// \brief Base interface for data handling classes
/// \author Nils Kohl (nils.kohl@fau.de)
///
/// The \ref PrimitiveDataHandling class shall be used as a base class for data handling classes.
///
/// When adding data to a \ref Primitive, it is necessary that the data structure has a corresponding
/// implementation of \ref PrimitiveDataHandling. This way it is assured, that the \ref Primitive is
/// able to initialize the data structure and to serialize and deserialize the data.
///
template< typename DataType, typename PrimitiveType >
class PrimitiveDataHandling
{
public:

  typedef DataType value_type;

  virtual ~PrimitiveDataHandling() {}

  /// Initializes the data of type \ref DataType and returns a pointer to the initialized data
  /// Must be thread-safe !
  /// \param primitive the primitive the data is initialized on
  virtual std::shared_ptr< DataType > initialize( const PrimitiveType * const primitive ) const = 0;

  /// Serializes the data of type \ref DataType to a \ref SendBuffer
  /// Must be thread-safe !
  /// \param primitive the primitive the data is taken from
  /// \param id the data index of the data that shall be serialized
  /// \param buffer the buffer it is serialized to
  virtual void serialize( const PrimitiveType * const primitive, const PrimitiveDataID< DataType, PrimitiveType > & id, SendBuffer & buffer ) const = 0;

  /// Deserializes the data of type \ref DataType from a \ref RecvBuffer
  /// Must be thread-safe !
  /// \param primitive the primitive the data shall be written to
  /// \param id the data index of the data that shall be deserialized
  /// \param buffer the buffer it is deserialized from
  virtual void deserialize( const PrimitiveType * const primitive, const PrimitiveDataID< DataType, PrimitiveType > & id, RecvBuffer & buffer ) const = 0;

};


/// \brief Base abstract class for data handling classes that only initialize
/// \author Nils Kohl (nils.kohl@fau.de)
///
/// Contains empty implementations for serialize and deserialize
///
template< typename DataType, typename PrimitiveType >
class OnlyInitializeDataHandling : public PrimitiveDataHandling< DataType, PrimitiveType >
{
public:

  ~OnlyInitializeDataHandling() {}

  /// Does nothing
  void serialize( const PrimitiveType * const primitive, const PrimitiveDataID< DataType, PrimitiveType > & id, SendBuffer & buffer ) const {};

  /// Does nothing
  void deserialize( const PrimitiveType * const primitive, const PrimitiveDataID< DataType, PrimitiveType > & id, RecvBuffer & buffer ) const {};

};




} // namespace hyteg

