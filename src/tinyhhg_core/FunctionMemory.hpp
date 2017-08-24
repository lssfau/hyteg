
#pragma once

#include "tinyhhg_core/support.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"

#include <core/DataTypes.h>
#include <core/logging/Logging.h>
#include <core/mpi/SendBuffer.h>
#include <core/mpi/RecvBuffer.h>

#include <map>
#include <memory>


namespace hhg {

using walberla::uint_t;
using walberla::real_t;
using walberla::mpi::SendBuffer;
using walberla::mpi::RecvBuffer;

class FunctionMemory
{
public:

  /// Constructs memory for a function
  FunctionMemory( const uint_t & numDependencies ) : numDependencies_( numDependencies ) {}

  virtual ~FunctionMemory(){}

  /// Returns the length of the allocated array for a certain level
  virtual uint_t getSize( uint_t level ) const = 0;

  /// Allocates an array for a certain level
  /// Uses the virtual member \ref getSize() to determine the length of the array
  inline std::unique_ptr< real_t[] > & addlevel( uint_t level )
  {
    if (data.count(level)>0)
    {
      WALBERLA_LOG_WARNING("Level already exists.");
    }
    else
    {
      data[level] = hhg::make_unique<real_t[]>(getSize(level));
    }
    return data[level];
  }

  /// Serializes the allocated data to a send buffer
  inline void serialize( SendBuffer & sendBuffer ) const
  {
    for ( const auto & it : data )
    {
      uint_t level = it.first;
      const real_t * dataPtr = it.second.get();
      sendBuffer << level;
      for ( uint_t idx = 0; idx < getSize( level ); idx++ )
      {
        sendBuffer << dataPtr[ idx ];
      }
    }
  }

  /// Deserializes the allocated data from a recv buffer
  inline void deserialize( RecvBuffer & recvBuffer )
  {
    data.clear();
    while ( !recvBuffer.isEmpty() )
    {
      uint_t level;
      recvBuffer >> level;
      data[ level ] = hhg::make_unique< real_t[] >( getSize( level ) );
      for ( uint_t idx = 0; idx < getSize( level ); idx++ )
      {
        recvBuffer >> data[ level ][ idx ];
      }
    }
  }

  /// Maps a level to the respective allocated data
  std::map< uint_t, std::unique_ptr< real_t[] > > data;

protected:

  uint_t numDependencies_;

};

}

namespace walberla {
namespace mpi {

template< typename T,  // Element type of SendBuffer
          typename G > // Growth policy of SendBuffer
inline mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const hhg::FunctionMemory & functionMemory )
{
   functionMemory.serialize( buffer );
   return buffer;
}

template< typename T > // Element type  of RecvBuffer
inline mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buffer, hhg::FunctionMemory & functionMemory )
{
   functionMemory.deserialize( buffer );
   return buffer;
}

template<>
struct BufferSizeTrait< hhg::FunctionMemory > { static const bool constantSize = false; };

} // namespace mpi
} // namespace walberla


namespace hhg {


template< typename DataType, typename PrimitiveType >
class FunctionMemoryDataHandling : public PrimitiveDataHandling< DataType, PrimitiveType >
{
public:
  static_assert( std::is_base_of< FunctionMemory, DataType >::value,
                 "FunctionMemoryDataHandling can only be used for subclasses of FunctionMemory!" );

  virtual ~FunctionMemoryDataHandling() {}

  virtual void serialize( const PrimitiveType * const primitive, const PrimitiveDataID< DataType, PrimitiveType > & id, SendBuffer & buffer ) const
  {
    DataType * data = primitive->getData( id );
    buffer << *data;
  }

  virtual void deserialize( const PrimitiveType * const primitive, const PrimitiveDataID< DataType, PrimitiveType > & id, RecvBuffer & buffer ) const
  {
    DataType * data = primitive->getData( id );
    buffer >> *data;
  }

};


}
