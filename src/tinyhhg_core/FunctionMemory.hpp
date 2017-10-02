
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"

#include <core/DataTypes.h>
#include <core/logging/Logging.h>
#include <core/mpi/SendBuffer.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/BufferSizeTrait.h>

#include <map>
#include <memory>


namespace hhg {

using walberla::uint_t;
using walberla::real_t;
using walberla::mpi::SendBuffer;
using walberla::mpi::RecvBuffer;

template< typename ValueType,
          typename = typename std::enable_if< std::is_arithmetic< ValueType >::value >::type >
class FunctionMemory
{
public:

  /// Constructs memory for a function
  FunctionMemory( const uint_t & numDependencies ) : numDependencies_( numDependencies ) {}

  virtual ~FunctionMemory(){}

  /// Returns the length of the allocated array for a certain level
  virtual uint_t getSize( uint_t level ) const = 0;

  // Returns a pointer to the first entry of the allocated array
  inline ValueType * getPointer( const uint_t & level ) const { return data_.at( level )->data(); }

  /// Allocates an array for a certain level
  /// Uses the virtual member \ref getSize() to determine the length of the array
  inline void addlevel( const uint_t & level )
  {
    WALBERLA_ASSERT_EQUAL( data_.count(level), 0, "Attempting to overwrite already existing level (level == " << level << ") in function memory!");
    data_[level] = std::unique_ptr< std::vector< ValueType > >( new std::vector< ValueType >( getSize( level ) ) );
  }

  /// Serializes the allocated data to a send buffer
  inline void serialize( SendBuffer & sendBuffer ) const
  {
    sendBuffer << numDependencies_;
    sendBuffer << data_.size();

    for ( const auto & it : data_ )
    {
      const uint_t level = it.first;

      sendBuffer << level;
      sendBuffer << getVector( level );
    }
  }

  /// Deserializes data from a recv buffer (clears all already allocated data and replaces it with the recv buffer's content)
  inline void deserialize( RecvBuffer & recvBuffer )
  {
    data_.clear();

    uint_t numLevels;

    recvBuffer >> numDependencies_;
    recvBuffer >> numLevels;

    for ( uint_t levelCnt = 0; levelCnt < numLevels; levelCnt++ )
    {
      uint_t level;
      recvBuffer >> level;

      addlevel( level );

      std::vector< ValueType > & dataVector = getVector( level );
      recvBuffer >> dataVector;
    }
  }

private:

  inline const std::vector< ValueType > & getVector( const uint_t & level ) const { return *( data_.at( level ) ); }
  inline       std::vector< ValueType > & getVector( const uint_t & level )       { return *( data_[ level ] ); }

  /// Maps a level to the respective allocated data
  std::map< uint_t, std::unique_ptr< std::vector< ValueType > > > data_;

protected:

  uint_t numDependencies_;

};

}

namespace walberla {
namespace mpi {

template< typename T,  // Element type of SendBuffer
          typename G,  // Growth policy of SendBuffer
          typename ValueType>
inline mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const hhg::FunctionMemory< ValueType > & functionMemory )
{
   functionMemory.serialize( buffer );
   return buffer;
}

template< typename T, // Element type  of RecvBuffer
          typename ValueType >
inline mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buffer, hhg::FunctionMemory< ValueType > & functionMemory )
{
   functionMemory.deserialize( buffer );
   return buffer;
}

template< typename ValueType >
struct BufferSizeTrait< hhg::FunctionMemory< ValueType > > { static const bool constantSize = false; };

} // namespace mpi
} // namespace walberla


namespace hhg {


template< typename DataType, typename PrimitiveType >
class FunctionMemoryDataHandling : public PrimitiveDataHandling< DataType, PrimitiveType >
{
public:

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
