
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
  FunctionMemory( const std::function< uint_t ( uint_t level, uint_t numDependencies ) > & sizeFunction,
                  const uint_t & numDependencies,
                  const uint_t & minLevel,
                  const uint_t & maxLevel,
                  const ValueType fillValue = ValueType() ) :
    fillValue_( fillValue )
  {
    WALBERLA_ASSERT_LESS_EQUAL( minLevel, maxLevel, "minLevel should be equal or less than maxLevel during FunctionMemory allocation." );
    for ( uint_t level = minLevel; level <= maxLevel; level++ )
    {
      addLevel( level, sizeFunction( level, numDependencies ), fillValue );
    }
  }

  uint_t getSize( const uint_t & level ) const { WALBERLA_ASSERT_GREATER( data_.count( level ), 0, "Requested level not allocated" ); return data_.at( level )->size(); }

  virtual ~FunctionMemory(){}

  // Returns a pointer to the first entry of the allocated array
  inline ValueType * getPointer( const uint_t & level ) const { return data_.at( level )->data(); }

  /// Serializes the allocated data to a send buffer
  inline void serialize( SendBuffer & sendBuffer ) const
  {
    const uint_t numLevels = data_.size();
    sendBuffer << numLevels;

    for ( const auto & it : data_ )
    {
      const uint_t level     = it.first;
      const uint_t levelSize = it.second->size();

      sendBuffer << level;
      sendBuffer << levelSize;
      sendBuffer << getVector( level );
    }
  }

  /// Deserializes data from a recv buffer (clears all already allocated data and replaces it with the recv buffer's content)
  inline void deserialize( RecvBuffer & recvBuffer )
  {
    data_.clear();

    uint_t numLevels;

    recvBuffer >> numLevels;

    for ( uint_t levelCnt = 0; levelCnt < numLevels; levelCnt++ )
    {
      uint_t level;
      uint_t levelSize;

      recvBuffer >> level;
      recvBuffer >> levelSize;

      addLevel( level, levelSize, fillValue_ );

      std::vector< ValueType > & dataVector = getVector( level );
      recvBuffer >> dataVector;
    }
  }

private:

  inline const std::vector< ValueType > & getVector( const uint_t & level ) const { return *( data_.at( level ) ); }
  inline       std::vector< ValueType > & getVector( const uint_t & level )       { return *( data_[ level ] ); }

  /// Allocates an array of size size for a certain level
  inline void addLevel( const uint_t & level, const uint_t & size, const ValueType & fillValue )
  {
    WALBERLA_ASSERT_EQUAL( data_.count(level), 0, "Attempting to overwrite already existing level (level == " << level << ") in function memory!");
    data_[level] = std::unique_ptr< std::vector< ValueType > >( new std::vector< ValueType >( size, fillValue ) );
  }

  /// Maps a level to the respective allocated data
  std::map< uint_t, std::unique_ptr< std::vector< ValueType > > > data_;

  const ValueType fillValue_;
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

template< typename DataType, typename PrimitiveType >
class MemoryDataHandling : public PrimitiveDataHandling< DataType, PrimitiveType >
{
public:

  MemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel, const std::function< uint_t ( uint_t level, uint_t numDependencies ) > & sizeFunction )
    : minLevel_( minLevel ),
      maxLevel_( maxLevel ),
      sizeFunction_( sizeFunction )
  {}

  virtual ~MemoryDataHandling() {}

    std::shared_ptr<FunctionMemory<real_t> > initialize(const PrimitiveType *const primitive) const {
      return std::make_shared< DataType >(sizeFunction_,
      primitive->getNumHigherDimNeighbors(),
      minLevel_,
      maxLevel_);
  }

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


private:

  const uint_t minLevel_;
  const uint_t maxLevel_;
  const std::function< uint_t ( uint_t level, uint_t numDependencies ) > & sizeFunction_;


};

}
