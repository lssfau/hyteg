
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "core/NonCopyable.h"
#include "tinyhhg_core/primitiveid.hpp"

#include <memory>
#include <vector>

namespace walberla {
namespace hhg {

namespace internal {

class PrimitiveData : private NonCopyable
{
public:

  template< typename DataType >
  PrimitiveData( DataType* ptr ) : ptr_( ptr )
  {
    WALBERLA_ASSERT_NOT_NULLPTR( ptr );
  }

  template< typename DataType >
  DataType* get()
  {
    return static_cast< DataType* >( ptr_ );
  }

private:

  void *ptr_;

};

}

class PrimitiveID;

class Primitive : private NonCopyable
{
public:

  friend class PrimitiveStorage;

  typedef internal::PrimitiveData PrimitiveData;

  template< typename DataType >
  DataType* getData( const PrimitiveDataID< DataType > & index )
  {
    return data_[ index ].first->template get< DataType >();
  }

  template< typename DataType >
  PrimitiveDataHandling< DataType >* getDataHandling( const PrimitiveDataID< DataType > & index )
  {
    return data_[ index ].second->template get< PrimitiveDataHandling< DataType > >();
  }

  uint_t getNumberOfPrimitiveDataEntries() const
  {
    return data_.size();
  }

private:

  /// Must stay private in order to guarantee that data is only added through the governing structure.
  /// This ensures valid DataIDs.
  template< typename DataType >
  void addData( const PrimitiveDataID< DataType > & index,
		      PrimitiveDataHandling< DataType > & dataHandling )
  {
    if( data_.size() <= index )
    {
      data_.resize( index + 1, std::pair< PrimitiveData*, PrimitiveData* >( NULL, NULL ) );
    }

    data_[index].first = new PrimitiveData( dataHandling.initialize( NULL ) );
    data_[index].second = new PrimitiveData( &dataHandling );
  }

  /// Holds a pointer to the actual data in the first entry and a pointer to the respective datahandling in the second entry.
  /// This way it is possible to loop over the data to for example serialize all registered data.
  std::vector< std::pair< PrimitiveData*, PrimitiveData* > > data_;

  PrimitiveID primitiveID_;

};



}
}
