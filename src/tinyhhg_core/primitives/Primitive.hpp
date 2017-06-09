
#pragma once

#include <tinyhhg_core/primitivedata/PrimitiveDataID.hpp>
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


class Primitive : private NonCopyable
{
public:

  typedef internal::PrimitiveData PrimitiveData;

  template< typename DataType >
  DataType* getData( const PrimitiveDataID< DataType > & index )
  {
    return data_[ index ]->template get< DataType >();
  }

  template< typename DataType >
  void addData( const PrimitiveDataID< DataType > & index, PrimitiveData * const data )
  {
    if( data_.size() <= index )
    {
      data_.resize( index + 1, NULL );
    }

    if( data != NULL )
    {
      WALBERLA_ASSERT_NULLPTR( data_[index] );
      data_[index] = data;
    }
  }

private:

  std::vector< PrimitiveData* > data_; ///< the data assigned to this primitive

  PrimitiveID primitiveID_;

};



}
}
