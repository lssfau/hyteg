
#include "IPrimitive.hpp"

namespace walberla {
namespace hhg {


bool IPrimitive::isBlockDataAllocated( const ConstBlockDataID & index ) const {

   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );

   return ( data_[index] != NULL );
}


template< typename T >
const T* IPrimitive::getData( const ConstBlockDataID & index ) const {

   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );

   if( data_[index] == NULL )
      return NULL;

   return data_[index]->template get< T >();
}


template< typename T >
const T* IPrimitive::getData( const BlockDataID & index ) const {

   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );

   if( data_[index] == NULL )
      return NULL;

   return data_[index]->template get< T >();
}


template< typename T >
T* IPrimitive::getData( const BlockDataID & index ) {

   return const_cast< T* >( static_cast< const IBlock* >( this )->getData<T>( index ) );
}




inline void IPrimitive::deleteData( const BlockDataID & index )
{
   WALBERLA_ASSERT_LESS( uint_t( index ), data_.size() );
   delete data_[index];
   data_[index] = NULL;
}

} // namespace hhg
} // namespace walberla
