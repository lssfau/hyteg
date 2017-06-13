
#pragma once

#include "core/logging/Logging.h"
#include "tinyhhg_core/primitives/Primitive.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"

#include <map>
#include <vector>

namespace hhg {

class PrimitiveStorage : private walberla::NonCopyable
{
public:

  PrimitiveStorage() : primitiveDataHandlers_( uint_c( 0 ) ) {}

        PrimitiveID addPrimitive();

  const Primitive* getPrimitive( const PrimitiveID & id ) const;
        Primitive* getPrimitive( const PrimitiveID & id );

  bool primitiveExistsLocally( const PrimitiveID & id ) const;

  template< typename DataType >
  inline PrimitiveDataID< DataType > addPrimitiveData( PrimitiveDataHandling< DataType > & dataHandling, const std::string & identifier );

private:

  std::map< PrimitiveID::IDType, Primitive* > primitives_;
  uint_t primitiveDataHandlers_;

};

template< typename DataType >
PrimitiveDataID< DataType > PrimitiveStorage::addPrimitiveData( PrimitiveDataHandling< DataType > & dataHandling, const std::string & identifier )
{
  WALBERLA_LOG_PROGRESS( "Adding block data (\"" << identifier << "\")" );

#ifndef NDEBUG
  for ( auto it = primitives_.begin(); it != primitives_.end(); it++ )
  {
    WALBERLA_ASSERT_EQUAL( primitiveDataHandlers_, it->second->getNumberOfPrimitiveDataEntries() );
  }
#endif

  PrimitiveDataID< DataType > dataID( primitiveDataHandlers_++ );

  for ( auto it = primitives_.begin(); it != primitives_.end(); it++ )
  {
    it->second->addData( dataID, dataHandling );
  }

  return dataID;
}


} // namespace hhg


