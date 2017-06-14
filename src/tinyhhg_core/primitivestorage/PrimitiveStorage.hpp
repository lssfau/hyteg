
#pragma once

#include "core/logging/Logging.h"
#include "tinyhhg_core/primitives/Primitive.hpp"
#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

#include <map>
#include <vector>

namespace hhg {

class PrimitiveStorage : private walberla::NonCopyable
{
public:

        PrimitiveStorage( const SetupPrimitiveStorage & setupStorage ) :
          primitiveDataHandlers_( 0 )
        {}

        PrimitiveID addVertex();

  const Vertex* getVertex( const PrimitiveID & id ) const;
        Vertex* getVertex( const PrimitiveID & id );

  bool primitiveExistsLocally( const PrimitiveID & id ) const;

  template< typename DataType >
  inline PrimitiveDataID< DataType > addPrimitiveData( const PrimitiveDataHandling< DataType > & dataHandling,
  						       const std::string & identifier );

  template< typename DataType >
  inline PrimitiveDataID< DataType > addVertexData( const PrimitiveDataHandling< DataType > & dataHandling,
						    const std::string & identifier );

  template< typename DataType >
  inline PrimitiveDataID< DataType > addEdgeData( const PrimitiveDataHandling< DataType > & dataHandling,
  						  const std::string & identifier );

  template< typename DataType >
  inline PrimitiveDataID< DataType > addFaceData( const PrimitiveDataHandling< DataType > & dataHandling,
  						  const std::string & identifier );

private:

  template< typename DataType >
  inline PrimitiveDataID< DataType > generateDataID();

  template< typename DataType,
            typename PrimitiveType,
	    typename = typename std::enable_if< std::is_base_of< Primitive, PrimitiveType >::value >::type >
  inline void addPrimitiveData( const PrimitiveDataHandling< DataType > & dataHandling,
				const std::string & identifier,
				const std::map< PrimitiveID::IDType, PrimitiveType* > & primitives,
				const PrimitiveDataID< DataType > & dataID );

  std::map< PrimitiveID::IDType, Vertex* > vertices_;
  std::map< PrimitiveID::IDType, Edge*  >  edges_;
  std::map< PrimitiveID::IDType, Face* >   faces_;

  uint_t primitiveDataHandlers_;

};


template< typename DataType >
PrimitiveDataID< DataType > PrimitiveStorage::addPrimitiveData( const PrimitiveDataHandling< DataType > & dataHandling,
						                const std::string & identifier )
{
  WALBERLA_LOG_PROGRESS( "Adding data to all primitives (\"" << identifier << "\")" );
  PrimitiveDataID< DataType > dataID = generateDataID< DataType >();
  addPrimitiveData( dataHandling, identifier, vertices_, dataID );
  addPrimitiveData( dataHandling, identifier, edges_, dataID );
  addPrimitiveData( dataHandling, identifier, faces_, dataID );
  return dataID;

}


template< typename DataType >
PrimitiveDataID< DataType > PrimitiveStorage::addVertexData( const PrimitiveDataHandling< DataType > & dataHandling,
						             const std::string & identifier )
{
  WALBERLA_LOG_PROGRESS( "Adding data to vertices (\"" << identifier << "\")" );
  PrimitiveDataID< DataType > dataID = generateDataID< DataType >();
  addPrimitiveData( dataHandling, identifier, vertices_, dataID );
  return dataID;
}


template< typename DataType >
PrimitiveDataID< DataType > PrimitiveStorage::addEdgeData( const PrimitiveDataHandling< DataType > & dataHandling,
						           const std::string & identifier )
{
  WALBERLA_LOG_PROGRESS( "Adding data to edges (\"" << identifier << "\")" );
  PrimitiveDataID< DataType > dataID = generateDataID< DataType >();
  addPrimitiveData( dataHandling, identifier, edges_, dataID );
  return dataID;
}


template< typename DataType >
PrimitiveDataID< DataType > PrimitiveStorage::addFaceData( const PrimitiveDataHandling< DataType > & dataHandling,
						             const std::string & identifier )
{
  WALBERLA_LOG_PROGRESS( "Adding data to faces (\"" << identifier << "\")" );
  PrimitiveDataID< DataType > dataID = generateDataID< DataType >();
  addPrimitiveData( dataHandling, identifier, faces_, dataID );
  return dataID;
}


template< typename DataType >
PrimitiveDataID< DataType > PrimitiveStorage::generateDataID()
{
#ifndef NDEBUG
  for ( auto it = vertices_.begin(); it != vertices_.end(); it++ )
  {
    WALBERLA_ASSERT_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
  }
  for ( auto it = edges_.begin(); it != edges_.end(); it++ )
  {
    WALBERLA_ASSERT_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
  }
  for ( auto it = faces_.begin(); it != faces_.end(); it++ )
  {
    WALBERLA_ASSERT_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
  }
#endif

  return PrimitiveDataID< DataType >( primitiveDataHandlers_++ );
}


template< typename DataType,
          typename PrimitiveType,
	  typename >
void PrimitiveStorage::addPrimitiveData( const PrimitiveDataHandling< DataType > & dataHandling,
					 const std::string & identifier,
					 const std::map< PrimitiveID::IDType, PrimitiveType* > & primitives,
					 const PrimitiveDataID< DataType > & dataID )
{
#ifndef NDEBUG
  for ( auto it = primitives.begin(); it != primitives.end(); it++ )
  {
    WALBERLA_ASSERT_EQUAL( primitiveDataHandlers_ - 1, it->second->getNumberOfDataEntries() );
  }
#endif

  for ( auto it = primitives.begin(); it != primitives.end(); it++ )
  {
    it->second->addData( dataID, dataHandling );
  }
}


} // namespace hhg


