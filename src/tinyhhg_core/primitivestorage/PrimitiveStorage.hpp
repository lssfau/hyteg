
#pragma once

#include "core/logging/Logging.h"
#include "tinyhhg_core/primitiveid.hpp"
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

  typedef std::map< PrimitiveID::IDType, Vertex* > VertexMap;
  typedef std::map< PrimitiveID::IDType, Edge* >   EdgeMap;
  typedef std::map< PrimitiveID::IDType, Face* >   FaceMap;

  PrimitiveStorage( const SetupPrimitiveStorage & setupStorage ) :
    primitiveDataHandlers_( 0 )
  {
    std::map< PrimitiveID::IDType, SetupVertex* > setupVertices;
    std::map< PrimitiveID::IDType, SetupEdge* >   setupEdges;
    std::map< PrimitiveID::IDType, SetupFace* >   setupFaces;

    setupStorage.getVertices( setupVertices );
    setupStorage.getEdges( setupEdges );
    setupStorage.getFaces( setupFaces );

    uint_t processRank = walberla::mpi::MPIManager::instance()->rank();

    for ( auto it = setupVertices.begin(); it != setupVertices.end(); it++  )
    {
      if ( processRank == it->second->getTargetRank() )
      {
	vertices_[ it->first ] = new Vertex( it->first, it->second->getCoordinates() );
      }
    }

    for ( auto it = setupEdges.begin(); it != setupEdges.end(); it++ )
    {
      if ( processRank == it->second->getTargetRank() )
      {
	PrimitiveID edgeID = it->first;
	DoFType edgeType = it->second->getDoFType();
	Vertex* edgeVertex0 = vertices_[ it->second->getVertexID0().getID() ];
	Vertex* edgeVertex1 = vertices_[ it->second->getVertexID1().getID() ];
	edges_[ edgeID.getID() ] = new Edge( edgeID.getID(), edgeType, edgeVertex0, edgeVertex1 );
      }
    }

    for ( auto it = setupFaces.begin(); it != setupFaces.end(); it++ )
    {
      if ( processRank == it->second->getTargetRank() )
      {
	PrimitiveID faceID = it->first;
	Edge* faceEdges[3];
	faceEdges[0] = edges_[ it->second->getEdgeID0().getID() ];
	faceEdges[1] = edges_[ it->second->getEdgeID1().getID() ];
	faceEdges[2] = edges_[ it->second->getEdgeID2().getID() ];
	faces_[ faceID.getID() ] = new Face( faceID.getID(), faceEdges );
      }
    }
  }

  PrimitiveID addVertex();

  ////////////////////
  // Access methods //
  ////////////////////

  /// Returns true, if the \ref Primitive that corresponds to the \ref PrimitiveID exists locally.
  bool primitiveExistsLocally( const PrimitiveID & id ) const { return vertexExistsLocally( id ) || edgeExistsLocally( id ) || faceExistsLocally( id ); }
  /// Returns true, if the \ref Vertex that corresponds to the \ref PrimitiveID exists locally.
  bool vertexExistsLocally( const PrimitiveID & id )    const { return vertices_.find( id.getID() ) != vertices_.end(); }
  /// Returns true, if the \ref Edge that corresponds to the \ref PrimitiveID exists locally.
  bool edgeExistsLocally( const PrimitiveID & id )      const { return edges_.find( id.getID() ) != edges_.end(); }
  /// Returns true, if the \ref Face that corresponds to the \ref PrimitiveID exists locally.
  bool faceExistsLocally( const PrimitiveID & id )      const { return faces_.find( id.getID() ) != faces_.end(); }

  /// Returns the \ref Vertex that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL otherwise.
  const Vertex* getVertex( const PrimitiveID & id ) const { return vertexExistsLocally( id ) ? vertices_.at( id.getID() ) : NULL; };
        Vertex* getVertex( const PrimitiveID & id )       { return vertexExistsLocally( id ) ? vertices_[ id.getID() ] : NULL; };

  /// Returns the \ref Edge that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL otherwise.
  const Edge* getEdge( const PrimitiveID & id ) const { return edgeExistsLocally( id ) ? edges_.at( id.getID() ) : NULL; };
	Edge* getEdge( const PrimitiveID & id )       { return edgeExistsLocally( id ) ? edges_[ id.getID() ] : NULL; };

  /// Returns the \ref Edge that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL otherwise.
  const Face* getFace( const PrimitiveID & id ) const { return faceExistsLocally( id ) ? faces_.at( id.getID() ) : NULL; };
	Face* getFace( const PrimitiveID & id )       { return faceExistsLocally( id ) ? faces_[ id.getID() ] : NULL; };

  VertexMap::iterator beginVertices() { return vertices_.begin(); }
  VertexMap::iterator endVertices()   { return vertices_.end(); }

  EdgeMap::iterator beginEdges()      { return edges_.begin(); }
  EdgeMap::iterator endEdges()        { return edges_.end(); }

  FaceMap::iterator beginFaces()      { return faces_.begin(); }
  FaceMap::iterator endFaces()        { return faces_.end(); }

  VertexMap::const_iterator beginVertices() const { return vertices_.begin(); }
  VertexMap::const_iterator endVertices()   const { return vertices_.end(); }

  EdgeMap::const_iterator beginEdges()      const { return edges_.begin(); }
  EdgeMap::const_iterator endEdges()        const { return edges_.end(); }

  FaceMap::const_iterator beginFaces()      const { return faces_.begin(); }
  FaceMap::const_iterator endFaces()        const { return faces_.end(); }


  ////////////////////////////
  // Primitive data methods //
  ////////////////////////////

  template< typename DataType >
  inline PrimitiveDataID< DataType, Primitive > addPrimitiveData( const PrimitiveDataHandling< DataType, Primitive > & dataHandling,
  						                  const std::string & identifier );

  template< typename DataType >
  inline PrimitiveDataID< DataType, Vertex > addVertexData( const PrimitiveDataHandling< DataType, Vertex > & dataHandling,
						            const std::string & identifier );

  template< typename DataType >
  inline PrimitiveDataID< DataType, Edge > addEdgeData( const PrimitiveDataHandling< DataType, Edge > & dataHandling,
  						        const std::string & identifier );

  template< typename DataType >
  inline PrimitiveDataID< DataType, Face > addFaceData( const PrimitiveDataHandling< DataType, Face > & dataHandling,
  						        const std::string & identifier );

private:

  template< typename DataType, typename PrimitiveType >
  inline PrimitiveDataID< DataType, PrimitiveType > generateDataID();

  template< typename DataType,
            typename PrimitiveType,
	    typename = typename std::enable_if< std::is_base_of< Primitive, PrimitiveType >::value >::type >
  inline void addPrimitiveData( const PrimitiveDataHandling< DataType, PrimitiveType > & dataHandling,
				const std::string & identifier,
				const std::map< PrimitiveID::IDType, PrimitiveType* > & primitives,
				const PrimitiveDataID< DataType, PrimitiveType > & dataID );

  VertexMap vertices_;
  EdgeMap   edges_;
  FaceMap   faces_;

  uint_t primitiveDataHandlers_;

};


template< typename DataType >
PrimitiveDataID< DataType, Primitive > PrimitiveStorage::addPrimitiveData( const PrimitiveDataHandling< DataType, Primitive > & dataHandling,
						                const std::string & identifier )
{
  WALBERLA_LOG_PROGRESS( "Adding data to all primitives (\"" << identifier << "\")" );
  PrimitiveDataID< DataType, Primitive > dataID = generateDataID< DataType, Primitive >();
  std::map< PrimitiveID::IDType, Primitive* > primitives;
  primitives.insert( vertices_.begin(), vertices_.end() );
  primitives.insert( edges_.begin(), edges_.end() );
  primitives.insert( faces_.begin(), faces_.end() );
  addPrimitiveData( dataHandling, identifier, primitives, dataID );
  return dataID;

}


template< typename DataType >
PrimitiveDataID< DataType, Vertex > PrimitiveStorage::addVertexData( const PrimitiveDataHandling< DataType, Vertex > & dataHandling,
						             const std::string & identifier )
{
  WALBERLA_LOG_PROGRESS( "Adding data to vertices (\"" << identifier << "\")" );
  PrimitiveDataID< DataType, Vertex > dataID = generateDataID< DataType, Vertex >();
  addPrimitiveData( dataHandling, identifier, vertices_, dataID );
  return dataID;
}


template< typename DataType >
PrimitiveDataID< DataType, Edge > PrimitiveStorage::addEdgeData( const PrimitiveDataHandling< DataType, Edge > & dataHandling,
						           const std::string & identifier )
{
  WALBERLA_LOG_PROGRESS( "Adding data to edges (\"" << identifier << "\")" );
  PrimitiveDataID< DataType, Edge > dataID = generateDataID< DataType, Edge >();
  addPrimitiveData( dataHandling, identifier, edges_, dataID );
  return dataID;
}


template< typename DataType >
PrimitiveDataID< DataType, Face > PrimitiveStorage::addFaceData( const PrimitiveDataHandling< DataType, Face > & dataHandling,
						             const std::string & identifier )
{
  WALBERLA_LOG_PROGRESS( "Adding data to faces (\"" << identifier << "\")" );
  PrimitiveDataID< DataType, Face > dataID = generateDataID< DataType, Face >();
  addPrimitiveData( dataHandling, identifier, faces_, dataID );
  return dataID;
}


template< typename DataType, typename PrimitiveType >
PrimitiveDataID< DataType, PrimitiveType > PrimitiveStorage::generateDataID()
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

  return PrimitiveDataID< DataType, PrimitiveType >( primitiveDataHandlers_++ );
}


template< typename DataType,
          typename PrimitiveType,
	  typename >
void PrimitiveStorage::addPrimitiveData( const PrimitiveDataHandling< DataType, PrimitiveType > & dataHandling,
					 const std::string & identifier,
					 const std::map< PrimitiveID::IDType, PrimitiveType* > & primitives,
					 const PrimitiveDataID< DataType, PrimitiveType > & dataID )
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


