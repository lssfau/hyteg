
#pragma once

#include "core/logging/Logging.h"
#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/primitives/Primitive.hpp"
#include "tinyhhg_core/primitives/face.hpp"
#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"

#include <map>
#include <vector>

namespace hhg {

class SetupPrimitiveStorage;

class PrimitiveStorage : private walberla::NonCopyable
{
public:

  typedef std::map< PrimitiveID::IDType, Primitive* > PrimitiveMap;
  typedef std::map< PrimitiveID::IDType, Vertex* >    VertexMap;
  typedef std::map< PrimitiveID::IDType, Edge* >      EdgeMap;
  typedef std::map< PrimitiveID::IDType, Face* >      FaceMap;

  PrimitiveStorage( const uint_t & rank, const SetupPrimitiveStorage & setupStorage );

  uint_t getRank() const { return rank_; }

  void checkConsistency();

  //////////////////////////////
  // Primitive access methods //
  //////////////////////////////

  uint_t getNumberOfLocalVertices() const { return vertices_.size(); }
  uint_t getNumberOfLocalEdges()    const { return edges_.size(); }
  uint_t getNumberOfLocalFaces()    const { return faces_.size(); }

  /// Returns true, if the \ref Primitive that corresponds to the \ref PrimitiveID exists locally.
  bool primitiveExistsLocally( const PrimitiveID & id ) const { return vertexExistsLocally( id ) || edgeExistsLocally( id ) || faceExistsLocally( id ); }
  /// Returns true, if the \ref Vertex that corresponds to the \ref PrimitiveID exists locally.
  bool vertexExistsLocally( const PrimitiveID & id )    const { return vertices_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Edge that corresponds to the \ref PrimitiveID exists locally.
  bool edgeExistsLocally( const PrimitiveID & id )      const { return edges_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Face that corresponds to the \ref PrimitiveID exists locally.
  bool faceExistsLocally( const PrimitiveID & id )      const { return faces_.count( id.getID() ) > 0; }

  /// Returns true, if the \ref Primitive that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  bool primitiveExistsInNeighborhood( const PrimitiveID & id ) const { return vertexExistsInNeighborhood( id ) || edgeExistsInNeighborhood( id ) || faceExistsInNeighborhood( id ); }
  /// Returns true, if the \ref Vertex that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  bool vertexExistsInNeighborhood( const PrimitiveID & id )    const { return neighborVertices_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Edge that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  bool edgeExistsInNeighborhood( const PrimitiveID & id )      const { return neighborEdges_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Face that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  bool faceExistsInNeighborhood( const PrimitiveID & id )      const { return neighborFaces_.count( id.getID() ) > 0; }

  template< typename PrimitiveType >
  inline bool primitiveExistsLocallyGenerically( const PrimitiveID & id ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  template< typename PrimitiveType >
  inline bool primitiveExistsInNeighborhoodGenerically( const PrimitiveID & id ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  /// Returns the \ref Primitive that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL if the \ref Primitive does not exist locally.
  const Primitive* getPrimitive( const PrimitiveID & id ) const;
        Primitive* getPrimitive( const PrimitiveID & id );

  /// Returns the \ref Vertex that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL if the \ref Vertex does not exist locally.
  const Vertex* getVertex( const PrimitiveID & id ) const { return vertexExistsLocally( id ) ? vertices_.at( id.getID() ) : nullptr; }
        Vertex* getVertex( const PrimitiveID & id )       { return vertexExistsLocally( id ) ? vertices_[ id.getID() ] : nullptr; }

  /// Returns the \ref Edge that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL if the \ref Edge does not exist locally.
  const Edge* getEdge( const PrimitiveID & id ) const { return edgeExistsLocally( id ) ? edges_.at( id.getID() ) : nullptr; }
        Edge* getEdge( const PrimitiveID & id )       { return edgeExistsLocally( id ) ? edges_[ id.getID() ] : nullptr; }

  /// Returns the \ref Face that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL if the \ref Face does not exist locally.
  const Face* getFace( const PrimitiveID & id ) const { return faceExistsLocally( id ) ? faces_.at( id.getID() ) : nullptr; }
        Face* getFace( const PrimitiveID & id )       { return faceExistsLocally( id ) ? faces_[ id.getID() ] : nullptr; }

  /// Returns the neighbor \ref Primitive that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL if the \ref Primitive does not exist in the direct neighborhood.
  const Primitive* getNeighborPrimitive( const PrimitiveID & id ) const;
        Primitive* getNeighborPrimitive( const PrimitiveID & id );

  /// Returns the neighbor \ref Vertex that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL if the \ref Vertex does not exist in the direct neighborhood.
  const Vertex* getNeighborVertex( const PrimitiveID & id ) const { return vertexExistsInNeighborhood( id ) ? neighborVertices_.at( id.getID() ) : nullptr; }
        Vertex* getNeighborVertex( const PrimitiveID & id )       { return vertexExistsInNeighborhood( id ) ? neighborVertices_[ id.getID() ] : nullptr; }

  /// Returns the neighbor \ref Edge that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL if the \ref Edge does not exist in the direct neighborhood.
  const Edge* getNeighborEdge( const PrimitiveID & id ) const { return edgeExistsInNeighborhood( id ) ? neighborEdges_.at( id.getID() ) : nullptr; }
        Edge* getNeighborEdge( const PrimitiveID & id )       { return edgeExistsInNeighborhood( id ) ? neighborEdges_[ id.getID() ] : nullptr; }

  /// Returns the neighbor \ref Face that is assigned to the passed \ref PrimitiveID.
  /// Returns NULL if the \ref Face does not exist in the direct neighborhood.
  const Face* getNeighborFace( const PrimitiveID & id ) const { return faceExistsInNeighborhood( id ) ? neighborFaces_.at( id.getID() ) : nullptr; }
        Face* getNeighborFace( const PrimitiveID & id )       { return faceExistsInNeighborhood( id ) ? neighborFaces_[ id.getID() ] : nullptr; }


  /// Generic versions of the getter methods.
  template< typename PrimitiveType >
  inline const PrimitiveType* getPrimitiveGenerically( const PrimitiveID & id ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  template< typename PrimitiveType >
  inline       PrimitiveType* getPrimitiveGenerically( const PrimitiveID & id )       { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  /// Generic versions of the getter methods.
  template< typename PrimitiveType >
  inline const PrimitiveType* getNeighborPrimitiveGenerically( const PrimitiveID & id ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  template< typename PrimitiveType >
  inline       PrimitiveType* getNeighborPrimitiveGenerically( const PrimitiveID & id )       { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  void getVertexIDs ( std::vector< PrimitiveID > & vertexIDs ) const;
  void getEdgeIDs   ( std::vector< PrimitiveID > & edgeIDs )   const;
  void getFaceIDs   ( std::vector< PrimitiveID > & faceIDs )   const;

  template< typename PrimitiveType >
  inline void getPrimitiveIDsGenerically( std::vector< PrimitiveID > & primitiveIDs ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  /// Fills the passed map with all PrimitiveIDs and the respective pointers to the primitives
  void getPrimitives( PrimitiveMap & primitiveMap ) const;

  const VertexMap & getVertices() const { return vertices_; }
  const EdgeMap   & getEdges()    const { return edges_;    }
  const FaceMap   & getFaces()    const { return faces_;    }

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

  /// Returns the rank of the process the primitive is located on.
  /// Returns the local MPI rank if it is a local primitive.
  /// Returns the correct rank if the primitive lies in the direct neighborhood.
  /// Should not be called for other primitives.
  uint_t getPrimitiveRank        ( const PrimitiveID & id ) const;

  /// Returns the correct rank if the primitive lies in the direct neighborhood.
  /// Should not be called for other primitives.
  uint_t getNeighborPrimitiveRank( const PrimitiveID & id ) const { WALBERLA_ASSERT( primitiveExistsInNeighborhood( id ) ); return neighborRanks_.at( id.getID() ); }

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

  VertexMap neighborVertices_;
  EdgeMap   neighborEdges_;
  FaceMap   neighborFaces_;

  uint_t rank_;
  uint_t primitiveDataHandlers_;

  std::map< PrimitiveID::IDType, uint_t > neighborRanks_;

};

template<>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Primitive >( const PrimitiveID & id ) const { return primitiveExistsLocally( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Vertex >( const PrimitiveID & id ) const { return vertexExistsLocally( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Edge >  ( const PrimitiveID & id ) const { return edgeExistsLocally( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Face >  ( const PrimitiveID & id ) const { return faceExistsLocally( id ); }


template<>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Primitive >( const PrimitiveID & id ) const { return primitiveExistsInNeighborhood( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Vertex >( const PrimitiveID & id ) const { return vertexExistsInNeighborhood( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Edge >  ( const PrimitiveID & id ) const { return edgeExistsInNeighborhood( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Face >  ( const PrimitiveID & id ) const { return faceExistsInNeighborhood( id ); }


template<>
inline const Primitive* PrimitiveStorage::getPrimitiveGenerically< Primitive >( const PrimitiveID & id ) const { return getPrimitive( id ); }

template<>
inline       Primitive* PrimitiveStorage::getPrimitiveGenerically< Primitive >( const PrimitiveID & id )       { return getPrimitive( id ); }

template<>
inline const Vertex* PrimitiveStorage::getPrimitiveGenerically< Vertex >( const PrimitiveID & id ) const { return getVertex( id ); }

template<>
inline       Vertex* PrimitiveStorage::getPrimitiveGenerically< Vertex >( const PrimitiveID & id )       { return getVertex( id ); }

template<>
inline const Edge*   PrimitiveStorage::getPrimitiveGenerically< Edge >  ( const PrimitiveID & id ) const { return getEdge( id ); }

template<>
inline       Edge*   PrimitiveStorage::getPrimitiveGenerically< Edge >  ( const PrimitiveID & id )       { return getEdge( id ); }

template<>
inline const Face*   PrimitiveStorage::getPrimitiveGenerically< Face >  ( const PrimitiveID & id ) const { return getFace( id ); }

template<>
inline       Face*   PrimitiveStorage::getPrimitiveGenerically< Face >  ( const PrimitiveID & id )       { return getFace( id ); }


template<>
inline const Primitive* PrimitiveStorage::getNeighborPrimitiveGenerically< Primitive >( const PrimitiveID & id ) const { return getNeighborPrimitive( id ); }

template<>
inline       Primitive* PrimitiveStorage::getNeighborPrimitiveGenerically< Primitive >( const PrimitiveID & id )       { return getNeighborPrimitive( id ); }

template<>
inline const Vertex* PrimitiveStorage::getNeighborPrimitiveGenerically< Vertex >( const PrimitiveID & id ) const { return getNeighborVertex( id ); }

template<>
inline       Vertex* PrimitiveStorage::getNeighborPrimitiveGenerically< Vertex >( const PrimitiveID & id )       { return getNeighborVertex( id ); }

template<>
inline const Edge*   PrimitiveStorage::getNeighborPrimitiveGenerically< Edge >  ( const PrimitiveID & id ) const { return getNeighborEdge( id ); }

template<>
inline       Edge*   PrimitiveStorage::getNeighborPrimitiveGenerically< Edge >  ( const PrimitiveID & id )       { return getNeighborEdge( id ); }

template<>
inline const Face*   PrimitiveStorage::getNeighborPrimitiveGenerically< Face >  ( const PrimitiveID & id ) const { return getNeighborFace( id ); }

template<>
inline       Face*   PrimitiveStorage::getNeighborPrimitiveGenerically< Face >  ( const PrimitiveID & id )       { return getNeighborFace( id ); }


template<>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Vertex >( std::vector< PrimitiveID > & primitiveIDs ) const { getVertexIDs( primitiveIDs ); }

template<>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Edge >( std::vector< PrimitiveID > & primitiveIDs ) const { getEdgeIDs( primitiveIDs ); }

template<>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Face >( std::vector< PrimitiveID > & primitiveIDs ) const { getFaceIDs( primitiveIDs ); }


template< typename DataType >
PrimitiveDataID< DataType, Primitive > PrimitiveStorage::addPrimitiveData( const PrimitiveDataHandling< DataType, Primitive > & dataHandling,
						                const std::string & identifier )
{
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
  PrimitiveDataID< DataType, Vertex > dataID = generateDataID< DataType, Vertex >();
  addPrimitiveData( dataHandling, identifier, vertices_, dataID );
  return dataID;
}


template< typename DataType >
PrimitiveDataID< DataType, Edge > PrimitiveStorage::addEdgeData( const PrimitiveDataHandling< DataType, Edge > & dataHandling,
						           const std::string & identifier )
{
  PrimitiveDataID< DataType, Edge > dataID = generateDataID< DataType, Edge >();
  addPrimitiveData( dataHandling, identifier, edges_, dataID );
  return dataID;
}


template< typename DataType >
PrimitiveDataID< DataType, Face > PrimitiveStorage::addFaceData( const PrimitiveDataHandling< DataType, Face > & dataHandling,
						             const std::string & identifier )
{
  PrimitiveDataID< DataType, Face > dataID = generateDataID< DataType, Face >();
  addPrimitiveData( dataHandling, identifier, faces_, dataID );
  return dataID;
}


template< typename DataType, typename PrimitiveType >
PrimitiveDataID< DataType, PrimitiveType > PrimitiveStorage::generateDataID()
{
#ifndef NDEBUG
  checkConsistency();
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
    WALBERLA_ASSERT_GREATER( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
  }
#endif

  for ( auto it = primitives.begin(); it != primitives.end(); it++ )
  {
    it->second->addData( dataID, dataHandling );
  }
}


} // namespace hhg


