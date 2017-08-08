
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

  typedef std::map< PrimitiveID::IDType, std::shared_ptr< Primitive > > PrimitiveMap;
  typedef std::map< PrimitiveID::IDType, std::shared_ptr< Vertex > >    VertexMap;
  typedef std::map< PrimitiveID::IDType, std::shared_ptr< Edge > >      EdgeMap;
  typedef std::map< PrimitiveID::IDType, std::shared_ptr< Face > >      FaceMap;

  PrimitiveStorage( const SetupPrimitiveStorage & setupStorage );

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
  /// The returned \ref Primitive is either local or lies in the direct neighborhood.
  /// Returns nullptr if the \ref Primitive does not exist locally nor in the direct neighborhood.
  const Primitive* getPrimitive( const PrimitiveID & id ) const;
        Primitive* getPrimitive( const PrimitiveID & id );

  /// Returns the \ref Vertex that is assigned to the passed \ref PrimitiveID.
  /// The returned \ref Vertex is either local or lies in the direct neighborhood.
  /// Returns nullptr if the \ref Vertex does not exist locally nor in the direct neighborhood.
  const Vertex* getVertex( const PrimitiveID & id ) const;
        Vertex* getVertex( const PrimitiveID & id );

  /// Returns the \ref Edge that is assigned to the passed \ref PrimitiveID.
  /// The returned \ref Edge is either local or lies in the direct neighborhood.
  /// Returns nullptr if the \ref Edge does not exist locally nor in the direct neighborhood.
  const Edge* getEdge( const PrimitiveID & id ) const;
        Edge* getEdge( const PrimitiveID & id );

  /// Returns the \ref Face that is assigned to the passed \ref PrimitiveID.
  /// The returned \ref Face is either local or lies in the direct neighborhood.
  /// Returns nullptr if the \ref Face does not exist locally nor in the direct neighborhood.
  const Face* getFace( const PrimitiveID & id ) const;
        Face* getFace( const PrimitiveID & id );

  /// Generic versions of the getter methods.
  template< typename PrimitiveType >
  inline const PrimitiveType* getPrimitiveGenerically( const PrimitiveID & id ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  template< typename PrimitiveType >
  inline       PrimitiveType* getPrimitiveGenerically( const PrimitiveID & id )       { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  /// Fills the passed vector with the IDs of the locally existing vertices
  void getVertexIDs ( std::vector< PrimitiveID > & vertexIDs ) const;

  /// Fills the passed vector with the IDs of the locally existing edges
  void getEdgeIDs   ( std::vector< PrimitiveID > & edgeIDs )   const;

  /// Fills the passed vector with the IDs of the locally existing faces
  void getFaceIDs   ( std::vector< PrimitiveID > & faceIDs )   const;

  template< typename PrimitiveType >
  inline void getPrimitiveIDsGenerically( std::vector< PrimitiveID > & primitiveIDs ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  /// Fills the passed map with all PrimitiveIDs and the respective pointers to the primitives
  void getPrimitives( PrimitiveMap & primitiveMap ) const;

  /// Returns a reference to a map of the locally existing \ref Vertex instances
  const VertexMap & getVertices() const { return vertices_; }

  /// Returns a reference to a map of the locally existing \ref Edge instances
  const EdgeMap   & getEdges()    const { return edges_;    }

  /// Returns a reference to a map of the locally existing \ref Face instances
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

  template< typename DataType, typename DataHandlingType >
  inline void addPrimitiveData(       PrimitiveDataID< DataType, Primitive > & dataID,
                                const std::shared_ptr< DataHandlingType > & dataHandling,
  						                  const std::string & identifier );

  template< typename DataType, typename DataHandlingType >
  inline void addVertexData(       PrimitiveDataID< DataType, Vertex > & dataID,
                             const std::shared_ptr< DataHandlingType > & dataHandling,
						                 const std::string & identifier );

  template< typename DataType, typename DataHandlingType >
  inline void addEdgeData(       PrimitiveDataID< DataType, Edge > & dataID,
                           const std::shared_ptr< DataHandlingType > & dataHandling,
  						                                          const std::string & identifier );

  template< typename DataType, typename DataHandlingType >
  inline void addFaceData(       PrimitiveDataID< DataType, Face > & dataID,
                           const std::shared_ptr< DataHandlingType > & dataHandling,
  						             const std::string & identifier );

private:

  template< typename DataType, typename PrimitiveType >
  inline PrimitiveDataID< DataType, PrimitiveType > generateDataID();

  template< typename DataType,
            typename PrimitiveType,
            typename DataHandlingType,
            typename = typename std::enable_if< std::is_base_of< Primitive, PrimitiveType >::value >::type >
  inline void addPrimitiveData( const std::shared_ptr< DataHandlingType > & dataHandling,
				                        const std::string & identifier,
				                        const std::map< PrimitiveID::IDType, std::shared_ptr< PrimitiveType > > & primitives,
				                        const PrimitiveDataID< DataType, PrimitiveType > & dataID );

  VertexMap vertices_;
  EdgeMap   edges_;
  FaceMap   faces_;

  VertexMap neighborVertices_;
  EdgeMap   neighborEdges_;
  FaceMap   neighborFaces_;

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
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Vertex >( std::vector< PrimitiveID > & primitiveIDs ) const { getVertexIDs( primitiveIDs ); }

template<>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Edge >( std::vector< PrimitiveID > & primitiveIDs ) const { getEdgeIDs( primitiveIDs ); }

template<>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Face >( std::vector< PrimitiveID > & primitiveIDs ) const { getFaceIDs( primitiveIDs ); }


template< typename DataType,
          typename DataHandlingType >
void PrimitiveStorage::addPrimitiveData(       PrimitiveDataID< DataType, Primitive > & dataID,
                                         const std::shared_ptr< DataHandlingType > & dataHandling,
                                         const std::string & identifier )
{
  dataID = generateDataID< DataType, Primitive >();
  PrimitiveMap primitives;
  primitives.insert( vertices_.begin(), vertices_.end() );
  primitives.insert( edges_.begin(), edges_.end() );
  primitives.insert( faces_.begin(), faces_.end() );
  addPrimitiveData( dataHandling, identifier, primitives, dataID );
}


template< typename DataType,
          typename DataHandlingType >
void PrimitiveStorage::addVertexData(       PrimitiveDataID< DataType, Vertex > & dataID,
                                      const std::shared_ptr< DataHandlingType > & dataHandling,
                                      const std::string & identifier )
{
  dataID = generateDataID< DataType, Vertex >();
  addPrimitiveData( dataHandling, identifier, vertices_, dataID );
}


template< typename DataType,
          typename DataHandlingType >
void PrimitiveStorage::addEdgeData(       PrimitiveDataID< DataType, Edge > & dataID,
                                    const std::shared_ptr< DataHandlingType > & dataHandling,
                                    const std::string & identifier )
{
  dataID = generateDataID< DataType, Edge >();
  addPrimitiveData( dataHandling, identifier, edges_, dataID );
}


template< typename DataType,
          typename DataHandlingType >
void PrimitiveStorage::addFaceData(       PrimitiveDataID< DataType, Face > & dataID,
                                    const std::shared_ptr< DataHandlingType > & dataHandling,
                                    const std::string & identifier )
{
  dataID = generateDataID< DataType, Face >();
  addPrimitiveData( dataHandling, identifier, faces_, dataID );
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
          typename DataHandlingType,
          typename >
void PrimitiveStorage::addPrimitiveData( const std::shared_ptr< DataHandlingType > & dataHandling,
                                         const std::string & identifier,
                                         const std::map< PrimitiveID::IDType, std::shared_ptr< PrimitiveType > > & primitives,
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


