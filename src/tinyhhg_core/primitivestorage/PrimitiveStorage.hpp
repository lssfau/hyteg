
#pragma once

#include "core/logging/Logging.h"
#include "core/timing/TimingTree.h"
#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/primitives/Primitive.hpp"
#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
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
  typedef std::map< PrimitiveID::IDType, std::shared_ptr< Cell > >      CellMap;

  PrimitiveStorage( const SetupPrimitiveStorage & setupStorage );

  void checkConsistency();

  /// @name \ref Primitive access methods
  /// Various methods to obtain primitives or IDs.
  ///@{
  uint_t getNumberOfLocalPrimitives() const { return getNumberOfLocalVertices() + getNumberOfLocalEdges() + getNumberOfLocalFaces() + getNumberOfLocalCells(); }
  uint_t getNumberOfLocalVertices() const { return vertices_.size(); }
  uint_t getNumberOfLocalEdges()    const { return edges_.size(); }
  uint_t getNumberOfLocalFaces()    const { return faces_.size(); }
  uint_t getNumberOfLocalCells()    const { return cells_.size(); }

  /// Returns true, if the \ref Primitive that corresponds to the \ref PrimitiveID exists locally.
  bool primitiveExistsLocally( const PrimitiveID & id ) const { return vertexExistsLocally( id ) || edgeExistsLocally( id ) || faceExistsLocally( id ) || cellExistsLocally( id ); }
  /// Returns true, if the \ref Vertex that corresponds to the \ref PrimitiveID exists locally.
  bool vertexExistsLocally( const PrimitiveID & id )    const { return vertices_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Edge that corresponds to the \ref PrimitiveID exists locally.
  bool edgeExistsLocally( const PrimitiveID & id )      const { return edges_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Face that corresponds to the \ref PrimitiveID exists locally.
  bool faceExistsLocally( const PrimitiveID & id )      const { return faces_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Cell that corresponds to the \ref PrimitiveID exists locally.
  bool cellExistsLocally( const PrimitiveID & id )      const { return cells_.count( id.getID() ) > 0; }

  /// Returns true, if the \ref Primitive that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  bool primitiveExistsInNeighborhood( const PrimitiveID & id ) const { return vertexExistsInNeighborhood( id ) || edgeExistsInNeighborhood( id ) || faceExistsInNeighborhood( id ) || cellExistsInNeighborhood( id ); }
  /// Returns true, if the \ref Vertex that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  bool vertexExistsInNeighborhood( const PrimitiveID & id )    const { return neighborVertices_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Edge that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  bool edgeExistsInNeighborhood( const PrimitiveID & id )      const { return neighborEdges_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Face that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  bool faceExistsInNeighborhood( const PrimitiveID & id )      const { return neighborFaces_.count( id.getID() ) > 0; }
  /// Returns true, if the \ref Cell that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  bool cellExistsInNeighborhood( const PrimitiveID & id )      const { return neighborCells_.count( id.getID() ) > 0; }

  /// Returns true, if the \ref Primitive of the generically passed type that corresponds to the \ref PrimitiveID exists locally.
  template< typename PrimitiveType >
  inline bool primitiveExistsLocallyGenerically( const PrimitiveID & id ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  /// Returns true, if the \ref Primitive of the generically passed type that corresponds to the \ref PrimitiveID exists in the direct neighborhood.
  template< typename PrimitiveType >
  inline bool primitiveExistsInNeighborhoodGenerically( const PrimitiveID & id ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  /// Returns the \ref Primitive that is assigned to the passed \ref PrimitiveID.
  /// The returned \ref Primitive is either local or lies in the direct neighborhood.
  /// Returns nullptr if the \ref Primitive does not exist locally nor in the direct neighborhood.
  ///@{
  const Primitive* getPrimitive( const PrimitiveID & id ) const;
        Primitive* getPrimitive( const PrimitiveID & id );
  ///@}

  /// Returns the \ref Vertex that is assigned to the passed \ref PrimitiveID.
  /// The returned \ref Vertex is either local or lies in the direct neighborhood.
  /// Returns nullptr if the \ref Vertex does not exist locally nor in the direct neighborhood.
  ///@{
  const Vertex* getVertex( const PrimitiveID & id ) const;
        Vertex* getVertex( const PrimitiveID & id );
  ///@}

  /// Returns the \ref Edge that is assigned to the passed \ref PrimitiveID.
  /// The returned \ref Edge is either local or lies in the direct neighborhood.
  /// Returns nullptr if the \ref Edge does not exist locally nor in the direct neighborhood.
  ///@{
  const Edge* getEdge( const PrimitiveID & id ) const;
        Edge* getEdge( const PrimitiveID & id );
  ///@}

  /// Returns the \ref Face that is assigned to the passed \ref PrimitiveID.
  /// The returned \ref Face is either local or lies in the direct neighborhood.
  /// Returns nullptr if the \ref Face does not exist locally nor in the direct neighborhood.
  ///@{
  const Face* getFace( const PrimitiveID & id ) const;
        Face* getFace( const PrimitiveID & id );
  ///@}

  /// Returns the \ref Cell that is assigned to the passed \ref PrimitiveID.
  /// The returned \ref Cell is either local or lies in the direct neighborhood.
  /// Returns nullptr if the \ref Cell does not exist locally nor in the direct neighborhood.
  ///@{
  const Cell* getCell( const PrimitiveID & id ) const;
        Cell* getCell( const PrimitiveID & id );
  ///@}

  /// @name Generic versions of the getter methods.
  ///@{
  template< typename PrimitiveType >
  inline const PrimitiveType* getPrimitiveGenerically( const PrimitiveID & id ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  template< typename PrimitiveType >
  inline       PrimitiveType* getPrimitiveGenerically( const PrimitiveID & id )       { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }
  ///@}

  /// Fills the passed vector with the IDs of the locally existing primitives
  void getPrimitiveIDs ( std::vector< PrimitiveID > & primitiveIDs ) const;

  /// Fills the passed vector with the IDs of the locally existing vertices
  void getVertexIDs ( std::vector< PrimitiveID > & vertexIDs ) const;

  /// Fills the passed vector with the IDs of the locally existing edges
  void getEdgeIDs   ( std::vector< PrimitiveID > & edgeIDs )   const;

  /// Fills the passed vector with the IDs of the locally existing faces
  void getFaceIDs   ( std::vector< PrimitiveID > & faceIDs )   const;

  /// Fills the passed vector with the IDs of the locally existing cells
  void getCellIDs   ( std::vector< PrimitiveID > & cellIDs )   const;

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

  /// Returns a reference to a map of the locally existing \ref Face instances
  const CellMap   & getCells()    const { return cells_;    }

  ///@}

  /// Returns the rank of the process the primitive is located on.
  /// Returns the local MPI rank if it is a local primitive.
  /// Returns the correct rank if the primitive lies in the direct neighborhood.
  /// Should not be called for other primitives.
  uint_t getPrimitiveRank        ( const PrimitiveID & id ) const;

  /// Returns the correct rank if the primitive lies in the direct neighborhood.
  /// Should not be called for other primitives.
  uint_t getNeighborPrimitiveRank( const PrimitiveID & id ) const { WALBERLA_ASSERT( primitiveExistsInNeighborhood( id ) ); return neighborRanks_.at( id.getID() ); }


  /// @name Primitive data methods
  /// Use these methods to add data to all primitives of a certain type using a respective \ref PrimitiveDataHandling implementation.
  /// \param dataID (out) the method creates a data ID and writes it to this parameter, the data can be obtained through a \ref Primitive using this ID
  /// \param dataHandling a pointer to the \ref PrimitiveDataHandling that shall be used to treat the data item
  /// \param identifier string that identifies the data that was added
  ///@{
  template< typename DataType, typename DataHandlingType >
  inline void addPrimitiveData(       PrimitiveDataID< DataType, Primitive > & dataID,
                                const std::shared_ptr< DataHandlingType >    & dataHandling,
  						                  const std::string                            & identifier );

  template< typename DataType, typename DataHandlingType >
  inline void addVertexData(       PrimitiveDataID< DataType, Vertex > & dataID,
                             const std::shared_ptr< DataHandlingType > & dataHandling,
						                 const std::string                         & identifier );

  template< typename DataType, typename DataHandlingType >
  inline void addEdgeData(       PrimitiveDataID< DataType, Edge >   & dataID,
                           const std::shared_ptr< DataHandlingType > & dataHandling,
                           const std::string                         & identifier );

  template< typename DataType, typename DataHandlingType >
  inline void addFaceData(       PrimitiveDataID< DataType, Face >   & dataID,
                           const std::shared_ptr< DataHandlingType > & dataHandling,
  						             const std::string                         & identifier );

  template< typename DataType, typename DataHandlingType >
  inline void addCellData(       PrimitiveDataID< DataType, Cell >   & dataID,
                           const std::shared_ptr< DataHandlingType > & dataHandling,
                           const std::string                         & identifier );
  ///@}

  /// Migrates the passed (local!) primitives to the respective target process.
  /// Must be called collectively, even if a processes does not send any primitives (pass empty map).
  /// Calls the serialization and deserialization methods of the data handling instances of all registered data items
  /// in order to transport the data over MPI.
  /// Automatically refreshes the neighborhood information.
  /// \param primitivesToMigrate key: primitive to migrate, value: target process
  void migratePrimitives( const std::map< PrimitiveID::IDType, uint_t > & primitivesToMigrate );

  /// Returns a stamp that is always increased when the topology of the storage somehow changes -
  /// e.g. after migration of primitives to other processes.
  uint_t getModificationStamp() const { return modificationStamp_; }

  /// Fills the passed set with all neighboring ranks (== all ranks from primitives that are located in the direct neighborhood)
  void getNeighboringRanks( std::set< uint_t >                 & neighboringRanks ) const;
  void getNeighboringRanks( std::set< walberla::mpi::MPIRank > & neighboringRanks ) const;

  inline void enableGlobalTiming(){
    if(!timingTree_) {
      timingTree_ = std::make_shared<walberla::WcTimingTree>();
    }
  }

  inline std::shared_ptr< walberla::WcTimingTree >& getTimingTree(){
    return timingTree_;
  }

private:

  // needed to differentiate when migrating primitives
  enum PrimitiveTypeEnum
  {
    VERTEX,
    EDGE,
    FACE,
    CELL,
    INVALID
  };

  /// Returns the primitive type of a local primitive.
  /// Returns invalid if the primitive is not locally availably.
  PrimitiveTypeEnum getPrimitiveType( const PrimitiveID & primitiveID ) const;

  /// Deserializes a primitive from the receive buffer and adds it to the storage.
  /// Reads:
  ///   - PrimitiveType
  ///   - Primitive
  /// \param recvBuffer the buffer the data shall be unpacked from
  /// \param isNeighborPrimitive if true, the primitive is added to the neighborhood instead of the local primitives
  /// \return the ID of the deserialized primitive
  PrimitiveID deserializeAndAddPrimitive( walberla::mpi::RecvBuffer & recvBuffer, const bool & isNeighborPrimitive );

  /// Serializes all data from a locally allocated primitive to the send buffer.
  void serializeAllPrimitiveData               ( walberla::mpi::SendBuffer & sendBuffer, const PrimitiveID & primitiveID );
  /// Initializes and then deserializes all data from the receive buffer to a locally allocated primitive.
  void initializeAndDeserializeAllPrimitiveData( walberla::mpi::RecvBuffer & recvBuffer, const PrimitiveID & primitiveID );


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
  CellMap   cells_;

  VertexMap neighborVertices_;
  EdgeMap   neighborEdges_;
  FaceMap   neighborFaces_;
  CellMap   neighborCells_;

  template< typename DataType >
  inline void addDataHandlingCallbacks( const PrimitiveDataID< DataType, Primitive > &                                                     dataID,
                                        const std::function< void( const std::shared_ptr< Primitive > & ) > &                              initializationFunction,
                                        const std::function< void( const std::shared_ptr< Primitive > &, walberla::mpi::SendBuffer & ) > & serializationFunction,
                                        const std::function< void( const std::shared_ptr< Primitive > &, walberla::mpi::RecvBuffer & ) > & deserializationFunction )
  {
    primitiveDataInitializationFunctions_[ dataID ]  = initializationFunction;
    primitiveDataSerializationFunctions_[ dataID ]   = serializationFunction;
    primitiveDataDeserializationFunctions_[ dataID ] = deserializationFunction;
  }

  template< typename DataType >
  inline void addDataHandlingCallbacks( const PrimitiveDataID< DataType, Vertex > &                                                     dataID,
                                        const std::function< void( const std::shared_ptr< Vertex > & ) > &                              initializationFunction,
                                        const std::function< void( const std::shared_ptr< Vertex > &, walberla::mpi::SendBuffer & ) > & serializationFunction,
                                        const std::function< void( const std::shared_ptr< Vertex > &, walberla::mpi::RecvBuffer & ) > & deserializationFunction )
  {
    vertexDataInitializationFunctions_[ dataID ]  = initializationFunction;
    vertexDataSerializationFunctions_[ dataID ]   = serializationFunction;
    vertexDataDeserializationFunctions_[ dataID ] = deserializationFunction;
  }

  template< typename DataType >
  inline void addDataHandlingCallbacks( const PrimitiveDataID< DataType, Edge > &                                                     dataID,
                                        const std::function< void( const std::shared_ptr< Edge > & ) > &                              initializationFunction,
                                        const std::function< void( const std::shared_ptr< Edge > &, walberla::mpi::SendBuffer & ) > & serializationFunction,
                                        const std::function< void( const std::shared_ptr< Edge > &, walberla::mpi::RecvBuffer & ) > & deserializationFunction )
  {
    edgeDataInitializationFunctions_[ dataID ]  = initializationFunction;
    edgeDataSerializationFunctions_[ dataID ]   = serializationFunction;
    edgeDataDeserializationFunctions_[ dataID ] = deserializationFunction;
  }

  template< typename DataType >
  inline void addDataHandlingCallbacks( const PrimitiveDataID< DataType, Face > &                                                     dataID,
                                        const std::function< void( const std::shared_ptr< Face > & ) > &                              initializationFunction,
                                        const std::function< void( const std::shared_ptr< Face > &, walberla::mpi::SendBuffer & ) > & serializationFunction,
                                        const std::function< void( const std::shared_ptr< Face > &, walberla::mpi::RecvBuffer & ) > & deserializationFunction )
  {
    faceDataInitializationFunctions_[ dataID ]  = initializationFunction;
    faceDataSerializationFunctions_[ dataID ]   = serializationFunction;
    faceDataDeserializationFunctions_[ dataID ] = deserializationFunction;
  }

  template< typename DataType >
  inline void addDataHandlingCallbacks( const PrimitiveDataID< DataType, Cell > &                                                     dataID,
                                        const std::function< void( const std::shared_ptr< Cell > & ) > &                              initializationFunction,
                                        const std::function< void( const std::shared_ptr< Cell > &, walberla::mpi::SendBuffer & ) > & serializationFunction,
                                        const std::function< void( const std::shared_ptr< Cell > &, walberla::mpi::RecvBuffer & ) > & deserializationFunction )
  {
    cellDataInitializationFunctions_[ dataID ]  = initializationFunction;
    cellDataSerializationFunctions_[ dataID ]   = serializationFunction;
    cellDataDeserializationFunctions_[ dataID ] = deserializationFunction;
  }

  // Maps from data ID to respective callback functions

  std::map< uint_t, std::function< void( const std::shared_ptr< Primitive > & ) > >                              primitiveDataInitializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Primitive > &, walberla::mpi::SendBuffer & ) > > primitiveDataSerializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Primitive > &, walberla::mpi::RecvBuffer & ) > > primitiveDataDeserializationFunctions_;

  std::map< uint_t, std::function< void( const std::shared_ptr< Vertex > & ) > >                              vertexDataInitializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Vertex > &, walberla::mpi::SendBuffer & ) > > vertexDataSerializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Vertex > &, walberla::mpi::RecvBuffer & ) > > vertexDataDeserializationFunctions_;

  std::map< uint_t, std::function< void( const std::shared_ptr< Edge   > & ) > >                              edgeDataInitializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Edge   > &, walberla::mpi::SendBuffer & ) > > edgeDataSerializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Edge   > &, walberla::mpi::RecvBuffer & ) > > edgeDataDeserializationFunctions_;

  std::map< uint_t, std::function< void( const std::shared_ptr< Face   > & ) > >                              faceDataInitializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Face   > &, walberla::mpi::SendBuffer & ) > > faceDataSerializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Face   > &, walberla::mpi::RecvBuffer & ) > > faceDataDeserializationFunctions_;

  std::map< uint_t, std::function< void( const std::shared_ptr< Cell   > & ) > >                              cellDataInitializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Cell   > &, walberla::mpi::SendBuffer & ) > > cellDataSerializationFunctions_;
  std::map< uint_t, std::function< void( const std::shared_ptr< Cell   > &, walberla::mpi::RecvBuffer & ) > > cellDataDeserializationFunctions_;

  uint_t primitiveDataHandlers_;

  std::map< PrimitiveID::IDType, uint_t > neighborRanks_;

  void wasModified() { modificationStamp_++; }
  uint_t modificationStamp_;

  std::shared_ptr< walberla::WcTimingTree > timingTree_;

};

////////////////////////////////////////////////
// Find various template specializations here //
////////////////////////////////////////////////
#include "PrimitiveStorage.tpp"
////////////////////////////////////////////////



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
  primitives.insert( cells_.begin(), cells_.end() );
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


template< typename DataType,
          typename DataHandlingType >
void PrimitiveStorage::addCellData(       PrimitiveDataID< DataType, Cell >   & dataID,
                                    const std::shared_ptr< DataHandlingType > & dataHandling,
                                    const std::string                         & identifier )
{
  dataID = generateDataID< DataType, Cell >();
  addPrimitiveData( dataHandling, identifier, cells_, dataID );
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
inline void PrimitiveStorage::addPrimitiveData( const std::shared_ptr< DataHandlingType > & dataHandling,
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

  // Set up initialization, serialization and deserialization callbacks
  auto initCallback = [ this, dataID, dataHandling ]( const std::shared_ptr< PrimitiveType > & primitive ) -> void
  {
    primitive->data_[ dataID ] = std::shared_ptr< internal::PrimitiveData >( new internal::PrimitiveData( dataHandling->initialize( primitive.get() ) ) );
  };

  std::function< void( const std::shared_ptr< PrimitiveType > &, walberla::mpi::SendBuffer & ) > serializationCallback =
      [ dataHandling, dataID ]( const std::shared_ptr< PrimitiveType > & primitive, walberla::mpi::SendBuffer & sendBuffer ) -> void
  {
    dataHandling->serialize( primitive.get(), dataID, sendBuffer );
  };

  std::function< void( const std::shared_ptr< PrimitiveType > &, walberla::mpi::RecvBuffer & ) > deserializationCallback =
      [ dataHandling, dataID ]( const std::shared_ptr< PrimitiveType > & primitive, walberla::mpi::RecvBuffer & recvBuffer ) -> void
  {
    dataHandling->deserialize( primitive.get(), dataID, recvBuffer );
  };

  addDataHandlingCallbacks( dataID, initCallback, serializationCallback, deserializationCallback );

  for ( const auto & primitive : primitives )
  {
    initCallback( primitive.second );
  }

  wasModified();
}


} // namespace hhg


