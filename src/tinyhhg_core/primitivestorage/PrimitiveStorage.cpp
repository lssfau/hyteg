
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/mpi/OpenMPBufferSystem.h"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitives/Primitive.hpp"
#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"


#include <map>
#include <vector>

namespace hhg {

using walberla::uint_t;

PrimitiveStorage::PrimitiveStorage( const SetupPrimitiveStorage & setupStorage ) :
  primitiveDataHandlers_( 0 )
{
  for ( auto it = setupStorage.beginVertices(); it != setupStorage.endVertices(); it++  )
  {
    if ( uint_c( walberla::mpi::MPIManager::instance()->rank() ) == setupStorage.getTargetRank( it->first ) )
    {
      vertices_[ it->first ] = std::make_shared< Vertex >( *it->second );
    }
  }

  for ( auto it = setupStorage.beginEdges(); it != setupStorage.endEdges(); it++ )
  {
    if ( uint_c( walberla::mpi::MPIManager::instance()->rank() )  == setupStorage.getTargetRank( it->first ) )
    {
      edges_[ it->first ] = std::make_shared< Edge >( *it->second );
    }
  }

  for ( auto it = setupStorage.beginFaces(); it != setupStorage.endFaces(); it++ )
  {
    if ( uint_c( walberla::mpi::MPIManager::instance()->rank() )  == setupStorage.getTargetRank( it->first ) )
    {
      faces_[ it->first ] = std::make_shared< Face >( *it->second );
    }
  }

  // Neighborhood

  for ( const auto & it : vertices_ )
  {
    auto vertex = it.second;

    for ( const auto & neighborVertexID : vertex->neighborVertices() )
    {
      const Vertex * neighborVertex = setupStorage.getVertex( neighborVertexID );
      if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
      {
        neighborVertices_[ neighborVertexID.getID() ] = std::make_shared< Vertex >( *neighborVertex );
        neighborRanks_[ neighborVertexID.getID() ] = setupStorage.getTargetRank( neighborVertexID.getID() );
      }
    }

    for ( const auto & neighborEdgeID : vertex->neighborEdges() )
    {
      const Edge * neighborEdge = setupStorage.getEdge( neighborEdgeID );
      if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
      {
        neighborEdges_[ neighborEdgeID.getID() ] = std::make_shared< Edge >( *neighborEdge );
        neighborRanks_[ neighborEdgeID.getID() ] = setupStorage.getTargetRank( neighborEdgeID.getID() );
      }
    }

    for ( const auto & neighborFaceID : vertex->neighborFaces() )
    {
      const Face * neighborFace = setupStorage.getFace( neighborFaceID );
      if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
      {
        neighborFaces_[ neighborFaceID.getID() ] = std::make_shared< Face >( *neighborFace );
        neighborRanks_[ neighborFaceID.getID() ] = setupStorage.getTargetRank( neighborFaceID.getID() );
      }
    }
  }

  for ( const auto & it : edges_ )
  {
    auto edge = it.second;

    for ( const auto & neighborVertexID : edge->neighborVertices() )
    {
      const Vertex * neighborVertex = setupStorage.getVertex( neighborVertexID );
      if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
      {
        neighborVertices_[ neighborVertexID.getID() ] = std::make_shared< Vertex >( *neighborVertex );
        neighborRanks_[ neighborVertexID.getID() ] = setupStorage.getTargetRank( neighborVertexID.getID() );
      }
    }

    for ( const auto & neighborEdgeID : edge->neighborEdges() )
    {
      const Edge * neighborEdge = setupStorage.getEdge( neighborEdgeID );
      if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
      {
        neighborEdges_[ neighborEdgeID.getID() ] = std::make_shared< Edge >( *neighborEdge );
        neighborRanks_[ neighborEdgeID.getID() ] = setupStorage.getTargetRank( neighborEdgeID.getID() );
      }
    }

    for ( const auto & neighborFaceID : edge->neighborFaces() )
    {
      const Face * neighborFace = setupStorage.getFace( neighborFaceID );
      if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
      {
        neighborFaces_[ neighborFaceID.getID() ] = std::make_shared< Face >( *neighborFace );
        neighborRanks_[ neighborFaceID.getID() ] = setupStorage.getTargetRank( neighborFaceID.getID() );
      }
    }
  }

  for ( const auto & it : faces_ )
  {
    auto face = it.second;

    for ( const auto & neighborVertexID : face->neighborVertices() )
    {
      const Vertex * neighborVertex = setupStorage.getVertex( neighborVertexID );
      if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
      {
        neighborVertices_[ neighborVertexID.getID() ] = std::make_shared< Vertex >( *neighborVertex );
        neighborRanks_[ neighborVertexID.getID() ] = setupStorage.getTargetRank( neighborVertexID.getID() );
      }
    }

    for ( const auto & neighborEdgeID : face->neighborEdges() )
    {
      const Edge * neighborEdge = setupStorage.getEdge( neighborEdgeID );
      if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
      {
        neighborEdges_[ neighborEdgeID.getID() ] = std::make_shared< Edge >( *neighborEdge );
        neighborRanks_[ neighborEdgeID.getID() ] = setupStorage.getTargetRank( neighborEdgeID.getID() );
      }
    }

    for ( const auto & neighborFaceID : face->neighborFaces() )
    {
      const Face * neighborFace = setupStorage.getFace( neighborFaceID );
      if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
      {
        neighborFaces_[ neighborFaceID.getID() ] = std::make_shared< Face >( *neighborFace );
        neighborRanks_[ neighborFaceID.getID() ] = setupStorage.getTargetRank( neighborFaceID.getID() );
      }
    }
  }


#ifndef NDEBUG
  checkConsistency();
#endif
}



void PrimitiveStorage::getPrimitives( PrimitiveMap & primitiveMap ) const
{
  primitiveMap.clear();

  primitiveMap.insert( beginVertices(), endVertices() );
  primitiveMap.insert( beginEdges(), endEdges() );
  primitiveMap.insert( beginFaces(), endFaces() );

  WALBERLA_ASSERT_EQUAL( primitiveMap.size(), vertices_.size() + edges_.size() + faces_.size() );
}


const Primitive* PrimitiveStorage::getPrimitive( const PrimitiveID & id ) const
{
  if ( vertexExistsLocally( id ) || vertexExistsInNeighborhood( id ) ) return getVertex( id );
  if (   edgeExistsLocally( id ) ||   edgeExistsInNeighborhood( id ) ) return   getEdge( id );
  if (   faceExistsLocally( id ) ||   faceExistsInNeighborhood( id ) ) return   getFace( id );
  return nullptr;
}

Primitive* PrimitiveStorage::getPrimitive( const PrimitiveID & id )
{
  if ( vertexExistsLocally( id ) || vertexExistsInNeighborhood( id ) ) return getVertex( id );
  if (   edgeExistsLocally( id ) ||   edgeExistsInNeighborhood( id ) ) return   getEdge( id );
  if (   faceExistsLocally( id ) ||   faceExistsInNeighborhood( id ) ) return   getFace( id );
  return nullptr;
}

const Vertex* PrimitiveStorage::getVertex( const PrimitiveID & id ) const
{
  if ( vertexExistsLocally( id ) )
  {
    return vertices_.at( id.getID() ).get();
  }
  else if ( vertexExistsInNeighborhood( id ) )
  {
    return neighborVertices_.at( id.getID() ).get();
  }
  else
  {
    return nullptr;
  }
}

Vertex* PrimitiveStorage::getVertex( const PrimitiveID & id )
{
  if ( vertexExistsLocally( id ) )
  {
    return vertices_[ id.getID() ].get();
  }
  else if ( vertexExistsInNeighborhood( id ) )
  {
    return neighborVertices_[ id.getID() ].get();
  }
  else
  {
    return nullptr;
  }
}

const Edge* PrimitiveStorage::getEdge( const PrimitiveID & id ) const
{
  if ( edgeExistsLocally( id ) )
  {
    return edges_.at( id.getID() ).get();
  }
  else if ( edgeExistsInNeighborhood( id ) )
  {
    return neighborEdges_.at( id.getID() ).get();
  }
  else
  {
    return nullptr;
  }
}

Edge* PrimitiveStorage::getEdge( const PrimitiveID & id )
{
  if ( edgeExistsLocally( id ) )
  {
    return edges_[ id.getID() ].get();
  }
  else if ( edgeExistsInNeighborhood( id ) )
  {
    return neighborEdges_[ id.getID() ].get();
  }
  else
  {
    return nullptr;
  }
}

const Face* PrimitiveStorage::getFace( const PrimitiveID & id ) const
{
  if ( faceExistsLocally( id ) )
  {
    return faces_.at( id.getID() ).get();
  }
  else if ( faceExistsInNeighborhood( id ) )
  {
    return neighborFaces_.at( id.getID() ).get();
  }
  else
  {
    return nullptr;
  }
}

Face* PrimitiveStorage::getFace( const PrimitiveID & id )
{
  if ( faceExistsLocally( id ) )
  {
    return faces_[ id.getID() ].get();
  }
  else if ( faceExistsInNeighborhood( id ) )
  {
    return neighborFaces_[ id.getID() ].get();
  }
  else
  {
    return nullptr;
  }
}

void PrimitiveStorage::getPrimitiveIDs ( std::vector< PrimitiveID > & primitiveIDs ) const
{
  primitiveIDs.clear();

  std::vector< PrimitiveID > someIDs;

  getVertexIDs( someIDs );
  primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );

  getEdgeIDs( someIDs );
  primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );

  getFaceIDs( someIDs );
  primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );
}

void PrimitiveStorage::getVertexIDs ( std::vector< PrimitiveID > & vertexIDs ) const
{
  vertexIDs.clear();
  for ( auto const & it : vertices_ )
  {
    vertexIDs.push_back( it.first );
  }
}

void PrimitiveStorage::getEdgeIDs ( std::vector< PrimitiveID > & edgeIDs ) const
{
  edgeIDs.clear();
  for ( auto const & it : edges_ )
  {
    edgeIDs.push_back( it.first );
  }
}

void PrimitiveStorage::getFaceIDs ( std::vector< PrimitiveID > & faceIDs ) const
{
  faceIDs.clear();
  for ( auto const & it : faces_ )
  {
    faceIDs.push_back( it.first );
  }
}

uint_t PrimitiveStorage::getPrimitiveRank ( const PrimitiveID & id ) const
{
  WALBERLA_ASSERT( primitiveExistsLocally( id ) || primitiveExistsInNeighborhood( id ) );
  if ( primitiveExistsLocally( id ) )
  {
    return uint_c( walberla::mpi::MPIManager::instance()->rank() );
  }
  else
  {
    return getNeighborPrimitiveRank( id );
  }
}

void PrimitiveStorage::migratePrimitives( const std::map< PrimitiveID::IDType, uint_t > & primitivesToMigrate )
{
  uint_t rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );
  uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  walberla::mpi::OpenMPBufferSystem bufferSystem( walberla::mpi::MPIManager::instance()->comm() );


  ///////////////////////////////////
  // Serialization and sender side //
  ///////////////////////////////////
  std::map< uint_t, std::vector< std::function< void( SendBuffer & ) > > > sendingFunctions;

  for ( const auto & it : primitivesToMigrate )
  {
    PrimitiveID primitiveID = it.first;
    uint_t      targetRank  = it.second;

    WALBERLA_CHECK( primitiveExistsLocally( primitiveID ), "Cannot migrate non-existent primitives." );
    WALBERLA_CHECK_LESS( targetRank, numProcesses );

    if ( targetRank == rank )
    {
      continue;
    }

    // serialize
    // - true (== hasContent)
    // - source rank
    // - primitive ID
    // - primitive type
    // - primitive (contains neighborhood IDs)
    // - primitive data
    // - for all neighbors (from lower to higher dimension):
    //   - neighbor primitive ID
    //   - neighbor primitive
    auto sendingFunction = [ = ]( SendBuffer & sendBuffer ) -> void
    {
      sendBuffer << true;
      WALBERLA_LOG_DEVEL( "Serializing source rank: " << rank );
      sendBuffer << rank;
      WALBERLA_LOG_DEVEL( "Serializing PrimitiveID: " << primitiveID.getID() );
      sendBuffer << primitiveID;

      sendBuffer << getPrimitiveType( primitiveID );

      Primitive * primitive = getPrimitive( primitiveID );

      sendBuffer << *primitive;

    };

    sendingFunctions[ targetRank ].push_back( sendingFunction );
  }


  auto emptySendingFunction = []( SendBuffer & sendBuffer ) -> void
  {
    // hasContent == false
    sendBuffer << false;
  };

  for ( uint_t receiverRank = 0; receiverRank < numProcesses; receiverRank++ )
  {
    if ( sendingFunctions.count( receiverRank ) == 0 )
    {
      sendingFunctions[ receiverRank ].push_back( emptySendingFunction );
    }
  }

  for ( const auto & sendFunctionVectors : sendingFunctions )
  {
    uint_t targetRank                                                       = sendFunctionVectors.first;
    std::vector< std::function< void( SendBuffer & ) > > sendFunctionVector = sendFunctionVectors.second;

    auto sendingFunctionExecuter = [ sendFunctionVector ]( SendBuffer & sendBuffer ) -> void
    {
      for ( const auto & sendingFunction : sendFunctionVector ) { sendingFunction( sendBuffer ); }
    };

    bufferSystem.addSendingFunction( static_cast< walberla::mpi::MPIRank >( targetRank ), sendingFunctionExecuter );
  }


  ///////////////////////////////////////
  // Deserialization and receiver side //
  ///////////////////////////////////////
  auto receivingFunction = [ = ]( RecvBuffer & recvBuffer ) -> void
  {
    while ( !recvBuffer.isEmpty() )
    {
      bool          hasContent;
      recvBuffer >> hasContent;

      if ( hasContent )
      {
        uint_t        sourceRank;
        PrimitiveID   primitiveID;

        PrimitiveTypeEnum primitiveType;

        recvBuffer >> sourceRank;
        WALBERLA_LOG_DEVEL( "Deserializing source rank: "<< sourceRank );
        recvBuffer >> primitiveID;
        WALBERLA_LOG_DEVEL( "Deserializing PrimitiveID: "<< primitiveID.getID() );
        recvBuffer >> primitiveType;

        switch ( primitiveType )
        {
        case VERTEX:
        {

          std::shared_ptr< Vertex > vertex = std::make_shared< Vertex >( primitiveID, Point3D() );
          recvBuffer >> *vertex;
          WALBERLA_LOG_INFO( *vertex );
          break;

        }
        default:
          break;

        }

      }
    }
  };

  for ( uint_t senderRank = 0; senderRank < numProcesses; senderRank++ )
  {
    bufferSystem.addReceivingFunction( static_cast< walberla::mpi::MPIRank >( senderRank ), receivingFunction );
  }


  //////////////////////////////
  // Performing communication //
  //////////////////////////////
  bufferSystem.startCommunication();
  bufferSystem.wait();
}


PrimitiveStorage::PrimitiveTypeEnum PrimitiveStorage::getPrimitiveType( const PrimitiveID & primitiveID ) const
{
  if ( vertexExistsLocally( primitiveID ) ) return VERTEX;
  if ( edgeExistsLocally  ( primitiveID ) ) return EDGE;
  if ( faceExistsLocally  ( primitiveID ) ) return FACE;
  return INVALID;
}


void PrimitiveStorage::checkConsistency()
{
  // 1. Number of data entries less than local counter
  // 2. PrimitiveIDs of maps match IDs of Primitives
  // 3. Neighborhood of Primitives
  for ( auto it = vertices_.begin(); it != vertices_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
  }
  for ( auto it = edges_.begin(); it != edges_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
  }
  for ( auto it = faces_.begin(); it != faces_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 3 );
  }

  // 4. Number of data entries of neighbor primitives is zero
  // 5. PrimitiveIDs of neighbor maps match IDs of neighbor Primitives
  // 6. Neighborhood of Primitives
  for ( auto it = neighborVertices_.begin(); it != neighborVertices_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
  }
  for ( auto it = neighborEdges_.begin(); it != neighborEdges_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
  }
  for ( auto it = neighborFaces_.begin(); it != neighborFaces_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 3 );
  }

  // 7. Number of callbacks is less or equal to the data handling counter
  WALBERLA_CHECK_LESS_EQUAL( primitiveDataInitializationFunctions_.size() , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( primitiveDataSerializationFunctions_.size()  , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( primitiveDataDeserializationFunctions_.size(), primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( vertexDataInitializationFunctions_.size()    , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( vertexDataSerializationFunctions_.size()     , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( vertexDataDeserializationFunctions_.size()   , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( edgeDataInitializationFunctions_.size()      , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( edgeDataSerializationFunctions_.size()       , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( edgeDataDeserializationFunctions_.size()     , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( faceDataInitializationFunctions_.size()      , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( faceDataSerializationFunctions_.size()       , primitiveDataHandlers_ );
  WALBERLA_CHECK_LESS_EQUAL( faceDataDeserializationFunctions_.size()     , primitiveDataHandlers_ );

}


} // namespace hhg

