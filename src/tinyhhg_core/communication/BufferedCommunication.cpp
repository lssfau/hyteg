
#include "tinyhhg_core/communication/BufferedCommunication.hpp"
#include "core/logging/Logging.h"

#include <functional>

namespace hhg {
namespace communication {

BufferedCommunicator::BufferedCommunicator( std::weak_ptr< PrimitiveStorage > primitiveStorage ) :
    primitiveStorage_( primitiveStorage )
{
  for ( auto & bufferSystem : bufferSystems_ )
  {
    bufferSystem = std::shared_ptr< walberla::mpi::OpenMPBufferSystem >( new walberla::mpi::OpenMPBufferSystem( walberla::mpi::MPIManager::instance()->comm() ) );
  }
}

void BufferedCommunicator::addPackInfo( const std::shared_ptr< PackInfo > & packInfo )
{
  packInfos_.push_back( packInfo );
}

void BufferedCommunicator::writeHeader( SendBuffer & sendBuffer, const PrimitiveID & senderID, const PrimitiveID & receiverID )
{
  sendBuffer << senderID << receiverID;
}

void BufferedCommunicator::readHeader ( RecvBuffer & recvBuffer,       PrimitiveID & senderID,       PrimitiveID & receiverID )
{
  recvBuffer >> senderID >> receiverID;
}

void BufferedCommunicator::receive( RecvBuffer & recvBuffer,
				                            const uint_t & numberOfMessages,
				                            const CommunicationDirection & communicationDirection )
{
  for ( uint_t message = 0; message < numberOfMessages; message++ )
  {
    PrimitiveID senderID;
    PrimitiveID receiverID;
    readHeader( recvBuffer, senderID, receiverID );

    std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();

    WALBERLA_ASSERT_NOT_NULLPTR( storage.get() );
    WALBERLA_ASSERT( storage->primitiveExistsLocally( receiverID ) );

    switch ( communicationDirection )
    {
    case VERTEX_TO_EDGE:
    {
      Edge * receivingEdge = storage->getEdge( receiverID );
      for ( auto & packInfo : packInfos_ )
      {
	      packInfo->unpackEdgeFromVertex( receivingEdge, senderID, recvBuffer );
      }
      break;
    }
    default:
      WALBERLA_ABORT( "Receive not implemented for this direction..." );
      break;
    }

  }
}

void BufferedCommunicator::startCommunication( const CommunicationDirection     & communicationDirection,
                                               const LocalCommunicationCallback & localCommunicationCallback,
                                               const PackCallback               & packCallback )
{
  if ( packInfos_.empty() )
  {
    return;
  }

  bool sendingToHigherDimension =    communicationDirection == VERTEX_TO_EDGE
                                  || communicationDirection == EDGE_TO_FACE;

  std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[ communicationDirection ];
  WALBERLA_CHECK_NOT_NULLPTR( bufferSystem.get() );
  std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();

  std::map< uint_t, std::vector< SendFunction > > sendFunctionsMap;   // rank -> sendFunctions
  std::map< uint_t, uint_t >                      ranksToReceiveFrom; // rank -> number of receives

  std::map< PrimitiveID::IDType, std::vector< Primitive::NeighborToProcessMap > > senderReceiverPairs;

  std::vector< PrimitiveID > senderIDs;

  switch ( communicationDirection )
  {
  case VERTEX_TO_EDGE:
    storage->getVertexIDs( senderIDs );
    break;
  default:
    WALBERLA_ABORT( "Not implemented" );
  }

  // Send functions
  for ( PrimitiveID senderID : senderIDs )
  {
    senderReceiverPairs[ senderID.getID() ];
    Primitive * sender = storage->getPrimitive( senderID );

    Primitive::NeighborToProcessMap neighborhood;
    if ( sendingToHigherDimension )
    {
      sender->getHigherDimNeighbors( neighborhood );
    }
    else
    {
      sender->getLowerDimNeighbors( neighborhood );
    }

    senderReceiverPairs[ senderID.getID() ].push_back( neighborhood );

    for ( auto const & neighbor : neighborhood )
    {
      PrimitiveID neighborID   = neighbor.first;
      uint_t      neighborRank = neighbor.second;

      if ( storage->primitiveExistsLocally( neighborID ) )
      {
        for ( auto & packInfo : packInfos_ )
        {
          localCommunicationCallback( senderID, neighborID, packInfo );
        }
      }
      else
      {
        if ( !packInfos_.empty() )
        {
          auto headerWriter = [ this, senderID, neighborID ]( SendBuffer & sendBuffer ) -> void { writeHeader( sendBuffer, senderID, neighborID ); };
          sendFunctionsMap[ neighborRank ].push_back( headerWriter );
        }

        for ( auto & packInfo : packInfos_ )
        {
          auto sendFunction = [ senderID, neighborID, packInfo, packCallback ]( SendBuffer & sendBuffer ) -> void { packCallback( senderID, neighborID, sendBuffer, packInfo ); };
          sendFunctionsMap[ neighborRank ].push_back( sendFunction );
        }
      }
    }
  }

  for ( auto it  = storage->beginEdges();
             it != storage->endEdges();
             it++ )
  {
    Edge * edge = it->second;
    for ( auto neighbor  = edge->beginLowerDimNeighbors();
               neighbor != edge->endLowerDimNeighbors();
               neighbor++ )
    {
      PrimitiveID neighborID   = neighbor->first;
      uint_t      neighborRank = neighbor->second;

      if ( !storage->vertexExistsLocally( neighborID ) )
      {
        if ( ranksToReceiveFrom.count( neighborRank ) == 0 )
        {
          ranksToReceiveFrom[ neighborRank ] = 1;
        }
        else
        {
          ranksToReceiveFrom[ neighborRank ] += 1;
        }
      }
    }
  }

  for ( auto it = sendFunctionsMap.begin(); it != sendFunctionsMap.end(); it++ )
  {
    uint_t                      receiverRank  = it->first;
    std::vector< SendFunction > sendFunctions = it->second;

    auto sendFunction = [ sendFunctions ]( SendBuffer & sendBuffer ) -> void { for ( auto & f : sendFunctions ) f( sendBuffer ); };

    bufferSystem->addSendingFunction( int_c( receiverRank ), sendFunction );
  }

  for ( const auto rankToReceiveFrom : ranksToReceiveFrom )
  {
    const uint_t senderRank       = rankToReceiveFrom.first;
    const uint_t numberOfMessages = rankToReceiveFrom.second;

    auto recvFunction = [ this, numberOfMessages ]( RecvBuffer & recvBuffer ) -> void { receive( recvBuffer, numberOfMessages, VERTEX_TO_EDGE ); };
    bufferSystem->addReceivingFunction( int_c( senderRank ), recvFunction );
  }

  bufferSystem->startCommunication();
}

void BufferedCommunicator::startCommunicationVertexToEdge()
{
  auto localCommunicationCallback = [ this ]( const PrimitiveID & senderID,
                                              const PrimitiveID & receiverID,
                                              const std::shared_ptr< PackInfo > & packInfo ) -> void
  {
    std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();
    WALBERLA_ASSERT_NOT_NULLPTR( storage.get() );
    WALBERLA_ASSERT( storage->vertexExistsLocally( senderID ) );
    WALBERLA_ASSERT( storage->edgeExistsLocally  ( receiverID ) );
    Vertex * sender   = storage->getVertex( senderID );
    Edge   * receiver = storage->getEdge  ( receiverID );
    packInfo->communicateLocalVertexToEdge( sender, receiver );
  };

  auto packCallback = [ this ]( const PrimitiveID & senderID,
                                const PrimitiveID & receiverID,
                                      SendBuffer  & sendBuffer,
                                const std::shared_ptr< PackInfo > & packInfo ) -> void
  {
    std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();
    WALBERLA_ASSERT_NOT_NULLPTR( storage.get() );
    WALBERLA_ASSERT( storage->vertexExistsLocally( senderID ) );
    Vertex * sender   = storage->getVertex( senderID );
    packInfo->packVertexForEdge( sender, receiverID, sendBuffer );
  };

  startCommunication( VERTEX_TO_EDGE, localCommunicationCallback, packCallback );
}

void BufferedCommunicator::endCommunicationVertexToEdge()
{
  if ( packInfos_.empty() )
  {
    return;
  }

  std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[ VERTEX_TO_EDGE ];
  bufferSystem->wait();
}

}
}
