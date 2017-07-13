
#pragma once

#include "tinyhhg_core/communication/PackInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/OpenMPBufferSystem.h"

namespace hhg {
namespace communication {

using walberla::int_c;

class BufferedCommunicator
{
public:

  typedef std::function<void ( SendBuffer & buf ) > SendFunction;
  typedef std::function<void ( RecvBuffer & buf ) > RecvFunction;

  typedef std::function< void ( const PrimitiveID & senderID,
                                const PrimitiveID & receiverID,
                                const std::shared_ptr< PrimitiveStorage > & storage,
                                const std::shared_ptr< PackInfo >         & packInfo ) > LocalCommunicationCallback;

  typedef std::function< void ( const PrimitiveID & senderID,
                                const PrimitiveID & receiverID,
                                const std::shared_ptr< PrimitiveStorage > & storage,
                                      SendBuffer                          & sendBuffer,
                                const std::shared_ptr< PackInfo >         & packInfo ) > PackCallback;

  typedef std::function< void ( const PrimitiveID & senderID,
                                const PrimitiveID & receiverID,
                                const std::shared_ptr< PrimitiveStorage > & storage,
                                      RecvBuffer                          & recvBuffer,
                                const std::shared_ptr< PackInfo >         & packInfo ) > UnpackCallback;

  BufferedCommunicator( std::weak_ptr< PrimitiveStorage > primitiveStorage );

  void addPackInfo( const std::shared_ptr< PackInfo > & packInfo );

  template< typename SenderType, typename ReceiverType >
  inline void startCommunication();

  void startCommunicationVertexToEdge();
  void endCommunicationVertexToEdge() { endCommunication( VERTEX_TO_EDGE ); }

private:

  enum CommunicationDirection
  {
    VERTEX_TO_EDGE,
    EDGE_TO_VERTEX,

    EDGE_TO_FACE,
    FACE_TO_EDGE,

    NUM_COMMUNICATION_DIRECTIONS
  };

  void writeHeader( SendBuffer & sendBuffer, const PrimitiveID & senderID, const PrimitiveID & receiverID );
  void readHeader ( RecvBuffer & recvBuffer,       PrimitiveID & senderID,       PrimitiveID & receiverID );

  void receive    (       RecvBuffer             & recvBuffer,
                    const uint_t                 & numberOfMessages,
                    const UnpackCallback         & unpackCallback );

  void startCommunication( const CommunicationDirection     & communicationDirection,
                           const LocalCommunicationCallback & localCommunicationCallback,
                           const PackCallback               & packCallback,
                           const UnpackCallback             & unpackCallback );

  void endCommunication( const CommunicationDirection & communicationDirection );

  std::weak_ptr< PrimitiveStorage > primitiveStorage_;
  std::vector< std::shared_ptr< PackInfo > > packInfos_;

  std::array< std::shared_ptr< walberla::mpi::OpenMPBufferSystem >, NUM_COMMUNICATION_DIRECTIONS > bufferSystems_;

};

template< typename SenderType, typename ReceiverType >
void BufferedCommunicator::startCommunication()
{
  if ( packInfos_.empty() )
  {
    return;
  }

  // ToDo:
  bool sendingToHigherDimension = true;
  CommunicationDirection communicationDirection = VERTEX_TO_EDGE;

  std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[ communicationDirection ];

  WALBERLA_CHECK_NOT_NULLPTR( bufferSystem.get() );
  std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();

  std::map< uint_t, std::vector< SendFunction > > sendFunctionsMap;   // rank -> sendFunctions
  std::map< uint_t, uint_t >                      ranksToReceiveFrom; // rank -> number of receives

  std::vector< PrimitiveID > senderIDs;
  std::vector< PrimitiveID > receiverIDs;

  switch ( communicationDirection )
  {
  case VERTEX_TO_EDGE:
    storage->getVertexIDs( senderIDs );
    storage->getEdgeIDs  ( receiverIDs );
    break;
  default:
    WALBERLA_ABORT( "Not implemented" );
  }

  // Send functions
  for ( const PrimitiveID & senderID : senderIDs )
  {
    WALBERLA_ASSERT( storage->primitiveExistsLocallyGenerically< SenderType >( senderID ) );
    SenderType * sender = storage->getPrimitiveGenerically< SenderType >( senderID );

    Primitive::NeighborToProcessMap receivingNeighborhood;
    if ( sendingToHigherDimension )
    {
      sender->getHigherDimNeighbors( receivingNeighborhood );
    }
    else
    {
      sender->getLowerDimNeighbors( receivingNeighborhood );
    }

    for ( const auto & neighbor : receivingNeighborhood )
    {
      PrimitiveID neighborID   = neighbor.first;
      uint_t      neighborRank = neighbor.second;

      if ( storage->primitiveExistsLocallyGenerically< ReceiverType >( neighborID ) )
      {
        ReceiverType * receiver = storage->getPrimitiveGenerically< ReceiverType >( neighborID );
        for ( auto & packInfo : packInfos_ )
        {
          packInfo->communicateLocal< SenderType, ReceiverType >( sender, receiver );
#if 0
          localCommunicationCallback( senderID, neighborID, storage, packInfo );
#endif
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
#if 0
          auto sendFunction = [ senderID, neighborID, storage, packInfo, packCallback ]( SendBuffer & sendBuffer ) -> void { packCallback( senderID, neighborID, storage, sendBuffer, packInfo ); };
#endif
          auto sendFunction = [ sender, neighborID, packInfo ]( SendBuffer & sendBuffer ) -> void { packInfo->pack< SenderType, ReceiverType >( sender, neighborID, sendBuffer ); };
          sendFunctionsMap[ neighborRank ].push_back( sendFunction );

        }
      }
    }
  }

  // Ranks to receive from
  for ( const PrimitiveID & receiverID : receiverIDs )
  {
    ReceiverType * receiver = storage->getPrimitiveGenerically< ReceiverType >( receiverID );

    Primitive::NeighborToProcessMap sendingNeighborhood;
    if ( sendingToHigherDimension )
    {
      receiver->getLowerDimNeighbors( sendingNeighborhood );
    }
    else
    {
      receiver->getHigherDimNeighbors( sendingNeighborhood );
    }

    for ( const auto & neighbor : sendingNeighborhood )
    {
      PrimitiveID neighborID   = neighbor.first;
      uint_t      neighborRank = neighbor.second;

      if ( !storage->primitiveExistsLocallyGenerically< SenderType >( neighborID ) )
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

#if 0
    auto recvFunction = [ this, numberOfMessages, unpackCallback ]( RecvBuffer & recvBuffer ) -> void { receive( recvBuffer, numberOfMessages, unpackCallback ); };
#endif

    auto recvFunction = [ this, numberOfMessages ]( RecvBuffer & recvBuffer ) -> void
    {
      for ( uint_t message = 0; message < numberOfMessages; message++ )
      {
        PrimitiveID senderID;
        PrimitiveID receiverID;
        readHeader( recvBuffer, senderID, receiverID );

        std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();

        WALBERLA_ASSERT_NOT_NULLPTR( storage.get() );
        WALBERLA_ASSERT( storage->primitiveExistsLocallyGenerically< ReceiverType >( receiverID ) );

        ReceiverType * receiver = storage->getPrimitiveGenerically< ReceiverType >( receiverID );

        for ( const auto & packInfo : packInfos_ )
        {
          packInfo->unpack< SenderType, ReceiverType >( receiver, senderID, recvBuffer);
        }
      }
    };

    bufferSystem->addReceivingFunction( int_c( senderRank ), recvFunction );
  }

  bufferSystem->startCommunication();
}



} // namespace communication
} // namespace hhg
