
#pragma once

#include "tinyhhg_core/communication/PackInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/OpenMPBufferSystem.h"
#include "core/timing/TimingTree.h"
#include "core/timing/TimingPool.h"

namespace hhg {
namespace communication {

using walberla::int_c;

/// \brief Executes communication between primitives
/// \author Nils Kohl (nils.kohl@fau.de)
///
/// The \ref BufferedCommunicator can be used to perform communication between primitives.
///
/// To communicate data, the respective \ref PackInfo objects have to be attached.
///
/// The communication is buffered since the data is packed into buffers and then sent to
/// the respective processes. If two primitives reside on the same process, the communication
/// can be performed directly (using the respective \ref PackInfo methods). There are, however
/// different methods that can be set for local communication.
///
/// Since the communication is non-blocking, other actions can be performed after the communication
/// was started. When the communicated data is required, the \ref BufferedCommunicator can be forced to
/// wait for the sends and receives to complete.
///
class BufferedCommunicator
{
public:

  enum LocalCommunicationMode 
  {
    /// Uses the direct communication callbacks of the respective PackInfos for local neighbors
    DIRECT, 
    /// Sends data to local neighbors over MPI
    BUFFERED_MPI,
    /// Number of differed modes
    NUM_LOCAL_COMMUNICATION_MODES
  };

  BufferedCommunicator( std::weak_ptr< PrimitiveStorage > primitiveStorage, const LocalCommunicationMode & localCommunicationMode = DIRECT );

  /// All data that are registered via respective \ref PackInfo objects are exchanged
  void addPackInfo( const std::shared_ptr< PackInfo > & packInfo );

  /// Performs blocking communication between two \ref Primitive types.
  /// The data of the sender AND the receiver can be modified after this method returns.
  /// \param SenderType type of the sending \ref Primitive (e.g. \ref Vertex or \ref Edge)
  /// \param ReceiverType type of the receiving \ref Primitive (e.g. \ref Vertex or \ref Edge)
  template< typename SenderType,
            typename ReceiverType,
            typename = typename std::enable_if<    ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Edge   >::value )
                                                || ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Face   >::value )
                                                || ( std::is_same< SenderType, Edge   >::value && std::is_same< ReceiverType, Vertex >::value )
                                                || ( std::is_same< SenderType, Edge   >::value && std::is_same< ReceiverType, Face   >::value )
                                                || ( std::is_same< SenderType, Face   >::value && std::is_same< ReceiverType, Vertex >::value )
                                                || ( std::is_same< SenderType, Face   >::value && std::is_same< ReceiverType, Edge   >::value ) >::type >
  inline void communicate() { startCommunication< SenderType, ReceiverType >(); endCommunication< SenderType, ReceiverType >(); }

  /// Starts the non-blocking communication between two \ref Primitive types.
  /// The data of the sender can be modified after this method returns.
  /// \param SenderType type of the sending \ref Primitive (e.g. \ref Vertex or \ref Edge)
  /// \param ReceiverType type of the receiving \ref Primitive (e.g. \ref Vertex or \ref Edge)
  template< typename SenderType,
            typename ReceiverType,
            typename = typename std::enable_if<    ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Edge   >::value )
                                                || ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Face   >::value )
                                                || ( std::is_same< SenderType, Edge   >::value && std::is_same< ReceiverType, Vertex >::value )
                                                || ( std::is_same< SenderType, Edge   >::value && std::is_same< ReceiverType, Face   >::value )
                                                || ( std::is_same< SenderType, Face   >::value && std::is_same< ReceiverType, Vertex >::value )
                                                || ( std::is_same< SenderType, Face   >::value && std::is_same< ReceiverType, Edge   >::value ) >::type >
  inline void startCommunication();

  /// Ends the non-blocking communication between two \ref Primitive types
  /// Waits for the started communication to be completed. It is only safe to modify the
  /// data of the receiver after this call returned.
  /// \param SenderType type of the sending \ref Primitive (e.g. \ref Vertex or \ref Edge)
  /// \param ReceiverType type of the receiving \ref Primitive (e.g. \ref Vertex or \ref Edge)
  template< typename SenderType,
            typename ReceiverType,
            typename = typename std::enable_if<    ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Edge   >::value )
                                                || ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Face   >::value )
                                                || ( std::is_same< SenderType, Edge   >::value && std::is_same< ReceiverType, Vertex >::value )
                                                || ( std::is_same< SenderType, Edge   >::value && std::is_same< ReceiverType, Face   >::value )
                                                || ( std::is_same< SenderType, Face   >::value && std::is_same< ReceiverType, Vertex >::value )
                                                || ( std::is_same< SenderType, Face   >::value && std::is_same< ReceiverType, Edge   >::value ) >::type >
  inline void endCommunication();

  /// @name Local communication mode
  /// Getter and setter to retrieve and set the local communication mode.
  /// Depending on the mode, data is transferred via different mechanisms if the sending
  /// and the receiving \ref Primitive are located on the same process.
  /// See also \ref LocalCommunicationMode
  ///@{
  LocalCommunicationMode getLocalCommunicationMode() const { return localCommunicationMode_; }
  void setLocalCommunicationMode( const LocalCommunicationMode & localCommunicationMode ) { setupBeforeNextCommunication(); localCommunicationMode_ = localCommunicationMode; }
  ///@}

  /// Writes timing data for the setup and for the wait phase to the passed \ref walberla::WcTimingTree
  void enableTiming( const std::shared_ptr< walberla::WcTimingTree > & timingTree ) { timingTree_ = timingTree; }

private:

  typedef std::function<void ( SendBuffer & buf ) > SendFunction;
  typedef std::function<void ( RecvBuffer & buf ) > RecvFunction;

  enum CommunicationDirection
  {
    VERTEX_TO_EDGE,
    VERTEX_TO_FACE,

    EDGE_TO_VERTEX,
    EDGE_TO_FACE,

    FACE_TO_VERTEX,
    FACE_TO_EDGE,

    NUM_COMMUNICATION_DIRECTIONS
  };

  static const std::array< std::string, CommunicationDirection::NUM_COMMUNICATION_DIRECTIONS >  COMMUNICATION_DIRECTION_STRINGS;
  static const std::array< std::string, LocalCommunicationMode::NUM_LOCAL_COMMUNICATION_MODES > LOCAL_COMMUNICATION_MODE_STRINGS;

  template< typename SenderType, typename ReceiverType >
  inline CommunicationDirection getCommunicationDirection() const;

  void writeHeader( SendBuffer & sendBuffer, const PrimitiveID & senderID, const PrimitiveID & receiverID );
  void readHeader ( RecvBuffer & recvBuffer,       PrimitiveID & senderID,       PrimitiveID & receiverID );

  void endCommunication( const CommunicationDirection & communicationDirection );

  void startTimer( const std::string & timerString );
  void stopTimer ( const std::string & timerString );

  void setupBeforeNextCommunication();

  std::weak_ptr< PrimitiveStorage > primitiveStorage_;

  uint_t primitiveStorageModificationStamp_;

  std::vector< std::shared_ptr< PackInfo > > packInfos_;

  std::array< std::shared_ptr< walberla::mpi::OpenMPBufferSystem >, NUM_COMMUNICATION_DIRECTIONS > bufferSystems_;
#ifndef NDEBUG
  std::array< bool,                                                 NUM_COMMUNICATION_DIRECTIONS > communicationInProgress_;
#endif

  LocalCommunicationMode localCommunicationMode_;

  // Cached communication setup
  std::array< bool,                                   NUM_COMMUNICATION_DIRECTIONS > setupBeforeNextCommunication_;
  std::array< std::vector< std::function< void() > >, NUM_COMMUNICATION_DIRECTIONS > directCommunicationFunctions_;

  std::shared_ptr< walberla::WcTimingTree > timingTree_;

};

template< typename SenderType, typename ReceiverType >
inline BufferedCommunicator::CommunicationDirection BufferedCommunicator::getCommunicationDirection() const
{
  if ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Edge   >::value ) return VERTEX_TO_EDGE;
  if ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Face   >::value ) return VERTEX_TO_FACE;

  if ( std::is_same< SenderType, Edge   >::value && std::is_same< ReceiverType, Vertex >::value ) return EDGE_TO_VERTEX;
  if ( std::is_same< SenderType, Edge   >::value && std::is_same< ReceiverType, Face   >::value ) return EDGE_TO_FACE;

  if ( std::is_same< SenderType, Face   >::value && std::is_same< ReceiverType, Vertex >::value ) return FACE_TO_VERTEX;
  if ( std::is_same< SenderType, Face   >::value && std::is_same< ReceiverType, Edge   >::value ) return FACE_TO_EDGE;

  WALBERLA_ABORT( "Sender and receiver types are invalid" );

  return VERTEX_TO_EDGE;
}

template< typename SenderType, typename ReceiverType, typename >
void BufferedCommunicator::startCommunication()
{
  CommunicationDirection communicationDirection = getCommunicationDirection< SenderType, ReceiverType >();

  const std::string      timerString            =   "Communication (setup)";

  startTimer( timerString );

  if ( packInfos_.empty() )
  {
    stopTimer( timerString );
    return;
  }

#ifndef NDEBUG
  WALBERLA_ASSERT( !communicationInProgress_[ communicationDirection ] );
  communicationInProgress_[ communicationDirection ] = true;
#endif

  std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[ communicationDirection ];
  WALBERLA_CHECK_NOT_NULLPTR( bufferSystem.get() );

  std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();
  WALBERLA_CHECK_NOT_NULLPTR( storage.get() );

  if ( storage->getModificationStamp() != primitiveStorageModificationStamp_ )
  {
    primitiveStorageModificationStamp_ = storage->getModificationStamp();
    setupBeforeNextCommunication();
  }

  if ( setupBeforeNextCommunication_[ communicationDirection ] )
  {
    bufferSystem->clearSendingFunctions();
    bufferSystem->clearReceivingFunctions();

    directCommunicationFunctions_[ communicationDirection ].clear();

    std::map< uint_t, std::vector< SendFunction > > sendFunctionsMap;   // rank -> sendFunctions
    std::map< uint_t, uint_t >                      ranksToReceiveFrom; // rank -> number of receives

    std::vector< PrimitiveID > senderIDs;
    std::vector< PrimitiveID > receiverIDs;

    storage->getPrimitiveIDsGenerically< SenderType >  ( senderIDs );
    storage->getPrimitiveIDsGenerically< ReceiverType >( receiverIDs );

    // Send functions
    for ( const PrimitiveID & senderID : senderIDs )
    {
      WALBERLA_ASSERT( storage->primitiveExistsLocallyGenerically< SenderType >( senderID ) );
      SenderType * sender = storage->getPrimitiveGenerically< SenderType >( senderID );

      std::vector< PrimitiveID > receivingNeighborhood;
      sender->template getNeighborPrimitivesGenerically< ReceiverType >( receivingNeighborhood );

      for ( const auto & neighborID : receivingNeighborhood )
      {
        WALBERLA_ASSERT(    storage->primitiveExistsLocallyGenerically< ReceiverType >( neighborID )
                         || storage->primitiveExistsInNeighborhoodGenerically< ReceiverType >( neighborID ) );

        if (    getLocalCommunicationMode() == DIRECT
             && storage->primitiveExistsLocallyGenerically< ReceiverType >( neighborID ) )
        {
          ReceiverType * receiver = storage->getPrimitiveGenerically< ReceiverType >( neighborID );
          for ( auto & packInfo : packInfos_ )
          {
            auto directCommunicationFunction = [ sender, receiver, packInfo ]() -> void { packInfo->communicateLocal< SenderType, ReceiverType >( sender, receiver ); };
            directCommunicationFunctions_[ communicationDirection ].push_back( directCommunicationFunction );
          }
        }
        else
        {
          uint_t neighborRank = storage->getPrimitiveRank( neighborID );

          if ( !packInfos_.empty() )
          {
            auto headerWriter = [ this, senderID, neighborID ]( SendBuffer & sendBuffer ) -> void { writeHeader( sendBuffer, senderID, neighborID ); };
            sendFunctionsMap[ neighborRank ].push_back( headerWriter );
          }

          for ( auto & packInfo : packInfos_ )
          {
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

      std::vector< PrimitiveID > sendingNeighborhood;
      receiver->template getNeighborPrimitivesGenerically< SenderType >( sendingNeighborhood );

      for ( const auto & neighborID : sendingNeighborhood )
      {
        WALBERLA_ASSERT(    storage->primitiveExistsLocallyGenerically< SenderType >( neighborID )
                         || storage->primitiveExistsInNeighborhoodGenerically< SenderType >( neighborID ) );

        if (    getLocalCommunicationMode() != DIRECT
             || !storage->primitiveExistsLocallyGenerically< SenderType >( neighborID ) )
        {
          uint_t neighborRank = storage->getPrimitiveRank( neighborID );

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

    setupBeforeNextCommunication_[ communicationDirection ] = false;

  } // setup

  // Local communication
  for ( auto & directCommunicationFunction : directCommunicationFunctions_[ communicationDirection ] )
  {
    directCommunicationFunction();
  }

  // Buffered communication
  bufferSystem->startCommunication();

  stopTimer( timerString );
}

template < typename SenderType, typename ReceiverType, typename >
void BufferedCommunicator::endCommunication()
{
  const CommunicationDirection communicationDirection = getCommunicationDirection< SenderType, ReceiverType >();
  const std::string            timerString            =   "Communication (wait )";

  startTimer( timerString );

  if ( packInfos_.empty() )
  {
    stopTimer( timerString );
    return;
  }


#ifndef NDEBUG
  WALBERLA_ASSERT( communicationInProgress_[ communicationDirection ] );
  communicationInProgress_[ communicationDirection ] = false;
#endif

  std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[ communicationDirection ];
  bufferSystem->wait();

  stopTimer( timerString );
}


} // namespace communication
} // namespace hhg
