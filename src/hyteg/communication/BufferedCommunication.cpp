
#include "hyteg/communication/BufferedCommunication.hpp"
#include "core/logging/Logging.h"

#include <functional>

namespace hyteg {
namespace communication {

std::atomic_uint BufferedCommunicator::bufferSystemTag_( 0 );

const uint_t BufferedCommunicator::SYNC_WORD( 1234 );

BufferedCommunicator::BufferedCommunicator( std::weak_ptr< PrimitiveStorage > primitiveStorage, const LocalCommunicationMode & localCommunicationMode ) :
    primitiveStorage_( primitiveStorage ),
    primitiveStorageModificationStamp_( primitiveStorage_.lock()->getModificationStamp() ),
    localCommunicationMode_( localCommunicationMode )
{
  const bool serialSends = true;
  const bool serialRecvs = true;
  for ( auto & bufferSystem : bufferSystems_ )
  {
    bufferSystem = std::shared_ptr< walberla::mpi::OpenMPBufferSystem >( new walberla::mpi::OpenMPBufferSystem( walberla::mpi::MPIManager::instance()->comm(), int_c( bufferSystemTag_++ ), serialSends, serialRecvs ) );
  }

  setupBeforeNextCommunication();

  for ( auto & communicationInProgress : communicationInProgress_ )
  {
    communicationInProgress = false;
  }
}

void BufferedCommunicator::addPackInfo( const std::shared_ptr< PackInfo > & packInfo )
{
  setupBeforeNextCommunication();
  packInfos_.push_back( packInfo );
}

void BufferedCommunicator::setLocalCommunicationMode( const LocalCommunicationMode & localCommunicationMode )
{
  for ( auto & communicationInProgress : communicationInProgress_ )
  {
    WALBERLA_CHECK( !communicationInProgress );
  }

  setupBeforeNextCommunication();
  localCommunicationMode_ = localCommunicationMode;
}

void BufferedCommunicator::writeHeader( SendBuffer & sendBuffer, const PrimitiveID & senderID, const PrimitiveID & receiverID )
{
  WALBERLA_DEBUG_SECTION()
  {
    sendBuffer << SYNC_WORD;
  }
  sendBuffer << senderID << receiverID;
}

void BufferedCommunicator::readHeader ( RecvBuffer & recvBuffer,       PrimitiveID & senderID,       PrimitiveID & receiverID )
{
  WALBERLA_DEBUG_SECTION()
  {
    uint_t syncWord;
    recvBuffer >> syncWord;
    WALBERLA_ASSERT_EQUAL( syncWord, SYNC_WORD, "Could not sync during unpacking. Chances are that the amount of data packed was not equal the amount of data unpacked." );
  }
  recvBuffer >> senderID >> receiverID;
}

void BufferedCommunicator::startTimer( const std::string & timerString )
{
  if ( timingTree_ )
  {
    timingTree_->start( timerString );
  }
}

void BufferedCommunicator::stopTimer( const std::string & timerString )
{
  if ( timingTree_ )
  {
    timingTree_->stop( timerString );
  }
}

void BufferedCommunicator::setupBeforeNextCommunication()
{
  setupBeforeNextCommunication_.fill( true );
}

}
}
