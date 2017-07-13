
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

void BufferedCommunicator::endCommunication( const CommunicationDirection & communicationDirection )
{
  if ( packInfos_.empty() )
  {
    return;
  }

  std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[ communicationDirection ];
  bufferSystem->wait();
}



}
}
