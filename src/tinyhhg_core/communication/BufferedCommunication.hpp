
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

  BufferedCommunicator( std::weak_ptr< PrimitiveStorage > primitiveStorage );

  void addPackInfo( const std::shared_ptr< PackInfo > & packInfo );

  void startCommunicationVertexToEdge();
  void endCommunicationVertexToEdge();

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
                    const CommunicationDirection & communicationDirection );

  void startCommunication( const CommunicationDirection & communicationDirection );

  std::weak_ptr< PrimitiveStorage > primitiveStorage_;
  std::vector< std::shared_ptr< PackInfo > > packInfos_;

  std::array< std::shared_ptr< walberla::mpi::OpenMPBufferSystem >, NUM_COMMUNICATION_DIRECTIONS > bufferSystems_;

};



} // namespace communication
} // namespace hhg
