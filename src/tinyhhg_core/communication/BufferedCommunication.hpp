
#pragma once

#include "tinyhhg_core/communication/PackInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/OpenMPBufferSystem.h"

namespace hhg {
namespace communication {

class BufferedCommunicator
{
public:

  typedef std::function<void ( SendBuffer & buf ) > SendFunction;
  typedef std::function<void ( RecvBuffer & buf ) > RecvFunction;

  enum COMMUNICATION_DIRECTIONS
  {
	VERTEX_TO_EDGE,
	EDGE_TO_VERTEX,

	EDGE_TO_FACE,
	FACE_TO_EDGE,

	NUM_COMMUNICATION_DIRECTIONS
  };

  BufferedCommunicator();

  void addPackInfo( const PackInfo & packInfo );

  void startCommunicationVertexToEdge();
  void endCommunicationVertexToEdge();

private:

  std::weak_ptr< PrimitiveStorage > primitiveStorage_;
  std::vector< PackInfo > packInfos_;

  std::array< std::shared_ptr< walberla::mpi::OpenMPBufferSystem >, NUM_COMMUNICATION_DIRECTIONS > bufferSystems_;

};

} // namespace communication
} // namespace hhg
