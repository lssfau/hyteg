
#include "tinyhhg_core/communication/BufferedCommunication.hpp"

namespace hhg {
namespace communication {

BufferedCommunicator::BufferedCommunicator()
{
  for ( auto bufferSystem : bufferSystems_ )
  {
	bufferSystem = std::shared_ptr< walberla::mpi::OpenMPBufferSystem >( new walberla::mpi::OpenMPBufferSystem( walberla::mpi::MPIManager::instance()->comm() ) );
  }
}

void BufferedCommunicator::startCommunicationVertexToEdge()
{
  if ( packInfos_.empty() )
  {
	return;
  }

  std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[ VERTEX_TO_EDGE ];
  std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();

  for ( auto it  = storage->beginVertices();
		     it != storage->endVertices();
		     it++ )
  {
	  Vertex * vertex = it->second;
	  for ( auto neighbor  = vertex->beginHigherDimNeighbors();
			     neighbor != vertex->endHigherDimNeighbors();
			     neighbor++ )
	  {
		PrimitiveID neighborID   = neighbor->first;
        uint_t      neighborRank = neighbor->second;

        for ( auto & packInfo : packInfos_ )
        {
          // Does not work since PackInfo pack functions are pure virtual.
          // Either we add empty implementations in the PackInfo class or we need a different solution here.
          // (Maybe this is screwed up anyway..)
          // SendFunction sendFunction = std::bind( &PackInfo::packVertexForEdge, packInfo, vertex, neighborID, _1 );
          // bufferSystem->addSendingFunction( sendFunction );
        }

	  }
  }

  // bufferSystem->sendAll();

}

void BufferedCommunicator::endCommunicationVertexToEdge()
{
}

}
}
