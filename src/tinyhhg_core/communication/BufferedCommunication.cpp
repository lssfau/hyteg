
#include "tinyhhg_core/communication/BufferedCommunication.hpp"

namespace hhg {
namespace communication {

BufferedCommunicator::BufferedCommunicator( std::weak_ptr< PrimitiveStorage > primitiveStorage ) :
    primitiveStorage_( primitiveStorage )
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

  // Send functions
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
        auto sendFunction = [ &packInfo, &vertex, &neighborID ]( SendBuffer & sendBuffer ) -> void { packInfo.packVertexForEdge( vertex, neighborID, sendBuffer ); };
        bufferSystem->addSendingFunction( int_c( neighborRank ), sendFunction );
      }
    }
  }

  // Recv functions
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

      for ( auto & packInfo : packInfos_ )
      {
	auto recvFunction = [ &packInfo, &edge, &neighborID ]( RecvBuffer & recvBuffer ) -> void { packInfo.unpackVertexFromEdge( edge, neighborID, recvBuffer ); };
	bufferSystem->addReceivingFunction( int_c( neighborRank ), recvFunction );
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
