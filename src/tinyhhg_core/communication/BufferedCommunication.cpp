
#include "tinyhhg_core/communication/BufferedCommunication.hpp"

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

void BufferedCommunicator::startCommunicationVertexToEdge()
{
  if ( packInfos_.empty() )
  {
    return;
  }

  std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[ VERTEX_TO_EDGE ];
  WALBERLA_CHECK_NOT_NULLPTR( bufferSystem.get() );
  std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();

  std::map< uint_t, std::vector< SendFunction > > sendFunctionsMap;

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
        auto sendFunction = [ packInfo, vertex, neighborID ]( SendBuffer & sendBuffer ) -> void { packInfo->packVertexForEdge( vertex, neighborID, sendBuffer ); };
        sendFunctionsMap[ neighborRank ].push_back( sendFunction );
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
	auto recvFunction = [ packInfo, edge, neighborID ]( RecvBuffer & recvBuffer ) -> void { packInfo->unpackEdgeFromVertex( edge, neighborID, recvBuffer ); };
	// bufferSystem->addReceivingFunction( int_c( neighborRank ), recvFunction );
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

  bufferSystem->startCommunication();

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
