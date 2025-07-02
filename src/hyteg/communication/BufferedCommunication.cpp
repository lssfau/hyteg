/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Nils Kohl, Marcus Mohr, Andreas Burkhart.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "hyteg/communication/BufferedCommunication.hpp"

#include <functional>

#include "core/logging/Logging.h"
#include "core/mpi/OpenMPBufferSystem.h"
#include "core/timing/TimingTree.h"

#include "hyteg/communication/MPITagProvider.hpp"
#include "hyteg/communication/PackInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {
namespace communication {
const uint_t BufferedCommunicator::SYNC_WORD( 1234 );

BufferedCommunicator::BufferedCommunicator( std::weak_ptr< PrimitiveStorage > primitiveStorage,
                                            const LocalCommunicationMode&     localCommunicationMode )
: primitiveStorage_( primitiveStorage )
, primitiveStorageModificationStamp_( primitiveStorage_.lock()->getModificationStamp() )
, localCommunicationMode_( localCommunicationMode )
{
   const bool serialSends = true;
   const bool serialRecvs = true;
   for ( auto& bufferSystem : bufferSystems_ )
   {
      auto newTag = MPITagProvider::getMPITag();
      claimedMPITags_.push_back( newTag );
      bufferSystem = std::shared_ptr< walberla::mpi::OpenMPBufferSystem >( new walberla::mpi::OpenMPBufferSystem(
          walberla::mpi::MPIManager::instance()->comm(), newTag, serialSends, serialRecvs ) );
   }

   setupBeforeNextCommunication();

   for ( auto& communicationInProgress : communicationInProgress_ )
   {
      communicationInProgress = false;
   }
}

BufferedCommunicator::~BufferedCommunicator()
{
   for ( auto& tag : claimedMPITags_ )
   {
      MPITagProvider::returnMPITag( tag );
   }
}

void BufferedCommunicator::addPackInfo( const std::shared_ptr< PackInfo >& packInfo )
{
   setupBeforeNextCommunication();
   packInfos_.push_back( packInfo );
}

void BufferedCommunicator::setLocalCommunicationMode( const LocalCommunicationMode& localCommunicationMode )
{
   for ( auto& communicationInProgress : communicationInProgress_ )
   {
      WALBERLA_CHECK( !communicationInProgress );
   }

   setupBeforeNextCommunication();
   localCommunicationMode_ = localCommunicationMode;
}

void BufferedCommunicator::enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree )
{
   timingTree_ = timingTree;
}

void BufferedCommunicator::writeHeader( SendBuffer& sendBuffer, const PrimitiveID& senderID, const PrimitiveID& receiverID )
{
   WALBERLA_DEBUG_SECTION()
   {
      sendBuffer << SYNC_WORD;
   }
   sendBuffer << senderID << receiverID;
}

void BufferedCommunicator::readHeader( RecvBuffer& recvBuffer, PrimitiveID& senderID, PrimitiveID& receiverID )
{
   WALBERLA_DEBUG_SECTION()
   {
      uint_t syncWord;
      recvBuffer >> syncWord;
      WALBERLA_ASSERT_EQUAL(
          syncWord,
          SYNC_WORD,
          "Could not sync during unpacking. Chances are that the amount of data packed was not equal the amount of data unpacked." );
   }
   recvBuffer >> senderID >> receiverID;
}

void BufferedCommunicator::startTimer( const std::string& timerString )
{
   if ( timingTree_ )
   {
      timingTree_->start( timerString );
   }
}

void BufferedCommunicator::stopTimer( const std::string& timerString )
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

template < typename SenderType, typename ReceiverType >
BufferedCommunicator::CommunicationDirection BufferedCommunicator::getCommunicationDirection() const
{
   staticAssertCommunicationDirections< SenderType, ReceiverType >();

   if ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Edge >::value )
      return VERTEX_TO_EDGE;
   if ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Face >::value )
      return VERTEX_TO_FACE;
   if ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Cell >::value )
      return VERTEX_TO_CELL;

   if ( std::is_same< SenderType, Edge >::value && std::is_same< ReceiverType, Vertex >::value )
      return EDGE_TO_VERTEX;
   if ( std::is_same< SenderType, Edge >::value && std::is_same< ReceiverType, Face >::value )
      return EDGE_TO_FACE;
   if ( std::is_same< SenderType, Edge >::value && std::is_same< ReceiverType, Cell >::value )
      return EDGE_TO_CELL;

   if ( std::is_same< SenderType, Face >::value && std::is_same< ReceiverType, Vertex >::value )
      return FACE_TO_VERTEX;
   if ( std::is_same< SenderType, Face >::value && std::is_same< ReceiverType, Edge >::value )
      return FACE_TO_EDGE;
   if ( std::is_same< SenderType, Face >::value && std::is_same< ReceiverType, Face >::value )
      return FACE_TO_FACE;
   if ( std::is_same< SenderType, Face >::value && std::is_same< ReceiverType, Cell >::value )
      return FACE_TO_CELL;

   if ( std::is_same< SenderType, Cell >::value && std::is_same< ReceiverType, Vertex >::value )
      return CELL_TO_VERTEX;
   if ( std::is_same< SenderType, Cell >::value && std::is_same< ReceiverType, Edge >::value )
      return CELL_TO_EDGE;
   if ( std::is_same< SenderType, Cell >::value && std::is_same< ReceiverType, Face >::value )
      return CELL_TO_FACE;
   if ( std::is_same< SenderType, Cell >::value && std::is_same< ReceiverType, Cell >::value )
      return CELL_TO_CELL;

   WALBERLA_ABORT( "Sender and receiver types are invalid" );

   return NUM_COMMUNICATION_DIRECTIONS;
}

template < typename SenderType, typename ReceiverType >
void BufferedCommunicator::startCommunication( std::vector< PrimitiveID > excludeReceivingIDs )
{
   staticAssertCommunicationDirections< SenderType, ReceiverType >();

   CommunicationDirection communicationDirection = getCommunicationDirection< SenderType, ReceiverType >();

   const std::string timerStringSetup    = "Communication (setup                   )";
   const std::string timerStringDirect   = "Communication (direct                  )";
   const std::string timerStringBuffered = "Communication (buffered / pack         )";

   startTimer( timerStringSetup );

   if ( packInfos_.empty() )
   {
      stopTimer( timerStringSetup );
      return;
   }

   WALBERLA_ASSERT( !communicationInProgress_[communicationDirection] );
   communicationInProgress_[communicationDirection] = true;

   std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[communicationDirection];
   WALBERLA_CHECK_NOT_NULLPTR( bufferSystem.get() );

   std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( storage.get() );

   if ( storage->getModificationStamp() != primitiveStorageModificationStamp_ )
   {
      primitiveStorageModificationStamp_ = storage->getModificationStamp();
      setupBeforeNextCommunication();
   }

   if ( setupBeforeNextCommunication_[communicationDirection] )
   {
      bufferSystem->clearSendingFunctions();
      bufferSystem->clearReceivingFunctions();

      directCommunicationFunctions_[communicationDirection].clear();

      std::map< uint_t, std::vector< SendFunction > > sendFunctionsMap;   // rank -> sendFunctions
      std::map< uint_t, uint_t >                      ranksToReceiveFrom; // rank -> number of receives

      std::vector< PrimitiveID > senderIDs;
      std::vector< PrimitiveID > receiverIDs;

      storage->getPrimitiveIDsGenerically< SenderType >( senderIDs );
      storage->getPrimitiveIDsGenerically< ReceiverType >( receiverIDs );

      for ( const PrimitiveID& excludeID : excludeReceivingIDs )
      {
         receiverIDs.erase( std::remove( receiverIDs.begin(), receiverIDs.end(), excludeID ), receiverIDs.end() );
      }

      // Send functions
      for ( const PrimitiveID& senderID : senderIDs )
      {
         WALBERLA_ASSERT( storage->primitiveExistsLocallyGenerically< SenderType >( senderID ) );
         SenderType* sender = storage->getPrimitiveGenerically< SenderType >( senderID );

         std::vector< PrimitiveID > receivingNeighborhood;
         if constexpr ( std::is_same_v< SenderType, Face > && std::is_same_v< ReceiverType, Face > )
         {
            for ( const auto& [_, neighbors] : sender->getIndirectTopLevelNeighborFaceIDsOverEdges() )
            {
               for ( const auto& nFacePID : neighbors )
               {
                  receivingNeighborhood.push_back( nFacePID );
               }
            }
         }
         else if constexpr ( std::is_same_v< SenderType, Cell > && std::is_same_v< ReceiverType, Cell > )
         {
            for ( const auto& [_, pid] : sender->getIndirectNeighborCellIDsOverFaces() )
            {
               receivingNeighborhood.push_back( pid );
            }
         }
         else
         {
            sender->template getNeighborPrimitivesGenerically< ReceiverType >( receivingNeighborhood );
         }

         for ( const PrimitiveID& excludeID : excludeReceivingIDs )
         {
            receivingNeighborhood.erase( std::remove( receivingNeighborhood.begin(), receivingNeighborhood.end(), excludeID ),
                                         receivingNeighborhood.end() );
         }

         for ( const auto& neighborID : receivingNeighborhood )
         {
            WALBERLA_ASSERT( storage->primitiveExistsLocallyGenerically< ReceiverType >( neighborID ) ||
                             storage->primitiveExistsInNeighborhoodGenerically< ReceiverType >( neighborID ) );

            if ( getLocalCommunicationMode() == DIRECT &&
                 storage->primitiveExistsLocallyGenerically< ReceiverType >( neighborID ) )
            {
               ReceiverType* receiver = storage->getPrimitiveGenerically< ReceiverType >( neighborID );
               for ( auto& packInfo : packInfos_ )
               {
                  auto directCommunicationFunction = [sender, receiver, packInfo]() -> void {
                     packInfo->communicateLocal< SenderType, ReceiverType >( sender, receiver );
                  };
                  directCommunicationFunctions_[communicationDirection].push_back( directCommunicationFunction );
               }
            }
            else
            {
               uint_t neighborRank = storage->getPrimitiveRank( neighborID );

               if ( !packInfos_.empty() )
               {
                  auto headerWriter = [this, senderID, neighborID]( SendBuffer& sendBuffer ) -> void {
                     writeHeader( sendBuffer, senderID, neighborID );
                  };
                  sendFunctionsMap[neighborRank].push_back( headerWriter );
               }

               for ( auto& packInfo : packInfos_ )
               {
                  auto sendFunction = [this, sender, neighborID, packInfo]( SendBuffer& sendBuffer ) -> void {
                     startTimer( "Packing" );
                     packInfo->pack< SenderType, ReceiverType >( sender, neighborID, sendBuffer );
                     stopTimer( "Packing" );
                  };
                  sendFunctionsMap[neighborRank].push_back( sendFunction );
               }
            }
         }
      }

      // Ranks to receive from
      for ( const PrimitiveID& receiverID : receiverIDs )
      {
         ReceiverType* receiver = storage->getPrimitiveGenerically< ReceiverType >( receiverID );

         std::vector< PrimitiveID > sendingNeighborhood;
         if constexpr ( std::is_same_v< SenderType, Face > && std::is_same_v< ReceiverType, Face > )
         {
            for ( const auto& [_, neighbors] : receiver->getIndirectTopLevelNeighborFaceIDsOverEdges() )
            {
               for ( const auto& nFacePID : neighbors )
               {
                  sendingNeighborhood.push_back( nFacePID );
               }
            }
         }
         else if constexpr ( std::is_same_v< SenderType, Cell > && std::is_same_v< ReceiverType, Cell > )
         {
            for ( const auto& [_, pid] : receiver->getIndirectNeighborCellIDsOverFaces() )
            {
               sendingNeighborhood.push_back( pid );
            }
         }
         else
         {
            receiver->template getNeighborPrimitivesGenerically< SenderType >( sendingNeighborhood );
         }

         for ( const auto& neighborID : sendingNeighborhood )
         {
            WALBERLA_ASSERT( storage->primitiveExistsLocallyGenerically< SenderType >( neighborID ) ||
                             storage->primitiveExistsInNeighborhoodGenerically< SenderType >( neighborID ) );

            if ( getLocalCommunicationMode() != DIRECT ||
                 !storage->primitiveExistsLocallyGenerically< SenderType >( neighborID ) )
            {
               uint_t neighborRank = storage->getPrimitiveRank( neighborID );

               if ( ranksToReceiveFrom.count( neighborRank ) == 0 )
               {
                  ranksToReceiveFrom[neighborRank] = 1;
               }
               else
               {
                  ranksToReceiveFrom[neighborRank] += 1;
               }
            }
         }
      }

      for ( auto it = sendFunctionsMap.begin(); it != sendFunctionsMap.end(); it++ )
      {
         uint_t                      receiverRank  = it->first;
         std::vector< SendFunction > sendFunctions = it->second;

         auto sendFunction = [sendFunctions]( SendBuffer& sendBuffer ) -> void {
            for ( auto& f : sendFunctions )
               f( sendBuffer );
         };

         bufferSystem->addSendingFunction( int_c( receiverRank ), sendFunction );
      }

      for ( const auto rankToReceiveFrom : ranksToReceiveFrom )
      {
         const uint_t senderRank       = rankToReceiveFrom.first;
         const uint_t numberOfMessages = rankToReceiveFrom.second;

         auto recvFunction = [this, numberOfMessages, storage]( RecvBuffer& recvBuffer ) -> void {
            for ( uint_t message = 0; message < numberOfMessages; message++ )
            {
               PrimitiveID senderID;
               PrimitiveID receiverID;
               readHeader( recvBuffer, senderID, receiverID );

               WALBERLA_ASSERT_NOT_NULLPTR( storage.get() );
               WALBERLA_ASSERT( storage->primitiveExistsLocallyGenerically< ReceiverType >( receiverID ) );

               ReceiverType* receiver = storage->getPrimitiveGenerically< ReceiverType >( receiverID );

               for ( const auto& packInfo : packInfos_ )
               {
                  startTimer( "Unpacking" );
                  packInfo->unpack< SenderType, ReceiverType >( receiver, senderID, recvBuffer );
                  stopTimer( "Unpacking" );
               }
            }
         };

         bufferSystem->addReceivingFunction( int_c( senderRank ), recvFunction );
      }

      setupBeforeNextCommunication_[communicationDirection] = false;
   } // setup

   stopTimer( timerStringSetup );

   // Buffered communication
   startTimer( timerStringBuffered );
   bufferSystem->startCommunication();
   stopTimer( timerStringBuffered );

   // Local communication
   startTimer( timerStringDirect );
   for ( auto& directCommunicationFunction : directCommunicationFunctions_[communicationDirection] )
   {
      directCommunicationFunction();
   }
   stopTimer( timerStringDirect );
}

template void BufferedCommunicator::startCommunication< Vertex, Edge >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Vertex, Cell >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Edge, Vertex >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Edge, Face >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Edge, Cell >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Face, Vertex >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Face, Edge >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Face, Face >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Face, Cell >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Cell, Vertex >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Cell, Edge >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Cell, Face >( std::vector< PrimitiveID > );

template void BufferedCommunicator::startCommunication< Cell, Cell >( std::vector< PrimitiveID > );

template < typename SenderType, typename ReceiverType >
void BufferedCommunicator::endCommunication()
{
   staticAssertCommunicationDirections< SenderType, ReceiverType >();

   const CommunicationDirection communicationDirection = getCommunicationDirection< SenderType, ReceiverType >();
   const std::string            timerString            = "Communication (buffered / wait + unpack)";

   startTimer( timerString );

   if ( packInfos_.empty() )
   {
      stopTimer( timerString );
      return;
   }

   WALBERLA_ASSERT( communicationInProgress_[communicationDirection] );
   communicationInProgress_[communicationDirection] = false;

   std::shared_ptr< walberla::mpi::OpenMPBufferSystem > bufferSystem = bufferSystems_[communicationDirection];
   bufferSystem->wait();

   stopTimer( timerString );
}

template void BufferedCommunicator::endCommunication< Vertex, Edge >();

template void BufferedCommunicator::endCommunication< Vertex, Cell >();

template void BufferedCommunicator::endCommunication< Edge, Vertex >();

template void BufferedCommunicator::endCommunication< Edge, Face >();

template void BufferedCommunicator::endCommunication< Edge, Cell >();

template void BufferedCommunicator::endCommunication< Face, Vertex >();

template void BufferedCommunicator::endCommunication< Face, Edge >();

template void BufferedCommunicator::endCommunication< Face, Face >();

template void BufferedCommunicator::endCommunication< Face, Cell >();

template void BufferedCommunicator::endCommunication< Cell, Vertex >();

template void BufferedCommunicator::endCommunication< Cell, Edge >();

template void BufferedCommunicator::endCommunication< Cell, Face >();

template void BufferedCommunicator::endCommunication< Cell, Cell >();

template < typename SenderType, typename ReceiverType >
void BufferedCommunicator::staticAssertCommunicationDirections() const
{
   static_assert( ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Edge >::value ) ||
                      ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Face >::value ) ||
                      ( std::is_same< SenderType, Vertex >::value && std::is_same< ReceiverType, Cell >::value )

                      || ( std::is_same< SenderType, Edge >::value && std::is_same< ReceiverType, Vertex >::value ) ||
                      ( std::is_same< SenderType, Edge >::value && std::is_same< ReceiverType, Face >::value ) ||
                      ( std::is_same< SenderType, Edge >::value && std::is_same< ReceiverType, Cell >::value )

                      || ( std::is_same< SenderType, Face >::value && std::is_same< ReceiverType, Vertex >::value ) ||
                      ( std::is_same< SenderType, Face >::value && std::is_same< ReceiverType, Edge >::value ) ||
                      ( std::is_same< SenderType, Face >::value && std::is_same< ReceiverType, Face >::value ) ||
                      ( std::is_same< SenderType, Face >::value && std::is_same< ReceiverType, Cell >::value )

                      || ( std::is_same< SenderType, Cell >::value && std::is_same< ReceiverType, Vertex >::value ) ||
                      ( std::is_same< SenderType, Cell >::value && std::is_same< ReceiverType, Edge >::value ) ||
                      ( std::is_same< SenderType, Cell >::value && std::is_same< ReceiverType, Face >::value ) ||
                      ( std::is_same< SenderType, Cell >::value && std::is_same< ReceiverType, Cell >::value ),

                  "BufferedCommunicator: illegal sender and receiver type combination." );
}
} // namespace communication
} // namespace hyteg
