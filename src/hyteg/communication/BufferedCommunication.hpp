/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#pragma once

#include <atomic>
#include <limits>

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/OpenMPBufferSystem.h"
#include "core/timing/TimingPool.h"
#include "core/timing/TimingTree.h"

#include "hyteg/communication/PackInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {
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
   /// Options for the communication mode that is used between primitives that belong to the same process
   enum LocalCommunicationMode
   {
      /// Uses the direct communication callbacks of the respective PackInfos for local neighbors
      DIRECT,
      /// Sends data to local neighbors over MPI
      BUFFERED_MPI,
      /// Number of differed modes
      NUM_LOCAL_COMMUNICATION_MODES
   };

   BufferedCommunicator( std::weak_ptr< PrimitiveStorage > primitiveStorage,
                         const LocalCommunicationMode&     localCommunicationMode = DIRECT );

   /// All data that are registered via respective \ref PackInfo objects are exchanged
   void addPackInfo( const std::shared_ptr< PackInfo >& packInfo );

   /// Starts the non-blocking communication between two \ref Primitive types.
   /// The data of the sender can be modified after this method returns.
   /// \tparam SenderType type of the sending \ref Primitive (e.g. \ref Vertex or \ref Edge)
   /// \tparam ReceiverType type of the receiving \ref Primitive (e.g. \ref Vertex or \ref Edge)
   /// \param excludeReceivingIDs exclude primtives with these IDs from receiving. The primitives will still send their data
   template < typename SenderType, typename ReceiverType >
   inline void startCommunication( std::vector< PrimitiveID > excludeReceivingIDs = {} );

   /// Ends the non-blocking communication between two \ref Primitive types
   /// Waits for the started communication to be completed. It is only safe to modify the
   /// data of the receiver after this call returned.
   /// \tparam SenderType type of the sending \ref Primitive (e.g. \ref Vertex or \ref Edge)
   /// \tparam ReceiverType type of the receiving \ref Primitive (e.g. \ref Vertex or \ref Edge)
   template < typename SenderType, typename ReceiverType >
   inline void endCommunication();

   /// @name Local communication mode
   /// Getter and setter to retrieve and set the local communication mode.
   /// Depending on the mode, data is transferred via different mechanisms if the sending
   /// and the receiving \ref Primitive are located on the same process.
   /// See also \ref LocalCommunicationMode
   ///@{
   LocalCommunicationMode getLocalCommunicationMode() const { return localCommunicationMode_; }
   void                   setLocalCommunicationMode( const LocalCommunicationMode& localCommunicationMode );
   ///@}

   /// Writes timing data for the setup and for the wait phase to \p timingTree
   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree ) { timingTree_ = timingTree; }

 private:
   typedef std::function< void( SendBuffer& buf ) > SendFunction;
   typedef std::function< void( RecvBuffer& buf ) > RecvFunction;

   enum CommunicationDirection
   {
      VERTEX_TO_EDGE,
      VERTEX_TO_FACE,
      VERTEX_TO_CELL,

      EDGE_TO_VERTEX,
      EDGE_TO_FACE,
      EDGE_TO_CELL,

      FACE_TO_VERTEX,
      FACE_TO_EDGE,
      FACE_TO_FACE,
      FACE_TO_CELL,

      CELL_TO_VERTEX,
      CELL_TO_EDGE,
      CELL_TO_FACE,
      CELL_TO_CELL,

      NUM_COMMUNICATION_DIRECTIONS
   };

   static const uint_t SYNC_WORD;

   static const std::array< std::string, CommunicationDirection::NUM_COMMUNICATION_DIRECTIONS >  COMMUNICATION_DIRECTION_STRINGS;
   static const std::array< std::string, LocalCommunicationMode::NUM_LOCAL_COMMUNICATION_MODES > LOCAL_COMMUNICATION_MODE_STRINGS;

   template < typename SenderType, typename ReceiverType >
   inline CommunicationDirection getCommunicationDirection() const;

   void writeHeader( SendBuffer& sendBuffer, const PrimitiveID& senderID, const PrimitiveID& receiverID );
   void readHeader( RecvBuffer& recvBuffer, PrimitiveID& senderID, PrimitiveID& receiverID );

   void endCommunication( const CommunicationDirection& communicationDirection );

   void startTimer( const std::string& timerString );
   void stopTimer( const std::string& timerString );

   void setupBeforeNextCommunication();

   template < typename SenderType, typename ReceiverType >
   inline void staticAssertCommunicationDirections() const;

   std::weak_ptr< PrimitiveStorage > primitiveStorage_;

   uint_t primitiveStorageModificationStamp_;

   std::vector< std::shared_ptr< PackInfo > > packInfos_;

   std::array< std::shared_ptr< walberla::mpi::OpenMPBufferSystem >, NUM_COMMUNICATION_DIRECTIONS > bufferSystems_;

   std::array< bool, NUM_COMMUNICATION_DIRECTIONS > communicationInProgress_;

   LocalCommunicationMode localCommunicationMode_;

   // Cached communication setup
   std::array< bool, NUM_COMMUNICATION_DIRECTIONS >                                   setupBeforeNextCommunication_;
   std::array< std::vector< std::function< void() > >, NUM_COMMUNICATION_DIRECTIONS > directCommunicationFunctions_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;
};

template < typename SenderType, typename ReceiverType >
inline BufferedCommunicator::CommunicationDirection BufferedCommunicator::getCommunicationDirection() const
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

         auto recvFunction = [this, numberOfMessages]( RecvBuffer& recvBuffer ) -> void {
            for ( uint_t message = 0; message < numberOfMessages; message++ )
            {
               PrimitiveID senderID;
               PrimitiveID receiverID;
               readHeader( recvBuffer, senderID, receiverID );

               std::shared_ptr< PrimitiveStorage > storage = primitiveStorage_.lock();

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

// Auxilliary class to provide MPI tag values for BufferedCommunicator
class MPITagProvider
{
   // Maximal tag value supported by the MPI library implementation used
   //
   // The MPI 4.1 standard requires the largest tag to be at least 2^15-1 = 32,767.
   // Larger values are possible, though. Since the value must be an int, it cannot
   // exceed 2,147,483,647 for the standard 32-bit signed int setting.
   static int maxMPITag_;

   // Stores the next tag that will be returned by getMPITag()
   static std::atomic_int nextMPITag_;

   // Marks whether class can still provide tag values
   static bool poolExhausted_;

 public:
   // Return the largest possible tag value supported by the MPI library in use
   static int getMaxMPITag()
   {
#ifdef WALBERLA_BUILD_WITH_MPI
      void* maxTag;
      int   status;
      MPI_Comm_get_attr( walberla::mpi::MPIManager::instance()->comm(), MPI_TAG_UB, &maxTag, &status );
      if ( status == 0 )
      {
         WALBERLA_ABORT( "Failed to query maximal tag value from MPI implementation!" );
      }
      return *static_cast< int* >( maxTag );
#else
      return std::numeric_limits< int >::max();
#endif
   }

   // Return another MPI tag value
   //
   // The current implementation is very simple. In order to return unique tag values it starts with
   // the smallest possbile value, i.e. 0, and then returns tags by incrementation until reaching the
   // limit. Thus, there is no re-use of values that are no longer needed, and the pool of tags might
   // get exhausted. In this case the class calls WALBERLA_ABORT().
   static int getMPITag()
   {
      // initialise largest available tag value (can only happen once MPI was activated)
      if ( maxMPITag_ == 0u )
      {
         maxMPITag_ = getMaxMPITag();
      }

      if ( poolExhausted_ )
      {
         WALBERLA_ABORT( "Your application exhausted the pool of available MPI tags.\n"
                         << "Your MPI implementation provides a maximum of " << maxMPITag_ << " tags." );
      }

      // this will store the return tag
      int freshTag{ nextMPITag_ };

      // check whether we can safely increase the tag
      if ( nextMPITag_ == maxMPITag_ )
      {
         poolExhausted_ = true;
      }
      else
      {
         ++nextMPITag_;
      }

      return freshTag;
   }
};

} // namespace communication
} // namespace hyteg
