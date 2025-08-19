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

#pragma once
#include <functional>
#include <vector>

#include "core/DataTypes.h"

#include "hyteg/types/BufferSystemForwardDeclare.hpp"
#include "hyteg/primitives/PrimitiveID.hpp"

namespace walberla {
namespace timing {
struct WcPolicy;
template < typename TP >
class TimingTree;
} // namespace timing
using WcTimingTree = timing::TimingTree< timing::WcPolicy >;
} // namespace walberla

namespace hyteg {
class PrimitiveStorage;
// class PrimitiveID;

namespace communication {

class PackInfo;
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

   ~BufferedCommunicator();

   /// All data that are registered via respective \ref PackInfo objects are exchanged
   void addPackInfo( const std::shared_ptr< PackInfo >& packInfo );

   /// Starts the non-blocking communication between two \ref Primitive types.
   /// The data of the sender can be modified after this method returns.
   /// \tparam SenderType type of the sending \ref Primitive (e.g. \ref Vertex or \ref Edge)
   /// \tparam ReceiverType type of the receiving \ref Primitive (e.g. \ref Vertex or \ref Edge)
   /// \param excludeReceivingIDs exclude primtives with these IDs from receiving. The primitives will still send their data
   template < typename SenderType, typename ReceiverType >
   void startCommunication( std::vector< PrimitiveID > excludeReceivingIDs = {} );

   /// Ends the non-blocking communication between two \ref Primitive types
   /// Waits for the started communication to be completed. It is only safe to modify the
   /// data of the receiver after this call returned.
   /// \tparam SenderType type of the sending \ref Primitive (e.g. \ref Vertex or \ref Edge)
   /// \tparam ReceiverType type of the receiving \ref Primitive (e.g. \ref Vertex or \ref Edge)
   template < typename SenderType, typename ReceiverType >
   void endCommunication();

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
   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree );

 private:
   typedef std::function< void( walberla::mpi::SendBuffer& buf ) > SendFunction;
   typedef std::function< void( walberla::mpi::RecvBuffer& buf ) > RecvFunction;

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

   static const walberla::uint_t SYNC_WORD;

   static const std::array< std::string, CommunicationDirection::NUM_COMMUNICATION_DIRECTIONS >  COMMUNICATION_DIRECTION_STRINGS;
   static const std::array< std::string, LocalCommunicationMode::NUM_LOCAL_COMMUNICATION_MODES > LOCAL_COMMUNICATION_MODE_STRINGS;

   template < typename SenderType, typename ReceiverType >
   CommunicationDirection getCommunicationDirection() const;

   void writeHeader( walberla::mpi::SendBuffer& sendBuffer, const PrimitiveID& senderID, const PrimitiveID& receiverID );

   void readHeader( walberla::mpi::RecvBuffer& recvBuffer, PrimitiveID& senderID, PrimitiveID& receiverID );

   void endCommunication( const CommunicationDirection& communicationDirection );

   void startTimer( const std::string& timerString );
   void stopTimer( const std::string& timerString );

   void setupBeforeNextCommunication();

   template < typename SenderType, typename ReceiverType >
   void staticAssertCommunicationDirections() const;

   std::weak_ptr< PrimitiveStorage > primitiveStorage_;

   walberla::uint_t primitiveStorageModificationStamp_;

   std::vector< std::shared_ptr< PackInfo > > packInfos_;

   std::array< std::shared_ptr< walberla::mpi::OpenMPBufferSystem >, NUM_COMMUNICATION_DIRECTIONS > bufferSystems_;

   std::array< bool, NUM_COMMUNICATION_DIRECTIONS > communicationInProgress_;

   LocalCommunicationMode localCommunicationMode_;

   // Cached communication setup
   std::array< bool, NUM_COMMUNICATION_DIRECTIONS >                                   setupBeforeNextCommunication_;
   std::array< std::vector< std::function< void() > >, NUM_COMMUNICATION_DIRECTIONS > directCommunicationFunctions_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

   // Claimed MPI Tags
   std::vector< int > claimedMPITags_;
};

} // namespace communication
} // namespace hyteg
