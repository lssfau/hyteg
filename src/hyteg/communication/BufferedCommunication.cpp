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

#include "hyteg/communication/BufferedCommunication.hpp"

#include <functional>

#include "core/logging/Logging.h"

#include "hyteg/communication/MPITagProvider.hpp"

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
      bufferSystem = std::shared_ptr< walberla::mpi::OpenMPBufferSystem >( new walberla::mpi::OpenMPBufferSystem(
          walberla::mpi::MPIManager::instance()->comm(), MPITagProvider::getMPITag(), serialSends, serialRecvs ) );
   }

   setupBeforeNextCommunication();

   for ( auto& communicationInProgress : communicationInProgress_ )
   {
      communicationInProgress = false;
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

} // namespace communication
} // namespace hyteg
