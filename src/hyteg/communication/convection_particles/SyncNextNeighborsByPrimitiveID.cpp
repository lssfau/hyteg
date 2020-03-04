//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file SyncNextNeighborsNoGhosts.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#include "SyncNextNeighborsByPrimitiveID.h"

#include <convection_particles/mpi/RemoveAndNotify.h>

namespace walberla {
namespace convection_particles {
namespace mpi {

void SyncNextNeighborsByPrimitiveID::operator()(data::ParticleStorage& ps,
                                                const hyteg::SetupPrimitiveStorage & setupStorage) const
{
   if (numProcesses_ == 1) return;

   for( uint_t nbProcessRank = 0; nbProcessRank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); nbProcessRank++ )
   {
      if (bs.sendBuffer(nbProcessRank).isEmpty())
      {
         // fill empty buffers with a dummy byte to force transmission
         bs.sendBuffer(nbProcessRank) << walberla::uint8_c(0);
      }
   }
   generateSynchronizationMessages(ps, setupStorage);

   // size of buffer is unknown and changes with each send
   bs.setReceiverInfoFromSendBufferState(false, true);
   bs.sendAll();

   // Receiving the updates for the remote rigid bodies from the connected processes
   WALBERLA_LOG_DETAIL( "Parsing of particle synchronization response starts..." );
   ParseMessagePrimitiveIDCommunication parseMessage;
   for( auto it = bs.begin(); it != bs.end(); ++it )
   {
      walberla::uint8_t tmp;
      it.buffer() >> tmp;
      while( !it.buffer().isEmpty() )
      {
         parseMessage(it.rank(), it.buffer(), ps);
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of particle synchronization response ended." );
}

void SyncNextNeighborsByPrimitiveID::generateSynchronizationMessages(data::ParticleStorage& ps,
                                                                     const hyteg::SetupPrimitiveStorage & setupStorage) const
{
   const uint_t ownRank = uint_c(rank_);

   WALBERLA_LOG_DETAIL( "Assembling of particle synchronization message starts..." );

   // position update
   for( auto pIt = ps.begin(); pIt != ps.end(); )
   {
      //skip all ghost particles
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::GHOST))
      {
         ++pIt;
         continue;
      }

      //skip all particles that do not communicate (create ghost particles) on other processes
      if (data::particle_flags::isSet( pIt->getFlags(), data::particle_flags::NON_COMMUNICATING))
      {
         ++pIt;
         continue;
      }

      // Note: At this point we know that the particle was locally owned before the position update.
      WALBERLA_CHECK_EQUAL(pIt->getOwner(), ownRank);

      WALBERLA_LOG_DETAIL( "Processing local particle " << pIt->getUid() );

      //particle has left subdomain?
      const auto ownerRank = setupStorage.getTargetRank( pIt->getContainingPrimitive() );
      if( ownerRank != ownRank )
      {
         WALBERLA_LOG_DETAIL( "Local particle " << pIt->getUid() << " is no longer on process " << ownRank << " but on process " << ownerRank );

         // create ghost on new owner process
         auto& buffer( bs.sendBuffer(ownerRank) );
         WALBERLA_LOG_DETAIL( "Sending ghost copy notification for particle " << pIt->getUid() << " to process " << ownerRank );
         packNotification(buffer, ParticleGhostCopyNotification( *pIt ));

         WALBERLA_LOG_DETAIL( "Sending migration notification for particle " << pIt->getUid() << " to process " << ownerRank << "." );
         //WALBERLA_LOG_DETAIL( "Process registration list before migration: " << pIt->getGhostOwners() );

         // Set new owner and transform to ghost particle
         pIt->setOwner(int_c( ownerRank ));
         data::particle_flags::set( pIt->getFlagsRef(), data::particle_flags::GHOST );

         // Send migration notification to new owner
         packNotification(buffer, ParticleMigrationNotification( *pIt ));

         //remove particle from local process
         pIt = ps.erase( pIt );

         continue;

      } else
      {
         // particle still is locally owned after position update.
         WALBERLA_LOG_DETAIL( "Owner of particle " << pIt->getUid() << " is still process " << pIt->getOwner() );
      }

      ++pIt;
   }

   WALBERLA_LOG_DETAIL( "Assembling of particle synchronization message ended." );
}

}  // namespace mpi
}  // namespace convection_particles
}  // namespace walberla