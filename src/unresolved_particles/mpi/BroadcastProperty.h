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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <unresolved_particles/data/DataTypes.h>
#include <unresolved_particles/data/Flags.h>
#include <unresolved_particles/data/ParticleStorage.h>

#include <core/mpi/BufferSystem.h>
#include <core/logging/Logging.h>

#include <type_traits>

namespace walberla {
namespace unresolved_particles {
namespace mpi {

/**
 * Broadcast a property from the master particle to all corresponding ghost particles.
 *
 * \par Usage:
 * The property which will be broadcasted can be selected by the Notification template
 * (see ForceTorqueNotification).
 * void update(data::Particle&& p, const Notification::Parameters& objparam)
 * will be called to update the ghost particles with the information transmitted.
 *
 * \post
 * - all ghost particles got updated
 *
 * \ingroup unresolved_particles_mpi
 */
class BroadcastProperty
{
public:
   template <typename Notification>
   void operator()(data::ParticleStorage& ps) const;

   int64_t getBytesSent() const { return bytesSent_; }
   int64_t getBytesReceived() const { return bytesReceived_; }

   int64_t getNumberOfSends() const { return numberOfSends_; }
   int64_t getNumberOfReceives() const { return numberOfReceives_; }
private:
   int numProcesses_ = walberla::mpi::MPIManager::instance()->numProcesses();

   mutable int64_t bytesSent_ = 0;
   mutable int64_t bytesReceived_ = 0;
   mutable int64_t numberOfSends_ = 0;
   mutable int64_t numberOfReceives_ = 0;
};

template <typename Notification>
void BroadcastProperty::operator()(data::ParticleStorage& ps) const
{
   if (numProcesses_ == 1) return;

   walberla::mpi::BufferSystem bs(walberla::mpi::MPIManager::instance()->comm());

   std::set<int> recvRanks; // potential message senders

   WALBERLA_LOG_DETAIL( "Assembling of property reduction message starts...");

   for( auto p : ps )
   {
      if (data::particle_flags::isSet( p.getFlags(), data::particle_flags::GHOST))
      {
         // Will receive message from the particles owner
         recvRanks.insert(p.getOwner());
      } else
      {
         //local particles should send the property to all ghost particles
         for (auto& ghostRank : p.getGhostOwners())
         {
            auto& sb = bs.sendBuffer(ghostRank);
            if (sb.isEmpty())
            {
               // fill empty buffers with a dummy byte to force transmission
               sb << walberla::uint8_c(0);
            }
            sb << Notification( p );
         }
      }
   }

   WALBERLA_LOG_DETAIL( "Assembling of property broadcasting message ended." );

   bs.setReceiverInfo(recvRanks, true);
   bs.sendAll();

   // Receiving the updates for the remote rigid bodies from the connected processes
   WALBERLA_LOG_DETAIL( "Parsing of property broadcasting message starts..." );
   for( auto it = bs.begin(); it != bs.end(); ++it )
   {
      walberla::uint8_t tmp;
      it.buffer() >> tmp;
      while( !it.buffer().isEmpty() )
      {
         typename Notification::Parameters objparam;
         it.buffer() >> objparam;

         WALBERLA_LOG_DETAIL( "Received reduction notification from neighboring process with rank " << it.rank() );

         auto pIt = ps.find( objparam.uid_ );
         WALBERLA_CHECK_UNEQUAL( pIt, ps.end() );

         update(*pIt, objparam);

         WALBERLA_LOG_DETAIL( "Processed broadcasting notification for particle " << objparam.uid_ << "."  );
      }
   }
   WALBERLA_LOG_DETAIL( "Parsing of property broadcasting message ended." );

   bytesSent_ = bs.getBytesSent();
   bytesReceived_ = bs.getBytesReceived();
   numberOfSends_ = bs.getNumberOfSends();
   numberOfReceives_ = bs.getNumberOfReceives();
}

}  // namespace mpi
}  // namespace unresolved_particles
}  // namespace walberla