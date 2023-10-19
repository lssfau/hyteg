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

#include <unresolved_particles/data/ParticleStorage.h>
#include <unresolved_particles/mpi/notifications/NotificationType.h>

#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace unresolved_particles {

/**
 * The ownership for one of the ghost particles has changed.
 */
class ParticleRemoteMigrationNotification {
public:
   struct Parameters {
      id_t uid_;
      int newOwner_;
   };

   inline ParticleRemoteMigrationNotification( const data::Particle& particle, const int& newOwner )
      : particle_(particle), newOwner_(newOwner) {}
   const data::Particle& particle_;
   const int newOwner_;
};

template<>
struct NotificationTrait<ParticleRemoteMigrationNotification>
{
   static const NotificationType id = PARTICLE_REMOTE_MIGRATION_NOTIFICATION;
};

}  // namespace unresolved_particles
}  // namespace walberla

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T,    // Element type of SendBuffer
          typename G>    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const unresolved_particles::ParticleRemoteMigrationNotification& obj )
{
   buf.addDebugMarker( "rm" );
   buf << obj.particle_.getUid();
   buf << obj.newOwner_;
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, unresolved_particles::ParticleRemoteMigrationNotification::Parameters& objparam )
{
   buf.readDebugMarker( "rm" );
   buf >> objparam.uid_;
   buf >> objparam.newOwner_;
   return buf;
}

template <>
struct BufferSizeTrait< unresolved_particles::ParticleRemoteMigrationNotification > {
   static const bool constantSize = true;
   static const uint_t size = BufferSizeTrait<id_t>::size + BufferSizeTrait<int>::size + mpi::BUFFER_DEBUG_OVERHEAD;
};

} // mpi
} // walberla