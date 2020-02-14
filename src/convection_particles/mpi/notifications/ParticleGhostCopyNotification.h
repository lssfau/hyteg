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
//! \file ParticleGhostCopyNotification.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <convection_particles/data/DataTypes.h>
#include <convection_particles/data/ParticleStorage.h>
#include <convection_particles/mpi/ShapePackUnpack.h>
#include <convection_particles/mpi/notifications/NotificationType.h>

#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace walberla {
namespace convection_particles {

/**
 * A complete particle copy for a new ghost particle.
 *
 * Copies all properties marked ON_GHOST_CREATION or ALWAYS.
 */
class ParticleGhostCopyNotification
{
public:
   struct Parameters
   {
      walberla::id_t uid {UniqueID<data::Particle>::invalidID()};
      walberla::convection_particles::Vec3 position {real_t(0)};
      walberla::real_t interactionRadius {real_t(0)};
      walberla::convection_particles::data::particle_flags::FlagT flags {};
      int owner {-1};
      walberla::convection_particles::Vec3 velocity {real_t(0)};
      walberla::convection_particles::Vec3 startPosition {real_t(0)};
      hyteg::indexing::Index startIndex {};
      uint_t startProcess {};
      hyteg::PrimitiveID startPrimitiveID {};
      uint_t startDoFType {};
      hyteg::edgedof::EdgeDoFOrientation startEdgeDoFOrientation {};
      std::vector< walberla::convection_particles::Vec3 > k {};
      real_t finalTemperature {};
   };

   inline explicit ParticleGhostCopyNotification( const data::Particle& particle ) : particle_(particle) {}
   const data::Particle& particle_;
};

inline data::ParticleStorage::iterator createNewParticle(data::ParticleStorage& ps, const ParticleGhostCopyNotification::Parameters& data)
{
   WALBERLA_ASSERT_EQUAL(ps.find(data.uid), ps.end(), "Particle with same uid already existent!");

   auto pIt = ps.create(data.uid);
   pIt->setUid(data.uid);
   pIt->setPosition(data.position);
   pIt->setInteractionRadius(data.interactionRadius);
   pIt->setFlags(data.flags);
   pIt->setOwner(data.owner);
   pIt->setVelocity(data.velocity);
   pIt->setStartPosition(data.startPosition);
   pIt->setStartIndex(data.startIndex);
   pIt->setStartProcess(data.startProcess);
   pIt->setStartPrimitiveID(data.startPrimitiveID);
   pIt->setStartDoFType(data.startDoFType);
   pIt->setStartEdgeDoFOrientation(data.startEdgeDoFOrientation);
   pIt->setK(data.k);
   pIt->setFinalTemperature(data.finalTemperature);
   return pIt;
}

template<>
struct NotificationTrait<ParticleGhostCopyNotification>
{
   static const NotificationType id = PARTICLE_GHOST_COPY_NOTIFICATION;
};

}  // namespace convection_particles
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const convection_particles::ParticleGhostCopyNotification& obj )
{
   buf.addDebugMarker( "cn" );
   buf << obj.particle_.getUid();
   buf << obj.particle_.getPosition();
   buf << obj.particle_.getInteractionRadius();
   buf << obj.particle_.getFlags();
   buf << obj.particle_.getOwner();
   buf << obj.particle_.getVelocity();
   buf << obj.particle_.getStartPosition();
   buf << obj.particle_.getStartIndex();
   buf << obj.particle_.getStartProcess();
   buf << obj.particle_.getStartPrimitiveID();
   buf << obj.particle_.getStartDoFType();
   buf << obj.particle_.getStartEdgeDoFOrientation();
   buf << obj.particle_.getK();
   buf << obj.particle_.getFinalTemperature();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, convection_particles::ParticleGhostCopyNotification::Parameters& objparam )
{
   buf.readDebugMarker( "cn" );
   buf >> objparam.uid;
   buf >> objparam.position;
   buf >> objparam.interactionRadius;
   buf >> objparam.flags;
   buf >> objparam.owner;
   buf >> objparam.velocity;
   buf >> objparam.startPosition;
   buf >> objparam.startIndex;
   buf >> objparam.startProcess;
   buf >> objparam.startPrimitiveID;
   buf >> objparam.startDoFType;
   buf >> objparam.startEdgeDoFOrientation;
   buf >> objparam.k;
   buf >> objparam.finalTemperature;
   return buf;
}

} // mpi
} // walberla