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
//! \file SyncNextNeighborsNoGhosts.h
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
#include <convection_particles/data/Flags.h>
#include <convection_particles/data/ParticleStorage.h>
#include <convection_particles/domain/IDomain.h>
#include <convection_particles/mpi/notifications/PackNotification.h>
#include <convection_particles/mpi/notifications/ParseMessage.h>
#include <convection_particles/mpi/notifications/ParticleGhostCopyNotification.h>
#include <convection_particles/mpi/notifications/ParticleMigrationNotification.h>
#include <convection_particles/mpi/notifications/ParticleRemoteMigrationNotification.h>
#include <convection_particles/mpi/notifications/ParticleRemovalNotification.h>
#include <convection_particles/mpi/notifications/ParticleUpdateNotification.h>

#include <core/mpi/BufferSystem.h>
#include <core/logging/Logging.h>

namespace walberla {
namespace convection_particles {
namespace mpi {

/**
 * Next neighbor synchronization kernel without ghost particles.
 *
 * This kernel conducts particle migrations if they move to another subdomain
 * but does not generate ghost particles on overlap. This can be useful for
 * particles without spatial extend like tracer particles.
 *
 * \ingroup convection_particles_mpi
 */
class SyncNextNeighborsNoGhosts
{
public:
   void operator()(data::ParticleStorage& ps,
                   const domain::IDomain& domain) const;

   int64_t getBytesSent() const { return bs.getBytesSent(); }
   int64_t getBytesReceived() const { return bs.getBytesReceived(); }

   int64_t getNumberOfSends() const { return bs.getNumberOfSends(); }
   int64_t getNumberOfReceives() const { return bs.getNumberOfReceives(); }
private:
   void generateSynchronizationMessages(data::ParticleStorage& ps,
                                        const domain::IDomain& domain) const;
   mutable std::vector<uint_t> neighborRanks_; ///cache for neighbor ranks -> will be updated in operator()

   mutable walberla::mpi::BufferSystem bs = walberla::mpi::BufferSystem( walberla::mpi::MPIManager::instance()->comm() );

   int numProcesses_ = walberla::mpi::MPIManager::instance()->numProcesses();
   int rank_         = walberla::mpi::MPIManager::instance()->rank();
};

}  // namespace mpi
}  // namespace convection_particles
}  // namespace walberla