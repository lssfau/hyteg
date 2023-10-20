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
#include <unresolved_particles/domain/IDomain.h>
#include <unresolved_particles/mpi/notifications/NewGhostParticleNotification.h>
#include <unresolved_particles/mpi/notifications/PackNotification.h>
#include <unresolved_particles/mpi/notifications/ParseMessage.h>
#include <unresolved_particles/mpi/notifications/ParticleGhostCopyNotification.h>
#include <unresolved_particles/mpi/notifications/ParticleMigrationNotification.h>
#include <unresolved_particles/mpi/notifications/ParticleRemoteMigrationNotification.h>
#include <unresolved_particles/mpi/notifications/ParticleRemovalInformationNotification.h>
#include <unresolved_particles/mpi/notifications/ParticleRemovalNotification.h>
#include <unresolved_particles/mpi/notifications/ParticleUpdateNotification.h>

#include <core/mpi/BufferSystem.h>
#include <core/logging/Logging.h>

namespace walberla {
namespace unresolved_particles {
namespace mpi {

/**
 * Kernel which updates all ghost particles.
 *
 * \ingroup unresolved_particles_mpi
 */
class SyncGhostOwners
{
public:
   void operator()( data::ParticleStorage& ps,
                    const domain::IDomain& domain,
                    const real_t dx = real_t(0),
                    const bool syncNonCommunicatingBodies = false ) const;

   int64_t getBytesSent() const { return bytesSent_; }
   int64_t getBytesReceived() const { return bytesReceived_; }

   int64_t getNumberOfSends() const { return numberOfSends_; }
   int64_t getNumberOfReceives() const { return numberOfReceives_; }
private:
   void updateAndMigrate( data::ParticleStorage& ps,
                          const domain::IDomain& domain,
                          const bool syncNonCommunicatingBodies ) const;

   void checkAndResolveOverlap( data::ParticleStorage& ps,
                                const domain::IDomain& domain,
                                const real_t dx,
                                const bool syncNonCommunicatingBodies ) const;

   mutable std::vector<uint_t> neighborRanks_; ///cache for neighbor ranks -> will be updated in operator()

   int numProcesses_ = walberla::mpi::MPIManager::instance()->numProcesses();
   int rank_         = walberla::mpi::MPIManager::instance()->rank();

   mutable int64_t bytesSent_ = 0;
   mutable int64_t bytesReceived_ = 0;
   mutable int64_t numberOfSends_ = 0;
   mutable int64_t numberOfReceives_ = 0;
};

}  // namespace mpi
}  // namespace unresolved_particles
}  // namespace walberla