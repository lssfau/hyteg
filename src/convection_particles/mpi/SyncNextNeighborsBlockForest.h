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
//! \file SyncNextNeighborsBlockForest.h
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
#include <convection_particles/domain/BlockForestDomain.h>
#include <convection_particles/mpi/notifications/PackNotification.h>
#include <convection_particles/mpi/notifications/ParseMessage.h>
#include <convection_particles/mpi/notifications/ParticleGhostCopyNotification.h>
#include <convection_particles/mpi/notifications/ParticleMigrationNotification.h>
#include <convection_particles/mpi/notifications/ParticleRemoteMigrationNotification.h>
#include <convection_particles/mpi/notifications/ParticleRemovalNotification.h>
#include <convection_particles/mpi/notifications/ParticleUpdateNotification.h>

#include <blockforest/BlockForest.h>

#include <core/mpi/BufferSystem.h>
#include <core/logging/Logging.h>

namespace walberla {
namespace convection_particles {
namespace mpi {

/**
 * Kernel which updates all ghost particles.
 *
 * \ingroup convection_particles_mpi
 */
class SyncNextNeighborsBlockForest
{
public:
   void operator()(data::ParticleStorage& ps,
                   const std::shared_ptr<blockforest::BlockForest>& blockforest,
                   const std::shared_ptr<domain::BlockForestDomain>& domain,
                   const real_t dx = real_t(0)) const;

   int64_t getBytesSent() const { return bs.getBytesSent(); }
   int64_t getBytesReceived() const { return bs.getBytesReceived(); }

   int64_t getNumberOfSends() const { return bs.getNumberOfSends(); }
   int64_t getNumberOfReceives() const { return bs.getNumberOfReceives(); }
private:
   void generateSynchronizationMessages(data::ParticleStorage& ps,
                                        const std::shared_ptr<blockforest::BlockForest>& blockforest,
                                        const real_t dx) const;

   mutable walberla::mpi::BufferSystem bs = walberla::mpi::BufferSystem( walberla::mpi::MPIManager::instance()->comm() );

   int numProcesses_ = walberla::mpi::MPIManager::instance()->numProcesses();
   int rank_         = walberla::mpi::MPIManager::instance()->rank();
};

}  // namespace mpi
}  // namespace convection_particles
}  // namespace walberla
