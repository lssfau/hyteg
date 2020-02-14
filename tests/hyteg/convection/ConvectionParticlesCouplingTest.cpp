/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include <convection_particles/data/ParticleStorage.h>
#include <convection_particles/domain/BlockForestDomain.h>
#include <convection_particles/mpi/SyncNextNeighborsNoGhosts.h>

#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/mpi/MPIManager.h>

using namespace walberla;
using namespace walberla::convection_particles;

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   //init domain partitioning
   auto forest = blockforest::createBlockForest( math::AABB(0,0,0,10,10,10), // simulation domain
                                                 Vector3<uint_t>(2,2,2), // blocks in each direction
                                                 Vector3<bool>(false, false, false) // periodicity
                                                 );
   domain::BlockForestDomain domain(forest);

   data::ParticleStorage particleStorage(100);

   Vec3 pt(2.5, 2.5, 2.5);
   if (forest->begin()->getAABB().contains(pt))
   {
      auto pIt = particleStorage.create();
      pIt->setPosition(pt);
      pIt->setInteractionRadius(real_t(0));
      pIt->setOwner(walberla::mpi::MPIManager::instance()->rank());
   }

   for (auto p : particleStorage)
   {
      WALBERLA_LOG_DEVEL(p.getPosition());
      p.setPosition(Vec3(7.5, 2.5, 2.5));
   }

   convection_particles::mpi::SyncNextNeighborsNoGhosts SNN;

   SNN(particleStorage, domain);

   for (auto p : particleStorage)
   {
      WALBERLA_LOG_DEVEL(p.getPosition());
      p.setPosition(Vec3(7.5, 2.5, 2.5));
   }

   return 0;
}
