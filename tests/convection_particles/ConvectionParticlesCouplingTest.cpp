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

#include <blockforest/Initialization.h>
#include <convection_particles/data/ParticleStorage.h>
#include <convection_particles/domain/BlockForestDomain.h>
#include <convection_particles/mpi/SyncNextNeighborsNoGhosts.h>
#include <core/Environment.h>
#include <core/mpi/MPIManager.h>
#include <hyteg/geometry/Intersection.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>

#include "coupling_hyteg_convection_particles/primitivestorage/PrimitiveStorageConvectionParticlesInterface.hpp"

using namespace walberla;
using namespace walberla::convection_particles;
using namespace hyteg;

void testBlockForestCoupling()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Block forest coupling ..." )
   //init domain partitioning
   auto forest = blockforest::createBlockForest( walberla::math::AABB(0,0,0,10,10,10), // simulation domain
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
}

void testSetupPrimitiveStorageCoupling2D()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Primitive storage coupling 2D ..." )

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( 0, 0), Point2D( 10, 10), MeshInfo::CRISS, 3, 3 );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   PrimitiveStorageConvectionParticlesInterface domain( storage );

   data::ParticleStorage particleStorage(100);

   Vec3 pt(2.5, 2.5, 0);

   // get any cell
   if ( domain.findContainingProcessRank( pt ) == walberla::mpi::MPIManager::instance()->rank() )
   {
      WALBERLA_LOG_INFO( "Creating particle ..." )
      auto pIt = particleStorage.create();
      pIt->setPosition(pt);
      pIt->setInteractionRadius(real_t(0));
      pIt->setOwner(walberla::mpi::MPIManager::instance()->rank());
   }

   for (auto p : particleStorage)
   {
      WALBERLA_LOG_DEVEL(p.getPosition());
      p.setPosition(Vec3(7.5, 2.5, 0));
   }

   convection_particles::mpi::SyncNextNeighborsNoGhosts SNN;

   SNN(particleStorage, domain);
   SNN(particleStorage, domain);

   for (auto p : particleStorage)
   {
      auto isGhost = p.getFlags().getRawDataRef() == walberla::convection_particles::data::particle_flags::GHOST;
      WALBERLA_LOG_INFO( "Received particle: pos " << p.getPosition() << ", ghost " << isGhost )
      p.setPosition(Vec3(7.5, 2.5, 0));
   }
}


void testSetupPrimitiveStorageCoupling3D()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Primitive storage coupling 3D ..." )

   MeshInfo meshInfo = MeshInfo::meshCuboid( Point3D( 0, 0, 0), Point3D( 10, 10, 10), 1, 1, 1 );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   PrimitiveStorageConvectionParticlesInterface domain( storage );

   data::ParticleStorage particleStorage(100);

   Vec3 pt(2.5, 2.5, 2.5);

   // get any cell
   if ( domain.findContainingProcessRank( pt ) == walberla::mpi::MPIManager::instance()->rank() )
   {
      WALBERLA_LOG_INFO( "Creating particle ..." )
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
   SNN(particleStorage, domain);

   for (auto p : particleStorage)
   {
      auto isGhost = p.getFlags().getRawDataRef() == walberla::convection_particles::data::particle_flags::GHOST;
      WALBERLA_LOG_INFO( "Received particle: pos " << p.getPosition() << ", ghost " << isGhost )
      p.setPosition(Vec3(7.5, 2.5, 2.5));
   }
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   testBlockForestCoupling();
   testSetupPrimitiveStorageCoupling2D();
   testSetupPrimitiveStorageCoupling3D();

   return 0;
}