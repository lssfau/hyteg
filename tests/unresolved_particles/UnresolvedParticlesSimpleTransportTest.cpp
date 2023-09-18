/*
* Copyright (c) 2017-2023 Nils Kohl.
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
#include <core/Environment.h>
#include <core/math/Random.h>
#include <core/mpi/MPIManager.h>
#include <hyteg/geometry/Intersection.hpp>
#include <hyteg/p1functionspace/P1VectorFunction.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include <hyteg/primitivestorage/Visualization.hpp>

#include "coupling_hyteg_unresolved_particles/UnresolvedParticles.hpp"

using namespace walberla;
using namespace walberla::unresolved_particles;
using walberla::unresolved_particles::Vec3;
using namespace hyteg;

void simpleTransportTest()
{
   real_t cubeSize             = 10;
   uint_t tunnelStretch        = 5;
   real_t particleInitCubeSize = cubeSize / 2;
   uint_t timesteps            = 1000;
   real_t dt                   = 1e-2;
   uint_t level                = 3;

   MeshInfo meshInfo = MeshInfo::meshCuboid(
       Point3D( 0, 0, 0 ), Point3D( real_c( tunnelStretch ) * cubeSize, cubeSize, cubeSize ), tunnelStretch, 1, 1 );
   auto setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   writeDomainPartitioningVTK( *storage, "vtk", "domain" );

   hyteg::unresolved_particles::UnresolvedParticles unresolvedParticles( storage );

   uint_t numParticles = 100;
   for ( uint_t i = 0; i < numParticles; i++ )
   {
      Point3D pos( ( walberla::math::realRandom() - 0.5 ) * particleInitCubeSize,
                   ( walberla::math::realRandom() - 0.5 ) * particleInitCubeSize,
                   ( walberla::math::realRandom() - 0.5 ) * particleInitCubeSize );
      pos += Point3D( 0.5 * cubeSize, 0.5 * cubeSize, 0.5 * cubeSize );

      auto particleCreated = unresolvedParticles.createParticle( pos );
      if ( particleCreated )
      {
         real_t massAndSize = std::clamp( walberla::math::realRandom(), 0.5, 1.0 );

         auto particle = particleCreated.value();

         particle->setInteractionRadius( massAndSize );
         particle->setInvMass( 1 / massAndSize );
      }
   }

   unresolvedParticles.initVTK();

   P1VectorFunction< real_t > forceField( "forceField", storage, 3, 3 );
   forceField[0].interpolate( 5.0, level );

   for ( uint_t ts = 0; ts < timesteps; ts++ )
   {
      WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep: " << ts );

      hyteg::unresolved_particles::applyForceField( unresolvedParticles, forceField, level );
      hyteg::unresolved_particles::explicitEulerStep( unresolvedParticles, dt );

      unresolvedParticles.writeVTK();
   }
}

int main( int argc, char** argv )
{
   Environment env( argc, argv );
   WALBERLA_UNUSED( env );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   simpleTransportTest();

   return 0;
}