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

/// This tests simply advects some particles through a "tube" and later counts if all of them arrived.
void simpleTransportTest()
{
   bool vtk = false;

   real_t cubeSize             = 10;
   uint_t tunnelStretch        = 2;
   real_t particleInitCubeSize = cubeSize / 3;

   uint_t timesteps =
       300; // this was determined experimentally (it`s just the number of time steps it takes such that all particles wander over)
   real_t dt    = 2e-2;
   uint_t level = 2;

   MeshInfo meshInfo = MeshInfo::meshCuboid(
       Point3D( 0, 0, 0 ), Point3D( real_c( tunnelStretch ) * cubeSize, cubeSize, cubeSize ), tunnelStretch, 1, 1 );
   auto setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   writeDomainPartitioningVTK( *storage, "vtk", "domain" );

   hyteg::unresolved_particles::UnresolvedParticles unresolvedParticles( storage );

   uint_t numParticles = 50;
   for ( uint_t i = 0; i < numParticles; i++ )
   {
      Point3D pos( ( walberla::math::realRandom() - 0.5 ) * particleInitCubeSize,
                   ( walberla::math::realRandom() - 0.5 ) * particleInitCubeSize,
                   ( walberla::math::realRandom() - 0.5 ) * particleInitCubeSize );
      pos += Point3D( 0.5 * cubeSize, 0.5 * cubeSize, 0.5 * cubeSize );

      auto particleCreated = unresolvedParticles.createParticle( pos );
      if ( particleCreated )
      {
         auto particle = particleCreated.value();
         particle->setInteractionRadius( 0.1 );
         particle->setInvMass( 0.1 );
      }
   }

   if ( vtk )
   {
      unresolvedParticles.initVTK();
   }

   P1VectorFunction< real_t > forceField( "forceField", storage, level, level );
   forceField[0].interpolate( 5.0, level );

   auto countParticles = [&]( real_t ballCenterX, real_t ballRadius ) -> uint_t {
      uint_t countedParticles = 0;

      for ( auto p : *unresolvedParticles.getParticleStorage() )
      {
         if ( ( p.getPosition() - Vec3( ballCenterX, cubeSize / 2, cubeSize / 2 ) ).length() < ballRadius )
         {
            countedParticles++;
         }
      }

      return countedParticles;
   };

   for ( uint_t ts = 0; ts < timesteps; ts++ )
   {
      auto particlesStart = countParticles( cubeSize / 2, cubeSize / 2 );
      auto particlesEnd   = countParticles( cubeSize + cubeSize / 2, cubeSize / 2 );

      particlesStart = walberla::mpi::reduce( particlesStart, walberla::mpi::SUM );
      particlesEnd   = walberla::mpi::reduce( particlesEnd, walberla::mpi::SUM );

      WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep: " << ts << " particles start: " << particlesStart
                                               << " particles end: " << particlesEnd );

      hyteg::unresolved_particles::applyField(
          unresolvedParticles, forceField, level, hyteg::unresolved_particles::BackgroundFieldType::FORCE );
      hyteg::unresolved_particles::explicitEulerStep( unresolvedParticles, dt );

      if ( vtk )
      {
         unresolvedParticles.writeVTK();
      }
   }

   auto particlesFinal = countParticles( cubeSize + cubeSize / 2, cubeSize / 2 );
   particlesFinal      = walberla::mpi::allReduce( particlesFinal, walberla::mpi::SUM );

   WALBERLA_CHECK_EQUAL( particlesFinal, numParticles );
}

int main( int argc, char** argv )
{
   Environment env( argc, argv );
   WALBERLA_UNUSED( env );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   simpleTransportTest();

   return 0;
}