/*
 * Copyright (c) 2023 Nils Kohl, Marcus Mohr.
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

/**
 * \page 12_UnresolvedParticles Tutorial for simple unresolved particulate flow.
 *
 * \dontinclude tutorials/12_UnresolvedParticles/12_UnresolvedParticles.cpp
 *
 * \brief In this tutorial, unresolved particles are transported subject to an underlying velocity field.
 *
 * This is a really basic case. The MESA-PD framework is quite powerful and allows much more than what is displayed here.
 * For instance, particles can be influenced by force fields (such that their mass is not neglected during integration),
 * and particle-particle-interactions can be implemented, too. For now we restrict ourselves to directly setting the particle
 * velocities, such that their mass and size does not influence the results.
 *
 * First we set up a simple cube shaped domain.
 * \snippet tutorials/12_UnresolvedParticles/12_UnresolvedParticles.cpp domain
 *
 * Then we define a velocity field that acts on our particles.
 * It resembles a convection cell in the x-y-plane. We will later scale it depending on the simulated time.
 * \snippet tutorials/12_UnresolvedParticles/12_UnresolvedParticles.cpp velocity
 *
 * We use two fields to avoid the usually slow interpolate function. Instead of interpolating the function every time step,
 * we just store the unscaled function, and scale it via the assign() method. This way, we replace the interpolate() call by an
 * assign() call which is generally much faster, at the cost of a little bit of memory. It clearly only works here as we simply
 * scale the convection cell.
 * \snippet tutorials/12_UnresolvedParticles/12_UnresolvedParticles.cpp fields
 *
 * Now we allocate our particles. The setup is commented inline below.
 * \snippet tutorials/12_UnresolvedParticles/12_UnresolvedParticles.cpp particles
 *
 * The VTK output needs to be initialized explicitly.
 * \snippet tutorials/12_UnresolvedParticles/12_UnresolvedParticles.cpp vtk
 *
 * Lets run our simulation. We alter the velocity field in each time step, apply it to the particles, and then integrate to
 * advect the particles.
 * \snippet tutorials/12_UnresolvedParticles/12_UnresolvedParticles.cpp simulation
 *
 * The particles are advected along the velocity field. Results are shown below.
 *
 * \htmlonly
   <center>
   <table>
   <tr><td><img src="12_UnresolvedParticles.0000.png" width="100%"/><center>initial</center></td></tr>
   <tr><td><img src="12_UnresolvedParticles.0050.png" width="100%"/><center>time step 50</center></td></tr>
   <tr><td><img src="12_UnresolvedParticles.0100.png" width="100%"/><center>time step 100</center></td></tr>
   <tr><td><img src="12_UnresolvedParticles.0500.png" width="100%"/><center>time step 500</center></td></tr>
   </table>
   </center>
   \endhtmlonly
 *
 * \section code Complete Program
 * \include tutorials/12_UnresolvedParticles/12_UnresolvedParticles.cpp
 *
 */

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "coupling_hyteg_unresolved_particles/UnresolvedParticles.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

struct Parameters
{
   uint_t timesteps    = 1000;
   real_t dt           = 1e-2;
   uint_t level        = 4;
   uint_t numParticles = 1000;
   real_t initRadius   = 0.3;

   real_t velFieldOscillation = 1;
};

void UnresolvedSpheres( Parameters parameters )
{
   /// [domain]
   MeshInfo meshInfo = MeshInfo::meshCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), 1, 1, 1 );
   auto     setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   writeDomainPartitioningVTK( *storage, "vtk", "domain" );
   /// [domain]

   /// [velocity]
   auto fx = []( const hyteg::Point3D& x ) { return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] ); };
   auto fy = []( const hyteg::Point3D& x ) { return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] ); };
   /// [velocity]

   /// [fields]
   P1VectorFunction< real_t > convcell( "convcell", storage, parameters.level, parameters.level );
   P1VectorFunction< real_t > vel( "vel", storage, parameters.level, parameters.level );

   convcell.interpolate( { fx, fy }, parameters.level );

   VTKOutput vtkOutput( "vtk", "fields", storage );
   vtkOutput.add( convcell );
   vtkOutput.add( vel );
   /// [fields]

   /// [particles]
   hyteg::unresolved_particles::UnresolvedParticles unresolvedParticles( storage );

   // Let`s randomly initialize particles in a ball at the center of our domain.
   for ( uint_t i = 0; i < parameters.numParticles; i++ )
   {
      real_t r     = std::sqrt( walberla::math::realRandom() ) * parameters.initRadius;
      real_t theta = walberla::math::realRandom() * pi;
      real_t phi   = walberla::math::realRandom() * 2 * pi;

      Point3D pos( r * std::sin( theta ) * std::cos( phi ) + 0.5, r * std::sin( theta ) * std::sin( phi ) + 0.5, r * std::cos( theta ) + 0.5 );

      // To create a new particle, we only specify the position and use the wrapper for convenience (it handles the parallel
      // process assignment automatically, and some further minor details). This needs to be called on all processes to make sure that
      // the particle is really created. It does not crash if not called collectively, but will only create the particle if the called
      // on the process whose local subdomain contains the passed position.
      //
      // The return value is an instance of std::optional. It is empty if the particle was not created on the _local_ process.
      // Otherwise, an iterator that points to the created particle can be accessed via std::optional's ::value() method.
      auto particleCreated = unresolvedParticles.createParticle( pos );

      // If the current process receives the particle, we can initialize its properties.
      if ( particleCreated )
      {
         // The particle was created locally. Let's access it via std::optional::value().
         auto particle = particleCreated.value();

         // The interaction radius is not physically relevant here, but may become in other applications if
         // particle-particle-interactions are implemented. It can be used for rendering later.
         particle->setInteractionRadius( 0.01 );

         // Since we do not apply forces, but set the particle velocity directly, their mass does not matter here.
         particle->setInvMass( 1 );

         // Custom properties can be added through the MESA-PD engine via code generation.
         // To aid simple, non-foreseeable use cases, it is possible to add real and integer scalars to vectors that are
         // properties by default.
         //
         // NOTE: If the vector is resized, this should be done for ALL particles. Especially the VTK output might crash
         // otherwise. Since we are creating all particles in this loop, this is not a problem here.
         //
         // Lets just add an integer to all particles, indicating whether it has been initialized
         // in the top or bottom half of the cube. This allows us to later see whether the particles are mixed.
         particle->getCustomRealRef().push_back( particle->getPosition()[1] < 0.5 ? -1 : 1 );
      }
   }
   /// [particles]

   /// [vtk]
   // Initializing VTK output, so that we can later just call the write function.
   unresolvedParticles.initVTK();
   // Adding our added custom property to the output.
   unresolvedParticles.addCustomRealElementToVTK( "original_cube_half", 0 );
   /// [vtk]

   /// [simulation]

   unresolvedParticles.writeVTK();
   vtkOutput.write( parameters.level, 0 );

   real_t simTime = 0;
   for ( uint_t ts = 0; ts < parameters.timesteps; ts++ )
   {
      if ( ts % 100 == 0 )
      {
         WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep: " << ts );
      }

      // We now alter our velocity field. This is just done here to produce some more interesting visuals.
      vel.assign( { real_c( 1 ) + real_c( 0.5 ) * std::cos( pi + parameters.velFieldOscillation * simTime ) },
                  { convcell },
                  parameters.level );

      // Applying the velocity field to the particles.
      hyteg::unresolved_particles::applyField(
          unresolvedParticles, vel, parameters.level, unresolved_particles::BackgroundFieldType::LINEAR_VELOCITY );

      // Particle time integration. Implementation in MESA-PD. Other integrators are available.
      hyteg::unresolved_particles::explicitEulerStep( unresolvedParticles, parameters.dt );

      unresolvedParticles.writeVTK();
      vtkOutput.write( parameters.level, ts + 1 );

      simTime += parameters.dt;
   }

   /// [simulation]
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   Parameters parameters;

   auto cfg = std::make_shared< walberla::config::Config >();
   WALBERLA_CHECK_NOT_NULLPTR( env.config(), "No parameter file given. Bye." );
   const walberla::Config::BlockHandle mainConf = env.config()->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION()
   {
      mainConf.listParameters();
   }

   parameters.timesteps           = mainConf.getParameter< uint_t >( "timesteps" );
   parameters.dt                  = mainConf.getParameter< real_t >( "dt" );
   parameters.level               = mainConf.getParameter< uint_t >( "level" );
   parameters.numParticles        = mainConf.getParameter< uint_t >( "numParticles" );
   parameters.initRadius          = mainConf.getParameter< real_t >( "initRadius" );
   parameters.velFieldOscillation = mainConf.getParameter< real_t >( "velFieldOscillation" );

   UnresolvedSpheres( parameters );

   return EXIT_SUCCESS;
}
