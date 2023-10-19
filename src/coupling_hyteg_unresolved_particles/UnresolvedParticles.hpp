/*
* Copyright (c) 2023 Nils Kohl.
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

#pragma once

#include <blockforest/Initialization.h>
#include <core/Environment.h>
#include <core/math/Random.h>
#include <core/mpi/MPIManager.h>
#include <hyteg/communication/Syncing.hpp>
#include <hyteg/geometry/Intersection.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include <hyteg/primitivestorage/Visualization.hpp>
#include <unresolved_particles/data/ParticleAccessor.h>
#include <unresolved_particles/data/ParticleStorage.h>
#include <unresolved_particles/domain/BlockForestDomain.h>
#include <unresolved_particles/kernel/ExplicitEuler.h>
#include <unresolved_particles/kernel/ParticleSelector.h>
#include <unresolved_particles/mpi/SyncNextNeighborsNoGhosts.h>
#include <unresolved_particles/vtk/ParticleVtkOutput.h>
#include <vtk/VTKOutput.h>

#include "hyteg/p1functionspace/P1VectorFunction.hpp"

#include "coupling_hyteg_unresolved_particles/data/CustomSelectors.hpp"
#include "coupling_hyteg_unresolved_particles/primitivestorage/PrimitiveStorageUnresolvedParticlesInterface.hpp"

namespace hyteg {
namespace unresolved_particles {

/// \brief Wrapper class around the unresolved particle implementation and the corresponding MESA-PD - HyTeG coupling to simplify
///        working with unresolved particles.
///
///
class UnresolvedParticles
{
 public:
   UnresolvedParticles( const std::shared_ptr< PrimitiveStorage >& storage )
   : storage_( storage )
   , domain_( std::make_shared< PrimitiveStorageUnresolvedParticlesInterface >( storage ) )
   , particleStorage_( std::make_shared< walberla::unresolved_particles::data::ParticleStorage >( 1 ) )
   , particleAccessor_( std::make_shared< walberla::unresolved_particles::data::ParticleAccessor >( particleStorage_ ) )
   {}

   /// Returns the underlying MESA-PD ParticleStorage object.
   /// This can be used to directly use MESA-PD features that are not implemented through this wrapper.
   std::shared_ptr< walberla::unresolved_particles::data::ParticleStorage > getParticleStorage() { return particleStorage_; }

   /// Returns the underlying MESA-PD ParticleAccessor object.
   /// This can be used to directly use MESA-PD features that are not implemented through this wrapper.
   std::shared_ptr< walberla::unresolved_particles::data::ParticleAccessor > getParticleAccessor() { return particleAccessor_; }

   /// Returns the underlying PrimitiveStorage.
   std::shared_ptr< PrimitiveStorage > getPrimitiveStorage() { return storage_; }

   /** Creates a particle and assigns it to the correct process.
    *
    * Should be called collectively to ensure that the particle is created (but won't crash otherwise - might just not create
    * the particle).
    *
    * This method returns an instance of std::optional. If the particle was _not_ created on this process, the std::optional
    * instance is empty. If the particle _was_ created, the particle properties can be accessed via the iterator object contained
    * inside the std::optional instance. This can simply be achieved via calling the ::value() method (see example below).
    *
    * All properties of the particle (apart from the initial position) need to be set manually via the returned iterator (see
    * example below).
    *
    * The creation and modification of particles can alternatively be realized directly through MESA-PD's ParticleStorage class.
    * An instance of that is wrapped by this class, and can also be accessed. This really is just a convenience method since the
    * initial position of a particle determines the process rank, and therefore the parallel particle creation requires some
    * additional boilerplate code.
    *
    * Example:
    *
        \code

            auto particleCreated = unresolvedParticles.createParticle( somePosition );

            if ( particleCreated )
            {
                auto particleIt = particleCreated.value();
                particleIt->setInteractionRadius( size );
                particleIt->setInvMass( 1 / mass );
            }

        \endcode
   **/
   std::optional< walberla::unresolved_particles::data::ParticleStorage::iterator > createParticle( const Point3D& pos )
   {
      auto posVec3 = toVec3( pos );
      if ( domain_->findContainingProcessRank( posVec3 ) == walberla::mpi::MPIManager::instance()->rank() )
      {
         auto particleIt = particleStorage_->create();
         particleIt->setPosition( posVec3 );
         particleIt->setOwner( walberla::mpi::MPIManager::instance()->rank() );
         return particleIt;
      }

      return {};
   }

   /// Initializes the VTK output. Subsequent calls to writeVTK() output data with the corresponding write interval.
   void initVTK( std::string identifier      = "particles",
                 uint_t      writeInterval   = 1,
                 std::string baseFolder      = "vtk",
                 std::string executionFolder = "simulation_step" )
   {
      particleVtkOutput_ = std::make_shared< walberla::unresolved_particles::vtk::ParticleVtkOutput >( particleStorage_ );
      particleVtkOutput_->addOutput< walberla::unresolved_particles::data::SelectParticleUid >( "uid" );
      particleVtkOutput_->addOutput< walberla::unresolved_particles::data::SelectParticleLinearVelocity >( "velocity" );
      particleVtkOutput_->addOutput< walberla::unresolved_particles::data::SelectParticleInteractionRadius >( "radius" );
      particleVtkOutput_->addOutput< walberla::unresolved_particles::data::SelectParticleOwner >( "owner" );

      particleVtkWriter_ = walberla::vtk::createVTKOutput_PointData(
          particleVtkOutput_, identifier, writeInterval, baseFolder, executionFolder, false, true, true, true );
   }

   /// The default particle implementation has properties to carry real scalars.
   /// Those scalars can of course be added through the VTK output.
   /// Since it may not be desirable to add all of them to the VTK output by default, this has to be explicitly asked for here.
   /// Also, at this point the VTK output cannot know how many scalar properties are added to the particles, so we need to do this
   /// explicitly anyway.
   ///
   /// \param name string to be displayed in VTK for this scalar
   /// \param component component in the custom vector property of the particle
   void addCustomRealElementToVTK( std::string name, uint_t component )
   {
      WALBERLA_CHECK_NOT_NULLPTR( particleVtkOutput_, "Particle VTK output has not been initialized." );

      particleVtkOutput_->addOutput( name,
                                     std::make_shared< walberla::unresolved_particles::vtk::OutputSelector<
                                         walberla::unresolved_particles::data::SelectParticleCustomRealElement > >(
                                         walberla::unresolved_particles::data::SelectParticleCustomRealElement( component ) ) );
   }

   /// The default particle implementation has properties to carry integer scalars.
   /// Those scalars can of course be added through the VTK output.
   /// Since it may not be desirable to add all of them to the VTK output by default, this has to be explicitly asked for here.
   /// Also, at this point the VTK output cannot know how many scalar properties are added to the particles, so we need to do this
   /// explicitly anyway.
   ///
   /// \param name string to be displayed in VTK for this scalar
   /// \param component component in the custom vector property of the particle
   void addCustomIntElementToVTK( std::string name, uint_t component )
   {
      WALBERLA_CHECK_NOT_NULLPTR( particleVtkOutput_, "Particle VTK output has not been initialized." );

      particleVtkOutput_->addOutput( name,
                                     std::make_shared< walberla::unresolved_particles::vtk::OutputSelector<
                                         walberla::unresolved_particles::data::SelectParticleCustomIntElement > >(
                                         walberla::unresolved_particles::data::SelectParticleCustomIntElement( component ) ) );
   }

   /// Writes VTK data. The VTK output must be initialized with initVTK() first.
   void writeVTK()
   {
      WALBERLA_CHECK_NOT_NULLPTR( particleVtkWriter_, "Particle VTK output has not been initialized." );
      sync();
      particleVtkWriter_->write();
   }

   /// Communicates particles after movement. Probably does not need to be called if wrapper functions are used.
   /// This also erases particles that are not found (i.e. when it leaves the domain).
   void sync()
   {
      walberla::unresolved_particles::mpi::SyncNextNeighborsNoGhosts SNN;
      SNN( *particleStorage_, *domain_ );

      // There is a small issue in the generated sync code.
      // Since it is skipped entirely in serial mode, particles that leave the domain are not erased when only one process is
      // active. This is handled below.
      if ( walberla::mpi::MPIManager::instance()->numProcesses() == 1 )
      {
         for ( auto pIt = particleStorage_->begin(); pIt != particleStorage_->end(); )
         {
            //skip all ghost particles
            if ( walberla::unresolved_particles::data::particle_flags::isSet(
                     pIt->getFlags(), walberla::unresolved_particles::data::particle_flags::GHOST ) )
            {
               ++pIt;
               continue;
            }

            //skip all particles that do not communicate (create ghost particles) on other processes
            if ( walberla::unresolved_particles::data::particle_flags::isSet(
                     pIt->getFlags(), walberla::unresolved_particles::data::particle_flags::NON_COMMUNICATING ) )
            {
               ++pIt;
               continue;
            }

            //correct position to make sure particle is always inside the domain!
            //everything is decided by the master particle therefore ghost particles are not touched
            if ( !walberla::unresolved_particles::data::particle_flags::isSet(
                     pIt->getFlags(), walberla::unresolved_particles::data::particle_flags::FIXED ) &&
                 !walberla::unresolved_particles::data::particle_flags::isSet(
                     pIt->getFlags(), walberla::unresolved_particles::data::particle_flags::GHOST ) )
            {
               domain_->periodicallyMapToDomain( pIt->getPositionRef() );
            }

            // Note: At this point we know that the particle was locally owned before the position update.
            WALBERLA_CHECK_EQUAL( pIt->getOwner(), 0 );

            WALBERLA_LOG_DETAIL( "Processing local particle " << pIt->getUid() );

            const auto ownerRank = domain_->findContainingProcessRank( pIt->getPosition() );
            if ( ownerRank < 0 )
            {
               // No owner found: Outflow condition.
               WALBERLA_LOG_DETAIL( "Sending deletion notifications for particle " << pIt->getUid() << " due to outflow." );

               // remove particle
               pIt = particleStorage_->erase( pIt );
               continue;
            }

            ++pIt;
         }
      }
   }

 private:
   std::shared_ptr< PrimitiveStorage >                                       storage_;
   std::shared_ptr< PrimitiveStorageUnresolvedParticlesInterface >           domain_;
   std::shared_ptr< walberla::unresolved_particles::data::ParticleStorage >  particleStorage_;
   std::shared_ptr< walberla::unresolved_particles::data::ParticleAccessor > particleAccessor_;

   std::shared_ptr< walberla::unresolved_particles::vtk::ParticleVtkOutput > particleVtkOutput_;
   std::shared_ptr< walberla::vtk::VTKOutput >                               particleVtkWriter_;
};

enum class BackgroundFieldType
{
   FORCE,
   LINEAR_VELOCITY
};

/// This function applies a field to all particles. The physical property that is applied can be set through a parameter.
/// Typically, force fields are applied - and this is probably the most physically sound approach.
/// However, alternatively, the linear velocity of particles can be set directly. This may be relevant for some applications.
void applyField( UnresolvedParticles&              unresolvedParticles,
                 const P1VectorFunction< real_t >& field,
                 uint_t                            level,
                 BackgroundFieldType               backgroundFieldType )
{
   unresolvedParticles.sync();
   hyteg::communication::syncVectorFunctionBetweenPrimitives( field, level );

   auto particleStorage = unresolvedParticles.getParticleStorage();

   bool threeD = unresolvedParticles.getPrimitiveStorage()->hasGlobalCells();

   for ( auto p : *particleStorage )
   {
      auto posVec3 = p.getPosition();
      auto pos     = hyteg::toPoint3D( posVec3 );

      real_t fx = 0;
      real_t fy = 0;
      real_t fz = 0;

      auto success = true;
      success &= field[0].evaluate( pos, level, fx );
      success &= field[1].evaluate( pos, level, fy );
      if ( threeD )
      {
         success &= field[2].evaluate( pos, level, fz );
      }

      WALBERLA_CHECK( success,
                      "Unresolved particles: force could not be evaluated on process "
                          << walberla::mpi::MPIManager::instance()->rank() << " at position " << pos );

      if ( success )
      {
         switch ( backgroundFieldType )
         {
         case BackgroundFieldType::FORCE:
            p.setForce( Vec3( fx, fy, fz ) );
            break;
         case BackgroundFieldType::LINEAR_VELOCITY:
            p.setLinearVelocity( Vec3( fx, fy, fz ) );
            break;
         }
      }
   }
}

/// Applies an explicit Euler step as implemented by MESA-PD.
/// Please refer to the corresponding documentation of the kernel.
void explicitEulerStep( UnresolvedParticles& unresolvedParticles, real_t dt )
{
   unresolvedParticles.sync();

   auto particleStorage  = unresolvedParticles.getParticleStorage();
   auto particleAccessor = unresolvedParticles.getParticleAccessor();

   walberla::unresolved_particles::kernel::ExplicitEuler explicitEuler( dt );

   particleStorage->forEachParticle(
       false, walberla::unresolved_particles::kernel::SelectLocal(), *particleAccessor, explicitEuler, *particleAccessor );
}

} // namespace unresolved_particles
} // namespace hyteg