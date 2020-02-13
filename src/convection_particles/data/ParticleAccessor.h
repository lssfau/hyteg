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
//! \file ParticleAccessor.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <convection_particles/data/IAccessor.h>
#include <convection_particles/data/ParticleStorage.h>

#include <core/UniqueID.h>

#include <limits>

namespace walberla {
namespace convection_particles {
namespace data {

/**
 * @brief Basic ParticleAccessor for the ParticleStorage
 *
 * Provides get, set and getRef for all members of the ParticleStorage.
 * Can be used as a basis class for a more advanced ParticleAccessor.
 */
class ParticleAccessor : public IAccessor
{
public:
   ParticleAccessor(const std::shared_ptr<data::ParticleStorage>& ps) : ps_(ps) {}
   virtual ~ParticleAccessor() = default;
   walberla::id_t const & getUid(const size_t p_idx) const {return ps_->getUid(p_idx);}
   walberla::id_t& getUidRef(const size_t p_idx) {return ps_->getUidRef(p_idx);}
   void setUid(const size_t p_idx, walberla::id_t const & v) { ps_->setUid(p_idx, v);}
   
   walberla::convection_particles::Vec3 const & getPosition(const size_t p_idx) const {return ps_->getPosition(p_idx);}
   walberla::convection_particles::Vec3& getPositionRef(const size_t p_idx) {return ps_->getPositionRef(p_idx);}
   void setPosition(const size_t p_idx, walberla::convection_particles::Vec3 const & v) { ps_->setPosition(p_idx, v);}
   
   walberla::real_t const & getInteractionRadius(const size_t p_idx) const {return ps_->getInteractionRadius(p_idx);}
   walberla::real_t& getInteractionRadiusRef(const size_t p_idx) {return ps_->getInteractionRadiusRef(p_idx);}
   void setInteractionRadius(const size_t p_idx, walberla::real_t const & v) { ps_->setInteractionRadius(p_idx, v);}
   
   walberla::convection_particles::data::particle_flags::FlagT const & getFlags(const size_t p_idx) const {return ps_->getFlags(p_idx);}
   walberla::convection_particles::data::particle_flags::FlagT& getFlagsRef(const size_t p_idx) {return ps_->getFlagsRef(p_idx);}
   void setFlags(const size_t p_idx, walberla::convection_particles::data::particle_flags::FlagT const & v) { ps_->setFlags(p_idx, v);}
   
   int const & getOwner(const size_t p_idx) const {return ps_->getOwner(p_idx);}
   int& getOwnerRef(const size_t p_idx) {return ps_->getOwnerRef(p_idx);}
   void setOwner(const size_t p_idx, int const & v) { ps_->setOwner(p_idx, v);}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getGhostOwners(const size_t p_idx) const {return ps_->getGhostOwners(p_idx);}
   std::unordered_set<walberla::mpi::MPIRank>& getGhostOwnersRef(const size_t p_idx) {return ps_->getGhostOwnersRef(p_idx);}
   void setGhostOwners(const size_t p_idx, std::unordered_set<walberla::mpi::MPIRank> const & v) { ps_->setGhostOwners(p_idx, v);}
   
   walberla::convection_particles::Vec3 const & getVelocity(const size_t p_idx) const {return ps_->getVelocity(p_idx);}
   walberla::convection_particles::Vec3& getVelocityRef(const size_t p_idx) {return ps_->getVelocityRef(p_idx);}
   void setVelocity(const size_t p_idx, walberla::convection_particles::Vec3 const & v) { ps_->setVelocity(p_idx, v);}
   
   blockforest::BlockID const & getCurrentBlock(const size_t p_idx) const {return ps_->getCurrentBlock(p_idx);}
   blockforest::BlockID& getCurrentBlockRef(const size_t p_idx) {return ps_->getCurrentBlockRef(p_idx);}
   void setCurrentBlock(const size_t p_idx, blockforest::BlockID const & v) { ps_->setCurrentBlock(p_idx, v);}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getNeighborState(const size_t p_idx) const {return ps_->getNeighborState(p_idx);}
   std::unordered_set<walberla::mpi::MPIRank>& getNeighborStateRef(const size_t p_idx) {return ps_->getNeighborStateRef(p_idx);}
   void setNeighborState(const size_t p_idx, std::unordered_set<walberla::mpi::MPIRank> const & v) { ps_->setNeighborState(p_idx, v);}
   

   id_t getInvalidUid() const {return UniqueID<data::Particle>::invalidID();}
   size_t getInvalidIdx() const {return std::numeric_limits<size_t>::max();}
   /**
   * @brief Returns the index of particle specified by uid.
   * @param uid unique id of the particle to be looked up
   * @return the index of the particle or std::numeric_limits<size_t>::max() if the particle is not found
   */
   size_t uidToIdx(const id_t& uid) const {auto it = ps_->find(uid); return it != ps_->end() ? it.getIdx() : std::numeric_limits<size_t>::max();}
   size_t size() const { return ps_->size(); }

   inline size_t create(const id_t& uid);
   inline size_t erase(const size_t& idx);
   inline size_t find(const id_t& uid);
protected:
   std::shared_ptr<data::ParticleStorage> ps_;
};

inline size_t ParticleAccessor::create(const id_t& uid)
{
   auto it = ps_->create(uid);
   return it.getIdx();
}
inline size_t ParticleAccessor::erase(const size_t& idx)
{
   data::ParticleStorage::iterator it(ps_.get(), idx);
   it = ps_->erase(it);
   return it.getIdx();
}
inline size_t ParticleAccessor::find(const id_t& uid)
{
   auto it = ps_->find(uid);
   return it.getIdx();
}

/**
 * @brief Basic ParticleAccessor which emulates a single particle in a ParticleStorage
 * without actually needing a ParticleStorage. This class is used mainly for testing purposes.
 *
 * Provides get, set and getRef.
 */
class SingleParticleAccessor : public IAccessor
{
public:
   virtual ~SingleParticleAccessor() = default;
   walberla::id_t const & getUid(const size_t /*p_idx*/) const {return uid_;}
   void setUid(const size_t /*p_idx*/, walberla::id_t const & v) { uid_ = v;}
   walberla::id_t& getUidRef(const size_t /*p_idx*/) {return uid_;}
   
   walberla::convection_particles::Vec3 const & getPosition(const size_t /*p_idx*/) const {return position_;}
   void setPosition(const size_t /*p_idx*/, walberla::convection_particles::Vec3 const & v) { position_ = v;}
   walberla::convection_particles::Vec3& getPositionRef(const size_t /*p_idx*/) {return position_;}
   
   walberla::real_t const & getInteractionRadius(const size_t /*p_idx*/) const {return interactionRadius_;}
   void setInteractionRadius(const size_t /*p_idx*/, walberla::real_t const & v) { interactionRadius_ = v;}
   walberla::real_t& getInteractionRadiusRef(const size_t /*p_idx*/) {return interactionRadius_;}
   
   walberla::convection_particles::data::particle_flags::FlagT const & getFlags(const size_t /*p_idx*/) const {return flags_;}
   void setFlags(const size_t /*p_idx*/, walberla::convection_particles::data::particle_flags::FlagT const & v) { flags_ = v;}
   walberla::convection_particles::data::particle_flags::FlagT& getFlagsRef(const size_t /*p_idx*/) {return flags_;}
   
   int const & getOwner(const size_t /*p_idx*/) const {return owner_;}
   void setOwner(const size_t /*p_idx*/, int const & v) { owner_ = v;}
   int& getOwnerRef(const size_t /*p_idx*/) {return owner_;}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getGhostOwners(const size_t /*p_idx*/) const {return ghostOwners_;}
   void setGhostOwners(const size_t /*p_idx*/, std::unordered_set<walberla::mpi::MPIRank> const & v) { ghostOwners_ = v;}
   std::unordered_set<walberla::mpi::MPIRank>& getGhostOwnersRef(const size_t /*p_idx*/) {return ghostOwners_;}
   
   walberla::convection_particles::Vec3 const & getVelocity(const size_t /*p_idx*/) const {return velocity_;}
   void setVelocity(const size_t /*p_idx*/, walberla::convection_particles::Vec3 const & v) { velocity_ = v;}
   walberla::convection_particles::Vec3& getVelocityRef(const size_t /*p_idx*/) {return velocity_;}
   
   blockforest::BlockID const & getCurrentBlock(const size_t /*p_idx*/) const {return currentBlock_;}
   void setCurrentBlock(const size_t /*p_idx*/, blockforest::BlockID const & v) { currentBlock_ = v;}
   blockforest::BlockID& getCurrentBlockRef(const size_t /*p_idx*/) {return currentBlock_;}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getNeighborState(const size_t /*p_idx*/) const {return neighborState_;}
   void setNeighborState(const size_t /*p_idx*/, std::unordered_set<walberla::mpi::MPIRank> const & v) { neighborState_ = v;}
   std::unordered_set<walberla::mpi::MPIRank>& getNeighborStateRef(const size_t /*p_idx*/) {return neighborState_;}
   

   id_t getInvalidUid() const {return UniqueID<data::Particle>::invalidID();}
   size_t getInvalidIdx() const {return std::numeric_limits<size_t>::max();}
   /**
   * @brief Returns the index of particle specified by uid.
   * @param uid unique id of the particle to be looked up
   * @return the index of the particle or std::numeric_limits<size_t>::max() if the particle is not found
   */
   size_t uidToIdx(const id_t& uid) const {return uid == uid_ ? 0 : std::numeric_limits<size_t>::max();}
   size_t size() const { return 1; }
private:
   walberla::id_t uid_;
   walberla::convection_particles::Vec3 position_;
   walberla::real_t interactionRadius_;
   walberla::convection_particles::data::particle_flags::FlagT flags_;
   int owner_;
   std::unordered_set<walberla::mpi::MPIRank> ghostOwners_;
   walberla::convection_particles::Vec3 velocity_;
   blockforest::BlockID currentBlock_;
   std::unordered_set<walberla::mpi::MPIRank> neighborState_;
};

} //namespace data
} //namespace convection_particles
} //namespace walberla