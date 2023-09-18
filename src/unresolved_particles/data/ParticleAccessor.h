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

#include <unresolved_particles/data/IAccessor.h>
#include <unresolved_particles/data/ParticleStorage.h>

#include <core/UniqueID.h>

#include <limits>

namespace walberla {
namespace unresolved_particles {
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
   ~ParticleAccessor() override = default;
   walberla::id_t const & getUid(const size_t p_idx) const {return ps_->getUid(p_idx);}
   walberla::id_t& getUidRef(const size_t p_idx) {return ps_->getUidRef(p_idx);}
   void setUid(const size_t p_idx, walberla::id_t const & v) { ps_->setUid(p_idx, v);}
   
   walberla::unresolved_particles::Vec3 const & getPosition(const size_t p_idx) const {return ps_->getPosition(p_idx);}
   walberla::unresolved_particles::Vec3& getPositionRef(const size_t p_idx) {return ps_->getPositionRef(p_idx);}
   void setPosition(const size_t p_idx, walberla::unresolved_particles::Vec3 const & v) { ps_->setPosition(p_idx, v);}
   
   walberla::real_t const & getInteractionRadius(const size_t p_idx) const {return ps_->getInteractionRadius(p_idx);}
   walberla::real_t& getInteractionRadiusRef(const size_t p_idx) {return ps_->getInteractionRadiusRef(p_idx);}
   void setInteractionRadius(const size_t p_idx, walberla::real_t const & v) { ps_->setInteractionRadius(p_idx, v);}
   
   walberla::unresolved_particles::data::particle_flags::FlagT const & getFlags(const size_t p_idx) const {return ps_->getFlags(p_idx);}
   walberla::unresolved_particles::data::particle_flags::FlagT& getFlagsRef(const size_t p_idx) {return ps_->getFlagsRef(p_idx);}
   void setFlags(const size_t p_idx, walberla::unresolved_particles::data::particle_flags::FlagT const & v) { ps_->setFlags(p_idx, v);}
   
   int const & getOwner(const size_t p_idx) const {return ps_->getOwner(p_idx);}
   int& getOwnerRef(const size_t p_idx) {return ps_->getOwnerRef(p_idx);}
   void setOwner(const size_t p_idx, int const & v) { ps_->setOwner(p_idx, v);}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getGhostOwners(const size_t p_idx) const {return ps_->getGhostOwners(p_idx);}
   std::unordered_set<walberla::mpi::MPIRank>& getGhostOwnersRef(const size_t p_idx) {return ps_->getGhostOwnersRef(p_idx);}
   void setGhostOwners(const size_t p_idx, std::unordered_set<walberla::mpi::MPIRank> const & v) { ps_->setGhostOwners(p_idx, v);}
   
   walberla::unresolved_particles::Vec3 const & getLinearVelocity(const size_t p_idx) const {return ps_->getLinearVelocity(p_idx);}
   walberla::unresolved_particles::Vec3& getLinearVelocityRef(const size_t p_idx) {return ps_->getLinearVelocityRef(p_idx);}
   void setLinearVelocity(const size_t p_idx, walberla::unresolved_particles::Vec3 const & v) { ps_->setLinearVelocity(p_idx, v);}
   
   walberla::real_t const & getInvMass(const size_t p_idx) const {return ps_->getInvMass(p_idx);}
   walberla::real_t& getInvMassRef(const size_t p_idx) {return ps_->getInvMassRef(p_idx);}
   void setInvMass(const size_t p_idx, walberla::real_t const & v) { ps_->setInvMass(p_idx, v);}
   
   walberla::unresolved_particles::Vec3 const & getForce(const size_t p_idx) const {return ps_->getForce(p_idx);}
   walberla::unresolved_particles::Vec3& getForceRef(const size_t p_idx) {return ps_->getForceRef(p_idx);}
   void setForce(const size_t p_idx, walberla::unresolved_particles::Vec3 const & v) { ps_->setForce(p_idx, v);}
   
   walberla::unresolved_particles::Vec3 const & getOldForce(const size_t p_idx) const {return ps_->getOldForce(p_idx);}
   walberla::unresolved_particles::Vec3& getOldForceRef(const size_t p_idx) {return ps_->getOldForceRef(p_idx);}
   void setOldForce(const size_t p_idx, walberla::unresolved_particles::Vec3 const & v) { ps_->setOldForce(p_idx, v);}
   
   walberla::unresolved_particles::Rot3 const & getRotation(const size_t p_idx) const {return ps_->getRotation(p_idx);}
   walberla::unresolved_particles::Rot3& getRotationRef(const size_t p_idx) {return ps_->getRotationRef(p_idx);}
   void setRotation(const size_t p_idx, walberla::unresolved_particles::Rot3 const & v) { ps_->setRotation(p_idx, v);}
   
   walberla::unresolved_particles::Vec3 const & getAngularVelocity(const size_t p_idx) const {return ps_->getAngularVelocity(p_idx);}
   walberla::unresolved_particles::Vec3& getAngularVelocityRef(const size_t p_idx) {return ps_->getAngularVelocityRef(p_idx);}
   void setAngularVelocity(const size_t p_idx, walberla::unresolved_particles::Vec3 const & v) { ps_->setAngularVelocity(p_idx, v);}
   
   walberla::unresolved_particles::Vec3 const & getTorque(const size_t p_idx) const {return ps_->getTorque(p_idx);}
   walberla::unresolved_particles::Vec3& getTorqueRef(const size_t p_idx) {return ps_->getTorqueRef(p_idx);}
   void setTorque(const size_t p_idx, walberla::unresolved_particles::Vec3 const & v) { ps_->setTorque(p_idx, v);}
   
   walberla::unresolved_particles::Vec3 const & getOldTorque(const size_t p_idx) const {return ps_->getOldTorque(p_idx);}
   walberla::unresolved_particles::Vec3& getOldTorqueRef(const size_t p_idx) {return ps_->getOldTorqueRef(p_idx);}
   void setOldTorque(const size_t p_idx, walberla::unresolved_particles::Vec3 const & v) { ps_->setOldTorque(p_idx, v);}
   
   walberla::unresolved_particles::Mat3 const & getInvInertiaBF(const size_t p_idx) const {return ps_->getInvInertiaBF(p_idx);}
   walberla::unresolved_particles::Mat3& getInvInertiaBFRef(const size_t p_idx) {return ps_->getInvInertiaBFRef(p_idx);}
   void setInvInertiaBF(const size_t p_idx, walberla::unresolved_particles::Mat3 const & v) { ps_->setInvInertiaBF(p_idx, v);}
   
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
   ~SingleParticleAccessor() override = default;
   walberla::id_t const & getUid(const size_t /*p_idx*/) const {return uid_;}
   void setUid(const size_t /*p_idx*/, walberla::id_t const & v) { uid_ = v;}
   walberla::id_t& getUidRef(const size_t /*p_idx*/) {return uid_;}
   
   walberla::unresolved_particles::Vec3 const & getPosition(const size_t /*p_idx*/) const {return position_;}
   void setPosition(const size_t /*p_idx*/, walberla::unresolved_particles::Vec3 const & v) { position_ = v;}
   walberla::unresolved_particles::Vec3& getPositionRef(const size_t /*p_idx*/) {return position_;}
   
   walberla::real_t const & getInteractionRadius(const size_t /*p_idx*/) const {return interactionRadius_;}
   void setInteractionRadius(const size_t /*p_idx*/, walberla::real_t const & v) { interactionRadius_ = v;}
   walberla::real_t& getInteractionRadiusRef(const size_t /*p_idx*/) {return interactionRadius_;}
   
   walberla::unresolved_particles::data::particle_flags::FlagT const & getFlags(const size_t /*p_idx*/) const {return flags_;}
   void setFlags(const size_t /*p_idx*/, walberla::unresolved_particles::data::particle_flags::FlagT const & v) { flags_ = v;}
   walberla::unresolved_particles::data::particle_flags::FlagT& getFlagsRef(const size_t /*p_idx*/) {return flags_;}
   
   int const & getOwner(const size_t /*p_idx*/) const {return owner_;}
   void setOwner(const size_t /*p_idx*/, int const & v) { owner_ = v;}
   int& getOwnerRef(const size_t /*p_idx*/) {return owner_;}
   
   std::unordered_set<walberla::mpi::MPIRank> const & getGhostOwners(const size_t /*p_idx*/) const {return ghostOwners_;}
   void setGhostOwners(const size_t /*p_idx*/, std::unordered_set<walberla::mpi::MPIRank> const & v) { ghostOwners_ = v;}
   std::unordered_set<walberla::mpi::MPIRank>& getGhostOwnersRef(const size_t /*p_idx*/) {return ghostOwners_;}
   
   walberla::unresolved_particles::Vec3 const & getLinearVelocity(const size_t /*p_idx*/) const {return linearVelocity_;}
   void setLinearVelocity(const size_t /*p_idx*/, walberla::unresolved_particles::Vec3 const & v) { linearVelocity_ = v;}
   walberla::unresolved_particles::Vec3& getLinearVelocityRef(const size_t /*p_idx*/) {return linearVelocity_;}
   
   walberla::real_t const & getInvMass(const size_t /*p_idx*/) const {return invMass_;}
   void setInvMass(const size_t /*p_idx*/, walberla::real_t const & v) { invMass_ = v;}
   walberla::real_t& getInvMassRef(const size_t /*p_idx*/) {return invMass_;}
   
   walberla::unresolved_particles::Vec3 const & getForce(const size_t /*p_idx*/) const {return force_;}
   void setForce(const size_t /*p_idx*/, walberla::unresolved_particles::Vec3 const & v) { force_ = v;}
   walberla::unresolved_particles::Vec3& getForceRef(const size_t /*p_idx*/) {return force_;}
   
   walberla::unresolved_particles::Vec3 const & getOldForce(const size_t /*p_idx*/) const {return oldForce_;}
   void setOldForce(const size_t /*p_idx*/, walberla::unresolved_particles::Vec3 const & v) { oldForce_ = v;}
   walberla::unresolved_particles::Vec3& getOldForceRef(const size_t /*p_idx*/) {return oldForce_;}
   
   walberla::unresolved_particles::Rot3 const & getRotation(const size_t /*p_idx*/) const {return rotation_;}
   void setRotation(const size_t /*p_idx*/, walberla::unresolved_particles::Rot3 const & v) { rotation_ = v;}
   walberla::unresolved_particles::Rot3& getRotationRef(const size_t /*p_idx*/) {return rotation_;}
   
   walberla::unresolved_particles::Vec3 const & getAngularVelocity(const size_t /*p_idx*/) const {return angularVelocity_;}
   void setAngularVelocity(const size_t /*p_idx*/, walberla::unresolved_particles::Vec3 const & v) { angularVelocity_ = v;}
   walberla::unresolved_particles::Vec3& getAngularVelocityRef(const size_t /*p_idx*/) {return angularVelocity_;}
   
   walberla::unresolved_particles::Vec3 const & getTorque(const size_t /*p_idx*/) const {return torque_;}
   void setTorque(const size_t /*p_idx*/, walberla::unresolved_particles::Vec3 const & v) { torque_ = v;}
   walberla::unresolved_particles::Vec3& getTorqueRef(const size_t /*p_idx*/) {return torque_;}
   
   walberla::unresolved_particles::Vec3 const & getOldTorque(const size_t /*p_idx*/) const {return oldTorque_;}
   void setOldTorque(const size_t /*p_idx*/, walberla::unresolved_particles::Vec3 const & v) { oldTorque_ = v;}
   walberla::unresolved_particles::Vec3& getOldTorqueRef(const size_t /*p_idx*/) {return oldTorque_;}
   
   walberla::unresolved_particles::Mat3 const & getInvInertiaBF(const size_t /*p_idx*/) const {return invInertiaBF_;}
   void setInvInertiaBF(const size_t /*p_idx*/, walberla::unresolved_particles::Mat3 const & v) { invInertiaBF_ = v;}
   walberla::unresolved_particles::Mat3& getInvInertiaBFRef(const size_t /*p_idx*/) {return invInertiaBF_;}
   
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
   walberla::unresolved_particles::Vec3 position_;
   walberla::real_t interactionRadius_;
   walberla::unresolved_particles::data::particle_flags::FlagT flags_;
   int owner_;
   std::unordered_set<walberla::mpi::MPIRank> ghostOwners_;
   walberla::unresolved_particles::Vec3 linearVelocity_;
   walberla::real_t invMass_;
   walberla::unresolved_particles::Vec3 force_;
   walberla::unresolved_particles::Vec3 oldForce_;
   walberla::unresolved_particles::Rot3 rotation_;
   walberla::unresolved_particles::Vec3 angularVelocity_;
   walberla::unresolved_particles::Vec3 torque_;
   walberla::unresolved_particles::Vec3 oldTorque_;
   walberla::unresolved_particles::Mat3 invInertiaBF_;
   blockforest::BlockID currentBlock_;
   std::unordered_set<walberla::mpi::MPIRank> neighborState_;
};

} //namespace data
} //namespace unresolved_particles
} //namespace walberla