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

#pragma once

#include <unresolved_particles/data/ParticleAccessor.h>
#include <unresolved_particles/data/ParticleStorage.h>
#include <unresolved_particles/data/ShapeStorage.h>

namespace walberla {
namespace unresolved_particles {
namespace data {

class ParticleAccessorWithShape : public data::ParticleAccessor
{
public:
   ParticleAccessorWithShape(std::shared_ptr<data::ParticleStorage>& ps, std::shared_ptr<data::ShapeStorage>& ss)
      : ParticleAccessor(ps)
      , ss_(ss)
   {}

   const auto& getInvMass(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvMass();}
   const auto& getMass(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getMass();}
   const auto& getInvInertiaBF(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInvInertiaBF();}
   const auto& getInertiaBF(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getInertiaBF();}
   auto getVolume(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)]->getVolume();}

   data::BaseShape* getShape(const size_t p_idx) const {return ss_->shapes[ps_->getShapeIDRef(p_idx)].get();}
private:
   std::shared_ptr<data::ShapeStorage> ss_;
};

} //namespace data
} //namespace unresolved_particles
} //namespace walberla
