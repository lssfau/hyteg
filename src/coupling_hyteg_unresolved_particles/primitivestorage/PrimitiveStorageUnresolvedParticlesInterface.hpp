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

#pragma once

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "unresolved_particles/domain/IDomain.h"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;
using walberla::unresolved_particles::Vec3;

class PrimitiveStorageUnresolvedParticlesInterface : public walberla::unresolved_particles::domain::IDomain
{
 public:
   PrimitiveStorageUnresolvedParticlesInterface( const std::shared_ptr< PrimitiveStorage >& primitiveStorage )
   : primitiveStorage_( primitiveStorage )
   {
      WALBERLA_CHECK_GREATER(
          primitiveStorage_->getAdditionalHaloDepth(),
          0,
          "You need to bump the additional halo depth of the PrimitiveStorage (via its constructor) to enable correct "
          "(unresolved) particle communication." );
   }

   virtual ~PrimitiveStorageUnresolvedParticlesInterface() {}

   bool isContainedInProcessSubdomain( const uint_t rank, const Vec3& pt ) const override;

   int findContainingProcessRank( const Vec3& pt ) const override;

   void periodicallyMapToDomain( Vec3& pt ) const override;

   std::vector< uint_t > getNeighborProcesses() const override;

   bool intersectsWithProcessSubdomain( const uint_t rank, const Vec3& pt, const real_t& radius ) const override;

   void correctParticlePosition( Vec3& pt ) const override;

 private:
   std::shared_ptr< PrimitiveStorage > primitiveStorage_;
   std::vector< uint_t >               neighborProcesses_;
};

} // namespace hyteg
