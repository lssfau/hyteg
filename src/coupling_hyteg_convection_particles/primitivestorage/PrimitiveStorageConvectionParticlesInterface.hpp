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

#pragma once

#pragma once

#include "convection_particles/domain/IDomain.h"

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;
using walberla::convection_particles::Vec3;

class PrimitiveStorageConvectionParticlesInterface : public walberla::convection_particles::domain::IDomain
{
 public:
   PrimitiveStorageConvectionParticlesInterface( const std::shared_ptr< PrimitiveStorage > & primitiveStorage ) :
       primitiveStorage_( primitiveStorage )
   {}

   virtual ~PrimitiveStorageConvectionParticlesInterface() {}

   bool isContainedInProcessSubdomain( const uint_t rank, const Vec3& pt ) const override;

   int findContainingProcessRank( const Vec3& pt ) const override;

   void periodicallyMapToDomain( Vec3 & pt ) const override;

   std::vector<uint_t> getNeighborProcesses() const override;

   bool intersectsWithProcessSubdomain( const uint_t rank, const Vec3& pt, const real_t& radius ) const override;

   void correctParticlePosition( Vec3 & pt ) const override;

 private:

   std::shared_ptr< PrimitiveStorage > primitiveStorage_;
   std::vector< uint_t > neighborProcesses_;

};

}

