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

#include "PrimitiveStorageUnresolvedParticlesInterface.hpp"

#include <algorithm>

#include "hyteg/geometry/BlendingHelpers.hpp"
#include "hyteg/geometry/Intersection.hpp"

#include "unresolved_particles/data/DataTypes.h"
#include "unresolved_particles/domain/IDomain.h"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;
using walberla::unresolved_particles::Vec3;

bool PrimitiveStorageUnresolvedParticlesInterface::isContainedInProcessSubdomain( const uint_t rank, const Vec3& pt ) const
{
   const int containingRank = findContainingProcessRank( pt );
   return (int) rank == containingRank;
}

int PrimitiveStorageUnresolvedParticlesInterface::findContainingProcessRank( const Vec3& pt ) const
{
   std::set< int > containingProcessRanks;
   if ( !primitiveStorage_->hasGlobalCells() )
   {
      Point3D    pointOfInterest( pt[0], pt[1], pt[2] );
      const auto target      = mapFromPhysicalToComputationalDomain2D( primitiveStorage_, pointOfInterest, 0, 0, false, true );
      const auto coordExists = std::get< 0 >( target );
      const auto id          = std::get< 1 >( target );
      const auto computationalCoords = std::get< 2 >( target );

      if ( coordExists )
      {
         containingProcessRanks.insert( static_cast< int >( primitiveStorage_->getPrimitiveRank( id ) ) );
      }
   }
   else
   {
      Point3D    pointOfInterest( pt[0], pt[1], pt[2] );
      const auto target      = mapFromPhysicalToComputationalDomain3D( primitiveStorage_, pointOfInterest, 0, 0, false, true );
      const auto coordExists = std::get< 0 >( target );
      const auto id          = std::get< 1 >( target );
      const auto computationalCoords = std::get< 2 >( target );

      if ( coordExists )
      {
         containingProcessRanks.insert( static_cast< int >( primitiveStorage_->getPrimitiveRank( id ) ) );
      }
   }

   if ( containingProcessRanks.empty() )
   {
      return -1;
   }
   else
   {
      // std::set is required to be sorted by the standard.
      // ::begin() returns the smallest element.
      return *containingProcessRanks.begin();
   }
}

void PrimitiveStorageUnresolvedParticlesInterface::periodicallyMapToDomain( Vec3& pt ) const
{
   WALBERLA_UNUSED( pt );
}

std::vector< uint_t > PrimitiveStorageUnresolvedParticlesInterface::getNeighborProcesses() const
{
   const std::set< uint_t > neighboringRanks = primitiveStorage_->getNeighboringRanks();
   const std::vector        neighboringRanksVector( neighboringRanks.begin(), neighboringRanks.end() );
   return neighboringRanksVector;
}

bool PrimitiveStorageUnresolvedParticlesInterface::intersectsWithProcessSubdomain( const uint_t  rank,
                                                                                   const Vec3&   pt,
                                                                                   const real_t& radius ) const
{
   if ( rank != uint_c( walberla::mpi::MPIManager::instance()->rank() ) )
   {
      return false;
   }

   if ( !primitiveStorage_->hasGlobalCells() )
   {
      Point3D    pointOfInterest( pt[0], pt[1], pt[2] );
      const auto target              = mapFromPhysicalToComputationalDomain2D( primitiveStorage_, pointOfInterest, radius );
      const auto coordExists         = std::get< 0 >( target );
      const auto id                  = std::get< 1 >( target );
      const auto computationalCoords = std::get< 2 >( target );

      return coordExists && primitiveStorage_->faceExistsLocally( id );
   }
   else
   {
      Point3D    pointOfInterest( pt[0], pt[1], pt[2] );
      const auto target              = mapFromPhysicalToComputationalDomain3D( primitiveStorage_, pointOfInterest, radius );
      const auto coordExists         = std::get< 0 >( target );
      const auto id                  = std::get< 1 >( target );
      const auto computationalCoords = std::get< 2 >( target );

      return coordExists && primitiveStorage_->cellExistsLocally( id );
   }
   return false;
}

void PrimitiveStorageUnresolvedParticlesInterface::correctParticlePosition( Vec3& pt ) const
{
   WALBERLA_UNUSED( pt );
}

} // namespace hyteg
