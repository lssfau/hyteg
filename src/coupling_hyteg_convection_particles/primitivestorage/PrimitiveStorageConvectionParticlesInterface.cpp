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

#include "PrimitiveStorageConvectionParticlesInterface.hpp"

#include <algorithm>

#include "hyteg/geometry/Intersection.hpp"

#include "convection_particles/data/DataTypes.h"
#include "convection_particles/domain/IDomain.h"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;
using walberla::convection_particles::Vec3;

bool PrimitiveStorageConvectionParticlesInterface::isContainedInProcessSubdomain( const uint_t rank, const Vec3& pt ) const
{
   const int containingRank = findContainingProcessRank( pt );
   return (int) rank == containingRank;
}

int PrimitiveStorageConvectionParticlesInterface::findContainingProcessRank( const Vec3& pt ) const
{
   std::set< int > containingProcessRanks;
   if ( !primitiveStorage_->hasGlobalCells() )
   {
      auto allFaces = primitiveStorage_->getFaces();
      auto neighborFaces = primitiveStorage_->getNeighborFaces();
      allFaces.insert( neighborFaces.begin(), neighborFaces.end() );

      for ( const auto& faceIt : allFaces )
      {
         auto faceID = faceIt.first;
         auto face   = faceIt.second;

         Point2D pointOfInterest( {pt[0], pt[1]} );

         if ( circleTriangleIntersection( pointOfInterest,
                                          1e-05,
                                          Point2D( {face->getCoordinates().at( 0 )[0], face->getCoordinates().at( 0 )[1]} ),
                                          Point2D( {face->getCoordinates().at( 1 )[0], face->getCoordinates().at( 1 )[1]} ),
                                          Point2D( {face->getCoordinates().at( 2 )[0], face->getCoordinates().at( 2 )[1]} ) ) )
         {
            containingProcessRanks.insert( static_cast< int >( primitiveStorage_->getPrimitiveRank( faceID ) ) );
         }
      }
   }
   else
   {
      auto allCells = primitiveStorage_->getCells();
      auto neighborCells = primitiveStorage_->getNeighborCells();
      allCells.insert( neighborCells.begin(), neighborCells.end() );

      for ( const auto& cellIt : allCells )
      {
         auto cellID = cellIt.first;
         auto cell   = cellIt.second;

         Point3D pointOfInterest( {pt[0], pt[1], pt[2]} );

         if ( sphereTetrahedronIntersection( pointOfInterest,
                                             1e-05,
                                             cell->getCoordinates().at( 0 ),
                                             cell->getCoordinates().at( 1 ),
                                             cell->getCoordinates().at( 2 ),
                                             cell->getCoordinates().at( 3 ) ) )
         {
            containingProcessRanks.insert( static_cast< int >( primitiveStorage_->getPrimitiveRank( cellID ) ) );
         }
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

void PrimitiveStorageConvectionParticlesInterface::periodicallyMapToDomain( Vec3& pt ) const
{
   WALBERLA_UNUSED( pt );
}

std::vector< uint_t > PrimitiveStorageConvectionParticlesInterface::getNeighborProcesses() const
{
   std::set< uint_t > neighborProcesses;
   if ( !primitiveStorage_->hasGlobalCells() )
   {
      for ( auto it : primitiveStorage_->getNeighborFaces() )
      {
         auto pID = it.first;
         neighborProcesses.insert( primitiveStorage_->getPrimitiveRank( pID ) );
      }
   }
   else
   {
      for ( auto it : primitiveStorage_->getNeighborCells() )
      {
         auto pID = it.first;
         neighborProcesses.insert( primitiveStorage_->getPrimitiveRank( pID ) );
      }
   }
   neighborProcesses.erase( uint_c( walberla::mpi::MPIManager::instance()->rank() ) );

   return std::vector< uint_t >( neighborProcesses.begin(), neighborProcesses.end() );
}

bool PrimitiveStorageConvectionParticlesInterface::intersectsWithProcessSubdomain( const uint_t  rank,
                                                                                   const Vec3&   pt,
                                                                                   const real_t& radius ) const
{
   if ( rank != uint_c( walberla::mpi::MPIManager::instance()->rank() ) )
   {
      return false;
   }

   if ( !primitiveStorage_->hasGlobalCells() )
   {
      for ( const auto& faceIt : primitiveStorage_->getFaces() )
      {
         auto face   = faceIt.second;

         Point2D pointOfInterest( {pt[0], pt[1]} );

         if ( circleTriangleIntersection( pointOfInterest,
                                          radius,
                                          Point2D( {face->getCoordinates().at( 0 )[0], face->getCoordinates().at( 0 )[1]} ),
                                          Point2D( {face->getCoordinates().at( 1 )[0], face->getCoordinates().at( 1 )[1]} ),
                                          Point2D( {face->getCoordinates().at( 2 )[0], face->getCoordinates().at( 2 )[1]} ) ) )
         {
            return true;
         }
      }
   }
   else
   {
      for ( const auto& cellIt : primitiveStorage_->getCells() )
      {
         auto cell   = cellIt.second;

         Point3D pointOfInterest( {pt[0], pt[1], pt[2]} );

         if ( sphereTetrahedronIntersection( pointOfInterest,
                                             radius,
                                             cell->getCoordinates().at( 0 ),
                                             cell->getCoordinates().at( 1 ),
                                             cell->getCoordinates().at( 2 ),
                                             cell->getCoordinates().at( 3 ) ) )
         {
            return true;
         }
      }
   }
   return false;
}

void PrimitiveStorageConvectionParticlesInterface::correctParticlePosition( Vec3& pt ) const
{
   WALBERLA_UNUSED( pt );
}

} // namespace hyteg
