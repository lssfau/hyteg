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

#include "SetupPrimitiveStorageConvectionParticlesInterface.hpp"

#include <algorithm>

#include "hyteg/geometry/Intersection.hpp"

#include "convection_particles/data/DataTypes.h"
#include "convection_particles/domain/IDomain.h"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;
using walberla::convection_particles::Vec3;

bool SetupPrimitiveStorageConvectionParticlesInterface::isContainedInProcessSubdomain( const uint_t rank, const Vec3& pt ) const
{
   const int containingRank = findContainingProcessRank( pt );
   return (int) rank == containingRank;
}

int SetupPrimitiveStorageConvectionParticlesInterface::findContainingProcessRank( const Vec3& pt ) const
{
   std::set< int > containingProcessRanks;
   if ( setupStorage_->getNumberOfCells() == 0 )
   {
      for ( const auto& faceIt : setupStorage_->getFaces() )
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
            containingProcessRanks.insert( static_cast< int >( setupStorage_->getTargetRank( faceID ) ) );
         }
      }
   }
   else
   {
      for ( const auto& cellIt : setupStorage_->getCells() )
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
            containingProcessRanks.insert( static_cast< int >( setupStorage_->getTargetRank( cellID ) ) );
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

void SetupPrimitiveStorageConvectionParticlesInterface::periodicallyMapToDomain( Vec3& pt ) const
{
   WALBERLA_UNUSED( pt );
}

std::vector< uint_t > SetupPrimitiveStorageConvectionParticlesInterface::getNeighborProcesses() const
{
   std::set< uint_t > neighborProcesses;
   if ( setupStorage_->getNumberOfCells() == 0 )
   {
      for ( const auto& faceIt : setupStorage_->getFaces() )
      {
         auto face = faceIt.second;

         for ( const auto& neighborVertexID : face->neighborVertices() )
         {
            auto neighborVertex = setupStorage_->getVertex( neighborVertexID );
            for ( const auto &neighborFaceID : neighborVertex->neighborFaces() )
            {
               neighborProcesses.insert( setupStorage_->getTargetRank( neighborFaceID ) );
            }
         }
      }
   }
   else
   {
      for ( const auto& cellIt : setupStorage_->getCells() )
      {
         auto cell = cellIt.second;

         for ( const auto& neighborVertexID : cell->neighborVertices() )
         {
            auto neighborVertex = setupStorage_->getVertex( neighborVertexID );
            for ( const auto &neighborCellID : neighborVertex->neighborCells() )
            {
               neighborProcesses.insert( setupStorage_->getTargetRank( neighborCellID ) );
            }
         }
      }
   }
   neighborProcesses.erase( uint_c( walberla::mpi::MPIManager::instance()->rank() ) );

   return std::vector< uint_t >( neighborProcesses.begin(), neighborProcesses.end() );
}

bool SetupPrimitiveStorageConvectionParticlesInterface::intersectsWithProcessSubdomain( const uint_t  rank,
                                                                                        const Vec3&   pt,
                                                                                        const real_t& radius ) const
{
   if ( setupStorage_->getNumberOfCells() == 0 )
   {
      for ( const auto& faceIt : setupStorage_->getFaces() )
      {
         auto faceID = faceIt.first;
         auto face   = faceIt.second;

         Point2D pointOfInterest( {pt[0], pt[1]} );

         if ( setupStorage_->getTargetRank( faceID ) == rank &&
              circleTriangleIntersection( pointOfInterest,
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
      for ( const auto& cellIt : setupStorage_->getCells() )
      {
         auto cellID = cellIt.first;
         auto cell   = cellIt.second;

         Point3D pointOfInterest( {pt[0], pt[1], pt[2]} );

         if ( setupStorage_->getTargetRank( cellID ) == rank && sphereTetrahedronIntersection( pointOfInterest,
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

void SetupPrimitiveStorageConvectionParticlesInterface::correctParticlePosition( Vec3& pt ) const
{
   WALBERLA_UNUSED( pt );
}

} // namespace hyteg
