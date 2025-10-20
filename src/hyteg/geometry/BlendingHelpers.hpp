/*
 * Copyright (c) 2022-2025 Nils Kohl, Marcus Mohr, Andreas Burkhart.
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

#include "core/DataTypes.h"
#include "core/math/Matrix3.h"

#include "hyteg/geometry/ClosestPoint.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;

/// Map point from 2D physical to 2D computational domain
///
/// Given the coordinates of a point in the physical domain try to figure out what its coordinates
/// were before blending, i.e. on the computational domain. Use this function in those cases were
/// you have no information to which macro-face the point belongs on the computational domain.
/// Figuring this out is then part of the task.
///
/// We apply three approaches to solve the problem:
///
///   1. For all face primitives of the local subdomain:
///      If a point-triangle inclusion test succeeds and the verifyPointPairing() method of
///      the corresponding face returns true, then we return true, the ID of the face and the
///      coordinates obtained from the inverse blending map of that face. While checking all
///      face primitives, we store the face with the smallest distance to the point that was mapped
///      back to the computational domain that also fulfills the verifyPointPairing() check.
///      If a face's distance is smaller than the given distanceTolerance parameter,
///      then the method can also return that face & respective computational domain point immediately.
///      Use distanceTolerance = 0 to disable that feature.
///
///   2. If approach #1 fails and searchToleranceRadius is positive, we perform for all face
///      primitives of the local subdomain a circle-triangle intersection test. If this
///      is successful and the verifyPointPairing() method of the corresponding face returns
///      true, then we return true, the ID of the face and the coordinates obtained from the
///      inverse blending map of that face.
///
///   3. If approach #2 fails and returnBestGuess is true, the face & respective computational domain point
///      fulfilling the verifyPointPairing() check with the smallest point face distance are returned,
///      provided we found one in step 1.
///
/// Additionally, if includeNeighboringFaces is true, then we do not only look for a process-local primitive, but also
/// search in the set of neighboring primitives.
///
inline std::tuple< bool, PrimitiveID, Point3D >
    mapFromPhysicalToComputationalDomain2D( std::shared_ptr< PrimitiveStorage > storage,
                                            const Point3D&                      physicalCoords,
                                            real_t                              searchToleranceRadius,
                                            real_t                              distanceTolerance       = real_c( 0 ),
                                            bool                                returnBestGuess         = false,
                                            bool                                includeNeighboringFaces = false )
{
   bool        foundCandidate = false;
   PrimitiveID faceID;
   Point3D     computationalCoords;

   real_t lowestTol = std::numeric_limits< real_t >::max();

   auto allFaces = storage->getFaces();

   if ( includeNeighboringFaces )
   {
      auto neighborFaces = storage->getNeighborFaces();
      allFaces.insert( neighborFaces.begin(), neighborFaces.end() );
   }

   for ( const auto& it : allFaces )
   {
      Face& face = *it.second;

      Point3D currentComputationalCoords;

      // map coordinates from physical to computational domain
      face.getGeometryMap()->evalFinv( physicalCoords, currentComputationalCoords );

      real_t dist = ( currentComputationalCoords - closestPointTriangle2D( currentComputationalCoords,
                                                                           face.getCoordinates()[0],
                                                                           face.getCoordinates()[1],
                                                                           face.getCoordinates()[2] ) )
                        .norm();

      if ( dist < lowestTol && face.getGeometryMap()->verifyPointPairing( currentComputationalCoords, physicalCoords ) )
      {
         foundCandidate      = true;
         lowestTol           = dist;
         faceID              = face.getID();
         computationalCoords = currentComputationalCoords;

         if ( std::fpclassify( dist ) == FP_ZERO || ( distanceTolerance > real_c( 0 ) && dist < distanceTolerance ) )
         {
            return { true, faceID, computationalCoords };
         }
      }
   }

   // No face found? Try different approach
   if ( searchToleranceRadius > real_c( 0 ) )
   {
      for ( const auto& it : allFaces )
      {
         Face& face = *it.second;

         // map coordinates from physical to computational domain
         face.getGeometryMap()->evalFinv( physicalCoords, computationalCoords );

         bool faceIsCandidate = circleTriangleIntersection( computationalCoords,
                                                            searchToleranceRadius,
                                                            face.getCoordinates()[0],
                                                            face.getCoordinates()[1],
                                                            face.getCoordinates()[2] );

         if ( faceIsCandidate && face.getGeometryMap()->verifyPointPairing( computationalCoords, physicalCoords ) )
         {
            faceID = face.getID();
            return { true, faceID, computationalCoords };
         }
      }
   }

   if ( foundCandidate && returnBestGuess )
   {
      return { true, faceID, computationalCoords };
   }

   return { false, faceID, computationalCoords };
}

/// Map point from 3D physical to 3D computational domain
///
/// Given the coordinates of a point in the physical domain try to figure out what its coordinates
/// were before blending, i.e. on the computational domain. Use this function in those cases were
/// you have no information to which macro-cell the point belongs on the computational domain.
/// Figuring this out is then part of the task.
///
/// We apply three approaches to solve the problem:
///
///   1. For all cell primitives of the local subdomain:
///      If a point-triangle inclusion test succeeds and the verifyPointPairing() method of
///      the corresponding cell returns true, then we return true, the ID of the cell and the
///      coordinates obtained from the inverse blending map of that cell. While checking all
///      cell primitives, we store the cell with the smallest distance to the point that was mapped
///      back to the computational domain that also fulfills the verifyPointPairing() check.
///      If a cell's distance is smaller than the given distanceTolerance parameter,
///      then the method can also return that cell & respective computational domain point immediately.
///      Use distanceTolerance = 0 to disable that feature.
///
///   2. If approach #1 fails and searchToleranceRadius is positive, we perform for all cell
///      primitives of the local subdomain a sphere-tetrahedron intersection test. If this
///      is successful and the verifyPointPairing() method of the corresponding cell returns
///      true, then we return true, the ID of the cell and the coordinates obtained from the
///      inverse blending map of that cell.
///
///   3. If approach #2 fails and returnBestGuess is true, the cell & respective computational domain point
///      fulfilling the verifyPointPairing() check with the smallest point cell distance are returned,
///      provided we found one in step 1.
///
/// Additionally, if includeNeighboringFaces is true, then we do not only look for a process-local primitive, but also
/// search in the set of neighboring primitives.
///
inline std::tuple< bool, PrimitiveID, Point3D >
    mapFromPhysicalToComputationalDomain3D( std::shared_ptr< PrimitiveStorage > storage,
                                            const Point3D&                      physicalCoords,
                                            real_t                              searchToleranceRadius,
                                            real_t                              distanceTolerance       = real_c( 0 ),
                                            bool                                returnBestGuess         = false,
                                            bool                                includeNeighboringCells = false )
{
   bool        foundCandidate = false;
   PrimitiveID cellID;
   Point3D     computationalCoords;

   real_t lowestTol = std::numeric_limits< real_t >::max();

   auto allCells = storage->getCells();

   if ( includeNeighboringCells )
   {
      auto neighborCells = storage->getNeighborCells();
      allCells.insert( neighborCells.begin(), neighborCells.end() );
   }

   for ( const auto& it : allCells )
   {
      Cell& cell = *it.second;

      Point3D currentComputationalCoords;

      // map coordinates from physical to computational domain
      cell.getGeometryMap()->evalFinv( physicalCoords, currentComputationalCoords );

      real_t dist = ( currentComputationalCoords - closestPointTetrahedron3D( currentComputationalCoords,
                                                                              cell.getCoordinates()[0],
                                                                              cell.getCoordinates()[1],
                                                                              cell.getCoordinates()[2],
                                                                              cell.getCoordinates()[3] ) )
                        .norm();

      if ( dist < lowestTol && cell.getGeometryMap()->verifyPointPairing( currentComputationalCoords, physicalCoords ) )
      {
         foundCandidate      = true;
         lowestTol           = dist;
         cellID              = cell.getID();
         computationalCoords = currentComputationalCoords;

         if ( std::fpclassify( dist ) == FP_ZERO || ( distanceTolerance > real_c( 0 ) && dist < distanceTolerance ) )
         {
            return { true, cellID, computationalCoords };
         }
      }
   }

   // No cell found? Try different approach
   if ( searchToleranceRadius > real_c( 0 ) )
   {
      for ( auto& it : allCells )
      {
         Cell& cell = *it.second;

         // map coordinates from physical to computational domain
         cell.getGeometryMap()->evalFinv( physicalCoords, computationalCoords );

         bool cellIsCandidate = sphereTetrahedronIntersection( computationalCoords,
                                                               searchToleranceRadius,
                                                               cell.getCoordinates()[0],
                                                               cell.getCoordinates()[1],
                                                               cell.getCoordinates()[2],
                                                               cell.getCoordinates()[3] );

         // if cell is candidate, check that this is actually the correct cell
         if ( cellIsCandidate && cell.getGeometryMap()->verifyPointPairing( computationalCoords, physicalCoords ) )
         {
            cellID = cell.getID();
            return { true, cellID, computationalCoords };
         }
      }
   }

   if ( foundCandidate && returnBestGuess )
   {
      return { true, cellID, computationalCoords };
   }

   return { false, cellID, computationalCoords };
}

} // namespace hyteg
