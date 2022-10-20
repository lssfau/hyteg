/*
 * Copyright (c) 2022 Nils Kohl, Marcus Mohr.
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

#include "hyteg/geometry/GeometryHelpers.hpp"

#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Face.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;

/// Find macro-face associated with a given point in 2D
std::tuple< bool, PrimitiveID > findFaceIDForPointIn2D( std::shared_ptr< PrimitiveStorage > storage,
                                                        const Point3D&                      computationalCoords,
                                                        real_t                              searchToleranceRadius )
{
   bool        found = false;
   PrimitiveID faceID;

   for ( const auto& it : storage->getFaces() )
   {
      Face& face = *it.second;

      found =
          isPointInTriangle( computationalCoords, face.getCoordinates()[0], face.getCoordinates()[1], face.getCoordinates()[2] );

      // leave on first hit
      if ( found )
      {
         faceID = face.getID();
         break;
      }
   }

   // No face found? Try different approach
   if ( !found && searchToleranceRadius > real_c( 0 ) )
   {
      for ( const auto& it : storage->getFaces() )
      {
         Face& face = *it.second;

         found = circleTriangleIntersection(
             computationalCoords, searchToleranceRadius, face.getCoordinates()[0], face.getCoordinates()[1], face.getCoordinates()[2] );

         if ( found )
         {
            faceID = face.getID();
            break;
         }
      }
   }

   return { found, faceID };
}

/// Find macro-cell associated with a given point in 3D
std::tuple< bool, PrimitiveID > findCellIDForPointIn3D( std::shared_ptr< PrimitiveStorage > storage,
                                                        const Point3D&                      computationalCoords,
                                                        real_t                              searchToleranceRadius )
{
   bool        found = false;
   PrimitiveID cellID;

   for ( auto& it : storage->getCells() )
   {
      Cell& cell = *it.second;

      found = isPointInTetrahedron( computationalCoords,
                                    cell.getCoordinates()[0],
                                    cell.getCoordinates()[1],
                                    cell.getCoordinates()[2],
                                    cell.getCoordinates()[3],
                                    cell.getFaceInwardNormal( 0 ),
                                    cell.getFaceInwardNormal( 1 ),
                                    cell.getFaceInwardNormal( 2 ),
                                    cell.getFaceInwardNormal( 3 ) );
      if ( found )
      {
         cellID = cell.getID();
         break;
      }
   }

   // No cell found? Try different approach
   if ( !found && searchToleranceRadius > real_c( 0 ) )
   {
      for ( auto& it : storage->getCells() )
      {
         Cell& cell = *it.second;

         found = sphereTetrahedronIntersection( computationalCoords,
                                                searchToleranceRadius,
                                                cell.getCoordinates()[0],
                                                cell.getCoordinates()[1],
                                                cell.getCoordinates()[2],
                                                cell.getCoordinates()[3] );
         if ( found )
         {
            cellID = cell.getID();
            break;
         }
      }
   }

   return { found, cellID };
}

} // namespace hyteg
