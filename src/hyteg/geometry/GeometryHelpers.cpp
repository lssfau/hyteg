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

   Point2D coordinates2D( { computationalCoords[0], computationalCoords[1] } );

   for ( const auto& it : storage->getFaces() )
   {
      Face& face = *it.second;

      Point2D faceCoodinates0( { face.getCoordinates()[0][0], face.getCoordinates()[0][1] } );
      Point2D faceCoodinates1( { face.getCoordinates()[1][0], face.getCoordinates()[1][1] } );
      Point2D faceCoodinates2( { face.getCoordinates()[2][0], face.getCoordinates()[2][1] } );

      found = isPointInTriangle( coordinates2D, faceCoodinates0, faceCoodinates1, faceCoodinates2 );

      // #define GEOMETRY_HELPERS_BE_VERBOSE
#ifdef GEOMETRY_HELPERS_BE_VERBOSE
      WALBERLA_LOG_INFO_ON_ROOT( " -----------------------------------------------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( " -> FaceID = " << face.getID() );
      WALBERLA_LOG_INFO_ON_ROOT( " -> Point found = " << ( found ? "TRUE" : "FALSE" ) );
      WALBERLA_LOG_INFO_ON_ROOT( " -> Computational Coordinates = " << computationalCoords );
      WALBERLA_LOG_INFO_ON_ROOT( " -> Face Vertices:" );
      WALBERLA_LOG_INFO_ON_ROOT( " -> " << faceCoodinates0 );
      WALBERLA_LOG_INFO_ON_ROOT( " -> " << faceCoodinates1 );
      WALBERLA_LOG_INFO_ON_ROOT( " -> " << faceCoodinates2 );
      WALBERLA_LOG_INFO_ON_ROOT( " -----------------------------------------------------------" );
#endif

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

         Point2D faceCoodinates0( { face.getCoordinates()[0][0], face.getCoordinates()[0][1] } );
         Point2D faceCoodinates1( { face.getCoordinates()[1][0], face.getCoordinates()[1][1] } );
         Point2D faceCoodinates2( { face.getCoordinates()[2][0], face.getCoordinates()[2][1] } );

         found = circleTriangleIntersection(
             coordinates2D, searchToleranceRadius, faceCoodinates0, faceCoodinates1, faceCoodinates2 );

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
