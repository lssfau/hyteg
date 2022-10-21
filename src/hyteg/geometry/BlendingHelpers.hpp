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

#pragma once

#include "core/DataTypes.h"
#include "core/math/Matrix3.h"

#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/types/pointnd.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;

/// Map point from 2D physical to 2D computational domain
///
/// Given the coordinates of a point in the physical domain try to figure out what its coordinates
/// were before blending, i.e. on the computational domain. Use this function in those cases were
/// you have no information to which macro-face the point belongs on the computational domain.
/// Figuring this out is then part of the task.
std::tuple< bool, PrimitiveID, Point3D > mapFromPhysicalToComputationalDomain2D( std::shared_ptr< PrimitiveStorage > storage,
                                                                                 Point3D& physicalCoords,
                                                                                 real_t   searchToleranceRadius )
{
   bool        found = false;
   PrimitiveID faceID;
   Point3D     computationalCoords;

   for ( const auto& it : storage->getFaces() )
   {
      Face& face = *it.second;

      // map coordinates from physical to computational domain
      face.getGeometryMap()->evalFinv( physicalCoords, computationalCoords );

      bool faceIsCandidate =
          isPointInTriangle( computationalCoords, face.getCoordinates()[0], face.getCoordinates()[1], face.getCoordinates()[2] );

      // #define BE_VERBOSE
#ifdef BE_VERBOSE
      WALBERLA_LOG_INFO_ON_ROOT( " -----------------------------------------------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( " -> FaceID = " << faceID );
      WALBERLA_LOG_INFO_ON_ROOT( " -> Face is candidate = " << ( faceIsCandidate ? "TRUE" : "FALSE" ) );
      WALBERLA_LOG_INFO_ON_ROOT( " -> Physical Coordinates = " << physicalCoords );
      WALBERLA_LOG_INFO_ON_ROOT( " -> Computational Coordinates = " << computationalCoords );
      WALBERLA_LOG_INFO_ON_ROOT( " -> Face Vertices:" );
      WALBERLA_LOG_INFO_ON_ROOT( " -> " << faceCoodinates0 );
      WALBERLA_LOG_INFO_ON_ROOT( " -> " << faceCoodinates1 );
      WALBERLA_LOG_INFO_ON_ROOT( " -> " << faceCoodinates2 );
      WALBERLA_LOG_INFO_ON_ROOT( " -----------------------------------------------------------" );
#endif

      // if face is candidate, check that this is actually the correct face
      if ( faceIsCandidate && face.getGeometryMap()->verifyPointPairing( computationalCoords, physicalCoords ) )
      {
         found  = true;
         faceID = face.getID();
         break;
      }
   }

   return { found, faceID, computationalCoords };
}

/// Map point from 3D physical to 3D computational domain
///
/// Given the coordinates of a point in the physical domain try to figure out what its coordinates
/// were before blending, i.e. on the computational domain. Use this function in those cases were
/// you have no information to which macro-cell the point belongs on the computational domain.
/// Figuring this out is then part of the task.
std::tuple< bool, PrimitiveID, Point3D > mapFromPhysicalToComputationalDomain3D( std::shared_ptr< PrimitiveStorage > storage,
                                                                                 Point3D& physicalCoords,
                                                                                 real_t   searchToleranceRadius )
{
   bool        found = false;
   PrimitiveID cellID;
   Point3D     computationalCoords;

   for ( const auto& it : storage->getCells() )
   {
      Cell& cell = *it.second;

      // map coordinates from physical to computational domain
      cell.getGeometryMap()->evalFinv( physicalCoords, computationalCoords );

      bool cellIsCandidate = isPointInTetrahedron( computationalCoords,
                                                   cell.getCoordinates()[0],
                                                   cell.getCoordinates()[1],
                                                   cell.getCoordinates()[2],
                                                   cell.getCoordinates()[3],
                                                   cell.getFaceInwardNormal( 0 ),
                                                   cell.getFaceInwardNormal( 1 ),
                                                   cell.getFaceInwardNormal( 2 ),
                                                   cell.getFaceInwardNormal( 3 ) );

      // if cell is candidate, check that this is actually the correct cell
      if ( cellIsCandidate && cell.getGeometryMap()->verifyPointPairing( computationalCoords, physicalCoords ) )
      {
         found  = true;
         cellID = cell.getID();
         break;
      }
   }

   return { found, cellID, computationalCoords };
}

} // namespace hyteg
