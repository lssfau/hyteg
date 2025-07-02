/*
 * Copyright (c) 2025 Andreas Burkhart.
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
#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;

/// \brief Given a point P, returns the closest point in a triangle in 2D ( z coordinate = 0 ) defined by the points A, B, C.
inline Point3D closestPointTriangle2D( const Point3D& P, const Point3D& A, const Point3D& B, const Point3D& C )
{
   // Vectors
   const Point3D AB = B - A;
   const Point3D AC = C - A;
   const Point3D AP = P - A;

   // First we check if one of the vertices of the triangle is closest, as this is most probable.

   // We are interpreting the segments AB and AC as normals here and check if P lies in a wedge
   // defined by having an angle >= 90° to both AB and AC. If this is the case P is closest to A.
   const real_t ABdotAP = AB.dot( AP );
   const real_t ACdotAP = AC.dot( AP );
   if ( ABdotAP <= 0 && ACdotAP <= 0 )
   {
      return A;
   }

   // Check the wedge defined by having an angle >= 90° to both BA and BC.
   // We can flip the sign on BA.dot(BP) and check AB.dot(BP) >= 0 instead.
   // Further we would like to check BC.dot(BP) <= 0 which we can rewrite as
   // (C - A + A - B).dot(BP) <= 0 or (AC).dot(BP) <= AB.dot(BP).
   // This saves computational effort.
   const Point3D BP      = P - B;
   const real_t  ABdotBP = AB.dot( BP );
   const real_t  ACdotBP = AC.dot( BP );
   if ( ABdotBP >= 0 && ACdotBP <= ABdotBP )
   {
      return B;
   }

   // Check the wedge defined by having an angle >= 90° to both CA and CB.
   // We can flip the sign on CA.dot(CP) and check AC.dot(CP) >= 0 instead.
   // Further we would like to check CB.dot(CP) <= 0 which we can rewrite as
   // (B - A + A - C).dot(CP) <= 0 or (AB).dot(CP) <= AC.dot(CP).
   // This saves computational effort.
   const Point3D CP      = P - C;
   const real_t  ACdotCP = AC.dot( CP );
   const real_t  ABdotCP = AB.dot( CP );
   if ( ACdotCP >= 0 && ABdotCP <= ACdotCP )
   {
      return C;
   }

   // Now we want to check if P is closest to one of the edges.

   // For example we can check wether a point lies on one side of the line
   // through AB by checking the signed area of the triangle ABP.
   // If this area is <= 0 then the point lies on the outward facing side
   // ( w.r.t. the outward facing normals of ABC ) or on the edge.
   // In addition we want the angle between AP and AB as well as the angle
   // between BP and BA to be <= 90°. Following the above this can be
   // checked via AB.dot(AP) >= 0 and AB.dot(BP) <= 0.
   //
   // For example:
   //   normalDir = AB.cross( AC )
   //   signedAreaSubTriangleCAP = AB.cross( P - A ).dot( normalDir )
   //   signedAreaSubTriangleBCP = BC.cross( P - B ).dot( normalDir )
   //   signedAreaSubTriangleABP = CA.cross( P - C ).dot( normalDir )
   //   signedAreaTriangle = signedAreaSubTriangleCAP + signedAreaSubTriangleBCP + signedAreaSubTriangleABP
   //
   // Barycentric coordinates:
   //   c1 = signedAreaSubTriangleBCP / signedAreaTriangle
   //   c2 = signedAreaSubTriangleCAP / signedAreaTriangle
   //   c3 = signedAreaSubTriangleABP / signedAreaTriangle
   //
   // Note that the area variables are scaled since we only need correct ratios!
   //
   // The cases AC and BC are treated in a similar fashion.
   const real_t signedAreaSubTriangleABP = ABdotAP * ACdotBP - ABdotBP * ACdotAP;
   if ( ABdotAP >= 0 and ABdotBP <= 0 and signedAreaSubTriangleABP <= 0 )
   {
      return A + ABdotAP / ( ABdotAP - ABdotBP ) * AB;
   }

   const real_t signedAreaSubTriangleCAP = ABdotCP * ACdotAP - ABdotAP * ACdotCP;
   if ( ACdotAP >= 0 and ACdotCP <= 0 and signedAreaSubTriangleCAP <= 0 )
   {
      return A + ACdotAP / ( ACdotAP - ACdotCP ) * AC;
   }

   // Here we want to check if ( AC.dot(BP) >= AB.dot(BP) ) and ( AB.dot(CP) >= AC.dot(CP) )
   // which we can rewrite as AC.dot(BP) >= AB.dot(BP)
   const real_t signedAreaSubTriangleBCP = ABdotBP * ACdotCP - ABdotCP * ACdotBP;
   if ( ACdotBP >= ABdotBP and ABdotCP >= ACdotCP and signedAreaSubTriangleBCP <= 0 )
   {
      return B + ( ABdotBP - ACdotBP ) / ( ACdotBP - ABdotBP + ABdotCP - ACdotCP ) * ( B - C );
   }

   // Now only the least likely case that P lies inside ABC remains.
   return P;
}

/// \brief Given a point P, returns the closest point in a triangle in 3D defined by the points A, B, C.
inline Point3D closestPointTriangle3D( const Point3D& P, const Point3D& A, const Point3D& B, const Point3D& C )
{
   // Follows the same logic as the 2D case ( see comments there ), but applies projection to the triangle plane if necessary.

   // Vectors
   const Point3D AB = B - A;
   const Point3D AC = C - A;
   const Point3D AP = P - A;

   // A is closest
   const real_t ABdotAP = AB.dot( AP );
   const real_t ACdotAP = AC.dot( AP );
   if ( ABdotAP <= 0 && ACdotAP <= 0 )
   {
      return A;
   }

   // B is closest
   const Point3D BP      = P - B;
   const real_t  ABdotBP = AB.dot( BP );
   const real_t  ACdotBP = AC.dot( BP );
   if ( ABdotBP >= 0 && ACdotBP <= ABdotBP )
   {
      return B;
   }

   // C is closest
   const Point3D CP      = P - C;
   const real_t  ACdotCP = AC.dot( CP );
   const real_t  ABdotCP = AB.dot( CP );
   if ( ACdotCP >= 0 && ABdotCP <= ACdotCP )
   {
      return C;
   }

   // AB is closest
   const real_t signedAreaSubTriangleABP = ABdotAP * ACdotBP - ABdotBP * ACdotAP;
   if ( ABdotAP >= 0 and ABdotBP <= 0 and signedAreaSubTriangleABP <= 0 )
   {
      return A + ABdotAP / ( ABdotAP - ABdotBP ) * AB;
   }

   // AC is closest
   const real_t signedAreaSubTriangleCAP = ABdotCP * ACdotAP - ABdotAP * ACdotCP;
   if ( ACdotAP >= 0 and ACdotCP <= 0 and signedAreaSubTriangleCAP <= 0 )
   {
      return A + ACdotAP / ( ACdotAP - ACdotCP ) * AC;
   }

   // BC is closest
   const real_t signedAreaSubTriangleBCP = ABdotBP * ACdotCP - ABdotCP * ACdotBP;
   if ( ACdotBP >= ABdotBP and ABdotCP >= ACdotCP and signedAreaSubTriangleBCP <= 0 )
   {
      return B + ( ABdotBP - ACdotBP ) / ( ACdotBP - ABdotBP + ABdotCP - ACdotCP ) * ( B - C );
   }

   const real_t signedAreaTriangle = signedAreaSubTriangleABP + signedAreaSubTriangleCAP + signedAreaSubTriangleBCP;

   // Now only the least likely case that the projection of P onto the triangle lies inside ABC remains.

   // Calculate projection to triangle plane via barycentric coordinates
   //   c1 * A + c2 * B + c3 * C
   // This now has to lie in ABC hence 0 <= c1,c2,c3 <= 1 and c1 + c2 + c3 = 1
   return ( signedAreaSubTriangleBCP * A + signedAreaSubTriangleCAP * B + signedAreaSubTriangleABP * C ) / signedAreaTriangle;
}

/// \brief Given a point P, returns the closest point in a tetrahedron in 3D defined by the points A, B, C, D.
inline Point3D
    closestPointTetrahedron3D( const Point3D& P, const Point3D& A, const Point3D& B, const Point3D& C, const Point3D& D )
{
   const auto normalABC = tetrahedronInwardNormal( A, B, C, D );
   const auto normalCAD = tetrahedronInwardNormal( C, A, D, B );
   const auto normalABD = tetrahedronInwardNormal( A, B, D, C );
   const auto normalBCD = tetrahedronInwardNormal( B, C, D, A );

   const auto distFace0 = distanceToPlane( P, A, normalABC );
   const auto distFace1 = distanceToPlane( P, A, normalCAD );
   const auto distFace2 = distanceToPlane( P, A, normalABD );
   const auto distFace3 = distanceToPlane( P, B, normalBCD );

   // Check if P is inside the tetrahedron
   if ( distFace0 >= 0 && distFace1 >= 0 && distFace2 >= 0 && distFace3 >= 0 )
   {
      return P;
   }

   // Now we check all faces that a visible from P and calculate the closest point on the face.
   // We are outside for at least one point, so now initialisation of closestPoint required.
   Point3D closestPoint;
   Point3D closestPointOnFace;
   real_t  minDist = std::numeric_limits< real_t >::max();
   real_t  distFace;

   // Check face ABC
   if ( distFace0 < 0 )
   {
      closestPointOnFace = closestPointTriangle3D( P, A, B, C );
      distFace           = ( P - closestPointOnFace ).norm();

      if ( distFace < minDist )
      {
         minDist      = distFace;
         closestPoint = closestPointOnFace;
      }
   }
   // Check face CAD
   if ( distFace1 < 0 )
   {
      closestPointOnFace = closestPointTriangle3D( P, C, A, D );
      distFace           = ( P - closestPointOnFace ).norm();

      if ( distFace < minDist )
      {
         minDist      = distFace;
         closestPoint = closestPointOnFace;
      }
   }
   // Check face ABD
   if ( distFace2 < 0 )
   {
      closestPointOnFace = closestPointTriangle3D( P, A, B, D );
      distFace           = ( P - closestPointOnFace ).norm();

      if ( distFace < minDist )
      {
         minDist      = distFace;
         closestPoint = closestPointOnFace;
      }
   }
   // Check face BCD
   if ( distFace3 < 0 )
   {
      closestPointOnFace = closestPointTriangle3D( P, B, C, D );
      distFace           = ( P - closestPointOnFace ).norm();

      if ( distFace < minDist )
      {
         minDist      = distFace;
         closestPoint = closestPointOnFace;
      }
   }

   // now return the closest point
   return closestPoint;
}

} // namespace hyteg
