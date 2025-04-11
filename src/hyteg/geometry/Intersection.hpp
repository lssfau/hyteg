/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;

/// \brief Returns the distance to a half plane in 3D.
/// If the returned distance is positive, the point of interest is located on the side of the plane the normal point to.
inline real_t distanceToPlane( const Point3D& pointOfInterest, const Point3D& pointOnPlane, const Point3D& planeNormal )
{
   return ( planeNormal / planeNormal.norm() ).dot( pointOfInterest - pointOnPlane );
}

/// \brief Returns the absolute distance to a half plane in 3D.
inline real_t distanceToPlane( const Point3D& planeVertex0,
                               const Point3D& planeVertex1,
                               const Point3D& planeVertex2,
                               const Point3D& opposingVertex )
{
   auto v      = opposingVertex - planeVertex0;
   auto normal = crossProduct( planeVertex1 - planeVertex0, planeVertex2 - planeVertex0 );
   normal /= normal.norm();
   auto dist           = normal.dot( v );
   auto projectedPoint = opposingVertex - dist * normal;
   auto inwardNormal   = opposingVertex - projectedPoint;
   return inwardNormal.norm();
}

/// \brief Returns the normal of the passed triangle in the direction of the opposite (4th) vertex of the tet.
///
/// This function may return a weird normal if the "opposing" vertex lies almost on the plane.
/// Use distanceToPlane() to be sure that this is not the case.
inline Point3D tetrahedronInwardNormal( const Point3D& planeVertex0,
                                        const Point3D& planeVertex1,
                                        const Point3D& planeVertex2,
                                        const Point3D& opposingVertex )
{
   auto v      = opposingVertex - planeVertex0;
   auto normal = crossProduct( planeVertex1 - planeVertex0, planeVertex2 - planeVertex0 );
   normal /= normal.norm();
   auto dist           = normal.dot( v );
   auto projectedPoint = opposingVertex - dist * normal;
   auto inwardNormal   = opposingVertex - projectedPoint;
   if ( inwardNormal.norm() < 1e-15 )
   {
      // The opposing vertex seems to be located on the plane of the triangle.
      // Let's return any normal..
      return normal;
   }
   else
   {
      return inwardNormal / inwardNormal.norm();
   }
}

/// Check whether in a planar geometry a point belongs to a triangle
/// \note If using Point3D as PointType the z-component will be ignored
template < typename PointType >
inline bool isPointInTriangle( const PointType& pointOfInterest, const PointType& v1, const PointType& v2, const PointType& v3 )
{
   if constexpr ( !std::is_same< PointType, Point2D >::value && !std::is_same< PointType, Point3D >::value )
   {
      WALBERLA_UNUSED( pointOfInterest );
      WALBERLA_UNUSED( v1 );
      WALBERLA_UNUSED( v2 );
      WALBERLA_UNUSED( v3 );
      return false;
   }
   else
   {
      const auto v1x     = v1[0];
      const auto v1y     = v1[1];
      const auto v2x     = v2[0];
      const auto v2y     = v2[1];
      const auto v3x     = v3[0];
      const auto v3y     = v3[1];
      const auto centrex = pointOfInterest[0];
      const auto centrey = pointOfInterest[1];

      const auto area = real_c( 0.5 ) * ( -v2y * v3x + v1y * ( -v2x + v3x ) + v1x * ( v2y - v3y ) + v2x * v3y );
      const auto s =
          real_c( 1 ) / ( real_c( 2 ) * area ) * ( v1y * v3x - v1x * v3y + ( v3y - v1y ) * centrex + ( v1x - v3x ) * centrey );
      const auto t =
          real_c( 1 ) / ( real_c( 2 ) * area ) * ( v1x * v2y - v1y * v2x + ( v1y - v2y ) * centrex + ( v2x - v1x ) * centrey );
      return ( s > real_c( 0 ) && t > real_c( 0 ) && real_c( 1 ) - s - t > real_c( 0 ) );
   }
}

/// \brief Returns true if the passed circle and triangle intersect.
/// \note If using Point3D as PointType the z-component will be ignored
template < typename PointType >
inline bool circleTriangleIntersection( const PointType& centre,
                                        const real_t&    radius,
                                        const PointType& v1,
                                        const PointType& v2,
                                        const PointType& v3 )
{
   if constexpr ( !std::is_same< PointType, Point2D >::value && !std::is_same< PointType, Point3D >::value )
   {
      WALBERLA_UNUSED( centre );
      WALBERLA_UNUSED( radius );
      WALBERLA_UNUSED( v1 );
      WALBERLA_UNUSED( v2 );
      WALBERLA_UNUSED( v3 );
      return false;
   }
   else
   {
      const auto v1x     = v1[0];
      const auto v1y     = v1[1];
      const auto v2x     = v2[0];
      const auto v2y     = v2[1];
      const auto v3x     = v3[0];
      const auto v3y     = v3[1];
      const auto centrex = centre[0];
      const auto centrey = centre[1];

      //
      // TEST 1: Vertex within circle
      //
      auto c1x = v1x - centrex;
      auto c1y = v1y - centrey;

      if ( std::sqrt( c1x * c1x + c1y * c1y ) <= radius )
         return true;

      auto c2x = v2x - centrex;
      auto c2y = v2y - centrey;

      if ( std::sqrt( c2x * c2x + c2y * c2y ) <= radius )
         return true;

      auto c3x = v3x - centrex;
      auto c3y = v3y - centrey;

      if ( std::sqrt( c3x * c3x + c3y * c3y ) <= radius )
         return true;

      //
      // TEST 2: Circle centre within triangle
      //
      if ( isPointInTriangle( centre, v1, v2, v3 ) )
         return true;

      //
      // TEST 3: Circle intersects edge
      //
      // Get the dot product...
      //
      c1x      = centrex - v1x;
      c1y      = centrey - v1y;
      auto e1x = v2x - v1x;
      auto e1y = v2y - v1y;

      auto k = c1x * e1x + c1y * e1y;

      if ( k > 0 )
      {
         auto len = std::sqrt( e1x * e1x + e1y * e1y );
         k        = k / len;

         if ( k < len )
         {
            if ( std::sqrt( c1x * c1x + c1y * c1y - k * k ) <= radius )
               return true;
         }
      }

      // Second edge
      c2x      = centrex - v2x;
      c2y      = centrey - v2y;
      auto e2x = v3x - v2x;
      auto e2y = v3y - v2y;

      k = c2x * e2x + c2y * e2y;

      if ( k > 0 )
      {
         auto len = std::sqrt( e2x * e2x + e2y * e2y );
         k        = k / len;

         if ( k < len )
         {
            if ( std::sqrt( c2x * c2x + c2y * c2y - k * k ) <= radius )
               return true;
         }
      }

      // Third edge
      c3x      = centrex - v3x;
      c3y      = centrey - v3y;
      auto e3x = v1x - v3x;
      auto e3y = v1y - v3y;

      k = c3x * e3x + c3y * e3y;

      if ( k > 0 )
      {
         auto len = std::sqrt( e3x * e3x + e3y * e3y );
         k        = k / len;

         if ( k < len )
         {
            if ( std::sqrt( c3x * c3x + c3y * c3y - k * k ) <= radius )
               return true;
         }
      }

      // We're done, no intersection
      return false;
   }
}

inline bool sphereTriangleIntersection( const Point3D& centre,
                                        const real_t&  radius,
                                        const Point3D& v1,
                                        const Point3D& v2,
                                        const Point3D& v3 )
{
   // If the distance of the sphere to the triangle plane is greater than the sphere's radius there is no intersection.
   // Otherwise, check intersection of the triangle with the circle that is common with the plane (where the plane cuts through the sphere).
   auto centreDistToPlane = distanceToPlane( v1, v2, v3, centre );
   if ( centreDistToPlane > radius )
      return false;

   auto planeNormal = tetrahedronInwardNormal( v1, v2, v3, centre );

   auto intersectionRadius = std::sqrt( radius * radius - centreDistToPlane * centreDistToPlane );

   // We need to rotate our coordinate system -> normal component in z direction.
   Point3D planeTangent0 = v2 - v1;
   planeTangent0 /= planeTangent0.norm();
   auto planeTangent1 = crossProduct( planeNormal, planeTangent0 );

   walberla::math::Matrix3< real_t > basisTrafo;

   basisTrafo( 0, 0 ) = planeTangent0[0];
   basisTrafo( 1, 0 ) = planeTangent0[1];
   basisTrafo( 2, 0 ) = planeTangent0[2];

   basisTrafo( 0, 1 ) = planeTangent1[0];
   basisTrafo( 1, 1 ) = planeTangent1[1];
   basisTrafo( 2, 1 ) = planeTangent1[2];

   basisTrafo( 0, 2 ) = planeNormal[0];
   basisTrafo( 1, 2 ) = planeNormal[1];
   basisTrafo( 2, 2 ) = planeNormal[2];

   basisTrafo.invert();

   walberla::math::Vector3< real_t > centreOldBasis( centre[0], centre[1], centre[2] );
   walberla::math::Vector3< real_t > v1OldBasis( v1[0], v1[1], v1[2] );
   walberla::math::Vector3< real_t > v2OldBasis( v2[0], v2[1], v2[2] );
   walberla::math::Vector3< real_t > v3OldBasis( v3[0], v3[1], v3[2] );

   auto centrePlaneBasis = basisTrafo * centreOldBasis;
   auto v1PlaneBasis     = basisTrafo * v1OldBasis;
   auto v2PlaneBasis     = basisTrafo * v2OldBasis;
   auto v3PlaneBasis     = basisTrafo * v3OldBasis;

   return circleTriangleIntersection( Point2D( centrePlaneBasis[0], centrePlaneBasis[1] ),
                                      intersectionRadius,
                                      Point2D( v1PlaneBasis[0], v1PlaneBasis[1] ),
                                      Point2D( v2PlaneBasis[0], v2PlaneBasis[1] ),
                                      Point2D( v3PlaneBasis[0], v3PlaneBasis[1] ) );
}

/// Returns true if the passed point is located in (or on) the passed tetrahedron.
/// Optimized version if the inward normals are already pre-computed.
inline bool isPointInTetrahedron( const Point3D& pointOfInterest,
                                  const Point3D& tetVertex0,
                                  const Point3D& tetVertex1,
                                  const Point3D& tetVertex2,
                                  const Point3D& tetVertex3,
                                  const Point3D& faceInwardNormalOpposingVertex0,
                                  const Point3D& faceInwardNormalOpposingVertex1,
                                  const Point3D& faceInwardNormalOpposingVertex2,
                                  const Point3D& faceInwardNormalOpposingVertex3 )
{
   const auto distFace0 = distanceToPlane( pointOfInterest, tetVertex0, faceInwardNormalOpposingVertex3 );
   const auto distFace1 = distanceToPlane( pointOfInterest, tetVertex1, faceInwardNormalOpposingVertex2 );
   const auto distFace2 = distanceToPlane( pointOfInterest, tetVertex2, faceInwardNormalOpposingVertex1 );
   const auto distFace3 = distanceToPlane( pointOfInterest, tetVertex3, faceInwardNormalOpposingVertex0 );

   return distFace0 >= 0 && distFace1 >= 0 && distFace2 >= 0 && distFace3 >= 0;
}

/// Returns true if the passed point is located in (or on) the passed tetrahedron.
inline bool isPointInTetrahedron( const Point3D& pointOfInterest,
                                  const Point3D& tetVertex0,
                                  const Point3D& tetVertex1,
                                  const Point3D& tetVertex2,
                                  const Point3D& tetVertex3 )
{
   const auto faceInwardNormalOpposingVertex0 = tetrahedronInwardNormal( tetVertex1, tetVertex2, tetVertex3, tetVertex0 );
   const auto faceInwardNormalOpposingVertex1 = tetrahedronInwardNormal( tetVertex0, tetVertex2, tetVertex3, tetVertex1 );
   const auto faceInwardNormalOpposingVertex2 = tetrahedronInwardNormal( tetVertex0, tetVertex1, tetVertex3, tetVertex2 );
   const auto faceInwardNormalOpposingVertex3 = tetrahedronInwardNormal( tetVertex0, tetVertex1, tetVertex2, tetVertex3 );

   return isPointInTetrahedron( pointOfInterest,
                                tetVertex0,
                                tetVertex1,
                                tetVertex2,
                                tetVertex3,
                                faceInwardNormalOpposingVertex0,
                                faceInwardNormalOpposingVertex1,
                                faceInwardNormalOpposingVertex2,
                                faceInwardNormalOpposingVertex3 );
}

/// Returns true if the passed sphere is completely located in the passed tetrahedron.
inline bool isSphereCompletelyInTetrahedron( const Point3D& sphereCenter,
                                             const real_t&  sphereRadius,
                                             const Point3D& tetVertex0,
                                             const Point3D& tetVertex1,
                                             const Point3D& tetVertex2,
                                             const Point3D& tetVertex3 )
{
   const auto normalFace0 = tetrahedronInwardNormal( tetVertex0, tetVertex1, tetVertex2, tetVertex3 );
   const auto normalFace1 = tetrahedronInwardNormal( tetVertex0, tetVertex1, tetVertex3, tetVertex2 );
   const auto normalFace2 = tetrahedronInwardNormal( tetVertex0, tetVertex2, tetVertex3, tetVertex1 );
   const auto normalFace3 = tetrahedronInwardNormal( tetVertex1, tetVertex2, tetVertex3, tetVertex0 );

   const auto distFace0 = distanceToPlane( sphereCenter, tetVertex0, normalFace0 );
   const auto distFace1 = distanceToPlane( sphereCenter, tetVertex0, normalFace1 );
   const auto distFace2 = distanceToPlane( sphereCenter, tetVertex0, normalFace2 );
   const auto distFace3 = distanceToPlane( sphereCenter, tetVertex1, normalFace3 );

   return distFace0 >= sphereRadius && distFace1 >= sphereRadius && distFace2 >= sphereRadius && distFace3 >= sphereRadius;
}

inline bool sphereTetrahedronIntersection( const Point3D& sphereCenter,
                                           const real_t&  sphereRadius,
                                           const Point3D& tetVertex0,
                                           const Point3D& tetVertex1,
                                           const Point3D& tetVertex2,
                                           const Point3D& tetVertex3 )
{
   const auto pointInTet = isPointInTetrahedron( sphereCenter, tetVertex0, tetVertex1, tetVertex2, tetVertex3 );
   if ( pointInTet )
      return true;

   const auto sphereInTet =
       isSphereCompletelyInTetrahedron( sphereCenter, sphereRadius, tetVertex0, tetVertex1, tetVertex2, tetVertex3 );
   if ( sphereInTet )
      return true;

   const auto sphereTriangleIntersection0 =
       sphereTriangleIntersection( sphereCenter, sphereRadius, tetVertex0, tetVertex1, tetVertex2 );
   if ( sphereTriangleIntersection0 )
      return true;

   const auto sphereTriangleIntersection1 =
       sphereTriangleIntersection( sphereCenter, sphereRadius, tetVertex0, tetVertex1, tetVertex3 );
   if ( sphereTriangleIntersection1 )
      return true;

   const auto sphereTriangleIntersection2 =
       sphereTriangleIntersection( sphereCenter, sphereRadius, tetVertex0, tetVertex2, tetVertex3 );
   if ( sphereTriangleIntersection2 )
      return true;

   const auto sphereTriangleIntersection3 =
       sphereTriangleIntersection( sphereCenter, sphereRadius, tetVertex1, tetVertex2, tetVertex3 );
   if ( sphereTriangleIntersection3 )
      return true;

   return false;
}

/// Returns a baricentric tolerance determining how close a point is to a face (triangle)
/// The baricentric tolerance measures the maximum deviation of the baricentric coordinates and their sum
/// from the interval [0,1] and 1 respectively. This criterion is a bit dodgy but it works in practice.
/// \note Returns 0 if the point lies inside the face
inline real_t pointFaceBaricentricTolerance2D( const Point3D& p, const Point3D& A, const Point3D& B, const Point3D& C )
{
   Point3D AB = B - A;
   Point3D BC = C - B;
   Point3D CA = A - C;

   Point3D normal = AB.cross( C - A );
   normal.normalize();

   real_t dist_                    = normal.dot( A );
   real_t doubleSignedAreaTriangle = AB.cross( C - A ).dot( normal );

   // signed area of the subtriangles ABP, BCP, CAP
   real_t doubleSignedAreaSubTriangle1 = AB.cross( p - A ).dot( normal );
   real_t doubleSignedAreaSubTriangle2 = BC.cross( p - B ).dot( normal );
   real_t doubleSignedAreaSubTriangle3 = CA.cross( p - C ).dot( normal );

   // barycentric coordinates
   real_t c1 = doubleSignedAreaSubTriangle2 / doubleSignedAreaTriangle;
   real_t c2 = doubleSignedAreaSubTriangle3 / doubleSignedAreaTriangle;
   real_t c3 = doubleSignedAreaSubTriangle1 / doubleSignedAreaTriangle;

   // check if all barycentric coordinates are between 0 and 1, their sum should be 1 by construction
   if ( c1 < 0 || c2 < 0 || c3 < 0 )
   {
      return std::max( { std::abs( std::min( { c1, c2, c3 } ) ),
                         std::abs( std::max( { real_c( 0 ), c1 - real_c( 1 ), c2 - real_c( 1 ), c3 - real_c( 1 ) } ) ),
                         std::abs( c1 + c2 + c3 - real_c( 1.0 ) ) } );
   }

   return real_c( 0 );
}

/// Returns a baricentric tolerance determining how close a point is to a cell (tetrahedron)
/// The baricentric tolerance measures the maximum deviation of the baricentric coordinates and their sum
/// from the interval [0,1] and 1 respectively. This criterion is a bit dodgy but it works in practice.
/// \note Returns 0 if the point lies inside the cell
inline real_t
    pointCellBaricentricTolerance3D( const Point3D& p, const Point3D& A, const Point3D& B, const Point3D& C, const Point3D& D )
{
   Point3D AP = p - A;
   Point3D BP = p - B;
   Point3D AB = B - A;
   Point3D AC = C - A;
   Point3D AD = D - A;
   Point3D BC = C - B;
   Point3D BD = D - B;

   // volumes of the four subtetrahedra
   real_t Va = BP.dot( BD.cross( BC ) );
   real_t Vb = AP.dot( AC.cross( AD ) );
   real_t Vc = AP.dot( AD.cross( AB ) );
   real_t Vd = AP.dot( AB.cross( AC ) );
   real_t V  = AB.dot( AC.cross( AD ) );

   // barycentric coordinates
   real_t c1 = Va / V;
   real_t c2 = Vb / V;
   real_t c3 = Vc / V;
   real_t c4 = Vd / V;

   if ( c1 < 0 || c2 < 0 || c3 < 0 || c4 < 0 )
   {
      return std::max(
          { std::abs( std::min( { c1, c2, c3, c4 } ) ),
            std::abs( std::max( { real_c( 0 ), c1 - real_c( 1 ), c2 - real_c( 1 ), c3 - real_c( 1 ), c4 - real_c( 1 ) } ) ),
            std::abs( c1 + c2 + c3 + c4 - real_c( 1.0 ) ) } );
   }

   return real_c( 0 );
}

} // namespace hyteg
