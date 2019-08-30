
#pragma once

#include "core/DataTypes.h"
#include "core/math/Matrix3.h"

#include "tinyhhg_core/types/pointnd.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;

/// \brief Returns the distance to a half plane in 3D.
/// If the returned distance is positive, the point of interest is located on the side of the plane the normal point to.
inline real_t distanceToPlane( const Point3D& pointOfInterest, const Point3D& pointOnPlane, const Point3D& planeNormal )
{
   return ( planeNormal / planeNormal.norm() ).dot( pointOfInterest - pointOnPlane );
}

/// \brief Returns the normal of the passed triangle in the direction of the opposite (4th) vertex of the tet.
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

/// \brief Returns true if the passed circle and triangle intersect.
inline bool circleTriangleIntersection( const Point2D& centre,
                                        const real_t&  radius,
                                        const Point2D& v1,
                                        const Point2D& v2,
                                        const Point2D& v3 )
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
   // NOTE: This works for clockwise ordered vertices!
   //
   const auto area = 0.5 * ( -v2y * v3x + v1y * ( -v2x + v3x ) + v1x * ( v2y - v3y ) + v2x * v3y );
   const auto s    = 1 / ( 2 * area ) * ( v1y * v3x - v1x * v3y + ( v3y - v1y ) * centrex + ( v1x - v3x ) * centrey );
   const auto t    = 1 / ( 2 * area ) * ( v1x * v2y - v1y * v2x + ( v1y - v2y ) * centrex + ( v2x - v1x ) * centrey );
   if ( s > 0 && t > 0 && 1 - s - t > 0 )
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

inline bool sphereTriangleIntersection( const Point3D& centre,
                                        const real_t&  radius,
                                        const Point3D& v1,
                                        const Point3D& v2,
                                        const Point3D& v3 )
{
   // If the distance of the sphere to the triangle plane is greater than the sphere's radius there is no intersection.
   // Otherwise, check intersection of the triangle with the circle that is common with the plane (where the plane cuts through the sphere).
   auto planeNormal       = tetrahedronInwardNormal( v1, v2, v3, centre );
   auto centreDistToPlane = distanceToPlane( centre, v1, planeNormal );
   if ( centreDistToPlane > radius )
      return false;

   auto intersectionRadius = std::sqrt( radius * radius - centreDistToPlane * centreDistToPlane );

   // We need to rotate our coordinate system -> normal component in z direction.
   auto planeTangent0 = v2 - v1;
   planeTangent0 /= planeTangent0.norm();
   auto planeTangent1 = crossProduct( planeNormal, planeTangent0 );

   walberla::math::Matrix3< real_t > basisTrafo;

   basisTrafo( 0, 0 ) = planeTangent0[0];
   basisTrafo( 0, 1 ) = planeTangent0[1];
   basisTrafo( 0, 2 ) = planeTangent0[2];

   basisTrafo( 1, 0 ) = planeTangent1[0];
   basisTrafo( 1, 1 ) = planeTangent1[1];
   basisTrafo( 1, 2 ) = planeTangent1[2];

   basisTrafo( 2, 0 ) = planeNormal[0];
   basisTrafo( 2, 1 ) = planeNormal[1];
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

   return circleTriangleIntersection( Point2D( {centrePlaneBasis[0], centrePlaneBasis[1]} ),
                                      intersectionRadius,
                                      Point2D( {v1PlaneBasis[0], v1PlaneBasis[1]} ),
                                      Point2D( {v2PlaneBasis[0], v2PlaneBasis[1]} ),
                                      Point2D( {v3PlaneBasis[0], v3PlaneBasis[1]} ) );
}

/// Returns true if the passed point is located in (or on) the passed tetrahedron.
inline bool isPointInTetrahedron( const Point3D& pointOfInterest,
                                  const Point3D& tetVertex0,
                                  const Point3D& tetVertex1,
                                  const Point3D& tetVertex2,
                                  const Point3D& tetVertex3 )
{
   const auto normalFace0 = tetrahedronInwardNormal( tetVertex0, tetVertex1, tetVertex2, tetVertex3 );
   const auto normalFace1 = tetrahedronInwardNormal( tetVertex0, tetVertex1, tetVertex3, tetVertex2 );
   const auto normalFace2 = tetrahedronInwardNormal( tetVertex0, tetVertex2, tetVertex3, tetVertex1 );
   const auto normalFace3 = tetrahedronInwardNormal( tetVertex1, tetVertex2, tetVertex3, tetVertex0 );

   const auto distFace0 = distanceToPlane( pointOfInterest, tetVertex0, normalFace0 );
   const auto distFace1 = distanceToPlane( pointOfInterest, tetVertex0, normalFace1 );
   const auto distFace2 = distanceToPlane( pointOfInterest, tetVertex0, normalFace2 );
   const auto distFace3 = distanceToPlane( pointOfInterest, tetVertex1, normalFace3 );

   return distFace0 >= -real_c( 1e-08 ) && distFace1 >= -real_c( 1e-08 ) && distFace2 >= -real_c( 1e-08 ) &&
          distFace3 >= -real_c( 1e-08 );
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
   return isSphereCompletelyInTetrahedron( sphereCenter, sphereRadius, tetVertex0, tetVertex1, tetVertex2, tetVertex3 ) ||
          sphereTriangleIntersection( sphereCenter, sphereRadius, tetVertex0, tetVertex1, tetVertex2 ) ||
          sphereTriangleIntersection( sphereCenter, sphereRadius, tetVertex0, tetVertex1, tetVertex3 ) ||
          sphereTriangleIntersection( sphereCenter, sphereRadius, tetVertex0, tetVertex2, tetVertex3 ) ||
          sphereTriangleIntersection( sphereCenter, sphereRadius, tetVertex1, tetVertex2, tetVertex3 );
}

} // namespace hyteg