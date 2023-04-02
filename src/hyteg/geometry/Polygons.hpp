/*
 * Copyright (c) 2017-2021 Nils Kohl.
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

#include <cmath>
#include <vector>

#include "core/math/Constants.h"
#include "core/Abort.h"

#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::math::pi;

/// \brief Given a convex polygon, an arbitrary "center" point in this polygon, and another point P in this polygon,
///        this function returns the fraction of the distance from the point P to the point P' on the boundary that intersects the
///        ray from the center through the point in [0, 1].
///
///        This can be interpreted as the radius to that point if the polygon is projected onto a circle around the chosen center.
///
/// \param P
/// \param center
/// \param polygonVertices
static void
    fractionalRadiusToPolygonBoundary( const Point3D& P, const Point3D& center, const std::vector< Point3D >& polygonVertices, real_t & r, real_t & angle )
{
   const real_t eps = 1e-12;

   auto   C        = center;
   real_t gammaSum = 0;
   auto   S        = polygonVertices[0];

   auto normal = crossProduct( polygonVertices[0], polygonVertices[1] );

   for ( uint_t i = 0; i < polygonVertices.size(); i++ )
   {
      auto A = polygonVertices[i];
      auto B = ( i == polygonVertices.size() - 1 ) ? polygonVertices[0] : polygonVertices[i + 1];

      auto a = ( B - C ).norm();
      auto b = ( A - C ).norm();
      auto c = ( A - B ).norm();

      auto gamma = std::acos( ( a * a + b * b - c * c ) / (2 * a * b) );
      gammaSum += gamma;

      auto alpha = std::acos( ( c * c + b * b - a * a ) / (2 * c * b) );

      //  s = length(cross_product(a,b))
      //  c = dot_product(a,b)
      //  angle = atan2(s, c)

      auto angleP_atC_full = std::atan2( crossProduct( ( S - C ), ( P - C ) ).norm(), ( S - C ).dot( ( P - C ) ) );
      auto angleP_atC      = std::atan2( crossProduct( ( A - C ), ( P - C ) ).norm(), ( A - C ).dot( ( P - C ) ) );

      if ( normal.dot( crossProduct( ( S - C ), ( P - C ) ) ) < 0 )
      {
         angleP_atC_full = - angleP_atC_full + 2 * pi;
      }

      if ( normal.dot( crossProduct( ( A - C ), ( P - C ) ) ) < 0 )
      {
         angleP_atC = - angleP_atC + 2 * pi;
      }

      if ( angleP_atC_full < 0 )
      {
         angleP_atC_full += 2 * pi;
      }

      if ( angleP_atC < 0 )
      {
         angleP_atC += 2 * pi;
      }

      if ( angleP_atC_full > gammaSum + eps )
      {
         continue;
      }

      auto angleP_Ptilde = pi - angleP_atC - alpha;
      auto Ptilde_length = std::sin( alpha ) * ( b / std::sin( angleP_Ptilde ) );

      r = std::abs( ( P - C ).norm() / Ptilde_length );
      angle = angleP_atC_full;

      WALBERLA_ASSERT_GREATER( r, -eps );
      WALBERLA_ASSERT_LESS( r, 1.0 + eps );

      return;
   }

   std::stringstream ss;
   for ( auto p : polygonVertices )
   {
      ss << p << "\n";
   }

   WALBERLA_ABORT( "Could not find point " << P << " in polygon with center " << center << ":\n" << ss.str() )
}

} // namespace hyteg