/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "hyteg/forms/N1E1Form.hpp"

using walberla::real_c;
using walberla::uint_c;

namespace hyteg {
namespace n1e1 {

class N1E1Form_curl_curl : public N1E1Form
{
 public:
   /// Edges are in FEniCS ordering:
   ///
   /// 3
   /// |\`\.
   /// | 0 `\.
   /// |  \   1
   /// 3  2 _  `\.
   /// |  /  `-2 `\.
   /// | 4      `\_`\
   /// 0------5------1
   void integrateAll( const std::array< Point3D, 4 >& coords,
                      const std::array< int, 6 >&     edgeDirections,
                      Matrix6r&                       elMat ) const final
   {
      elMat.setZero();

      // F maps from reference tet K' to affine tet K
      // K = F(K') = Bx' + b
      // B is the Jacobian of the transformation
      hyteg::Matrix3r B;
      B.col( 0 ) = coords[1] - coords[0];
      B.col( 1 ) = coords[2] - coords[0];
      B.col( 2 ) = coords[3] - coords[0];

      const real_t absDetB = std::abs( B.determinant() );

      // for first order elements, the curl of the basis functions φ is constant
      // clang-format off
      std::array< Point3D, 6 > curlPhi = { Point3D{  2,  0,  0 } * real_c( edgeDirections[0] ),
                                           Point3D{  0, -2,  0 } * real_c( edgeDirections[1] ),
                                           Point3D{  0,  0,  2 } * real_c( edgeDirections[2] ),
                                           Point3D{ -2,  2,  0 } * real_c( edgeDirections[3] ),
                                           Point3D{  2,  0, -2 } * real_c( edgeDirections[4] ),
                                           Point3D{  0, -2,  2 } * real_c( edgeDirections[5] ) };
      // clang-format on

      for ( int i = 0; i < 6; i++ )
      {
         for ( int j = i; j < 6; j++ )
         {
            const real_t val = 1.0 / ( 6.0 * absDetB ) * ( B * curlPhi[uint_c( i )] ).dot( B * curlPhi[uint_c( j )] );
            elMat( i, j )    = val;
            elMat( j, i )    = val;
         }
      }
   };
};

} // namespace n1e1
} // namespace hyteg
