/*
 * Copyright (c) 2023 Daniel Bauer.
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

#include "N1E1Form_mass_LU_fp32.hpp"

namespace hyteg {
namespace n1e1 {

void N1E1Form_mass_LU_fp32::integrateAll( const std::array< PointND< float, 3 >, 4 >& coords,
                                          const std::array< walberla::int16_t, 6 >&   edgeDirections,
                                          Matrix< float, 6, 6 >&                      elMat ) const
{
   elMat.setAll( 0.0f );

   // F maps from reference tet K' to affine tet K
   // K = F(K') = Bx' + b
   // B is the Jacobian of the transformation
   Eigen::Matrix< float, 3, 3 > BT;
   BT.row( 0 ) = coords[1].vector_ - coords[0].vector_;
   BT.row( 1 ) = coords[2].vector_ - coords[0].vector_;
   BT.row( 2 ) = coords[3].vector_ - coords[0].vector_;

   const float                                      absDetB  = std::abs( BT.determinant() );
   Eigen::FullPivLU< Eigen::Matrix< float, 3, 3 > > luDecomp = BT.fullPivLu();
   luDecomp.setThreshold( 1e-20f );

   // φᵢ = phi[i] * (x, y, z, 1)ᵀ
   std::array< Eigen::Matrix< float, 3, 4 >, 6 > phi;
   // clang-format off
      //         x   y   z   1
      phi[0] <<  0,  0,  0,  0, // 0
                 0,  0, -1,  0, // -z
                 0,  1,  0,  0; // y
      phi[1] <<  0,  0, -1,  0, // -z
                 0,  0,  0,  0, // 0
                 1,  0,  0,  0; // x
      phi[2] <<  0, -1,  0,  0, // -y
                 1,  0,  0,  0, // x
                 0,  0,  0,  0; // 0
      phi[3] <<  0,  0,  1,  0, // z
                 0,  0,  1,  0, // z
                -1, -1,  0,  1; // -x-y+1
      phi[4] <<  0,  1,  0,  0, // y
                -1,  0, -1,  1, // -x-z+1
                 0,  1,  0,  0; // y
      phi[5] <<  0, -1, -1,  1, // -y-z+1
                 1,  0,  0,  0, // x
                 1,  0,  0,  0; // x
   // clang-format on

   // B⁻ᵀ φᵢ = BInvTphi[i] * (x, y, z, 1)ᵀ
   std::array< Eigen::Matrix< float, 3, 4 >, 6 > BInvTphi;
   for ( uint_t i = 0; i < 6; i++ )
   {
      phi[i] *= edgeDirections[i];
      BInvTphi[i] = luDecomp.solve( phi[i] );
   }

   // integral = ∫_K' (x, y, z, 1) (x, y, z, 1)ᵀ
   Eigen::Matrix< float, 4, 4 > integral;
   // clang-format off
      integral << 1.0f/ 60.0f, 1.0f/120.0f, 1.0f/120.0f, 1.0f/ 24.0f,
                  1.0f/120.0f, 1.0f/ 60.0f, 1.0f/120.0f, 1.0f/ 24.0f,
                  1.0f/120.0f, 1.0f/120.0f, 1.0f/ 60.0f, 1.0f/ 24.0f,
                  1.0f/ 24.0f, 1.0f/ 24.0f, 1.0f/ 24.0f, 1.0f/  6.0f;
   // clang-format on

   for ( uint_t i = 0; i < 6; i++ )
   {
      for ( uint_t j = i; j < 6; j++ )
      {
         const float val = absDetB * ( BInvTphi[i].array() * ( BInvTphi[j] * integral ).array() ).sum();
         elMat( walberla::numeric_cast< Eigen::Index >( i ), walberla::numeric_cast< Eigen::Index >( j ) ) = val;
         elMat( walberla::numeric_cast< Eigen::Index >( j ), walberla::numeric_cast< Eigen::Index >( i ) ) = val;
      }
   }
}

} // namespace n1e1
} // namespace hyteg
