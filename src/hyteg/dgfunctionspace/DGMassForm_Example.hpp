/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/dgfunctionspace/DGBasisInfo.hpp"
#include "hyteg/dgfunctionspace/DGForm.hpp"
#include "hyteg/dgfunctionspace/DGFormVolume.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

class DGMassForm_Example : public DGFormVolume
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                     elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_1 = ( 1.0 / 12.0 ) * tmp_0;
      real_t tmp_2 = ( 1.0 / 24.0 ) * tmp_0;
      real_t a_0_0 = tmp_1;
      real_t a_0_1 = tmp_2;
      real_t a_0_2 = tmp_2;
      real_t a_1_0 = tmp_2;
      real_t a_1_1 = tmp_1;
      real_t a_1_2 = tmp_2;
      real_t a_2_0 = tmp_2;
      real_t a_2_1 = tmp_2;
      real_t a_2_2 = tmp_1;

      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;

      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;

      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
   }

   virtual void integrateVolume3D( const std::vector< Point3D >& coords,
                                   const DGBasisInfo&                                       trialBasis,
                                   const DGBasisInfo&                                       testBasis,
                                   int                                                      trialDegree,
                                   int                                                      testDegree,
                                   MatrixXr&                     elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 3, testDegree ), trialBasis.numDoFsPerElement( 3, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );
      const auto p_affine_0_2 = coords[0]( 2 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );
      const auto p_affine_1_2 = coords[1]( 2 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );
      const auto p_affine_2_2 = coords[2]( 2 );

      const auto p_affine_3_0 = coords[3]( 0 );
      const auto p_affine_3_1 = coords[3]( 1 );
      const auto p_affine_3_2 = coords[3]( 2 );

      real_t tmp_0  = p_affine_0_0 * p_affine_1_1;
      real_t tmp_1  = p_affine_0_0 * p_affine_1_2;
      real_t tmp_2  = p_affine_2_1 * p_affine_3_2;
      real_t tmp_3  = p_affine_0_1 * p_affine_1_0;
      real_t tmp_4  = p_affine_0_1 * p_affine_1_2;
      real_t tmp_5  = p_affine_2_2 * p_affine_3_0;
      real_t tmp_6  = p_affine_0_2 * p_affine_1_0;
      real_t tmp_7  = p_affine_0_2 * p_affine_1_1;
      real_t tmp_8  = p_affine_2_0 * p_affine_3_1;
      real_t tmp_9  = p_affine_2_2 * p_affine_3_1;
      real_t tmp_10 = p_affine_2_0 * p_affine_3_2;
      real_t tmp_11 = p_affine_2_1 * p_affine_3_0;
      real_t tmp_12 = std::abs( p_affine_0_0 * tmp_2 - p_affine_0_0 * tmp_9 - p_affine_0_1 * tmp_10 + p_affine_0_1 * tmp_5 -
                                p_affine_0_2 * tmp_11 + p_affine_0_2 * tmp_8 - p_affine_1_0 * tmp_2 + p_affine_1_0 * tmp_9 +
                                p_affine_1_1 * tmp_10 - p_affine_1_1 * tmp_5 + p_affine_1_2 * tmp_11 - p_affine_1_2 * tmp_8 +
                                p_affine_2_0 * tmp_4 - p_affine_2_0 * tmp_7 - p_affine_2_1 * tmp_1 + p_affine_2_1 * tmp_6 +
                                p_affine_2_2 * tmp_0 - p_affine_2_2 * tmp_3 - p_affine_3_0 * tmp_4 + p_affine_3_0 * tmp_7 +
                                p_affine_3_1 * tmp_1 - p_affine_3_1 * tmp_6 - p_affine_3_2 * tmp_0 + p_affine_3_2 * tmp_3 );
      real_t tmp_13 = ( 1.0 / 60.0 ) * tmp_12;
      real_t tmp_14 = ( 1.0 / 120.0 ) * tmp_12;
      real_t a_0_0  = tmp_13;
      real_t a_0_1  = tmp_14;
      real_t a_0_2  = tmp_14;
      real_t a_0_3  = tmp_14;
      real_t a_1_0  = tmp_14;
      real_t a_1_1  = tmp_13;
      real_t a_1_2  = tmp_14;
      real_t a_1_3  = tmp_14;
      real_t a_2_0  = tmp_14;
      real_t a_2_1  = tmp_14;
      real_t a_2_2  = tmp_13;
      real_t a_2_3  = tmp_14;
      real_t a_3_0  = tmp_14;
      real_t a_3_1  = tmp_14;
      real_t a_3_2  = tmp_14;
      real_t a_3_3  = tmp_13;

      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;

      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 1, 3 ) = a_1_3;

      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
      elMat( 2, 3 ) = a_2_3;

      elMat( 3, 0 ) = a_3_0;
      elMat( 3, 1 ) = a_3_1;
      elMat( 3, 2 ) = a_3_2;
      elMat( 3, 3 ) = a_3_3;
   }
};

} // namespace dg
} // namespace hyteg