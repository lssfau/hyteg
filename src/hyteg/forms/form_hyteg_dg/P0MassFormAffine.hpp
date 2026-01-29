/*
* Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr, Daniel Bauer, Fabian Böhm.
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
#include "hyteg/forms/DGForm.hpp"
#include "hyteg/forms/form_hyteg_dg/DGFormVolume.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

class P0MassFormAffine : public DGFormVolume
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
   {
      WALBERLA_ASSERT( trialDegree == 0 );
      WALBERLA_ASSERT( testDegree == 0 );

      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      real_t p_affine_0_0       = coords[0][0];
      real_t p_affine_0_1       = coords[0][1];
      real_t p_affine_1_0       = coords[1][0];
      real_t p_affine_1_1       = coords[1][1];
      real_t p_affine_2_0       = coords[2][0];
      real_t p_affine_2_1       = coords[2][1];
      real_t jac_affine_0_0     = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1     = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_1_0     = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1     = -p_affine_0_1 + p_affine_2_1;
      real_t abs_det_jac_affine = std::abs( jac_affine_0_0 * jac_affine_1_1 - jac_affine_0_1 * jac_affine_1_0 );
      real_t a_0_0              = 0.5 * abs_det_jac_affine;
      elMat( 0, 0 )             = a_0_0;
   }

   virtual void integrateVolume3D( const std::vector< Point3D >& coords,
                                   const DGBasisInfo&            trialBasis,
                                   const DGBasisInfo&            testBasis,
                                   int                           trialDegree,
                                   int                           testDegree,
                                   MatrixXr&                     elMat ) const override
   {
      WALBERLA_ASSERT( trialDegree == 0 );
      WALBERLA_ASSERT( testDegree == 0 );

      elMat.resize( testBasis.numDoFsPerElement( 3, testDegree ), trialBasis.numDoFsPerElement( 3, trialDegree ) );

      real_t p_affine_0_0   = coords[0][0];
      real_t p_affine_0_1   = coords[0][1];
      real_t p_affine_0_2   = coords[0][2];
      real_t p_affine_1_0   = coords[1][0];
      real_t p_affine_1_1   = coords[1][1];
      real_t p_affine_1_2   = coords[1][2];
      real_t p_affine_2_0   = coords[2][0];
      real_t p_affine_2_1   = coords[2][1];
      real_t p_affine_2_2   = coords[2][2];
      real_t p_affine_3_0   = coords[3][0];
      real_t p_affine_3_1   = coords[3][1];
      real_t p_affine_3_2   = coords[3][2];
      real_t jac_affine_0_0 = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1 = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_0_2 = -p_affine_0_0 + p_affine_3_0;
      real_t jac_affine_1_0 = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1 = -p_affine_0_1 + p_affine_2_1;
      real_t jac_affine_1_2 = -p_affine_0_1 + p_affine_3_1;
      real_t jac_affine_2_0 = -p_affine_0_2 + p_affine_1_2;
      real_t jac_affine_2_1 = -p_affine_0_2 + p_affine_2_2;
      real_t jac_affine_2_2 = -p_affine_0_2 + p_affine_3_2;
      real_t abs_det_jac_affine =
          std::abs( jac_affine_0_0 * jac_affine_1_1 * jac_affine_2_2 - jac_affine_0_0 * jac_affine_1_2 * jac_affine_2_1 -
                    jac_affine_0_1 * jac_affine_1_0 * jac_affine_2_2 + jac_affine_0_1 * jac_affine_1_2 * jac_affine_2_0 +
                    jac_affine_0_2 * jac_affine_1_0 * jac_affine_2_1 - jac_affine_0_2 * jac_affine_1_1 * jac_affine_2_0 );
      real_t a_0_0  = ( 1.0 / 6.0 ) * abs_det_jac_affine;
      elMat( 0, 0 ) = a_0_0;
   }
};

} // namespace dg
} // namespace hyteg
