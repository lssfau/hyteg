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
#include "hyteg/types/matrix.hpp"
#include "hyteg/types/pointnd.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

class DGMassForm_Example : public DGForm
{
 public:
   void integrateVolume( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coords,
                         const DGBasisInfo&                                       trialBasis,
                         const DGBasisInfo&                                       testBasis,
                         int                                                      trialDegree,
                         int                                                      testDegree,
                         Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( testDegree ), trialBasis.numDoFsPerElement( trialDegree ) );

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

   virtual void integrateFacetInner( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElement,
                                     const std::array< Eigen::Matrix< real_t, 2, 1 >, 2 >&    coordsFacet,
                                     const DGBasisInfo&                                       trialBasis,
                                     const DGBasisInfo&                                       testBasis,
                                     int                                                      trialDegree,
                                     int                                                      testDegree,
                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( trialBasis );
      WALBERLA_UNUSED( testBasis );
      WALBERLA_UNUSED( trialDegree );
      WALBERLA_UNUSED( testDegree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }

   virtual void integrateFacetCoupling( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElementInner,
                                        const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElementOuter,
                                        const std::array< Eigen::Matrix< real_t, 2, 1 >, 2 >&    coordsFacet,
                                        const DGBasisInfo&                                       trialBasis,
                                        const DGBasisInfo&                                       testBasis,
                                        int                                                      trialDegree,
                                        int                                                      testDegree,
                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      WALBERLA_UNUSED( coordsElementInner );
      WALBERLA_UNUSED( coordsElementOuter );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( trialBasis );
      WALBERLA_UNUSED( testBasis );
      WALBERLA_UNUSED( trialDegree );
      WALBERLA_UNUSED( testDegree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   };
};

} // namespace dg
} // namespace hyteg