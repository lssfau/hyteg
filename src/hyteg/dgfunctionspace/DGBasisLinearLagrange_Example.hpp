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
#include "hyteg/types/matrix.hpp"
#include "hyteg/types/pointnd.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

/// \brief Quick reference implementation of linear DG elements with Lagrange basis functions.
class DGBasisLinearLagrange_Example : public DGBasisInfo
{
 public:
   virtual int minPolynomialDegree() const { return 1; }

   virtual int maxPolynomialDegree() const { return 1; }

   virtual int numDoFsPerElement( int degree ) const
   {
      WALBERLA_CHECK_EQUAL( degree, 1, "Only linear elements supported." );
      return 3;
   }

   virtual int quadratureDegreeForLinearFunctional() const { return 4; }

   virtual void
       evaluate( int degree, const Eigen::Matrix< real_t, 2, 1 >& pos, const std::vector< real_t >& dofs, real_t& value ) const
   {
      WALBERLA_CHECK_EQUAL( degree, 1, "Only degree 1 supported." );
      WALBERLA_CHECK_EQUAL(
          dofs.size(), numDoFsPerElement( degree ), "Linear Lagrange basis requires exactly 3 DoFs per element." );

      const auto x_ref_0 = pos( 0 );
      const auto x_ref_1 = pos( 1 );

      const auto dof_0 = dofs[0];
      const auto dof_1 = dofs[1];
      const auto dof_2 = dofs[2];

      real_t a_0_0 = dof_0 * ( -x_ref_0 - x_ref_1 + 1 ) + dof_1 * x_ref_0 + dof_2 * x_ref_1;

      value = a_0_0;
   };

   virtual void integrateBasisFunction( int                                                   degree,
                                        const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >& coords,
                                        const std::function< real_t( const Point3D& ) >&      f,
                                        int                                                   basisFunctionIdx,
                                        real_t&                                               value )
   {
      WALBERLA_CHECK_EQUAL( degree, 1, "Only degree 1 supported." );
      WALBERLA_CHECK_GREATER_EQUAL( basisFunctionIdx, 0, "Basis function index should be larger equal 0." );
      WALBERLA_CHECK_LESS( basisFunctionIdx,
                           numDoFsPerElement( degree ),
                           "Basis function index too large: linear Lagrange basis has exactly 3 basis functions per element." );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      callback_Scalar_Variable_Coefficient_2D_k = f;

      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_k(
          0.091576213509770743 * p_affine_0_0 + 0.091576213509770743 * p_affine_1_0 + 0.81684757298045851 * p_affine_2_0,
          0.091576213509770743 * p_affine_0_1 + 0.091576213509770743 * p_affine_1_1 + 0.81684757298045851 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      Scalar_Variable_Coefficient_2D_k(
          0.44594849091596489 * p_affine_0_0 + 0.44594849091596489 * p_affine_1_0 + 0.10810301816807022 * p_affine_2_0,
          0.44594849091596489 * p_affine_0_1 + 0.44594849091596489 * p_affine_1_1 + 0.10810301816807022 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_k_out0_id1 );
      Scalar_Variable_Coefficient_2D_k(
          0.091576213509770743 * p_affine_0_0 + 0.81684757298045851 * p_affine_1_0 + 0.091576213509770743 * p_affine_2_0,
          0.091576213509770743 * p_affine_0_1 + 0.81684757298045851 * p_affine_1_1 + 0.091576213509770743 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_k_out0_id2 );
      Scalar_Variable_Coefficient_2D_k(
          0.44594849091596489 * p_affine_0_0 + 0.10810301816807022 * p_affine_1_0 + 0.44594849091596489 * p_affine_2_0,
          0.44594849091596489 * p_affine_0_1 + 0.10810301816807022 * p_affine_1_1 + 0.44594849091596489 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_k_out0_id3 );
      Scalar_Variable_Coefficient_2D_k(
          0.81684757298045851 * p_affine_0_0 + 0.091576213509770743 * p_affine_1_0 + 0.091576213509770743 * p_affine_2_0,
          0.81684757298045851 * p_affine_0_1 + 0.091576213509770743 * p_affine_1_1 + 0.091576213509770743 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_k_out0_id4 );
      Scalar_Variable_Coefficient_2D_k(
          0.10810301816807022 * p_affine_0_0 + 0.44594849091596489 * p_affine_1_0 + 0.44594849091596489 * p_affine_2_0,
          0.10810301816807022 * p_affine_0_1 + 0.44594849091596489 * p_affine_1_1 + 0.44594849091596489 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_k_out0_id5 );
      real_t tmp_0 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_1 = 0.054975871827660928 * Scalar_Variable_Coefficient_2D_k_out0_id0 * tmp_0;
      real_t tmp_2 = 0.11169079483900572 * Scalar_Variable_Coefficient_2D_k_out0_id1 * tmp_0;
      real_t tmp_3 = 0.054975871827660928 * Scalar_Variable_Coefficient_2D_k_out0_id2 * tmp_0;
      real_t tmp_4 = 0.11169079483900572 * Scalar_Variable_Coefficient_2D_k_out0_id3 * tmp_0;
      real_t tmp_5 = 0.054975871827660928 * Scalar_Variable_Coefficient_2D_k_out0_id4 * tmp_0;
      real_t tmp_6 = 0.11169079483900572 * Scalar_Variable_Coefficient_2D_k_out0_id5 * tmp_0;
      real_t a_0_0 = 0.091576213509770743 * tmp_1 + 0.44594849091596489 * tmp_2 + 0.091576213509770743 * tmp_3 +
                     0.44594849091596489 * tmp_4 + 0.81684757298045851 * tmp_5 + 0.10810301816807022 * tmp_6;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0.091576213509770743 * tmp_1 + 0.44594849091596489 * tmp_2 + 0.81684757298045851 * tmp_3 +
                     0.10810301816807022 * tmp_4 + 0.091576213509770743 * tmp_5 + 0.44594849091596489 * tmp_6;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0.81684757298045851 * tmp_1 + 0.10810301816807022 * tmp_2 + 0.091576213509770743 * tmp_3 +
                     0.44594849091596489 * tmp_4 + 0.091576213509770743 * tmp_5 + 0.44594849091596489 * tmp_6;

      Eigen::Matrix< real_t, 3, 1 > result( { a_0_0, a_1_1, a_2_2 } );
      value = result( basisFunctionIdx );
   }

 private:
   void Scalar_Variable_Coefficient_2D_k( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_k( Point3D( { in_0, in_1, 0 } ) );
   }

   std::function< real_t( const Point3D& ) > callback_Scalar_Variable_Coefficient_2D_k;
};

} // namespace dg
} // namespace hyteg