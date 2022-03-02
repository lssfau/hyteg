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

/// \brief Quick reference implementation of const and linear DG elements with Lagrange basis functions.
class DGBasisLinearLagrange_Example : public DGBasisInfo
{
 public:
   virtual uint_t minPolynomialDegree() const { return 0; }

   virtual uint_t maxPolynomialDegree() const { return 1; }

   virtual uint_t numDoFsPerElement( uint_t degree ) const
   {
      WALBERLA_CHECK_LESS_EQUAL( degree, 1, "Only const or linear elements supported." );
      const std::map< uint_t, uint_t > ndofs = { { 0, 1 }, { 1, 3 } };
      return ndofs.at( degree );
   }

   virtual uint_t quadratureDegreeForLinearFunctional() const { return 4; }

   virtual void
       evaluate( uint_t degree, const Eigen::Matrix< real_t, 2, 1 >& pos, const std::vector< real_t >& dofs, real_t& value ) const
   {
      WALBERLA_CHECK_EQUAL( dofs.size(), numDoFsPerElement( degree ), "Number of DoFs does not match degree." );

      switch ( degree )
      {
      case 0: {
         value = 1;
         break;
      }

      case 1: {
         const auto dof_0 = dofs[0];
         const auto dof_1 = dofs[1];
         const auto dof_2 = dofs[2];
#if 0
      // affine space

      const auto x_affine_eval_0 = pos( 0 );
      const auto x_affine_eval_1 = pos( 1 );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_2_1 + tmp_0;
      real_t tmp_2  = -p_affine_0_0;
      real_t tmp_3  = p_affine_1_0 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_0 ) * ( p_affine_2_0 + tmp_2 ) );
      real_t tmp_5  = tmp_4 * ( tmp_2 + x_affine_eval_0 );
      real_t tmp_6  = tmp_1 * tmp_5;
      real_t tmp_7  = tmp_4 * ( tmp_0 + x_affine_eval_1 );
      real_t tmp_8  = tmp_7 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_9  = tmp_3 * tmp_7;
      real_t tmp_10 = tmp_5 * ( p_affine_0_1 - p_affine_1_1 );
      real_t a_0_0  = dof_0 * ( -tmp_10 - tmp_6 - tmp_8 - tmp_9 + 1 ) + dof_1 * ( tmp_6 + tmp_8 ) + dof_2 * ( tmp_10 + tmp_9 );
#endif
         // ref space

         const auto x_ref_0 = pos( 0 );
         const auto x_ref_1 = pos( 1 );

         real_t a_0_0 = dof_0 * ( -x_ref_0 - x_ref_1 + 1 ) + dof_1 * x_ref_0 + dof_2 * x_ref_1;

         value = a_0_0;
      }
      }
   };

   virtual void integrateBasisFunction( uint_t                                                degree,
                                        const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >& coords,
                                        const std::function< real_t( const Point3D& ) >&      f,
                                        std::vector< real_t >&                                values )
   {
      values.resize( numDoFsPerElement( degree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      callback_Scalar_Variable_Coefficient_2D_k = f;

      switch ( degree )
      {
      case 0: {
         WALBERLA_ABORT( "Not implemented." );
         break;
      }
      case 1: {
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

         values[0] = a_0_0;
         values[1] = a_1_1;
         values[2] = a_2_2;
         break;
      }
      }
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