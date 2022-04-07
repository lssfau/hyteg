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

   virtual uint_t numDoFsPerElement( uint_t dim, uint_t degree ) const
   {
      if ( dim == 2 )
      {
         WALBERLA_CHECK_LESS_EQUAL( degree, 1, "Only const or linear elements supported." );
         const std::map< uint_t, uint_t > ndofs = { { 0, 1 }, { 1, 3 } };
         return ndofs.at( degree );
      }
      else if ( dim == 3 )
      {
         WALBERLA_CHECK_LESS_EQUAL( degree, 1, "Only const or linear elements supported." );
         const std::map< uint_t, uint_t > ndofs = { { 0, 1 }, { 1, 4 } };
         return ndofs.at( degree );
      }
      else
      {
         WALBERLA_ABORT( "Dimensionality invalid." );
      }
   }

   virtual uint_t quadratureDegreeForLinearFunctional() const { return 4; }

   virtual void
       evaluate( uint_t degree, const Eigen::Matrix< real_t, 2, 1 >& pos, const std::vector< real_t >& dofs, real_t& value ) const
   {
      WALBERLA_CHECK_EQUAL( dofs.size(), numDoFsPerElement( 2, degree ), "Number of DoFs does not match degree." );

      switch ( degree )
      {
      case 0: {
         value = dofs[0];
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

   virtual void
       evaluate( uint_t degree, const Eigen::Matrix< real_t, 3, 1 >& pos, const std::vector< real_t >& dofs, real_t& value ) const
   {
      WALBERLA_CHECK_EQUAL( dofs.size(), numDoFsPerElement( 3, degree ), "Number of DoFs does not match degree." );

      switch ( degree )
      {
      case 0: {
         value = dofs[0];
         break;
      }

      case 1: {
         const auto dof_0 = dofs[0];
         const auto dof_1 = dofs[1];
         const auto dof_2 = dofs[2];
         const auto dof_3 = dofs[3];

         // ref space

         const auto x_ref_0 = pos( 0 );
         const auto x_ref_1 = pos( 1 );
         const auto x_ref_2 = pos( 2 );

         real_t a_0_0 = dof_0 * ( -x_ref_0 - x_ref_1 - x_ref_2 + 1 ) + dof_1 * x_ref_0 + dof_2 * x_ref_1 + dof_3 * x_ref_2;

         value = a_0_0;
      }
      }
   };

   virtual void integrateBasisFunction( uint_t                                                degree,
                                        const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >& coords,
                                        const std::function< real_t( const Point3D& ) >&      f,
                                        std::vector< real_t >&                                values )
   {
      values.resize( numDoFsPerElement( 2, degree ) );

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

   virtual void integrateBasisFunction( uint_t                                                degree,
                                        const std::array< Eigen::Matrix< real_t, 3, 1 >, 4 >& coords,
                                        const std::function< real_t( const Point3D& ) >&      f,
                                        std::vector< real_t >&                                values )
   {
      values.resize( numDoFsPerElement( 3, degree ) );

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

      callback_Scalar_Variable_Coefficient_3D_k = f;

      switch ( degree )
      {
      case 0: {
         WALBERLA_ABORT( "Not implemented." );
         break;
      }
      case 1: {
         real_t Scalar_Variable_Coefficient_3D_k_out0_id0  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id1  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id2  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id3  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id4  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id5  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id6  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id7  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id8  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id9  = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id10 = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id11 = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id12 = 0;
         real_t Scalar_Variable_Coefficient_3D_k_out0_id13 = 0;
         Scalar_Variable_Coefficient_3D_k( 0.3108859192633005 * p_affine_0_0 + 0.31088591926330061 * p_affine_1_0 +
                                               0.31088591926330061 * p_affine_2_0 + 0.067342242210098213 * p_affine_3_0,
                                           0.3108859192633005 * p_affine_0_1 + 0.31088591926330061 * p_affine_1_1 +
                                               0.31088591926330061 * p_affine_2_1 + 0.067342242210098213 * p_affine_3_1,
                                           0.3108859192633005 * p_affine_0_2 + 0.31088591926330061 * p_affine_1_2 +
                                               0.31088591926330061 * p_affine_2_2 + 0.067342242210098213 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id0 );
         Scalar_Variable_Coefficient_3D_k( 0.092735250310891248 * p_affine_0_0 + 0.092735250310891248 * p_affine_1_0 +
                                               0.092735250310891248 * p_affine_2_0 + 0.72179424906732625 * p_affine_3_0,
                                           0.092735250310891248 * p_affine_0_1 + 0.092735250310891248 * p_affine_1_1 +
                                               0.092735250310891248 * p_affine_2_1 + 0.72179424906732625 * p_affine_3_1,
                                           0.092735250310891248 * p_affine_0_2 + 0.092735250310891248 * p_affine_1_2 +
                                               0.092735250310891248 * p_affine_2_2 + 0.72179424906732625 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id1 );
         Scalar_Variable_Coefficient_3D_k( 0.45449629587435036 * p_affine_0_0 + 0.045503704125649642 * p_affine_1_0 +
                                               0.045503704125649642 * p_affine_2_0 + 0.45449629587435036 * p_affine_3_0,
                                           0.45449629587435036 * p_affine_0_1 + 0.045503704125649642 * p_affine_1_1 +
                                               0.045503704125649642 * p_affine_2_1 + 0.45449629587435036 * p_affine_3_1,
                                           0.45449629587435036 * p_affine_0_2 + 0.045503704125649642 * p_affine_1_2 +
                                               0.045503704125649642 * p_affine_2_2 + 0.45449629587435036 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id2 );
         Scalar_Variable_Coefficient_3D_k( 0.045503704125649629 * p_affine_0_0 + 0.45449629587435036 * p_affine_1_0 +
                                               0.45449629587435036 * p_affine_2_0 + 0.045503704125649642 * p_affine_3_0,
                                           0.045503704125649629 * p_affine_0_1 + 0.45449629587435036 * p_affine_1_1 +
                                               0.45449629587435036 * p_affine_2_1 + 0.045503704125649642 * p_affine_3_1,
                                           0.045503704125649629 * p_affine_0_2 + 0.45449629587435036 * p_affine_1_2 +
                                               0.45449629587435036 * p_affine_2_2 + 0.045503704125649642 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id3 );
         Scalar_Variable_Coefficient_3D_k( 0.45449629587435036 * p_affine_0_0 + 0.045503704125649642 * p_affine_1_0 +
                                               0.45449629587435036 * p_affine_2_0 + 0.045503704125649642 * p_affine_3_0,
                                           0.45449629587435036 * p_affine_0_1 + 0.045503704125649642 * p_affine_1_1 +
                                               0.45449629587435036 * p_affine_2_1 + 0.045503704125649642 * p_affine_3_1,
                                           0.45449629587435036 * p_affine_0_2 + 0.045503704125649642 * p_affine_1_2 +
                                               0.45449629587435036 * p_affine_2_2 + 0.045503704125649642 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id4 );
         Scalar_Variable_Coefficient_3D_k( 0.45449629587435036 * p_affine_0_0 + 0.45449629587435036 * p_affine_1_0 +
                                               0.045503704125649642 * p_affine_2_0 + 0.045503704125649642 * p_affine_3_0,
                                           0.45449629587435036 * p_affine_0_1 + 0.45449629587435036 * p_affine_1_1 +
                                               0.045503704125649642 * p_affine_2_1 + 0.045503704125649642 * p_affine_3_1,
                                           0.45449629587435036 * p_affine_0_2 + 0.45449629587435036 * p_affine_1_2 +
                                               0.045503704125649642 * p_affine_2_2 + 0.045503704125649642 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id5 );
         Scalar_Variable_Coefficient_3D_k( 0.3108859192633005 * p_affine_0_0 + 0.31088591926330061 * p_affine_1_0 +
                                               0.067342242210098213 * p_affine_2_0 + 0.31088591926330061 * p_affine_3_0,
                                           0.3108859192633005 * p_affine_0_1 + 0.31088591926330061 * p_affine_1_1 +
                                               0.067342242210098213 * p_affine_2_1 + 0.31088591926330061 * p_affine_3_1,
                                           0.3108859192633005 * p_affine_0_2 + 0.31088591926330061 * p_affine_1_2 +
                                               0.067342242210098213 * p_affine_2_2 + 0.31088591926330061 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id6 );
         Scalar_Variable_Coefficient_3D_k( 0.092735250310891248 * p_affine_0_0 + 0.092735250310891248 * p_affine_1_0 +
                                               0.72179424906732625 * p_affine_2_0 + 0.092735250310891248 * p_affine_3_0,
                                           0.092735250310891248 * p_affine_0_1 + 0.092735250310891248 * p_affine_1_1 +
                                               0.72179424906732625 * p_affine_2_1 + 0.092735250310891248 * p_affine_3_1,
                                           0.092735250310891248 * p_affine_0_2 + 0.092735250310891248 * p_affine_1_2 +
                                               0.72179424906732625 * p_affine_2_2 + 0.092735250310891248 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id7 );
         Scalar_Variable_Coefficient_3D_k( 0.3108859192633005 * p_affine_0_0 + 0.067342242210098213 * p_affine_1_0 +
                                               0.31088591926330061 * p_affine_2_0 + 0.31088591926330061 * p_affine_3_0,
                                           0.3108859192633005 * p_affine_0_1 + 0.067342242210098213 * p_affine_1_1 +
                                               0.31088591926330061 * p_affine_2_1 + 0.31088591926330061 * p_affine_3_1,
                                           0.3108859192633005 * p_affine_0_2 + 0.067342242210098213 * p_affine_1_2 +
                                               0.31088591926330061 * p_affine_2_2 + 0.31088591926330061 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id8 );
         Scalar_Variable_Coefficient_3D_k( 0.092735250310891248 * p_affine_0_0 + 0.72179424906732625 * p_affine_1_0 +
                                               0.092735250310891248 * p_affine_2_0 + 0.092735250310891248 * p_affine_3_0,
                                           0.092735250310891248 * p_affine_0_1 + 0.72179424906732625 * p_affine_1_1 +
                                               0.092735250310891248 * p_affine_2_1 + 0.092735250310891248 * p_affine_3_1,
                                           0.092735250310891248 * p_affine_0_2 + 0.72179424906732625 * p_affine_1_2 +
                                               0.092735250310891248 * p_affine_2_2 + 0.092735250310891248 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id9 );
         Scalar_Variable_Coefficient_3D_k( 0.067342242210098102 * p_affine_0_0 + 0.31088591926330061 * p_affine_1_0 +
                                               0.31088591926330061 * p_affine_2_0 + 0.31088591926330061 * p_affine_3_0,
                                           0.067342242210098102 * p_affine_0_1 + 0.31088591926330061 * p_affine_1_1 +
                                               0.31088591926330061 * p_affine_2_1 + 0.31088591926330061 * p_affine_3_1,
                                           0.067342242210098102 * p_affine_0_2 + 0.31088591926330061 * p_affine_1_2 +
                                               0.31088591926330061 * p_affine_2_2 + 0.31088591926330061 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id10 );
         Scalar_Variable_Coefficient_3D_k( 0.72179424906732625 * p_affine_0_0 + 0.092735250310891248 * p_affine_1_0 +
                                               0.092735250310891248 * p_affine_2_0 + 0.092735250310891248 * p_affine_3_0,
                                           0.72179424906732625 * p_affine_0_1 + 0.092735250310891248 * p_affine_1_1 +
                                               0.092735250310891248 * p_affine_2_1 + 0.092735250310891248 * p_affine_3_1,
                                           0.72179424906732625 * p_affine_0_2 + 0.092735250310891248 * p_affine_1_2 +
                                               0.092735250310891248 * p_affine_2_2 + 0.092735250310891248 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id11 );
         Scalar_Variable_Coefficient_3D_k( 0.045503704125649636 * p_affine_0_0 + 0.045503704125649642 * p_affine_1_0 +
                                               0.45449629587435036 * p_affine_2_0 + 0.45449629587435036 * p_affine_3_0,
                                           0.045503704125649636 * p_affine_0_1 + 0.045503704125649642 * p_affine_1_1 +
                                               0.45449629587435036 * p_affine_2_1 + 0.45449629587435036 * p_affine_3_1,
                                           0.045503704125649636 * p_affine_0_2 + 0.045503704125649642 * p_affine_1_2 +
                                               0.45449629587435036 * p_affine_2_2 + 0.45449629587435036 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id12 );
         Scalar_Variable_Coefficient_3D_k( 0.045503704125649636 * p_affine_0_0 + 0.45449629587435036 * p_affine_1_0 +
                                               0.045503704125649642 * p_affine_2_0 + 0.45449629587435036 * p_affine_3_0,
                                           0.045503704125649636 * p_affine_0_1 + 0.45449629587435036 * p_affine_1_1 +
                                               0.045503704125649642 * p_affine_2_1 + 0.45449629587435036 * p_affine_3_1,
                                           0.045503704125649636 * p_affine_0_2 + 0.45449629587435036 * p_affine_1_2 +
                                               0.045503704125649642 * p_affine_2_2 + 0.45449629587435036 * p_affine_3_2,
                                           &Scalar_Variable_Coefficient_3D_k_out0_id13 );
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
         real_t tmp_13 = 0.018781320953002646 * Scalar_Variable_Coefficient_3D_k_out0_id0 * tmp_12;
         real_t tmp_14 = 0.012248840519393657 * Scalar_Variable_Coefficient_3D_k_out0_id1 * tmp_12;
         real_t tmp_15 = 0.018781320953002646 * Scalar_Variable_Coefficient_3D_k_out0_id10 * tmp_12;
         real_t tmp_16 = 0.012248840519393657 * Scalar_Variable_Coefficient_3D_k_out0_id11 * tmp_12;
         real_t tmp_17 = 0.0070910034628469103 * Scalar_Variable_Coefficient_3D_k_out0_id12 * tmp_12;
         real_t tmp_18 = 0.0070910034628469103 * Scalar_Variable_Coefficient_3D_k_out0_id13 * tmp_12;
         real_t tmp_19 = 0.0070910034628469103 * Scalar_Variable_Coefficient_3D_k_out0_id2 * tmp_12;
         real_t tmp_20 = 0.0070910034628469103 * Scalar_Variable_Coefficient_3D_k_out0_id3 * tmp_12;
         real_t tmp_21 = 0.0070910034628469103 * Scalar_Variable_Coefficient_3D_k_out0_id4 * tmp_12;
         real_t tmp_22 = 0.0070910034628469103 * Scalar_Variable_Coefficient_3D_k_out0_id5 * tmp_12;
         real_t tmp_23 = 0.018781320953002646 * Scalar_Variable_Coefficient_3D_k_out0_id6 * tmp_12;
         real_t tmp_24 = 0.012248840519393657 * Scalar_Variable_Coefficient_3D_k_out0_id7 * tmp_12;
         real_t tmp_25 = 0.018781320953002646 * Scalar_Variable_Coefficient_3D_k_out0_id8 * tmp_12;
         real_t tmp_26 = 0.012248840519393657 * Scalar_Variable_Coefficient_3D_k_out0_id9 * tmp_12;
         real_t a_0_0  = 0.3108859192633005 * tmp_13 + 0.092735250310891248 * tmp_14 + 0.067342242210098102 * tmp_15 +
                        0.72179424906732625 * tmp_16 + 0.045503704125649636 * tmp_17 + 0.045503704125649636 * tmp_18 +
                        0.45449629587435036 * tmp_19 + 0.045503704125649629 * tmp_20 + 0.45449629587435036 * tmp_21 +
                        0.45449629587435036 * tmp_22 + 0.3108859192633005 * tmp_23 + 0.092735250310891248 * tmp_24 +
                        0.3108859192633005 * tmp_25 + 0.092735250310891248 * tmp_26;
         real_t a_0_1 = 0;
         real_t a_0_2 = 0;
         real_t a_0_3 = 0;
         real_t a_1_0 = 0;
         real_t a_1_1 = 0.31088591926330061 * tmp_13 + 0.092735250310891248 * tmp_14 + 0.31088591926330061 * tmp_15 +
                        0.092735250310891248 * tmp_16 + 0.045503704125649642 * tmp_17 + 0.45449629587435036 * tmp_18 +
                        0.045503704125649642 * tmp_19 + 0.45449629587435036 * tmp_20 + 0.045503704125649642 * tmp_21 +
                        0.45449629587435036 * tmp_22 + 0.31088591926330061 * tmp_23 + 0.092735250310891248 * tmp_24 +
                        0.067342242210098213 * tmp_25 + 0.72179424906732625 * tmp_26;
         real_t a_1_2 = 0;
         real_t a_1_3 = 0;
         real_t a_2_0 = 0;
         real_t a_2_1 = 0;
         real_t a_2_2 = 0.31088591926330061 * tmp_13 + 0.092735250310891248 * tmp_14 + 0.31088591926330061 * tmp_15 +
                        0.092735250310891248 * tmp_16 + 0.45449629587435036 * tmp_17 + 0.045503704125649642 * tmp_18 +
                        0.045503704125649642 * tmp_19 + 0.45449629587435036 * tmp_20 + 0.45449629587435036 * tmp_21 +
                        0.045503704125649642 * tmp_22 + 0.067342242210098213 * tmp_23 + 0.72179424906732625 * tmp_24 +
                        0.31088591926330061 * tmp_25 + 0.092735250310891248 * tmp_26;
         real_t a_2_3 = 0;
         real_t a_3_0 = 0;
         real_t a_3_1 = 0;
         real_t a_3_2 = 0;
         real_t a_3_3 = 0.067342242210098213 * tmp_13 + 0.72179424906732625 * tmp_14 + 0.31088591926330061 * tmp_15 +
                        0.092735250310891248 * tmp_16 + 0.45449629587435036 * tmp_17 + 0.45449629587435036 * tmp_18 +
                        0.45449629587435036 * tmp_19 + 0.045503704125649642 * tmp_20 + 0.045503704125649642 * tmp_21 +
                        0.045503704125649642 * tmp_22 + 0.31088591926330061 * tmp_23 + 0.092735250310891248 * tmp_24 +
                        0.31088591926330061 * tmp_25 + 0.092735250310891248 * tmp_26;

         values[0] = a_0_0;
         values[1] = a_1_1;
         values[2] = a_2_2;
         values[3] = a_3_3;
         break;
      }
      }
   }

 private:
   void Scalar_Variable_Coefficient_2D_k( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_k( Point3D( { in_0, in_1, 0 } ) );
   }

   void Scalar_Variable_Coefficient_3D_k( real_t in_0, real_t in_1, real_t in_2, real_t* out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_k( Point3D( { in_0, in_1, in_2 } ) );
   }

   std::function< real_t( const Point3D& ) > callback_Scalar_Variable_Coefficient_2D_k;
   std::function< real_t( const Point3D& ) > callback_Scalar_Variable_Coefficient_3D_k;
};

} // namespace dg
} // namespace hyteg