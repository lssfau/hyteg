
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

# pragma once

# include "core/DataTypes.h"

# include "hyteg/dgfunctionspace/DGBasisInfo.hpp"
# include "hyteg/dgfunctionspace/DGForm.hpp"
# include "hyteg/dgfunctionspace/DGForm2D.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"


# include "Eigen/Eigen"

namespace hyteg {
namespace dg{
namespace eg{

class EGEpsilonFormEP1_0 : public hyteg::dg::DGForm
{

 public:
    EGEpsilonFormEP1_0(std::function< real_t ( const Point3D & ) > mu)
: callback_Scalar_Variable_Coefficient_3D_mu (mu)
, callback_Scalar_Variable_Coefficient_2D_mu (mu)
    {}

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_mu;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;

void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
}

      
 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = -tmp_3;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_2 - tmp_4*tmp_6);
      real_t tmp_8 = 1.0*tmp_7;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = tmp_2*tmp_8 + tmp_4*tmp_9;
      real_t tmp_11 = Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_10;
      real_t tmp_12 = -2*tmp_1*tmp_8 - 2*tmp_9;
      real_t tmp_13 = 0.5*tmp_7;
      real_t tmp_14 = tmp_0*tmp_13;
      real_t tmp_15 = tmp_13*tmp_3;
      real_t tmp_16 = tmp_1*tmp_13;
      real_t tmp_17 = tmp_0*tmp_15 + tmp_14*tmp_4 + tmp_16*tmp_5 + tmp_16*tmp_6;
      real_t tmp_18 = Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_17;
      real_t tmp_19 = -4*tmp_14 - 4*tmp_15;
      real_t tmp_20 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_21 = 0.054975871827660928*tmp_20;
      real_t tmp_22 = tmp_10*tmp_12;
      real_t tmp_23 = tmp_17*tmp_19;
      real_t tmp_24 = 0.11169079483900572*tmp_20;
      real_t tmp_25 = 0.054975871827660928*tmp_20;
      real_t tmp_26 = 0.11169079483900572*tmp_20;
      real_t tmp_27 = 0.054975871827660928*tmp_20;
      real_t tmp_28 = 0.11169079483900572*tmp_20;
      real_t tmp_29 = 2.0*tmp_7;
      real_t tmp_30 = tmp_1*tmp_29;
      real_t tmp_31 = tmp_29*tmp_3;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_17*tmp_31;
      real_t tmp_34 = tmp_29*tmp_5;
      real_t tmp_35 = tmp_0*tmp_29;
      real_t tmp_36 = tmp_10*tmp_34;
      real_t tmp_37 = tmp_17*tmp_35;
      real_t a_0_0 = tmp_21*(tmp_11*tmp_12 + tmp_18*tmp_19) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_23) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_23) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_23) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_23) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_23);
      real_t a_0_1 = tmp_21*(tmp_11*tmp_30 + tmp_18*tmp_31) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_33) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_33) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_33) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_33) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_33);
      real_t a_0_2 = tmp_21*(tmp_11*tmp_34 + tmp_18*tmp_35) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_37) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_37) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_37) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_37) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_37);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -tmp_2;
      real_t tmp_7 = 1.0 / (-tmp_1*tmp_6 + tmp_5);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_10 = tmp_7*(0.21132486540518713*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = tmp_7*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_10*tmp_2 + tmp_13*tmp_4;
      real_t tmp_15 = tmp_14 - 1.0/3.0;
      real_t tmp_16 = tmp_0*tmp_13 + tmp_10*tmp_3;
      real_t tmp_17 = tmp_16 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_0*(tmp_1*tmp_15 + tmp_17*tmp_4);
      real_t tmp_19 = 0.5*tmp_7;
      real_t tmp_20 = tmp_19*tmp_3;
      real_t tmp_21 = tmp_19*tmp_2;
      real_t tmp_22 = -tmp_20 - tmp_21;
      real_t tmp_23 = 0.5*tmp_22;
      real_t tmp_24 = tmp_15*tmp_3 + tmp_17*tmp_6;
      real_t tmp_25 = 1.0*tmp_7;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_0*tmp_25;
      real_t tmp_28 = 0.5*p_affine_10_0*(-tmp_26 - tmp_27) + 0.5*p_affine_10_1*tmp_22;
      real_t tmp_29 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_30 = 1.0 / (tmp_29);
      real_t tmp_31 = -tmp_14 - tmp_16 + 1;
      real_t tmp_32 = tmp_19*tmp_4;
      real_t tmp_33 = 0.5*p_affine_10_0*(tmp_25*tmp_5 + tmp_27*tmp_6) + 0.5*p_affine_10_1*(tmp_0*tmp_32 + tmp_1*tmp_32 + tmp_20*tmp_6 + tmp_21*tmp_3);
      real_t tmp_34 = 2*tmp_29;
      real_t tmp_35 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_34;
      real_t tmp_36 = tmp_7*(0.78867513459481287*tmp_8 + tmp_9);
      real_t tmp_37 = tmp_7*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_38 = tmp_2*tmp_36 + tmp_37*tmp_4;
      real_t tmp_39 = tmp_38 - 1.0/3.0;
      real_t tmp_40 = tmp_0*tmp_37 + tmp_3*tmp_36;
      real_t tmp_41 = tmp_40 - 1.0/3.0;
      real_t tmp_42 = p_affine_10_0*(tmp_1*tmp_39 + tmp_4*tmp_41);
      real_t tmp_43 = tmp_3*tmp_39 + tmp_41*tmp_6;
      real_t tmp_44 = -tmp_38 - tmp_40 + 1;
      real_t tmp_45 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_34;
      real_t tmp_46 = 0.25*tmp_7;
      real_t tmp_47 = tmp_2*tmp_46;
      real_t tmp_48 = 0.5*p_affine_10_0*tmp_26 + 0.5*p_affine_10_1*tmp_21;
      real_t tmp_49 = tmp_3*tmp_46;
      real_t tmp_50 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*tmp_20;
      real_t a_0_0 = tmp_35*(-tmp_18*tmp_23 - tmp_24*tmp_28 + tmp_24*tmp_30*tmp_31 - tmp_31*tmp_33) + tmp_45*(-tmp_23*tmp_42 - tmp_28*tmp_43 + tmp_30*tmp_43*tmp_44 - tmp_33*tmp_44);
      real_t a_0_1 = tmp_35*(tmp_14*tmp_24*tmp_30 - tmp_14*tmp_33 - tmp_18*tmp_47 - tmp_24*tmp_48) + tmp_45*(tmp_30*tmp_38*tmp_43 - tmp_33*tmp_38 - tmp_42*tmp_47 - tmp_43*tmp_48);
      real_t a_0_2 = tmp_35*(tmp_16*tmp_24*tmp_30 - tmp_16*tmp_33 - tmp_18*tmp_49 - tmp_24*tmp_50) + tmp_45*(tmp_30*tmp_40*tmp_43 - tmp_33*tmp_40 - tmp_42*tmp_49 - tmp_43*tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElementInner[0]( 0 );
      const auto p_affine_0_1 = coordsElementInner[0]( 1 );

      const auto p_affine_1_0 = coordsElementInner[1]( 0 );
      const auto p_affine_1_1 = coordsElementInner[1]( 1 );

      const auto p_affine_2_0 = coordsElementInner[2]( 0 );
      const auto p_affine_2_1 = coordsElementInner[2]( 1 );

      const auto p_affine_3_0 = coordsElementOuter[0]( 0 );
      const auto p_affine_3_1 = coordsElementOuter[0]( 1 );

      const auto p_affine_4_0 = coordsElementOuter[1]( 0 );
      const auto p_affine_4_1 = coordsElementOuter[1]( 1 );

      const auto p_affine_5_0 = coordsElementOuter[2]( 0 );
      const auto p_affine_5_1 = coordsElementOuter[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -tmp_2;
      real_t tmp_7 = 1.0 / (-tmp_1*tmp_6 + tmp_5);
      real_t tmp_8 = -p_affine_0_1;
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.21132486540518713*tmp_9;
      real_t tmp_11 = tmp_7*(tmp_10 + tmp_8);
      real_t tmp_12 = -p_affine_0_0;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + 0.21132486540518713*tmp_13;
      real_t tmp_15 = tmp_7*(tmp_12 + tmp_14);
      real_t tmp_16 = tmp_11*tmp_2 + tmp_15*tmp_4 - 1.0/3.0;
      real_t tmp_17 = tmp_0*tmp_15 + tmp_11*tmp_3 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_0*(tmp_1*tmp_16 + tmp_17*tmp_4);
      real_t tmp_19 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_20 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_21 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_22 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_23 = 1.0 / (tmp_19*tmp_20 - tmp_21*tmp_22);
      real_t tmp_24 = 0.5*tmp_23;
      real_t tmp_25 = tmp_19*tmp_24;
      real_t tmp_26 = tmp_21*tmp_24;
      real_t tmp_27 = -tmp_25 - tmp_26;
      real_t tmp_28 = 0.5*tmp_27;
      real_t tmp_29 = tmp_16*tmp_3 + tmp_17*tmp_6;
      real_t tmp_30 = 1.0*tmp_23;
      real_t tmp_31 = tmp_20*tmp_30;
      real_t tmp_32 = tmp_22*tmp_30;
      real_t tmp_33 = 0.5*p_affine_10_0*(-tmp_31 - tmp_32) + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_34 = -p_affine_3_1;
      real_t tmp_35 = tmp_23*(tmp_10 + tmp_34);
      real_t tmp_36 = -p_affine_3_0;
      real_t tmp_37 = tmp_23*(tmp_14 + tmp_36);
      real_t tmp_38 = tmp_20*tmp_37 + tmp_21*tmp_35;
      real_t tmp_39 = tmp_19*tmp_35 + tmp_22*tmp_37;
      real_t tmp_40 = -tmp_38 - tmp_39 + 1;
      real_t tmp_41 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_42 = 1.0 / (tmp_41);
      real_t tmp_43 = tmp_29*tmp_42;
      real_t tmp_44 = 1.0*tmp_7;
      real_t tmp_45 = 0.5*tmp_7;
      real_t tmp_46 = tmp_3*tmp_45;
      real_t tmp_47 = tmp_4*tmp_45;
      real_t tmp_48 = p_affine_10_0*(tmp_0*tmp_44*tmp_6 + tmp_44*tmp_5) + p_affine_10_1*(tmp_0*tmp_47 + tmp_1*tmp_47 + tmp_2*tmp_46 + tmp_46*tmp_6);
      real_t tmp_49 = 2*tmp_41;
      real_t tmp_50 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_49;
      real_t tmp_51 = p_affine_6_1 + 0.78867513459481287*tmp_9;
      real_t tmp_52 = tmp_7*(tmp_51 + tmp_8);
      real_t tmp_53 = p_affine_6_0 + 0.78867513459481287*tmp_13;
      real_t tmp_54 = tmp_7*(tmp_12 + tmp_53);
      real_t tmp_55 = tmp_2*tmp_52 + tmp_4*tmp_54 - 1.0/3.0;
      real_t tmp_56 = tmp_0*tmp_54 + tmp_3*tmp_52 - 1.0/3.0;
      real_t tmp_57 = p_affine_10_0*(tmp_1*tmp_55 + tmp_4*tmp_56);
      real_t tmp_58 = tmp_3*tmp_55 + tmp_56*tmp_6;
      real_t tmp_59 = tmp_23*(tmp_34 + tmp_51);
      real_t tmp_60 = tmp_23*(tmp_36 + tmp_53);
      real_t tmp_61 = tmp_20*tmp_60 + tmp_21*tmp_59;
      real_t tmp_62 = tmp_19*tmp_59 + tmp_22*tmp_60;
      real_t tmp_63 = -tmp_61 - tmp_62 + 1;
      real_t tmp_64 = tmp_42*tmp_58;
      real_t tmp_65 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_49;
      real_t tmp_66 = 0.25*tmp_23;
      real_t tmp_67 = tmp_21*tmp_66;
      real_t tmp_68 = 0.5*p_affine_10_0*tmp_31 + 0.5*p_affine_10_1*tmp_26;
      real_t tmp_69 = tmp_19*tmp_66;
      real_t tmp_70 = 0.5*p_affine_10_0*tmp_32 + 0.5*p_affine_10_1*tmp_25;
      real_t a_0_0 = tmp_50*(-tmp_18*tmp_28 - tmp_29*tmp_33 - tmp_40*tmp_43 + 0.5*tmp_40*tmp_48) + tmp_65*(-tmp_28*tmp_57 - tmp_33*tmp_58 + 0.5*tmp_48*tmp_63 - tmp_63*tmp_64);
      real_t a_0_1 = tmp_50*(-tmp_18*tmp_67 - tmp_29*tmp_68 - tmp_38*tmp_43 + 0.5*tmp_38*tmp_48) + tmp_65*(0.5*tmp_48*tmp_61 - tmp_57*tmp_67 - tmp_58*tmp_68 - tmp_61*tmp_64);
      real_t a_0_2 = tmp_50*(-tmp_18*tmp_69 - tmp_29*tmp_70 - tmp_39*tmp_43 + 0.5*tmp_39*tmp_48) + tmp_65*(0.5*tmp_48*tmp_62 - tmp_57*tmp_69 - tmp_58*tmp_70 - tmp_62*tmp_64);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t tmp_0 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = -tmp_2;
      real_t tmp_6 = 1.0 / (-tmp_1*tmp_5 + tmp_3*tmp_4);
      real_t tmp_7 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_8 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_9 = tmp_6*(0.21132486540518713*tmp_7 + tmp_8);
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_12 = tmp_6*(0.21132486540518713*tmp_10 + tmp_11);
      real_t tmp_13 = tmp_12*tmp_4 + tmp_2*tmp_9 - 1.0/3.0;
      real_t tmp_14 = tmp_0*tmp_12 + tmp_3*tmp_9 - 1.0/3.0;
      real_t tmp_15 = tmp_1*tmp_13 + tmp_14*tmp_4;
      real_t tmp_16 = 0.5*tmp_6;
      real_t tmp_17 = tmp_16*tmp_3;
      real_t tmp_18 = tmp_16*tmp_2;
      real_t tmp_19 = -tmp_17 - tmp_18;
      real_t tmp_20 = p_affine_10_0*tmp_19;
      real_t tmp_21 = 1.0*tmp_6;
      real_t tmp_22 = tmp_21*tmp_4;
      real_t tmp_23 = tmp_0*tmp_21;
      real_t tmp_24 = p_affine_10_0*(-tmp_22 - tmp_23) + p_affine_10_1*tmp_19;
      real_t tmp_25 = tmp_13*tmp_3 + tmp_14*tmp_5;
      real_t tmp_26 = tmp_6*(0.78867513459481287*tmp_7 + tmp_8);
      real_t tmp_27 = tmp_6*(0.78867513459481287*tmp_10 + tmp_11);
      real_t tmp_28 = tmp_2*tmp_26 + tmp_27*tmp_4 - 1.0/3.0;
      real_t tmp_29 = tmp_0*tmp_27 + tmp_26*tmp_3 - 1.0/3.0;
      real_t tmp_30 = tmp_1*tmp_28 + tmp_29*tmp_4;
      real_t tmp_31 = tmp_28*tmp_3 + tmp_29*tmp_5;
      real_t tmp_32 = p_affine_10_0*tmp_18;
      real_t tmp_33 = p_affine_10_0*tmp_22 + p_affine_10_1*tmp_18;
      real_t tmp_34 = p_affine_10_0*tmp_17;
      real_t tmp_35 = p_affine_10_0*tmp_23 + p_affine_10_1*tmp_17;
      real_t a_0_0 = -0.5*tmp_15*tmp_20 - 0.5*tmp_20*tmp_30 - 0.5*tmp_24*tmp_25 - 0.5*tmp_24*tmp_31;
      real_t a_0_1 = -0.5*tmp_15*tmp_32 - 0.5*tmp_25*tmp_33 - 0.5*tmp_30*tmp_32 - 0.5*tmp_31*tmp_33;
      real_t a_0_2 = -0.5*tmp_15*tmp_34 - 0.5*tmp_25*tmp_35 - 0.5*tmp_30*tmp_34 - 0.5*tmp_31*tmp_35;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

  void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      elMat( 0, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

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

      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id6 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id7 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id8 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id9 );
      Scalar_Variable_Coefficient_3D_mu( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id10 );
      Scalar_Variable_Coefficient_3D_mu( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id11 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id12 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id13 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_5 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_9 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_10 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_11 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = tmp_11*tmp_2;
      real_t tmp_14 = tmp_1*tmp_9;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_3 - tmp_0*tmp_6 + tmp_10*tmp_12 - tmp_10*tmp_14 - tmp_13*tmp_8 + tmp_4*tmp_8*tmp_9);
      real_t tmp_16 = 1.0*tmp_15;
      real_t tmp_17 = tmp_16*tmp_7;
      real_t tmp_18 = -tmp_13 + tmp_4*tmp_9;
      real_t tmp_19 = tmp_16*tmp_18;
      real_t tmp_20 = tmp_12 - tmp_14;
      real_t tmp_21 = tmp_16*tmp_20;
      real_t tmp_22 = tmp_0*tmp_17 + tmp_10*tmp_21 + tmp_19*tmp_8;
      real_t tmp_23 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_22;
      real_t tmp_24 = -2*tmp_17 - 2*tmp_19 - 2*tmp_21;
      real_t tmp_25 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_26 = tmp_0*tmp_1 - tmp_11*tmp_8;
      real_t tmp_27 = 0.5*tmp_15;
      real_t tmp_28 = tmp_26*tmp_27;
      real_t tmp_29 = -tmp_0*tmp_4 + tmp_10*tmp_11;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = -tmp_1*tmp_10 + tmp_4*tmp_8;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_27*tmp_7;
      real_t tmp_34 = tmp_18*tmp_27;
      real_t tmp_35 = tmp_20*tmp_27;
      real_t tmp_36 = tmp_0*tmp_32 + tmp_10*tmp_28 + tmp_2*tmp_35 + tmp_30*tmp_8 + tmp_33*tmp_9 + tmp_34*tmp_5;
      real_t tmp_37 = tmp_36*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_38 = -tmp_0*tmp_5 + tmp_8*tmp_9;
      real_t tmp_39 = tmp_27*tmp_38;
      real_t tmp_40 = tmp_0*tmp_2 - tmp_10*tmp_9;
      real_t tmp_41 = tmp_27*tmp_40;
      real_t tmp_42 = tmp_10*tmp_5 - tmp_2*tmp_8;
      real_t tmp_43 = tmp_27*tmp_42;
      real_t tmp_44 = tmp_0*tmp_43 + tmp_1*tmp_34 + tmp_10*tmp_39 + tmp_11*tmp_33 + tmp_35*tmp_4 + tmp_41*tmp_8;
      real_t tmp_45 = tmp_44*(-tmp_39 - tmp_41 - tmp_43);
      real_t tmp_46 = p_affine_0_0*p_affine_1_1;
      real_t tmp_47 = p_affine_0_0*p_affine_1_2;
      real_t tmp_48 = p_affine_2_1*p_affine_3_2;
      real_t tmp_49 = p_affine_0_1*p_affine_1_0;
      real_t tmp_50 = p_affine_0_1*p_affine_1_2;
      real_t tmp_51 = p_affine_2_2*p_affine_3_0;
      real_t tmp_52 = p_affine_0_2*p_affine_1_0;
      real_t tmp_53 = p_affine_0_2*p_affine_1_1;
      real_t tmp_54 = p_affine_2_0*p_affine_3_1;
      real_t tmp_55 = p_affine_2_2*p_affine_3_1;
      real_t tmp_56 = p_affine_2_0*p_affine_3_2;
      real_t tmp_57 = p_affine_2_1*p_affine_3_0;
      real_t tmp_58 = std::abs(p_affine_0_0*tmp_48 - p_affine_0_0*tmp_55 + p_affine_0_1*tmp_51 - p_affine_0_1*tmp_56 + p_affine_0_2*tmp_54 - p_affine_0_2*tmp_57 - p_affine_1_0*tmp_48 + p_affine_1_0*tmp_55 - p_affine_1_1*tmp_51 + p_affine_1_1*tmp_56 - p_affine_1_2*tmp_54 + p_affine_1_2*tmp_57 + p_affine_2_0*tmp_50 - p_affine_2_0*tmp_53 - p_affine_2_1*tmp_47 + p_affine_2_1*tmp_52 + p_affine_2_2*tmp_46 - p_affine_2_2*tmp_49 - p_affine_3_0*tmp_50 + p_affine_3_0*tmp_53 + p_affine_3_1*tmp_47 - p_affine_3_1*tmp_52 - p_affine_3_2*tmp_46 + p_affine_3_2*tmp_49);
      real_t tmp_59 = 0.018781320953002646*tmp_58;
      real_t tmp_60 = tmp_22*tmp_24;
      real_t tmp_61 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id1;
      real_t tmp_62 = 0.012248840519393657*tmp_58;
      real_t tmp_63 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id2;
      real_t tmp_64 = 0.0070910034628469103*tmp_58;
      real_t tmp_65 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id3;
      real_t tmp_66 = 0.0070910034628469103*tmp_58;
      real_t tmp_67 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id4;
      real_t tmp_68 = 0.0070910034628469103*tmp_58;
      real_t tmp_69 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id5;
      real_t tmp_70 = 0.0070910034628469103*tmp_58;
      real_t tmp_71 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id6;
      real_t tmp_72 = 0.018781320953002646*tmp_58;
      real_t tmp_73 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id7;
      real_t tmp_74 = 0.012248840519393657*tmp_58;
      real_t tmp_75 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id8;
      real_t tmp_76 = 0.018781320953002646*tmp_58;
      real_t tmp_77 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id9;
      real_t tmp_78 = 0.012248840519393657*tmp_58;
      real_t tmp_79 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id10;
      real_t tmp_80 = 0.018781320953002646*tmp_58;
      real_t tmp_81 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id11;
      real_t tmp_82 = 0.012248840519393657*tmp_58;
      real_t tmp_83 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id12;
      real_t tmp_84 = 0.0070910034628469103*tmp_58;
      real_t tmp_85 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id13;
      real_t tmp_86 = 0.0070910034628469103*tmp_58;
      real_t tmp_87 = 2.0*tmp_15;
      real_t tmp_88 = tmp_7*tmp_87;
      real_t tmp_89 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_87;
      real_t tmp_90 = tmp_31*tmp_36;
      real_t tmp_91 = tmp_42*tmp_44;
      real_t tmp_92 = tmp_22*tmp_88;
      real_t tmp_93 = Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_87;
      real_t tmp_94 = Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_87;
      real_t tmp_95 = Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_87;
      real_t tmp_96 = Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_87;
      real_t tmp_97 = Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_87;
      real_t tmp_98 = Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_87;
      real_t tmp_99 = Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_87;
      real_t tmp_100 = Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_87;
      real_t tmp_101 = Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_87;
      real_t tmp_102 = Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_87;
      real_t tmp_103 = Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_87;
      real_t tmp_104 = Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_87;
      real_t tmp_105 = Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_87;
      real_t tmp_106 = tmp_23*tmp_87;
      real_t tmp_107 = tmp_29*tmp_36;
      real_t tmp_108 = tmp_40*tmp_44;
      real_t tmp_109 = tmp_18*tmp_22;
      real_t tmp_110 = tmp_26*tmp_36;
      real_t tmp_111 = tmp_38*tmp_44;
      real_t tmp_112 = tmp_20*tmp_22;
      real_t a_0_0 = tmp_59*(tmp_23*tmp_24 + tmp_25*tmp_37 + tmp_25*tmp_45) + tmp_62*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_60 + tmp_37*tmp_61 + tmp_45*tmp_61) + tmp_64*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_60 + tmp_37*tmp_63 + tmp_45*tmp_63) + tmp_66*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_60 + tmp_37*tmp_65 + tmp_45*tmp_65) + tmp_68*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_60 + tmp_37*tmp_67 + tmp_45*tmp_67) + tmp_70*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_60 + tmp_37*tmp_69 + tmp_45*tmp_69) + tmp_72*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_60 + tmp_37*tmp_71 + tmp_45*tmp_71) + tmp_74*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_60 + tmp_37*tmp_73 + tmp_45*tmp_73) + tmp_76*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_60 + tmp_37*tmp_75 + tmp_45*tmp_75) + tmp_78*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_60 + tmp_37*tmp_77 + tmp_45*tmp_77) + tmp_80*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_60 + tmp_37*tmp_79 + tmp_45*tmp_79) + tmp_82*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_60 + tmp_37*tmp_81 + tmp_45*tmp_81) + tmp_84*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_60 + tmp_37*tmp_83 + tmp_45*tmp_83) + tmp_86*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_60 + tmp_37*tmp_85 + tmp_45*tmp_85);
      real_t a_0_1 = tmp_59*(tmp_23*tmp_88 + tmp_89*tmp_90 + tmp_89*tmp_91) + tmp_62*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_92 + tmp_90*tmp_93 + tmp_91*tmp_93) + tmp_64*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_92 + tmp_90*tmp_94 + tmp_91*tmp_94) + tmp_66*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_92 + tmp_90*tmp_95 + tmp_91*tmp_95) + tmp_68*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_92 + tmp_90*tmp_96 + tmp_91*tmp_96) + tmp_70*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_92 + tmp_90*tmp_97 + tmp_91*tmp_97) + tmp_72*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_92 + tmp_90*tmp_98 + tmp_91*tmp_98) + tmp_74*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_92 + tmp_90*tmp_99 + tmp_91*tmp_99) + tmp_76*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_92 + tmp_100*tmp_90 + tmp_100*tmp_91) + tmp_78*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_92 + tmp_101*tmp_90 + tmp_101*tmp_91) + tmp_80*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_92 + tmp_102*tmp_90 + tmp_102*tmp_91) + tmp_82*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_92 + tmp_103*tmp_90 + tmp_103*tmp_91) + tmp_84*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_92 + tmp_104*tmp_90 + tmp_104*tmp_91) + tmp_86*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_92 + tmp_105*tmp_90 + tmp_105*tmp_91);
      real_t a_0_2 = tmp_59*(tmp_106*tmp_18 + tmp_107*tmp_89 + tmp_108*tmp_89) + tmp_62*(tmp_107*tmp_93 + tmp_108*tmp_93 + tmp_109*tmp_93) + tmp_64*(tmp_107*tmp_94 + tmp_108*tmp_94 + tmp_109*tmp_94) + tmp_66*(tmp_107*tmp_95 + tmp_108*tmp_95 + tmp_109*tmp_95) + tmp_68*(tmp_107*tmp_96 + tmp_108*tmp_96 + tmp_109*tmp_96) + tmp_70*(tmp_107*tmp_97 + tmp_108*tmp_97 + tmp_109*tmp_97) + tmp_72*(tmp_107*tmp_98 + tmp_108*tmp_98 + tmp_109*tmp_98) + tmp_74*(tmp_107*tmp_99 + tmp_108*tmp_99 + tmp_109*tmp_99) + tmp_76*(tmp_100*tmp_107 + tmp_100*tmp_108 + tmp_100*tmp_109) + tmp_78*(tmp_101*tmp_107 + tmp_101*tmp_108 + tmp_101*tmp_109) + tmp_80*(tmp_102*tmp_107 + tmp_102*tmp_108 + tmp_102*tmp_109) + tmp_82*(tmp_103*tmp_107 + tmp_103*tmp_108 + tmp_103*tmp_109) + tmp_84*(tmp_104*tmp_107 + tmp_104*tmp_108 + tmp_104*tmp_109) + tmp_86*(tmp_105*tmp_107 + tmp_105*tmp_108 + tmp_105*tmp_109);
      real_t a_0_3 = tmp_59*(tmp_106*tmp_20 + tmp_110*tmp_89 + tmp_111*tmp_89) + tmp_62*(tmp_110*tmp_93 + tmp_111*tmp_93 + tmp_112*tmp_93) + tmp_64*(tmp_110*tmp_94 + tmp_111*tmp_94 + tmp_112*tmp_94) + tmp_66*(tmp_110*tmp_95 + tmp_111*tmp_95 + tmp_112*tmp_95) + tmp_68*(tmp_110*tmp_96 + tmp_111*tmp_96 + tmp_112*tmp_96) + tmp_70*(tmp_110*tmp_97 + tmp_111*tmp_97 + tmp_112*tmp_97) + tmp_72*(tmp_110*tmp_98 + tmp_111*tmp_98 + tmp_112*tmp_98) + tmp_74*(tmp_110*tmp_99 + tmp_111*tmp_99 + tmp_112*tmp_99) + tmp_76*(tmp_100*tmp_110 + tmp_100*tmp_111 + tmp_100*tmp_112) + tmp_78*(tmp_101*tmp_110 + tmp_101*tmp_111 + tmp_101*tmp_112) + tmp_80*(tmp_102*tmp_110 + tmp_102*tmp_111 + tmp_102*tmp_112) + tmp_82*(tmp_103*tmp_110 + tmp_103*tmp_111 + tmp_103*tmp_112) + tmp_84*(tmp_104*tmp_110 + tmp_104*tmp_111 + tmp_104*tmp_112) + tmp_86*(tmp_105*tmp_110 + tmp_105*tmp_111 + tmp_105*tmp_112);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_0*tmp_10;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = tmp_4*tmp_9;
      real_t tmp_15 = 1.0 / (-tmp_1*tmp_13 + tmp_1*tmp_2*tmp_9 + tmp_11*tmp_3 - tmp_12*tmp_6 - tmp_14*tmp_3 + tmp_6*tmp_8);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_8;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = tmp_33 - 1.0/4.0;
      real_t tmp_35 = tmp_0*tmp_3 - tmp_2*tmp_6;
      real_t tmp_36 = -tmp_3*tmp_9 + tmp_6*tmp_7;
      real_t tmp_37 = -tmp_13 + tmp_2*tmp_9;
      real_t tmp_38 = tmp_20*tmp_35 + tmp_26*tmp_36 + tmp_32*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/4.0;
      real_t tmp_40 = -tmp_0*tmp_1 + tmp_4*tmp_6;
      real_t tmp_41 = tmp_1*tmp_9 - tmp_10*tmp_6;
      real_t tmp_42 = tmp_11 - tmp_14;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_26*tmp_41 + tmp_32*tmp_42;
      real_t tmp_44 = tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_34 + tmp_2*tmp_44 + tmp_39*tmp_4;
      real_t tmp_46 = 0.5*tmp_15;
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_36*tmp_46;
      real_t tmp_49 = tmp_21*tmp_46;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = 0.5*p_affine_13_0;
      real_t tmp_52 = tmp_50*tmp_51;
      real_t tmp_53 = tmp_10*tmp_39 + tmp_34*tmp_9 + tmp_44*tmp_7;
      real_t tmp_54 = tmp_40*tmp_46;
      real_t tmp_55 = tmp_35*tmp_46;
      real_t tmp_56 = tmp_46*tmp_5;
      real_t tmp_57 = -tmp_54 - tmp_55 - tmp_56;
      real_t tmp_58 = tmp_51*tmp_57;
      real_t tmp_59 = tmp_1*tmp_39 + tmp_3*tmp_44 + tmp_34*tmp_6;
      real_t tmp_60 = 1.0*tmp_15;
      real_t tmp_61 = tmp_42*tmp_60;
      real_t tmp_62 = tmp_37*tmp_60;
      real_t tmp_63 = tmp_27*tmp_60;
      real_t tmp_64 = 0.5*p_affine_13_0*(-tmp_61 - tmp_62 - tmp_63) + 0.5*p_affine_13_1*tmp_50 + 0.5*p_affine_13_2*tmp_57;
      real_t tmp_65 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_66 = std::pow(tmp_65, -0.25);
      real_t tmp_67 = -tmp_33 - tmp_38 - tmp_43 + 1;
      real_t tmp_68 = tmp_27*tmp_46;
      real_t tmp_69 = tmp_37*tmp_46;
      real_t tmp_70 = tmp_42*tmp_46;
      real_t tmp_71 = 0.5*p_affine_13_0*(tmp_1*tmp_62 + tmp_3*tmp_61 + tmp_6*tmp_63) + 0.5*p_affine_13_1*(tmp_0*tmp_68 + tmp_1*tmp_48 + tmp_2*tmp_70 + tmp_3*tmp_47 + tmp_4*tmp_69 + tmp_49*tmp_6) + 0.5*p_affine_13_2*(tmp_1*tmp_55 + tmp_10*tmp_69 + tmp_3*tmp_54 + tmp_56*tmp_6 + tmp_68*tmp_9 + tmp_7*tmp_70);
      real_t tmp_72 = 2.0*std::pow(tmp_65, 1.0/2.0);
      real_t tmp_73 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_72;
      real_t tmp_74 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_75 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_76 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_77 = tmp_21*tmp_75 + tmp_27*tmp_76 + tmp_5*tmp_74;
      real_t tmp_78 = tmp_77 - 1.0/4.0;
      real_t tmp_79 = tmp_35*tmp_74 + tmp_36*tmp_75 + tmp_37*tmp_76;
      real_t tmp_80 = tmp_79 - 1.0/4.0;
      real_t tmp_81 = tmp_40*tmp_74 + tmp_41*tmp_75 + tmp_42*tmp_76;
      real_t tmp_82 = tmp_81 - 1.0/4.0;
      real_t tmp_83 = tmp_0*tmp_78 + tmp_2*tmp_82 + tmp_4*tmp_80;
      real_t tmp_84 = tmp_10*tmp_80 + tmp_7*tmp_82 + tmp_78*tmp_9;
      real_t tmp_85 = tmp_1*tmp_80 + tmp_3*tmp_82 + tmp_6*tmp_78;
      real_t tmp_86 = -tmp_77 - tmp_79 - tmp_81 + 1;
      real_t tmp_87 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_72;
      real_t tmp_88 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_89 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_90 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_91 = tmp_21*tmp_89 + tmp_27*tmp_90 + tmp_5*tmp_88;
      real_t tmp_92 = tmp_91 - 1.0/4.0;
      real_t tmp_93 = tmp_35*tmp_88 + tmp_36*tmp_89 + tmp_37*tmp_90;
      real_t tmp_94 = tmp_93 - 1.0/4.0;
      real_t tmp_95 = tmp_40*tmp_88 + tmp_41*tmp_89 + tmp_42*tmp_90;
      real_t tmp_96 = tmp_95 - 1.0/4.0;
      real_t tmp_97 = tmp_0*tmp_92 + tmp_2*tmp_96 + tmp_4*tmp_94;
      real_t tmp_98 = tmp_10*tmp_94 + tmp_7*tmp_96 + tmp_9*tmp_92;
      real_t tmp_99 = tmp_1*tmp_94 + tmp_3*tmp_96 + tmp_6*tmp_92;
      real_t tmp_100 = -tmp_91 - tmp_93 - tmp_95 + 1;
      real_t tmp_101 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_72;
      real_t tmp_102 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_103 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_104 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_105 = tmp_102*tmp_5 + tmp_103*tmp_21 + tmp_104*tmp_27;
      real_t tmp_106 = tmp_105 - 1.0/4.0;
      real_t tmp_107 = tmp_102*tmp_35 + tmp_103*tmp_36 + tmp_104*tmp_37;
      real_t tmp_108 = tmp_107 - 1.0/4.0;
      real_t tmp_109 = tmp_102*tmp_40 + tmp_103*tmp_41 + tmp_104*tmp_42;
      real_t tmp_110 = tmp_109 - 1.0/4.0;
      real_t tmp_111 = tmp_0*tmp_106 + tmp_108*tmp_4 + tmp_110*tmp_2;
      real_t tmp_112 = tmp_10*tmp_108 + tmp_106*tmp_9 + tmp_110*tmp_7;
      real_t tmp_113 = tmp_1*tmp_108 + tmp_106*tmp_6 + tmp_110*tmp_3;
      real_t tmp_114 = -tmp_105 - tmp_107 - tmp_109 + 1;
      real_t tmp_115 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_72;
      real_t tmp_116 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_117 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_118 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_119 = tmp_116*tmp_5 + tmp_117*tmp_21 + tmp_118*tmp_27;
      real_t tmp_120 = tmp_119 - 1.0/4.0;
      real_t tmp_121 = tmp_116*tmp_35 + tmp_117*tmp_36 + tmp_118*tmp_37;
      real_t tmp_122 = tmp_121 - 1.0/4.0;
      real_t tmp_123 = tmp_116*tmp_40 + tmp_117*tmp_41 + tmp_118*tmp_42;
      real_t tmp_124 = tmp_123 - 1.0/4.0;
      real_t tmp_125 = tmp_0*tmp_120 + tmp_122*tmp_4 + tmp_124*tmp_2;
      real_t tmp_126 = tmp_10*tmp_122 + tmp_120*tmp_9 + tmp_124*tmp_7;
      real_t tmp_127 = tmp_1*tmp_122 + tmp_120*tmp_6 + tmp_124*tmp_3;
      real_t tmp_128 = -tmp_119 - tmp_121 - tmp_123 + 1;
      real_t tmp_129 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_72;
      real_t tmp_130 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_131 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_132 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_133 = tmp_130*tmp_5 + tmp_131*tmp_21 + tmp_132*tmp_27;
      real_t tmp_134 = tmp_133 - 1.0/4.0;
      real_t tmp_135 = tmp_130*tmp_35 + tmp_131*tmp_36 + tmp_132*tmp_37;
      real_t tmp_136 = tmp_135 - 1.0/4.0;
      real_t tmp_137 = tmp_130*tmp_40 + tmp_131*tmp_41 + tmp_132*tmp_42;
      real_t tmp_138 = tmp_137 - 1.0/4.0;
      real_t tmp_139 = tmp_0*tmp_134 + tmp_136*tmp_4 + tmp_138*tmp_2;
      real_t tmp_140 = tmp_10*tmp_136 + tmp_134*tmp_9 + tmp_138*tmp_7;
      real_t tmp_141 = tmp_1*tmp_136 + tmp_134*tmp_6 + tmp_138*tmp_3;
      real_t tmp_142 = -tmp_133 - tmp_135 - tmp_137 + 1;
      real_t tmp_143 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_72;
      real_t tmp_144 = 0.25*p_affine_13_0*tmp_15;
      real_t tmp_145 = tmp_144*tmp_5;
      real_t tmp_146 = tmp_144*tmp_21;
      real_t tmp_147 = 0.5*p_affine_13_0*tmp_63 + 0.5*p_affine_13_1*tmp_49 + 0.5*p_affine_13_2*tmp_56;
      real_t tmp_148 = tmp_144*tmp_35;
      real_t tmp_149 = tmp_144*tmp_36;
      real_t tmp_150 = 0.5*p_affine_13_0*tmp_62 + 0.5*p_affine_13_1*tmp_48 + 0.5*p_affine_13_2*tmp_55;
      real_t tmp_151 = tmp_144*tmp_40;
      real_t tmp_152 = tmp_144*tmp_41;
      real_t tmp_153 = 0.5*p_affine_13_0*tmp_61 + 0.5*p_affine_13_1*tmp_47 + 0.5*p_affine_13_2*tmp_54;
      real_t a_0_0 = tmp_101*(1.0*tmp_100*tmp_66*tmp_99 - tmp_100*tmp_71 - tmp_52*tmp_97 - tmp_58*tmp_98 - tmp_64*tmp_99) + tmp_115*(-tmp_111*tmp_52 - tmp_112*tmp_58 + 1.0*tmp_113*tmp_114*tmp_66 - tmp_113*tmp_64 - tmp_114*tmp_71) + tmp_129*(-tmp_125*tmp_52 - tmp_126*tmp_58 + 1.0*tmp_127*tmp_128*tmp_66 - tmp_127*tmp_64 - tmp_128*tmp_71) + tmp_143*(-tmp_139*tmp_52 - tmp_140*tmp_58 + 1.0*tmp_141*tmp_142*tmp_66 - tmp_141*tmp_64 - tmp_142*tmp_71) + tmp_73*(-tmp_45*tmp_52 - tmp_53*tmp_58 - tmp_59*tmp_64 + 1.0*tmp_59*tmp_66*tmp_67 - tmp_67*tmp_71) + tmp_87*(-tmp_52*tmp_83 - tmp_58*tmp_84 - tmp_64*tmp_85 + 1.0*tmp_66*tmp_85*tmp_86 - tmp_71*tmp_86);
      real_t a_0_1 = tmp_101*(-tmp_145*tmp_98 - tmp_146*tmp_97 - tmp_147*tmp_99 + 1.0*tmp_66*tmp_91*tmp_99 - tmp_71*tmp_91) + tmp_115*(1.0*tmp_105*tmp_113*tmp_66 - tmp_105*tmp_71 - tmp_111*tmp_146 - tmp_112*tmp_145 - tmp_113*tmp_147) + tmp_129*(1.0*tmp_119*tmp_127*tmp_66 - tmp_119*tmp_71 - tmp_125*tmp_146 - tmp_126*tmp_145 - tmp_127*tmp_147) + tmp_143*(1.0*tmp_133*tmp_141*tmp_66 - tmp_133*tmp_71 - tmp_139*tmp_146 - tmp_140*tmp_145 - tmp_141*tmp_147) + tmp_73*(-tmp_145*tmp_53 - tmp_146*tmp_45 - tmp_147*tmp_59 + 1.0*tmp_33*tmp_59*tmp_66 - tmp_33*tmp_71) + tmp_87*(-tmp_145*tmp_84 - tmp_146*tmp_83 - tmp_147*tmp_85 + 1.0*tmp_66*tmp_77*tmp_85 - tmp_71*tmp_77);
      real_t a_0_2 = tmp_101*(-tmp_148*tmp_98 - tmp_149*tmp_97 - tmp_150*tmp_99 + 1.0*tmp_66*tmp_93*tmp_99 - tmp_71*tmp_93) + tmp_115*(1.0*tmp_107*tmp_113*tmp_66 - tmp_107*tmp_71 - tmp_111*tmp_149 - tmp_112*tmp_148 - tmp_113*tmp_150) + tmp_129*(1.0*tmp_121*tmp_127*tmp_66 - tmp_121*tmp_71 - tmp_125*tmp_149 - tmp_126*tmp_148 - tmp_127*tmp_150) + tmp_143*(1.0*tmp_135*tmp_141*tmp_66 - tmp_135*tmp_71 - tmp_139*tmp_149 - tmp_140*tmp_148 - tmp_141*tmp_150) + tmp_73*(-tmp_148*tmp_53 - tmp_149*tmp_45 - tmp_150*tmp_59 + 1.0*tmp_38*tmp_59*tmp_66 - tmp_38*tmp_71) + tmp_87*(-tmp_148*tmp_84 - tmp_149*tmp_83 - tmp_150*tmp_85 + 1.0*tmp_66*tmp_79*tmp_85 - tmp_71*tmp_79);
      real_t a_0_3 = tmp_101*(-tmp_151*tmp_98 - tmp_152*tmp_97 - tmp_153*tmp_99 + 1.0*tmp_66*tmp_95*tmp_99 - tmp_71*tmp_95) + tmp_115*(1.0*tmp_109*tmp_113*tmp_66 - tmp_109*tmp_71 - tmp_111*tmp_152 - tmp_112*tmp_151 - tmp_113*tmp_153) + tmp_129*(1.0*tmp_123*tmp_127*tmp_66 - tmp_123*tmp_71 - tmp_125*tmp_152 - tmp_126*tmp_151 - tmp_127*tmp_153) + tmp_143*(1.0*tmp_137*tmp_141*tmp_66 - tmp_137*tmp_71 - tmp_139*tmp_152 - tmp_140*tmp_151 - tmp_141*tmp_153) + tmp_73*(-tmp_151*tmp_53 - tmp_152*tmp_45 - tmp_153*tmp_59 + 1.0*tmp_43*tmp_59*tmp_66 - tmp_43*tmp_71) + tmp_87*(-tmp_151*tmp_84 - tmp_152*tmp_83 - tmp_153*tmp_85 + 1.0*tmp_66*tmp_81*tmp_85 - tmp_71*tmp_81);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_0*tmp_10;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = tmp_4*tmp_9;
      real_t tmp_15 = 1.0 / (-tmp_1*tmp_13 + tmp_1*tmp_2*tmp_9 + tmp_11*tmp_3 - tmp_12*tmp_6 - tmp_14*tmp_3 + tmp_6*tmp_8);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = 0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18;
      real_t tmp_20 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_21 = tmp_15*(tmp_19 + tmp_20);
      real_t tmp_22 = -tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_23 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_24 = -tmp_23;
      real_t tmp_25 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_26 = 0.091576213509770743*tmp_24 + 0.81684757298045851*tmp_25;
      real_t tmp_27 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_28 = tmp_15*(tmp_26 + tmp_27);
      real_t tmp_29 = -tmp_12 + tmp_8;
      real_t tmp_30 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_31 = -tmp_30;
      real_t tmp_32 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_33 = 0.091576213509770743*tmp_31 + 0.81684757298045851*tmp_32;
      real_t tmp_34 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_35 = tmp_15*(tmp_33 + tmp_34);
      real_t tmp_36 = tmp_21*tmp_5 + tmp_22*tmp_28 + tmp_29*tmp_35 - 1.0/4.0;
      real_t tmp_37 = tmp_0*tmp_3 - tmp_2*tmp_6;
      real_t tmp_38 = -tmp_3*tmp_9 + tmp_6*tmp_7;
      real_t tmp_39 = -tmp_13 + tmp_2*tmp_9;
      real_t tmp_40 = tmp_21*tmp_37 + tmp_28*tmp_38 + tmp_35*tmp_39 - 1.0/4.0;
      real_t tmp_41 = -tmp_0*tmp_1 + tmp_4*tmp_6;
      real_t tmp_42 = tmp_1*tmp_9 - tmp_10*tmp_6;
      real_t tmp_43 = tmp_11 - tmp_14;
      real_t tmp_44 = tmp_21*tmp_41 + tmp_28*tmp_42 + tmp_35*tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_36 + tmp_2*tmp_44 + tmp_4*tmp_40;
      real_t tmp_46 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_47 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_48 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_49 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_50 = -tmp_46*tmp_47 + tmp_48*tmp_49;
      real_t tmp_51 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_52 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_53 = tmp_46*tmp_52;
      real_t tmp_54 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_55 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_56 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_57 = tmp_47*tmp_56;
      real_t tmp_58 = tmp_46*tmp_54;
      real_t tmp_59 = tmp_48*tmp_56;
      real_t tmp_60 = tmp_49*tmp_55;
      real_t tmp_61 = 1.0 / (-tmp_47*tmp_58 + tmp_48*tmp_49*tmp_54 + tmp_51*tmp_53 - tmp_51*tmp_59 - tmp_52*tmp_60 + tmp_55*tmp_57);
      real_t tmp_62 = 0.5*tmp_61;
      real_t tmp_63 = tmp_50*tmp_62;
      real_t tmp_64 = tmp_46*tmp_51 - tmp_60;
      real_t tmp_65 = tmp_62*tmp_64;
      real_t tmp_66 = tmp_47*tmp_55 - tmp_48*tmp_51;
      real_t tmp_67 = tmp_62*tmp_66;
      real_t tmp_68 = -tmp_63 - tmp_65 - tmp_67;
      real_t tmp_69 = 0.5*p_affine_13_0;
      real_t tmp_70 = tmp_68*tmp_69;
      real_t tmp_71 = tmp_10*tmp_40 + tmp_36*tmp_9 + tmp_44*tmp_7;
      real_t tmp_72 = tmp_53 - tmp_59;
      real_t tmp_73 = tmp_62*tmp_72;
      real_t tmp_74 = tmp_55*tmp_56 - tmp_58;
      real_t tmp_75 = tmp_62*tmp_74;
      real_t tmp_76 = tmp_48*tmp_54 - tmp_52*tmp_55;
      real_t tmp_77 = tmp_62*tmp_76;
      real_t tmp_78 = -tmp_73 - tmp_75 - tmp_77;
      real_t tmp_79 = tmp_69*tmp_78;
      real_t tmp_80 = tmp_1*tmp_40 + tmp_3*tmp_44 + tmp_36*tmp_6;
      real_t tmp_81 = -tmp_49*tmp_52 + tmp_57;
      real_t tmp_82 = 1.0*tmp_61;
      real_t tmp_83 = tmp_81*tmp_82;
      real_t tmp_84 = tmp_49*tmp_54 - tmp_51*tmp_56;
      real_t tmp_85 = tmp_82*tmp_84;
      real_t tmp_86 = -tmp_47*tmp_54 + tmp_51*tmp_52;
      real_t tmp_87 = tmp_82*tmp_86;
      real_t tmp_88 = 0.5*p_affine_13_0*(-tmp_83 - tmp_85 - tmp_87) + 0.5*p_affine_13_1*tmp_68 + 0.5*p_affine_13_2*tmp_78;
      real_t tmp_89 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_90 = tmp_61*(tmp_19 + tmp_89);
      real_t tmp_91 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_92 = tmp_61*(tmp_26 + tmp_91);
      real_t tmp_93 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_94 = tmp_61*(tmp_33 + tmp_93);
      real_t tmp_95 = tmp_66*tmp_92 + tmp_76*tmp_90 + tmp_86*tmp_94;
      real_t tmp_96 = tmp_64*tmp_92 + tmp_74*tmp_90 + tmp_84*tmp_94;
      real_t tmp_97 = tmp_50*tmp_92 + tmp_72*tmp_90 + tmp_81*tmp_94;
      real_t tmp_98 = -tmp_95 - tmp_96 - tmp_97 + 1;
      real_t tmp_99 = (std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)*std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)) + (std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)*std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)) + (std::abs(tmp_23*tmp_32 - tmp_25*tmp_30)*std::abs(tmp_23*tmp_32 - tmp_25*tmp_30));
      real_t tmp_100 = 1.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_80;
      real_t tmp_102 = 1.0*tmp_15;
      real_t tmp_103 = 0.5*tmp_15;
      real_t tmp_104 = tmp_103*tmp_6;
      real_t tmp_105 = tmp_1*tmp_103;
      real_t tmp_106 = tmp_103*tmp_3;
      real_t tmp_107 = tmp_103*tmp_29;
      real_t tmp_108 = tmp_103*tmp_39;
      real_t tmp_109 = tmp_103*tmp_43;
      real_t tmp_110 = p_affine_13_0*(tmp_1*tmp_102*tmp_39 + tmp_102*tmp_29*tmp_6 + tmp_102*tmp_3*tmp_43) + p_affine_13_1*(tmp_0*tmp_107 + tmp_104*tmp_22 + tmp_105*tmp_38 + tmp_106*tmp_42 + tmp_108*tmp_4 + tmp_109*tmp_2) + p_affine_13_2*(tmp_10*tmp_108 + tmp_104*tmp_5 + tmp_105*tmp_37 + tmp_106*tmp_41 + tmp_107*tmp_9 + tmp_109*tmp_7);
      real_t tmp_111 = 2.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_112 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_111;
      real_t tmp_113 = 0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18;
      real_t tmp_114 = tmp_15*(tmp_113 + tmp_20);
      real_t tmp_115 = 0.44594849091596489*tmp_24 + 0.10810301816807022*tmp_25;
      real_t tmp_116 = tmp_15*(tmp_115 + tmp_27);
      real_t tmp_117 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32;
      real_t tmp_118 = tmp_15*(tmp_117 + tmp_34);
      real_t tmp_119 = tmp_114*tmp_5 + tmp_116*tmp_22 + tmp_118*tmp_29 - 1.0/4.0;
      real_t tmp_120 = tmp_114*tmp_37 + tmp_116*tmp_38 + tmp_118*tmp_39 - 1.0/4.0;
      real_t tmp_121 = tmp_114*tmp_41 + tmp_116*tmp_42 + tmp_118*tmp_43 - 1.0/4.0;
      real_t tmp_122 = tmp_0*tmp_119 + tmp_120*tmp_4 + tmp_121*tmp_2;
      real_t tmp_123 = tmp_10*tmp_120 + tmp_119*tmp_9 + tmp_121*tmp_7;
      real_t tmp_124 = tmp_1*tmp_120 + tmp_119*tmp_6 + tmp_121*tmp_3;
      real_t tmp_125 = tmp_61*(tmp_113 + tmp_89);
      real_t tmp_126 = tmp_61*(tmp_115 + tmp_91);
      real_t tmp_127 = tmp_61*(tmp_117 + tmp_93);
      real_t tmp_128 = tmp_125*tmp_76 + tmp_126*tmp_66 + tmp_127*tmp_86;
      real_t tmp_129 = tmp_125*tmp_74 + tmp_126*tmp_64 + tmp_127*tmp_84;
      real_t tmp_130 = tmp_125*tmp_72 + tmp_126*tmp_50 + tmp_127*tmp_81;
      real_t tmp_131 = -tmp_128 - tmp_129 - tmp_130 + 1;
      real_t tmp_132 = tmp_100*tmp_124;
      real_t tmp_133 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_111;
      real_t tmp_134 = 0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_135 = tmp_15*(tmp_134 + tmp_20);
      real_t tmp_136 = 0.81684757298045851*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_137 = tmp_15*(tmp_136 + tmp_27);
      real_t tmp_138 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_139 = tmp_15*(tmp_138 + tmp_34);
      real_t tmp_140 = tmp_135*tmp_5 + tmp_137*tmp_22 + tmp_139*tmp_29 - 1.0/4.0;
      real_t tmp_141 = tmp_135*tmp_37 + tmp_137*tmp_38 + tmp_139*tmp_39 - 1.0/4.0;
      real_t tmp_142 = tmp_135*tmp_41 + tmp_137*tmp_42 + tmp_139*tmp_43 - 1.0/4.0;
      real_t tmp_143 = tmp_0*tmp_140 + tmp_141*tmp_4 + tmp_142*tmp_2;
      real_t tmp_144 = tmp_10*tmp_141 + tmp_140*tmp_9 + tmp_142*tmp_7;
      real_t tmp_145 = tmp_1*tmp_141 + tmp_140*tmp_6 + tmp_142*tmp_3;
      real_t tmp_146 = tmp_61*(tmp_134 + tmp_89);
      real_t tmp_147 = tmp_61*(tmp_136 + tmp_91);
      real_t tmp_148 = tmp_61*(tmp_138 + tmp_93);
      real_t tmp_149 = tmp_146*tmp_76 + tmp_147*tmp_66 + tmp_148*tmp_86;
      real_t tmp_150 = tmp_146*tmp_74 + tmp_147*tmp_64 + tmp_148*tmp_84;
      real_t tmp_151 = tmp_146*tmp_72 + tmp_147*tmp_50 + tmp_148*tmp_81;
      real_t tmp_152 = -tmp_149 - tmp_150 - tmp_151 + 1;
      real_t tmp_153 = tmp_100*tmp_145;
      real_t tmp_154 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_111;
      real_t tmp_155 = 0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_156 = tmp_15*(tmp_155 + tmp_20);
      real_t tmp_157 = 0.10810301816807022*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_158 = tmp_15*(tmp_157 + tmp_27);
      real_t tmp_159 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_160 = tmp_15*(tmp_159 + tmp_34);
      real_t tmp_161 = tmp_156*tmp_5 + tmp_158*tmp_22 + tmp_160*tmp_29 - 1.0/4.0;
      real_t tmp_162 = tmp_156*tmp_37 + tmp_158*tmp_38 + tmp_160*tmp_39 - 1.0/4.0;
      real_t tmp_163 = tmp_156*tmp_41 + tmp_158*tmp_42 + tmp_160*tmp_43 - 1.0/4.0;
      real_t tmp_164 = tmp_0*tmp_161 + tmp_162*tmp_4 + tmp_163*tmp_2;
      real_t tmp_165 = tmp_10*tmp_162 + tmp_161*tmp_9 + tmp_163*tmp_7;
      real_t tmp_166 = tmp_1*tmp_162 + tmp_161*tmp_6 + tmp_163*tmp_3;
      real_t tmp_167 = tmp_61*(tmp_155 + tmp_89);
      real_t tmp_168 = tmp_61*(tmp_157 + tmp_91);
      real_t tmp_169 = tmp_61*(tmp_159 + tmp_93);
      real_t tmp_170 = tmp_167*tmp_76 + tmp_168*tmp_66 + tmp_169*tmp_86;
      real_t tmp_171 = tmp_167*tmp_74 + tmp_168*tmp_64 + tmp_169*tmp_84;
      real_t tmp_172 = tmp_167*tmp_72 + tmp_168*tmp_50 + tmp_169*tmp_81;
      real_t tmp_173 = -tmp_170 - tmp_171 - tmp_172 + 1;
      real_t tmp_174 = tmp_100*tmp_166;
      real_t tmp_175 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_111;
      real_t tmp_176 = 0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_177 = tmp_15*(tmp_176 + tmp_20);
      real_t tmp_178 = 0.091576213509770743*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_179 = tmp_15*(tmp_178 + tmp_27);
      real_t tmp_180 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_181 = tmp_15*(tmp_180 + tmp_34);
      real_t tmp_182 = tmp_177*tmp_5 + tmp_179*tmp_22 + tmp_181*tmp_29 - 1.0/4.0;
      real_t tmp_183 = tmp_177*tmp_37 + tmp_179*tmp_38 + tmp_181*tmp_39 - 1.0/4.0;
      real_t tmp_184 = tmp_177*tmp_41 + tmp_179*tmp_42 + tmp_181*tmp_43 - 1.0/4.0;
      real_t tmp_185 = tmp_0*tmp_182 + tmp_183*tmp_4 + tmp_184*tmp_2;
      real_t tmp_186 = tmp_10*tmp_183 + tmp_182*tmp_9 + tmp_184*tmp_7;
      real_t tmp_187 = tmp_1*tmp_183 + tmp_182*tmp_6 + tmp_184*tmp_3;
      real_t tmp_188 = tmp_61*(tmp_176 + tmp_89);
      real_t tmp_189 = tmp_61*(tmp_178 + tmp_91);
      real_t tmp_190 = tmp_61*(tmp_180 + tmp_93);
      real_t tmp_191 = tmp_188*tmp_76 + tmp_189*tmp_66 + tmp_190*tmp_86;
      real_t tmp_192 = tmp_188*tmp_74 + tmp_189*tmp_64 + tmp_190*tmp_84;
      real_t tmp_193 = tmp_188*tmp_72 + tmp_189*tmp_50 + tmp_190*tmp_81;
      real_t tmp_194 = -tmp_191 - tmp_192 - tmp_193 + 1;
      real_t tmp_195 = tmp_100*tmp_187;
      real_t tmp_196 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_111;
      real_t tmp_197 = 0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_198 = tmp_15*(tmp_197 + tmp_20);
      real_t tmp_199 = 0.44594849091596489*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_200 = tmp_15*(tmp_199 + tmp_27);
      real_t tmp_201 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_202 = tmp_15*(tmp_201 + tmp_34);
      real_t tmp_203 = tmp_198*tmp_5 + tmp_200*tmp_22 + tmp_202*tmp_29 - 1.0/4.0;
      real_t tmp_204 = tmp_198*tmp_37 + tmp_200*tmp_38 + tmp_202*tmp_39 - 1.0/4.0;
      real_t tmp_205 = tmp_198*tmp_41 + tmp_200*tmp_42 + tmp_202*tmp_43 - 1.0/4.0;
      real_t tmp_206 = tmp_0*tmp_203 + tmp_2*tmp_205 + tmp_204*tmp_4;
      real_t tmp_207 = tmp_10*tmp_204 + tmp_203*tmp_9 + tmp_205*tmp_7;
      real_t tmp_208 = tmp_1*tmp_204 + tmp_203*tmp_6 + tmp_205*tmp_3;
      real_t tmp_209 = tmp_61*(tmp_197 + tmp_89);
      real_t tmp_210 = tmp_61*(tmp_199 + tmp_91);
      real_t tmp_211 = tmp_61*(tmp_201 + tmp_93);
      real_t tmp_212 = tmp_209*tmp_76 + tmp_210*tmp_66 + tmp_211*tmp_86;
      real_t tmp_213 = tmp_209*tmp_74 + tmp_210*tmp_64 + tmp_211*tmp_84;
      real_t tmp_214 = tmp_209*tmp_72 + tmp_210*tmp_50 + tmp_211*tmp_81;
      real_t tmp_215 = -tmp_212 - tmp_213 - tmp_214 + 1;
      real_t tmp_216 = tmp_100*tmp_208;
      real_t tmp_217 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_111;
      real_t tmp_218 = 0.25*p_affine_13_0*tmp_61;
      real_t tmp_219 = tmp_218*tmp_76;
      real_t tmp_220 = tmp_218*tmp_66;
      real_t tmp_221 = 0.5*p_affine_13_0*tmp_87 + 0.5*p_affine_13_1*tmp_67 + 0.5*p_affine_13_2*tmp_77;
      real_t tmp_222 = tmp_218*tmp_74;
      real_t tmp_223 = tmp_218*tmp_64;
      real_t tmp_224 = 0.5*p_affine_13_0*tmp_85 + 0.5*p_affine_13_1*tmp_65 + 0.5*p_affine_13_2*tmp_75;
      real_t tmp_225 = tmp_218*tmp_72;
      real_t tmp_226 = tmp_218*tmp_50;
      real_t tmp_227 = 0.5*p_affine_13_0*tmp_83 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_73;
      real_t a_0_0 = tmp_112*(-tmp_101*tmp_98 + 0.5*tmp_110*tmp_98 - tmp_45*tmp_70 - tmp_71*tmp_79 - tmp_80*tmp_88) + tmp_133*(0.5*tmp_110*tmp_131 - tmp_122*tmp_70 - tmp_123*tmp_79 - tmp_124*tmp_88 - tmp_131*tmp_132) + tmp_154*(0.5*tmp_110*tmp_152 - tmp_143*tmp_70 - tmp_144*tmp_79 - tmp_145*tmp_88 - tmp_152*tmp_153) + tmp_175*(0.5*tmp_110*tmp_173 - tmp_164*tmp_70 - tmp_165*tmp_79 - tmp_166*tmp_88 - tmp_173*tmp_174) + tmp_196*(0.5*tmp_110*tmp_194 - tmp_185*tmp_70 - tmp_186*tmp_79 - tmp_187*tmp_88 - tmp_194*tmp_195) + tmp_217*(0.5*tmp_110*tmp_215 - tmp_206*tmp_70 - tmp_207*tmp_79 - tmp_208*tmp_88 - tmp_215*tmp_216);
      real_t a_0_1 = tmp_112*(-tmp_101*tmp_95 + 0.5*tmp_110*tmp_95 - tmp_219*tmp_71 - tmp_220*tmp_45 - tmp_221*tmp_80) + tmp_133*(0.5*tmp_110*tmp_128 - tmp_122*tmp_220 - tmp_123*tmp_219 - tmp_124*tmp_221 - tmp_128*tmp_132) + tmp_154*(0.5*tmp_110*tmp_149 - tmp_143*tmp_220 - tmp_144*tmp_219 - tmp_145*tmp_221 - tmp_149*tmp_153) + tmp_175*(0.5*tmp_110*tmp_170 - tmp_164*tmp_220 - tmp_165*tmp_219 - tmp_166*tmp_221 - tmp_170*tmp_174) + tmp_196*(0.5*tmp_110*tmp_191 - tmp_185*tmp_220 - tmp_186*tmp_219 - tmp_187*tmp_221 - tmp_191*tmp_195) + tmp_217*(0.5*tmp_110*tmp_212 - tmp_206*tmp_220 - tmp_207*tmp_219 - tmp_208*tmp_221 - tmp_212*tmp_216);
      real_t a_0_2 = tmp_112*(-tmp_101*tmp_96 + 0.5*tmp_110*tmp_96 - tmp_222*tmp_71 - tmp_223*tmp_45 - tmp_224*tmp_80) + tmp_133*(0.5*tmp_110*tmp_129 - tmp_122*tmp_223 - tmp_123*tmp_222 - tmp_124*tmp_224 - tmp_129*tmp_132) + tmp_154*(0.5*tmp_110*tmp_150 - tmp_143*tmp_223 - tmp_144*tmp_222 - tmp_145*tmp_224 - tmp_150*tmp_153) + tmp_175*(0.5*tmp_110*tmp_171 - tmp_164*tmp_223 - tmp_165*tmp_222 - tmp_166*tmp_224 - tmp_171*tmp_174) + tmp_196*(0.5*tmp_110*tmp_192 - tmp_185*tmp_223 - tmp_186*tmp_222 - tmp_187*tmp_224 - tmp_192*tmp_195) + tmp_217*(0.5*tmp_110*tmp_213 - tmp_206*tmp_223 - tmp_207*tmp_222 - tmp_208*tmp_224 - tmp_213*tmp_216);
      real_t a_0_3 = tmp_112*(-tmp_101*tmp_97 + 0.5*tmp_110*tmp_97 - tmp_225*tmp_71 - tmp_226*tmp_45 - tmp_227*tmp_80) + tmp_133*(0.5*tmp_110*tmp_130 - tmp_122*tmp_226 - tmp_123*tmp_225 - tmp_124*tmp_227 - tmp_130*tmp_132) + tmp_154*(0.5*tmp_110*tmp_151 - tmp_143*tmp_226 - tmp_144*tmp_225 - tmp_145*tmp_227 - tmp_151*tmp_153) + tmp_175*(0.5*tmp_110*tmp_172 - tmp_164*tmp_226 - tmp_165*tmp_225 - tmp_166*tmp_227 - tmp_172*tmp_174) + tmp_196*(0.5*tmp_110*tmp_193 - tmp_185*tmp_226 - tmp_186*tmp_225 - tmp_187*tmp_227 - tmp_193*tmp_195) + tmp_217*(0.5*tmp_110*tmp_214 - tmp_206*tmp_226 - tmp_207*tmp_225 - tmp_208*tmp_227 - tmp_214*tmp_216);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_0*tmp_10;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = tmp_4*tmp_9;
      real_t tmp_15 = 1.0 / (-tmp_1*tmp_13 + tmp_1*tmp_2*tmp_9 + tmp_11*tmp_3 - tmp_12*tmp_6 - tmp_14*tmp_3 + tmp_6*tmp_8);
      real_t tmp_16 = -p_affine_8_2 + p_affine_9_2;
      real_t tmp_17 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_18 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_19 = tmp_15*(0.091576213509770743*tmp_16 + 0.81684757298045851*tmp_17 + tmp_18);
      real_t tmp_20 = -tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_21 = -p_affine_8_1 + p_affine_9_1;
      real_t tmp_22 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_23 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_24 = tmp_15*(0.091576213509770743*tmp_21 + 0.81684757298045851*tmp_22 + tmp_23);
      real_t tmp_25 = -tmp_12 + tmp_8;
      real_t tmp_26 = -p_affine_8_0 + p_affine_9_0;
      real_t tmp_27 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_28 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_29 = tmp_15*(0.091576213509770743*tmp_26 + 0.81684757298045851*tmp_27 + tmp_28);
      real_t tmp_30 = tmp_19*tmp_5 + tmp_20*tmp_24 + tmp_25*tmp_29 - 1.0/4.0;
      real_t tmp_31 = tmp_0*tmp_3 - tmp_2*tmp_6;
      real_t tmp_32 = -tmp_3*tmp_9 + tmp_6*tmp_7;
      real_t tmp_33 = -tmp_13 + tmp_2*tmp_9;
      real_t tmp_34 = tmp_19*tmp_31 + tmp_24*tmp_32 + tmp_29*tmp_33 - 1.0/4.0;
      real_t tmp_35 = -tmp_0*tmp_1 + tmp_4*tmp_6;
      real_t tmp_36 = tmp_1*tmp_9 - tmp_10*tmp_6;
      real_t tmp_37 = tmp_11 - tmp_14;
      real_t tmp_38 = tmp_19*tmp_35 + tmp_24*tmp_36 + tmp_29*tmp_37 - 1.0/4.0;
      real_t tmp_39 = tmp_0*tmp_30 + tmp_2*tmp_38 + tmp_34*tmp_4;
      real_t tmp_40 = 0.5*tmp_15;
      real_t tmp_41 = tmp_36*tmp_40;
      real_t tmp_42 = tmp_32*tmp_40;
      real_t tmp_43 = tmp_20*tmp_40;
      real_t tmp_44 = -tmp_41 - tmp_42 - tmp_43;
      real_t tmp_45 = p_affine_13_0*tmp_44;
      real_t tmp_46 = tmp_10*tmp_34 + tmp_30*tmp_9 + tmp_38*tmp_7;
      real_t tmp_47 = tmp_35*tmp_40;
      real_t tmp_48 = tmp_31*tmp_40;
      real_t tmp_49 = tmp_40*tmp_5;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = p_affine_13_0*tmp_50;
      real_t tmp_52 = 1.0*tmp_15;
      real_t tmp_53 = tmp_37*tmp_52;
      real_t tmp_54 = tmp_33*tmp_52;
      real_t tmp_55 = tmp_25*tmp_52;
      real_t tmp_56 = p_affine_13_0*(-tmp_53 - tmp_54 - tmp_55) + p_affine_13_1*tmp_44 + p_affine_13_2*tmp_50;
      real_t tmp_57 = tmp_1*tmp_34 + tmp_3*tmp_38 + tmp_30*tmp_6;
      real_t tmp_58 = tmp_15*(0.44594849091596489*tmp_16 + 0.10810301816807022*tmp_17 + tmp_18);
      real_t tmp_59 = tmp_15*(0.44594849091596489*tmp_21 + 0.10810301816807022*tmp_22 + tmp_23);
      real_t tmp_60 = tmp_15*(0.44594849091596489*tmp_26 + 0.10810301816807022*tmp_27 + tmp_28);
      real_t tmp_61 = tmp_20*tmp_59 + tmp_25*tmp_60 + tmp_5*tmp_58 - 1.0/4.0;
      real_t tmp_62 = tmp_31*tmp_58 + tmp_32*tmp_59 + tmp_33*tmp_60 - 1.0/4.0;
      real_t tmp_63 = tmp_35*tmp_58 + tmp_36*tmp_59 + tmp_37*tmp_60 - 1.0/4.0;
      real_t tmp_64 = tmp_0*tmp_61 + tmp_2*tmp_63 + tmp_4*tmp_62;
      real_t tmp_65 = tmp_10*tmp_62 + tmp_61*tmp_9 + tmp_63*tmp_7;
      real_t tmp_66 = tmp_1*tmp_62 + tmp_3*tmp_63 + tmp_6*tmp_61;
      real_t tmp_67 = tmp_15*(0.81684757298045851*tmp_16 + 0.091576213509770743*tmp_17 + tmp_18);
      real_t tmp_68 = tmp_15*(0.81684757298045851*tmp_21 + 0.091576213509770743*tmp_22 + tmp_23);
      real_t tmp_69 = tmp_15*(0.81684757298045851*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28);
      real_t tmp_70 = tmp_20*tmp_68 + tmp_25*tmp_69 + tmp_5*tmp_67 - 1.0/4.0;
      real_t tmp_71 = tmp_31*tmp_67 + tmp_32*tmp_68 + tmp_33*tmp_69 - 1.0/4.0;
      real_t tmp_72 = tmp_35*tmp_67 + tmp_36*tmp_68 + tmp_37*tmp_69 - 1.0/4.0;
      real_t tmp_73 = tmp_0*tmp_70 + tmp_2*tmp_72 + tmp_4*tmp_71;
      real_t tmp_74 = tmp_10*tmp_71 + tmp_7*tmp_72 + tmp_70*tmp_9;
      real_t tmp_75 = tmp_1*tmp_71 + tmp_3*tmp_72 + tmp_6*tmp_70;
      real_t tmp_76 = tmp_15*(0.10810301816807022*tmp_16 + 0.44594849091596489*tmp_17 + tmp_18);
      real_t tmp_77 = tmp_15*(0.10810301816807022*tmp_21 + 0.44594849091596489*tmp_22 + tmp_23);
      real_t tmp_78 = tmp_15*(0.10810301816807022*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28);
      real_t tmp_79 = tmp_20*tmp_77 + tmp_25*tmp_78 + tmp_5*tmp_76 - 1.0/4.0;
      real_t tmp_80 = tmp_31*tmp_76 + tmp_32*tmp_77 + tmp_33*tmp_78 - 1.0/4.0;
      real_t tmp_81 = tmp_35*tmp_76 + tmp_36*tmp_77 + tmp_37*tmp_78 - 1.0/4.0;
      real_t tmp_82 = tmp_0*tmp_79 + tmp_2*tmp_81 + tmp_4*tmp_80;
      real_t tmp_83 = tmp_10*tmp_80 + tmp_7*tmp_81 + tmp_79*tmp_9;
      real_t tmp_84 = tmp_1*tmp_80 + tmp_3*tmp_81 + tmp_6*tmp_79;
      real_t tmp_85 = tmp_15*(0.091576213509770743*tmp_16 + 0.091576213509770743*tmp_17 + tmp_18);
      real_t tmp_86 = tmp_15*(0.091576213509770743*tmp_21 + 0.091576213509770743*tmp_22 + tmp_23);
      real_t tmp_87 = tmp_15*(0.091576213509770743*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28);
      real_t tmp_88 = tmp_20*tmp_86 + tmp_25*tmp_87 + tmp_5*tmp_85 - 1.0/4.0;
      real_t tmp_89 = tmp_31*tmp_85 + tmp_32*tmp_86 + tmp_33*tmp_87 - 1.0/4.0;
      real_t tmp_90 = tmp_35*tmp_85 + tmp_36*tmp_86 + tmp_37*tmp_87 - 1.0/4.0;
      real_t tmp_91 = tmp_0*tmp_88 + tmp_2*tmp_90 + tmp_4*tmp_89;
      real_t tmp_92 = tmp_10*tmp_89 + tmp_7*tmp_90 + tmp_88*tmp_9;
      real_t tmp_93 = tmp_1*tmp_89 + tmp_3*tmp_90 + tmp_6*tmp_88;
      real_t tmp_94 = tmp_15*(0.44594849091596489*tmp_16 + 0.44594849091596489*tmp_17 + tmp_18);
      real_t tmp_95 = tmp_15*(0.44594849091596489*tmp_21 + 0.44594849091596489*tmp_22 + tmp_23);
      real_t tmp_96 = tmp_15*(0.44594849091596489*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28);
      real_t tmp_97 = tmp_20*tmp_95 + tmp_25*tmp_96 + tmp_5*tmp_94 - 1.0/4.0;
      real_t tmp_98 = tmp_31*tmp_94 + tmp_32*tmp_95 + tmp_33*tmp_96 - 1.0/4.0;
      real_t tmp_99 = tmp_35*tmp_94 + tmp_36*tmp_95 + tmp_37*tmp_96 - 1.0/4.0;
      real_t tmp_100 = tmp_0*tmp_97 + tmp_2*tmp_99 + tmp_4*tmp_98;
      real_t tmp_101 = tmp_10*tmp_98 + tmp_7*tmp_99 + tmp_9*tmp_97;
      real_t tmp_102 = tmp_1*tmp_98 + tmp_3*tmp_99 + tmp_6*tmp_97;
      real_t tmp_103 = p_affine_13_0*tmp_49;
      real_t tmp_104 = p_affine_13_0*tmp_43;
      real_t tmp_105 = p_affine_13_0*tmp_55 + p_affine_13_1*tmp_43 + p_affine_13_2*tmp_49;
      real_t tmp_106 = p_affine_13_0*tmp_48;
      real_t tmp_107 = p_affine_13_0*tmp_42;
      real_t tmp_108 = p_affine_13_0*tmp_54 + p_affine_13_1*tmp_42 + p_affine_13_2*tmp_48;
      real_t tmp_109 = p_affine_13_0*tmp_47;
      real_t tmp_110 = p_affine_13_0*tmp_41;
      real_t tmp_111 = p_affine_13_0*tmp_53 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_47;
      real_t a_0_0 = -0.11169079483900572*tmp_100*tmp_45 - 0.11169079483900572*tmp_101*tmp_51 - 0.11169079483900572*tmp_102*tmp_56 - 0.054975871827660928*tmp_39*tmp_45 - 0.11169079483900572*tmp_45*tmp_64 - 0.054975871827660928*tmp_45*tmp_73 - 0.11169079483900572*tmp_45*tmp_82 - 0.054975871827660928*tmp_45*tmp_91 - 0.054975871827660928*tmp_46*tmp_51 - 0.11169079483900572*tmp_51*tmp_65 - 0.054975871827660928*tmp_51*tmp_74 - 0.11169079483900572*tmp_51*tmp_83 - 0.054975871827660928*tmp_51*tmp_92 - 0.054975871827660928*tmp_56*tmp_57 - 0.11169079483900572*tmp_56*tmp_66 - 0.054975871827660928*tmp_56*tmp_75 - 0.11169079483900572*tmp_56*tmp_84 - 0.054975871827660928*tmp_56*tmp_93;
      real_t a_0_1 = -0.11169079483900572*tmp_100*tmp_104 - 0.11169079483900572*tmp_101*tmp_103 - 0.11169079483900572*tmp_102*tmp_105 - 0.054975871827660928*tmp_103*tmp_46 - 0.11169079483900572*tmp_103*tmp_65 - 0.054975871827660928*tmp_103*tmp_74 - 0.11169079483900572*tmp_103*tmp_83 - 0.054975871827660928*tmp_103*tmp_92 - 0.054975871827660928*tmp_104*tmp_39 - 0.11169079483900572*tmp_104*tmp_64 - 0.054975871827660928*tmp_104*tmp_73 - 0.11169079483900572*tmp_104*tmp_82 - 0.054975871827660928*tmp_104*tmp_91 - 0.054975871827660928*tmp_105*tmp_57 - 0.11169079483900572*tmp_105*tmp_66 - 0.054975871827660928*tmp_105*tmp_75 - 0.11169079483900572*tmp_105*tmp_84 - 0.054975871827660928*tmp_105*tmp_93;
      real_t a_0_2 = -0.11169079483900572*tmp_100*tmp_107 - 0.11169079483900572*tmp_101*tmp_106 - 0.11169079483900572*tmp_102*tmp_108 - 0.054975871827660928*tmp_106*tmp_46 - 0.11169079483900572*tmp_106*tmp_65 - 0.054975871827660928*tmp_106*tmp_74 - 0.11169079483900572*tmp_106*tmp_83 - 0.054975871827660928*tmp_106*tmp_92 - 0.054975871827660928*tmp_107*tmp_39 - 0.11169079483900572*tmp_107*tmp_64 - 0.054975871827660928*tmp_107*tmp_73 - 0.11169079483900572*tmp_107*tmp_82 - 0.054975871827660928*tmp_107*tmp_91 - 0.054975871827660928*tmp_108*tmp_57 - 0.11169079483900572*tmp_108*tmp_66 - 0.054975871827660928*tmp_108*tmp_75 - 0.11169079483900572*tmp_108*tmp_84 - 0.054975871827660928*tmp_108*tmp_93;
      real_t a_0_3 = -0.11169079483900572*tmp_100*tmp_110 - 0.11169079483900572*tmp_101*tmp_109 - 0.11169079483900572*tmp_102*tmp_111 - 0.054975871827660928*tmp_109*tmp_46 - 0.11169079483900572*tmp_109*tmp_65 - 0.054975871827660928*tmp_109*tmp_74 - 0.11169079483900572*tmp_109*tmp_83 - 0.054975871827660928*tmp_109*tmp_92 - 0.054975871827660928*tmp_110*tmp_39 - 0.11169079483900572*tmp_110*tmp_64 - 0.054975871827660928*tmp_110*tmp_73 - 0.11169079483900572*tmp_110*tmp_82 - 0.054975871827660928*tmp_110*tmp_91 - 0.054975871827660928*tmp_111*tmp_57 - 0.11169079483900572*tmp_111*tmp_66 - 0.054975871827660928*tmp_111*tmp_75 - 0.11169079483900572*tmp_111*tmp_84 - 0.054975871827660928*tmp_111*tmp_93;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }

};




class EGEpsilonFormP1E_0 : public hyteg::dg::DGForm
{

 public:
    EGEpsilonFormP1E_0(std::function< real_t ( const Point3D & ) > mu)
: callback_Scalar_Variable_Coefficient_3D_mu (mu)
, callback_Scalar_Variable_Coefficient_2D_mu (mu)
    {}

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_mu;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;

void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
}

      
 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = -tmp_3;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_2 - tmp_4*tmp_6);
      real_t tmp_8 = 1.0*tmp_7;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = tmp_2*tmp_8 + tmp_4*tmp_9;
      real_t tmp_11 = Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_10;
      real_t tmp_12 = -2*tmp_1*tmp_8 - 2*tmp_9;
      real_t tmp_13 = 0.5*tmp_7;
      real_t tmp_14 = tmp_0*tmp_13;
      real_t tmp_15 = tmp_13*tmp_3;
      real_t tmp_16 = tmp_1*tmp_13;
      real_t tmp_17 = tmp_0*tmp_15 + tmp_14*tmp_4 + tmp_16*tmp_5 + tmp_16*tmp_6;
      real_t tmp_18 = Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_17;
      real_t tmp_19 = -4*tmp_14 - 4*tmp_15;
      real_t tmp_20 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_21 = 0.054975871827660928*tmp_20;
      real_t tmp_22 = tmp_10*tmp_12;
      real_t tmp_23 = tmp_17*tmp_19;
      real_t tmp_24 = 0.11169079483900572*tmp_20;
      real_t tmp_25 = 0.054975871827660928*tmp_20;
      real_t tmp_26 = 0.11169079483900572*tmp_20;
      real_t tmp_27 = 0.054975871827660928*tmp_20;
      real_t tmp_28 = 0.11169079483900572*tmp_20;
      real_t tmp_29 = 2.0*tmp_7;
      real_t tmp_30 = tmp_1*tmp_29;
      real_t tmp_31 = tmp_29*tmp_3;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_17*tmp_31;
      real_t tmp_34 = tmp_29*tmp_5;
      real_t tmp_35 = tmp_0*tmp_29;
      real_t tmp_36 = tmp_10*tmp_34;
      real_t tmp_37 = tmp_17*tmp_35;
      real_t a_0_0 = tmp_21*(tmp_11*tmp_12 + tmp_18*tmp_19) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_23) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_23) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_23) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_23) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_23);
      real_t a_1_0 = tmp_21*(tmp_11*tmp_30 + tmp_18*tmp_31) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_33) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_33) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_33) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_33) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_33);
      real_t a_2_0 = tmp_21*(tmp_11*tmp_34 + tmp_18*tmp_35) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_37) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_37) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_37) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_37) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_37);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -tmp_2;
      real_t tmp_7 = 1.0 / (-tmp_1*tmp_6 + tmp_5);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_10 = tmp_7*(0.21132486540518713*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = tmp_7*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_10*tmp_2 + tmp_13*tmp_4;
      real_t tmp_15 = tmp_14 - 1.0/3.0;
      real_t tmp_16 = tmp_0*tmp_13 + tmp_10*tmp_3;
      real_t tmp_17 = tmp_16 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_0*(tmp_1*tmp_15 + tmp_17*tmp_4);
      real_t tmp_19 = 0.5*tmp_7;
      real_t tmp_20 = tmp_19*tmp_3;
      real_t tmp_21 = tmp_19*tmp_2;
      real_t tmp_22 = -tmp_20 - tmp_21;
      real_t tmp_23 = 0.5*tmp_22;
      real_t tmp_24 = tmp_15*tmp_3 + tmp_17*tmp_6;
      real_t tmp_25 = 1.0*tmp_7;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_0*tmp_25;
      real_t tmp_28 = 0.5*p_affine_10_0*(-tmp_26 - tmp_27) + 0.5*p_affine_10_1*tmp_22;
      real_t tmp_29 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_30 = 1.0 / (tmp_29);
      real_t tmp_31 = -tmp_14 - tmp_16 + 1;
      real_t tmp_32 = tmp_19*tmp_4;
      real_t tmp_33 = 0.5*p_affine_10_0*(tmp_25*tmp_5 + tmp_27*tmp_6) + 0.5*p_affine_10_1*(tmp_0*tmp_32 + tmp_1*tmp_32 + tmp_20*tmp_6 + tmp_21*tmp_3);
      real_t tmp_34 = 2*tmp_29;
      real_t tmp_35 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_34;
      real_t tmp_36 = tmp_7*(0.78867513459481287*tmp_8 + tmp_9);
      real_t tmp_37 = tmp_7*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_38 = tmp_2*tmp_36 + tmp_37*tmp_4;
      real_t tmp_39 = tmp_38 - 1.0/3.0;
      real_t tmp_40 = tmp_0*tmp_37 + tmp_3*tmp_36;
      real_t tmp_41 = tmp_40 - 1.0/3.0;
      real_t tmp_42 = p_affine_10_0*(tmp_1*tmp_39 + tmp_4*tmp_41);
      real_t tmp_43 = tmp_3*tmp_39 + tmp_41*tmp_6;
      real_t tmp_44 = -tmp_38 - tmp_40 + 1;
      real_t tmp_45 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_34;
      real_t tmp_46 = 0.25*tmp_7;
      real_t tmp_47 = tmp_2*tmp_46;
      real_t tmp_48 = 0.5*p_affine_10_0*tmp_26 + 0.5*p_affine_10_1*tmp_21;
      real_t tmp_49 = tmp_3*tmp_46;
      real_t tmp_50 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*tmp_20;
      real_t a_0_0 = tmp_35*(-tmp_18*tmp_23 - tmp_24*tmp_28 + tmp_24*tmp_30*tmp_31 - tmp_31*tmp_33) + tmp_45*(-tmp_23*tmp_42 - tmp_28*tmp_43 + tmp_30*tmp_43*tmp_44 - tmp_33*tmp_44);
      real_t a_1_0 = tmp_35*(tmp_14*tmp_24*tmp_30 - tmp_14*tmp_33 - tmp_18*tmp_47 - tmp_24*tmp_48) + tmp_45*(tmp_30*tmp_38*tmp_43 - tmp_33*tmp_38 - tmp_42*tmp_47 - tmp_43*tmp_48);
      real_t a_2_0 = tmp_35*(tmp_16*tmp_24*tmp_30 - tmp_16*tmp_33 - tmp_18*tmp_49 - tmp_24*tmp_50) + tmp_45*(tmp_30*tmp_40*tmp_43 - tmp_33*tmp_40 - tmp_42*tmp_49 - tmp_43*tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElementInner[0]( 0 );
      const auto p_affine_0_1 = coordsElementInner[0]( 1 );

      const auto p_affine_1_0 = coordsElementInner[1]( 0 );
      const auto p_affine_1_1 = coordsElementInner[1]( 1 );

      const auto p_affine_2_0 = coordsElementInner[2]( 0 );
      const auto p_affine_2_1 = coordsElementInner[2]( 1 );

      const auto p_affine_3_0 = coordsElementOuter[0]( 0 );
      const auto p_affine_3_1 = coordsElementOuter[0]( 1 );

      const auto p_affine_4_0 = coordsElementOuter[1]( 0 );
      const auto p_affine_4_1 = coordsElementOuter[1]( 1 );

      const auto p_affine_5_0 = coordsElementOuter[2]( 0 );
      const auto p_affine_5_1 = coordsElementOuter[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_3 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_4 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -tmp_2;
      real_t tmp_7 = 1.0 / (-tmp_1*tmp_6 + tmp_5);
      real_t tmp_8 = -p_affine_3_1;
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.21132486540518713*tmp_9;
      real_t tmp_11 = tmp_7*(tmp_10 + tmp_8);
      real_t tmp_12 = -p_affine_3_0;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + 0.21132486540518713*tmp_13;
      real_t tmp_15 = tmp_7*(tmp_12 + tmp_14);
      real_t tmp_16 = tmp_11*tmp_2 + tmp_15*tmp_4 - 1.0/3.0;
      real_t tmp_17 = tmp_0*tmp_15 + tmp_11*tmp_3 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_0*(tmp_1*tmp_16 + tmp_17*tmp_4);
      real_t tmp_19 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_20 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_21 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_22 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_23 = 1.0 / (tmp_19*tmp_20 - tmp_21*tmp_22);
      real_t tmp_24 = 0.5*tmp_23;
      real_t tmp_25 = tmp_19*tmp_24;
      real_t tmp_26 = tmp_21*tmp_24;
      real_t tmp_27 = -tmp_25 - tmp_26;
      real_t tmp_28 = 0.5*tmp_27;
      real_t tmp_29 = tmp_16*tmp_3 + tmp_17*tmp_6;
      real_t tmp_30 = 1.0*tmp_23;
      real_t tmp_31 = tmp_20*tmp_30;
      real_t tmp_32 = tmp_22*tmp_30;
      real_t tmp_33 = 0.5*p_affine_10_0*(-tmp_31 - tmp_32) + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_34 = -p_affine_0_1;
      real_t tmp_35 = tmp_23*(tmp_10 + tmp_34);
      real_t tmp_36 = -p_affine_0_0;
      real_t tmp_37 = tmp_23*(tmp_14 + tmp_36);
      real_t tmp_38 = tmp_20*tmp_37 + tmp_21*tmp_35;
      real_t tmp_39 = tmp_19*tmp_35 + tmp_22*tmp_37;
      real_t tmp_40 = -tmp_38 - tmp_39 + 1;
      real_t tmp_41 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_42 = 1.0 / (tmp_41);
      real_t tmp_43 = tmp_29*tmp_42;
      real_t tmp_44 = 1.0*tmp_7;
      real_t tmp_45 = 0.5*tmp_7;
      real_t tmp_46 = tmp_3*tmp_45;
      real_t tmp_47 = tmp_4*tmp_45;
      real_t tmp_48 = 0.5*p_affine_10_0*(tmp_0*tmp_44*tmp_6 + tmp_44*tmp_5) + 0.5*p_affine_10_1*(tmp_0*tmp_47 + tmp_1*tmp_47 + tmp_2*tmp_46 + tmp_46*tmp_6);
      real_t tmp_49 = 2*tmp_41;
      real_t tmp_50 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_49;
      real_t tmp_51 = p_affine_6_1 + 0.78867513459481287*tmp_9;
      real_t tmp_52 = tmp_7*(tmp_51 + tmp_8);
      real_t tmp_53 = p_affine_6_0 + 0.78867513459481287*tmp_13;
      real_t tmp_54 = tmp_7*(tmp_12 + tmp_53);
      real_t tmp_55 = tmp_2*tmp_52 + tmp_4*tmp_54 - 1.0/3.0;
      real_t tmp_56 = tmp_0*tmp_54 + tmp_3*tmp_52 - 1.0/3.0;
      real_t tmp_57 = p_affine_10_0*(tmp_1*tmp_55 + tmp_4*tmp_56);
      real_t tmp_58 = tmp_3*tmp_55 + tmp_56*tmp_6;
      real_t tmp_59 = tmp_23*(tmp_34 + tmp_51);
      real_t tmp_60 = tmp_23*(tmp_36 + tmp_53);
      real_t tmp_61 = tmp_20*tmp_60 + tmp_21*tmp_59;
      real_t tmp_62 = tmp_19*tmp_59 + tmp_22*tmp_60;
      real_t tmp_63 = -tmp_61 - tmp_62 + 1;
      real_t tmp_64 = tmp_42*tmp_58;
      real_t tmp_65 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_49;
      real_t tmp_66 = 0.25*tmp_23;
      real_t tmp_67 = tmp_21*tmp_66;
      real_t tmp_68 = 0.5*p_affine_10_0*tmp_31 + 0.5*p_affine_10_1*tmp_26;
      real_t tmp_69 = tmp_19*tmp_66;
      real_t tmp_70 = 0.5*p_affine_10_0*tmp_32 + 0.5*p_affine_10_1*tmp_25;
      real_t a_0_0 = tmp_50*(tmp_18*tmp_28 + tmp_29*tmp_33 - tmp_40*tmp_43 - tmp_40*tmp_48) + tmp_65*(tmp_28*tmp_57 + tmp_33*tmp_58 - tmp_48*tmp_63 - tmp_63*tmp_64);
      real_t a_1_0 = tmp_50*(tmp_18*tmp_67 + tmp_29*tmp_68 - tmp_38*tmp_43 - tmp_38*tmp_48) + tmp_65*(-tmp_48*tmp_61 + tmp_57*tmp_67 + tmp_58*tmp_68 - tmp_61*tmp_64);
      real_t a_2_0 = tmp_50*(tmp_18*tmp_69 + tmp_29*tmp_70 - tmp_39*tmp_43 - tmp_39*tmp_48) + tmp_65*(-tmp_48*tmp_62 + tmp_57*tmp_69 + tmp_58*tmp_70 - tmp_62*tmp_64);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -tmp_2;
      real_t tmp_7 = 1.0 / (-tmp_1*tmp_6 + tmp_5);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_10 = tmp_7*(0.21132486540518713*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = tmp_7*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_10*tmp_2 + tmp_13*tmp_4;
      real_t tmp_15 = tmp_14 - 1.0/3.0;
      real_t tmp_16 = tmp_0*tmp_13 + tmp_10*tmp_3;
      real_t tmp_17 = tmp_16 - 1.0/3.0;
      real_t tmp_18 = tmp_1*tmp_15 + tmp_17*tmp_4;
      real_t tmp_19 = 0.5*tmp_7;
      real_t tmp_20 = tmp_19*tmp_3;
      real_t tmp_21 = tmp_19*tmp_2;
      real_t tmp_22 = -tmp_20 - tmp_21;
      real_t tmp_23 = p_affine_10_0*tmp_22;
      real_t tmp_24 = 1.0*tmp_7;
      real_t tmp_25 = tmp_24*tmp_4;
      real_t tmp_26 = tmp_0*tmp_24;
      real_t tmp_27 = p_affine_10_0*(-tmp_25 - tmp_26) + p_affine_10_1*tmp_22;
      real_t tmp_28 = tmp_15*tmp_3 + tmp_17*tmp_6;
      real_t tmp_29 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_30 = 1.0 / (tmp_29);
      real_t tmp_31 = -tmp_14 - tmp_16 + 1;
      real_t tmp_32 = tmp_19*tmp_4;
      real_t tmp_33 = p_affine_10_0*(tmp_24*tmp_5 + tmp_26*tmp_6) + p_affine_10_1*(tmp_0*tmp_32 + tmp_1*tmp_32 + tmp_20*tmp_6 + tmp_21*tmp_3);
      real_t tmp_34 = 2*tmp_29;
      real_t tmp_35 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_34;
      real_t tmp_36 = tmp_7*(0.78867513459481287*tmp_8 + tmp_9);
      real_t tmp_37 = tmp_7*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_38 = tmp_2*tmp_36 + tmp_37*tmp_4;
      real_t tmp_39 = tmp_38 - 1.0/3.0;
      real_t tmp_40 = tmp_0*tmp_37 + tmp_3*tmp_36;
      real_t tmp_41 = tmp_40 - 1.0/3.0;
      real_t tmp_42 = tmp_1*tmp_39 + tmp_4*tmp_41;
      real_t tmp_43 = tmp_3*tmp_39 + tmp_41*tmp_6;
      real_t tmp_44 = -tmp_38 - tmp_40 + 1;
      real_t tmp_45 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_34;
      real_t tmp_46 = p_affine_10_0*tmp_21;
      real_t tmp_47 = p_affine_10_0*tmp_25 + p_affine_10_1*tmp_21;
      real_t tmp_48 = p_affine_10_0*tmp_20;
      real_t tmp_49 = p_affine_10_0*tmp_26 + p_affine_10_1*tmp_20;
      real_t a_0_0 = tmp_35*(-tmp_18*tmp_23 - tmp_27*tmp_28 + tmp_28*tmp_30*tmp_31 - tmp_31*tmp_33) + tmp_45*(-tmp_23*tmp_42 - tmp_27*tmp_43 + tmp_30*tmp_43*tmp_44 - tmp_33*tmp_44);
      real_t a_1_0 = tmp_35*(tmp_14*tmp_28*tmp_30 - tmp_14*tmp_33 - tmp_18*tmp_46 - tmp_28*tmp_47) + tmp_45*(tmp_30*tmp_38*tmp_43 - tmp_33*tmp_38 - tmp_42*tmp_46 - tmp_43*tmp_47);
      real_t a_2_0 = tmp_35*(tmp_16*tmp_28*tmp_30 - tmp_16*tmp_33 - tmp_18*tmp_48 - tmp_28*tmp_49) + tmp_45*(tmp_30*tmp_40*tmp_43 - tmp_33*tmp_40 - tmp_42*tmp_48 - tmp_43*tmp_49);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

  void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

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

      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id6 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id7 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id8 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id9 );
      Scalar_Variable_Coefficient_3D_mu( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id10 );
      Scalar_Variable_Coefficient_3D_mu( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id11 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id12 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id13 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_5 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_9 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_10 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_11 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = tmp_11*tmp_2;
      real_t tmp_14 = tmp_1*tmp_9;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_3 - tmp_0*tmp_6 + tmp_10*tmp_12 - tmp_10*tmp_14 - tmp_13*tmp_8 + tmp_4*tmp_8*tmp_9);
      real_t tmp_16 = 1.0*tmp_15;
      real_t tmp_17 = tmp_16*tmp_7;
      real_t tmp_18 = -tmp_13 + tmp_4*tmp_9;
      real_t tmp_19 = tmp_16*tmp_18;
      real_t tmp_20 = tmp_12 - tmp_14;
      real_t tmp_21 = tmp_16*tmp_20;
      real_t tmp_22 = tmp_0*tmp_17 + tmp_10*tmp_21 + tmp_19*tmp_8;
      real_t tmp_23 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_22;
      real_t tmp_24 = -2*tmp_17 - 2*tmp_19 - 2*tmp_21;
      real_t tmp_25 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_26 = tmp_0*tmp_1 - tmp_11*tmp_8;
      real_t tmp_27 = 0.5*tmp_15;
      real_t tmp_28 = tmp_26*tmp_27;
      real_t tmp_29 = -tmp_0*tmp_4 + tmp_10*tmp_11;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = -tmp_1*tmp_10 + tmp_4*tmp_8;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_27*tmp_7;
      real_t tmp_34 = tmp_18*tmp_27;
      real_t tmp_35 = tmp_20*tmp_27;
      real_t tmp_36 = tmp_0*tmp_32 + tmp_10*tmp_28 + tmp_2*tmp_35 + tmp_30*tmp_8 + tmp_33*tmp_9 + tmp_34*tmp_5;
      real_t tmp_37 = tmp_36*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_38 = -tmp_0*tmp_5 + tmp_8*tmp_9;
      real_t tmp_39 = tmp_27*tmp_38;
      real_t tmp_40 = tmp_0*tmp_2 - tmp_10*tmp_9;
      real_t tmp_41 = tmp_27*tmp_40;
      real_t tmp_42 = tmp_10*tmp_5 - tmp_2*tmp_8;
      real_t tmp_43 = tmp_27*tmp_42;
      real_t tmp_44 = tmp_0*tmp_43 + tmp_1*tmp_34 + tmp_10*tmp_39 + tmp_11*tmp_33 + tmp_35*tmp_4 + tmp_41*tmp_8;
      real_t tmp_45 = tmp_44*(-tmp_39 - tmp_41 - tmp_43);
      real_t tmp_46 = p_affine_0_0*p_affine_1_1;
      real_t tmp_47 = p_affine_0_0*p_affine_1_2;
      real_t tmp_48 = p_affine_2_1*p_affine_3_2;
      real_t tmp_49 = p_affine_0_1*p_affine_1_0;
      real_t tmp_50 = p_affine_0_1*p_affine_1_2;
      real_t tmp_51 = p_affine_2_2*p_affine_3_0;
      real_t tmp_52 = p_affine_0_2*p_affine_1_0;
      real_t tmp_53 = p_affine_0_2*p_affine_1_1;
      real_t tmp_54 = p_affine_2_0*p_affine_3_1;
      real_t tmp_55 = p_affine_2_2*p_affine_3_1;
      real_t tmp_56 = p_affine_2_0*p_affine_3_2;
      real_t tmp_57 = p_affine_2_1*p_affine_3_0;
      real_t tmp_58 = std::abs(p_affine_0_0*tmp_48 - p_affine_0_0*tmp_55 + p_affine_0_1*tmp_51 - p_affine_0_1*tmp_56 + p_affine_0_2*tmp_54 - p_affine_0_2*tmp_57 - p_affine_1_0*tmp_48 + p_affine_1_0*tmp_55 - p_affine_1_1*tmp_51 + p_affine_1_1*tmp_56 - p_affine_1_2*tmp_54 + p_affine_1_2*tmp_57 + p_affine_2_0*tmp_50 - p_affine_2_0*tmp_53 - p_affine_2_1*tmp_47 + p_affine_2_1*tmp_52 + p_affine_2_2*tmp_46 - p_affine_2_2*tmp_49 - p_affine_3_0*tmp_50 + p_affine_3_0*tmp_53 + p_affine_3_1*tmp_47 - p_affine_3_1*tmp_52 - p_affine_3_2*tmp_46 + p_affine_3_2*tmp_49);
      real_t tmp_59 = 0.018781320953002646*tmp_58;
      real_t tmp_60 = tmp_22*tmp_24;
      real_t tmp_61 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id1;
      real_t tmp_62 = 0.012248840519393657*tmp_58;
      real_t tmp_63 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id2;
      real_t tmp_64 = 0.0070910034628469103*tmp_58;
      real_t tmp_65 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id3;
      real_t tmp_66 = 0.0070910034628469103*tmp_58;
      real_t tmp_67 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id4;
      real_t tmp_68 = 0.0070910034628469103*tmp_58;
      real_t tmp_69 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id5;
      real_t tmp_70 = 0.0070910034628469103*tmp_58;
      real_t tmp_71 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id6;
      real_t tmp_72 = 0.018781320953002646*tmp_58;
      real_t tmp_73 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id7;
      real_t tmp_74 = 0.012248840519393657*tmp_58;
      real_t tmp_75 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id8;
      real_t tmp_76 = 0.018781320953002646*tmp_58;
      real_t tmp_77 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id9;
      real_t tmp_78 = 0.012248840519393657*tmp_58;
      real_t tmp_79 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id10;
      real_t tmp_80 = 0.018781320953002646*tmp_58;
      real_t tmp_81 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id11;
      real_t tmp_82 = 0.012248840519393657*tmp_58;
      real_t tmp_83 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id12;
      real_t tmp_84 = 0.0070910034628469103*tmp_58;
      real_t tmp_85 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id13;
      real_t tmp_86 = 0.0070910034628469103*tmp_58;
      real_t tmp_87 = 2.0*tmp_15;
      real_t tmp_88 = tmp_7*tmp_87;
      real_t tmp_89 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_87;
      real_t tmp_90 = tmp_31*tmp_36;
      real_t tmp_91 = tmp_42*tmp_44;
      real_t tmp_92 = tmp_22*tmp_88;
      real_t tmp_93 = Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_87;
      real_t tmp_94 = Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_87;
      real_t tmp_95 = Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_87;
      real_t tmp_96 = Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_87;
      real_t tmp_97 = Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_87;
      real_t tmp_98 = Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_87;
      real_t tmp_99 = Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_87;
      real_t tmp_100 = Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_87;
      real_t tmp_101 = Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_87;
      real_t tmp_102 = Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_87;
      real_t tmp_103 = Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_87;
      real_t tmp_104 = Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_87;
      real_t tmp_105 = Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_87;
      real_t tmp_106 = tmp_23*tmp_87;
      real_t tmp_107 = tmp_29*tmp_36;
      real_t tmp_108 = tmp_40*tmp_44;
      real_t tmp_109 = tmp_18*tmp_22;
      real_t tmp_110 = tmp_26*tmp_36;
      real_t tmp_111 = tmp_38*tmp_44;
      real_t tmp_112 = tmp_20*tmp_22;
      real_t a_0_0 = tmp_59*(tmp_23*tmp_24 + tmp_25*tmp_37 + tmp_25*tmp_45) + tmp_62*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_60 + tmp_37*tmp_61 + tmp_45*tmp_61) + tmp_64*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_60 + tmp_37*tmp_63 + tmp_45*tmp_63) + tmp_66*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_60 + tmp_37*tmp_65 + tmp_45*tmp_65) + tmp_68*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_60 + tmp_37*tmp_67 + tmp_45*tmp_67) + tmp_70*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_60 + tmp_37*tmp_69 + tmp_45*tmp_69) + tmp_72*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_60 + tmp_37*tmp_71 + tmp_45*tmp_71) + tmp_74*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_60 + tmp_37*tmp_73 + tmp_45*tmp_73) + tmp_76*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_60 + tmp_37*tmp_75 + tmp_45*tmp_75) + tmp_78*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_60 + tmp_37*tmp_77 + tmp_45*tmp_77) + tmp_80*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_60 + tmp_37*tmp_79 + tmp_45*tmp_79) + tmp_82*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_60 + tmp_37*tmp_81 + tmp_45*tmp_81) + tmp_84*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_60 + tmp_37*tmp_83 + tmp_45*tmp_83) + tmp_86*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_60 + tmp_37*tmp_85 + tmp_45*tmp_85);
      real_t a_1_0 = tmp_59*(tmp_23*tmp_88 + tmp_89*tmp_90 + tmp_89*tmp_91) + tmp_62*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_92 + tmp_90*tmp_93 + tmp_91*tmp_93) + tmp_64*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_92 + tmp_90*tmp_94 + tmp_91*tmp_94) + tmp_66*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_92 + tmp_90*tmp_95 + tmp_91*tmp_95) + tmp_68*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_92 + tmp_90*tmp_96 + tmp_91*tmp_96) + tmp_70*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_92 + tmp_90*tmp_97 + tmp_91*tmp_97) + tmp_72*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_92 + tmp_90*tmp_98 + tmp_91*tmp_98) + tmp_74*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_92 + tmp_90*tmp_99 + tmp_91*tmp_99) + tmp_76*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_92 + tmp_100*tmp_90 + tmp_100*tmp_91) + tmp_78*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_92 + tmp_101*tmp_90 + tmp_101*tmp_91) + tmp_80*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_92 + tmp_102*tmp_90 + tmp_102*tmp_91) + tmp_82*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_92 + tmp_103*tmp_90 + tmp_103*tmp_91) + tmp_84*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_92 + tmp_104*tmp_90 + tmp_104*tmp_91) + tmp_86*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_92 + tmp_105*tmp_90 + tmp_105*tmp_91);
      real_t a_2_0 = tmp_59*(tmp_106*tmp_18 + tmp_107*tmp_89 + tmp_108*tmp_89) + tmp_62*(tmp_107*tmp_93 + tmp_108*tmp_93 + tmp_109*tmp_93) + tmp_64*(tmp_107*tmp_94 + tmp_108*tmp_94 + tmp_109*tmp_94) + tmp_66*(tmp_107*tmp_95 + tmp_108*tmp_95 + tmp_109*tmp_95) + tmp_68*(tmp_107*tmp_96 + tmp_108*tmp_96 + tmp_109*tmp_96) + tmp_70*(tmp_107*tmp_97 + tmp_108*tmp_97 + tmp_109*tmp_97) + tmp_72*(tmp_107*tmp_98 + tmp_108*tmp_98 + tmp_109*tmp_98) + tmp_74*(tmp_107*tmp_99 + tmp_108*tmp_99 + tmp_109*tmp_99) + tmp_76*(tmp_100*tmp_107 + tmp_100*tmp_108 + tmp_100*tmp_109) + tmp_78*(tmp_101*tmp_107 + tmp_101*tmp_108 + tmp_101*tmp_109) + tmp_80*(tmp_102*tmp_107 + tmp_102*tmp_108 + tmp_102*tmp_109) + tmp_82*(tmp_103*tmp_107 + tmp_103*tmp_108 + tmp_103*tmp_109) + tmp_84*(tmp_104*tmp_107 + tmp_104*tmp_108 + tmp_104*tmp_109) + tmp_86*(tmp_105*tmp_107 + tmp_105*tmp_108 + tmp_105*tmp_109);
      real_t a_3_0 = tmp_59*(tmp_106*tmp_20 + tmp_110*tmp_89 + tmp_111*tmp_89) + tmp_62*(tmp_110*tmp_93 + tmp_111*tmp_93 + tmp_112*tmp_93) + tmp_64*(tmp_110*tmp_94 + tmp_111*tmp_94 + tmp_112*tmp_94) + tmp_66*(tmp_110*tmp_95 + tmp_111*tmp_95 + tmp_112*tmp_95) + tmp_68*(tmp_110*tmp_96 + tmp_111*tmp_96 + tmp_112*tmp_96) + tmp_70*(tmp_110*tmp_97 + tmp_111*tmp_97 + tmp_112*tmp_97) + tmp_72*(tmp_110*tmp_98 + tmp_111*tmp_98 + tmp_112*tmp_98) + tmp_74*(tmp_110*tmp_99 + tmp_111*tmp_99 + tmp_112*tmp_99) + tmp_76*(tmp_100*tmp_110 + tmp_100*tmp_111 + tmp_100*tmp_112) + tmp_78*(tmp_101*tmp_110 + tmp_101*tmp_111 + tmp_101*tmp_112) + tmp_80*(tmp_102*tmp_110 + tmp_102*tmp_111 + tmp_102*tmp_112) + tmp_82*(tmp_103*tmp_110 + tmp_103*tmp_111 + tmp_103*tmp_112) + tmp_84*(tmp_104*tmp_110 + tmp_104*tmp_111 + tmp_104*tmp_112) + tmp_86*(tmp_105*tmp_110 + tmp_105*tmp_111 + tmp_105*tmp_112);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_0*tmp_10;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = tmp_4*tmp_9;
      real_t tmp_15 = 1.0 / (-tmp_1*tmp_13 + tmp_1*tmp_2*tmp_9 + tmp_11*tmp_3 - tmp_12*tmp_6 - tmp_14*tmp_3 + tmp_6*tmp_8);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_8;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = tmp_33 - 1.0/4.0;
      real_t tmp_35 = tmp_0*tmp_3 - tmp_2*tmp_6;
      real_t tmp_36 = -tmp_3*tmp_9 + tmp_6*tmp_7;
      real_t tmp_37 = -tmp_13 + tmp_2*tmp_9;
      real_t tmp_38 = tmp_20*tmp_35 + tmp_26*tmp_36 + tmp_32*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/4.0;
      real_t tmp_40 = -tmp_0*tmp_1 + tmp_4*tmp_6;
      real_t tmp_41 = tmp_1*tmp_9 - tmp_10*tmp_6;
      real_t tmp_42 = tmp_11 - tmp_14;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_26*tmp_41 + tmp_32*tmp_42;
      real_t tmp_44 = tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_34 + tmp_2*tmp_44 + tmp_39*tmp_4;
      real_t tmp_46 = 0.5*tmp_15;
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_36*tmp_46;
      real_t tmp_49 = tmp_21*tmp_46;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = 0.5*p_affine_13_0;
      real_t tmp_52 = tmp_50*tmp_51;
      real_t tmp_53 = tmp_10*tmp_39 + tmp_34*tmp_9 + tmp_44*tmp_7;
      real_t tmp_54 = tmp_40*tmp_46;
      real_t tmp_55 = tmp_35*tmp_46;
      real_t tmp_56 = tmp_46*tmp_5;
      real_t tmp_57 = -tmp_54 - tmp_55 - tmp_56;
      real_t tmp_58 = tmp_51*tmp_57;
      real_t tmp_59 = tmp_1*tmp_39 + tmp_3*tmp_44 + tmp_34*tmp_6;
      real_t tmp_60 = 1.0*tmp_15;
      real_t tmp_61 = tmp_42*tmp_60;
      real_t tmp_62 = tmp_37*tmp_60;
      real_t tmp_63 = tmp_27*tmp_60;
      real_t tmp_64 = 0.5*p_affine_13_0*(-tmp_61 - tmp_62 - tmp_63) + 0.5*p_affine_13_1*tmp_50 + 0.5*p_affine_13_2*tmp_57;
      real_t tmp_65 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_66 = std::pow(tmp_65, -0.25);
      real_t tmp_67 = -tmp_33 - tmp_38 - tmp_43 + 1;
      real_t tmp_68 = tmp_27*tmp_46;
      real_t tmp_69 = tmp_37*tmp_46;
      real_t tmp_70 = tmp_42*tmp_46;
      real_t tmp_71 = 0.5*p_affine_13_0*(tmp_1*tmp_62 + tmp_3*tmp_61 + tmp_6*tmp_63) + 0.5*p_affine_13_1*(tmp_0*tmp_68 + tmp_1*tmp_48 + tmp_2*tmp_70 + tmp_3*tmp_47 + tmp_4*tmp_69 + tmp_49*tmp_6) + 0.5*p_affine_13_2*(tmp_1*tmp_55 + tmp_10*tmp_69 + tmp_3*tmp_54 + tmp_56*tmp_6 + tmp_68*tmp_9 + tmp_7*tmp_70);
      real_t tmp_72 = 2.0*std::pow(tmp_65, 1.0/2.0);
      real_t tmp_73 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_72;
      real_t tmp_74 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_75 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_76 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_77 = tmp_21*tmp_75 + tmp_27*tmp_76 + tmp_5*tmp_74;
      real_t tmp_78 = tmp_77 - 1.0/4.0;
      real_t tmp_79 = tmp_35*tmp_74 + tmp_36*tmp_75 + tmp_37*tmp_76;
      real_t tmp_80 = tmp_79 - 1.0/4.0;
      real_t tmp_81 = tmp_40*tmp_74 + tmp_41*tmp_75 + tmp_42*tmp_76;
      real_t tmp_82 = tmp_81 - 1.0/4.0;
      real_t tmp_83 = tmp_0*tmp_78 + tmp_2*tmp_82 + tmp_4*tmp_80;
      real_t tmp_84 = tmp_10*tmp_80 + tmp_7*tmp_82 + tmp_78*tmp_9;
      real_t tmp_85 = tmp_1*tmp_80 + tmp_3*tmp_82 + tmp_6*tmp_78;
      real_t tmp_86 = -tmp_77 - tmp_79 - tmp_81 + 1;
      real_t tmp_87 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_72;
      real_t tmp_88 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_89 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_90 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_91 = tmp_21*tmp_89 + tmp_27*tmp_90 + tmp_5*tmp_88;
      real_t tmp_92 = tmp_91 - 1.0/4.0;
      real_t tmp_93 = tmp_35*tmp_88 + tmp_36*tmp_89 + tmp_37*tmp_90;
      real_t tmp_94 = tmp_93 - 1.0/4.0;
      real_t tmp_95 = tmp_40*tmp_88 + tmp_41*tmp_89 + tmp_42*tmp_90;
      real_t tmp_96 = tmp_95 - 1.0/4.0;
      real_t tmp_97 = tmp_0*tmp_92 + tmp_2*tmp_96 + tmp_4*tmp_94;
      real_t tmp_98 = tmp_10*tmp_94 + tmp_7*tmp_96 + tmp_9*tmp_92;
      real_t tmp_99 = tmp_1*tmp_94 + tmp_3*tmp_96 + tmp_6*tmp_92;
      real_t tmp_100 = -tmp_91 - tmp_93 - tmp_95 + 1;
      real_t tmp_101 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_72;
      real_t tmp_102 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_103 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_104 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_105 = tmp_102*tmp_5 + tmp_103*tmp_21 + tmp_104*tmp_27;
      real_t tmp_106 = tmp_105 - 1.0/4.0;
      real_t tmp_107 = tmp_102*tmp_35 + tmp_103*tmp_36 + tmp_104*tmp_37;
      real_t tmp_108 = tmp_107 - 1.0/4.0;
      real_t tmp_109 = tmp_102*tmp_40 + tmp_103*tmp_41 + tmp_104*tmp_42;
      real_t tmp_110 = tmp_109 - 1.0/4.0;
      real_t tmp_111 = tmp_0*tmp_106 + tmp_108*tmp_4 + tmp_110*tmp_2;
      real_t tmp_112 = tmp_10*tmp_108 + tmp_106*tmp_9 + tmp_110*tmp_7;
      real_t tmp_113 = tmp_1*tmp_108 + tmp_106*tmp_6 + tmp_110*tmp_3;
      real_t tmp_114 = -tmp_105 - tmp_107 - tmp_109 + 1;
      real_t tmp_115 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_72;
      real_t tmp_116 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_117 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_118 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_119 = tmp_116*tmp_5 + tmp_117*tmp_21 + tmp_118*tmp_27;
      real_t tmp_120 = tmp_119 - 1.0/4.0;
      real_t tmp_121 = tmp_116*tmp_35 + tmp_117*tmp_36 + tmp_118*tmp_37;
      real_t tmp_122 = tmp_121 - 1.0/4.0;
      real_t tmp_123 = tmp_116*tmp_40 + tmp_117*tmp_41 + tmp_118*tmp_42;
      real_t tmp_124 = tmp_123 - 1.0/4.0;
      real_t tmp_125 = tmp_0*tmp_120 + tmp_122*tmp_4 + tmp_124*tmp_2;
      real_t tmp_126 = tmp_10*tmp_122 + tmp_120*tmp_9 + tmp_124*tmp_7;
      real_t tmp_127 = tmp_1*tmp_122 + tmp_120*tmp_6 + tmp_124*tmp_3;
      real_t tmp_128 = -tmp_119 - tmp_121 - tmp_123 + 1;
      real_t tmp_129 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_72;
      real_t tmp_130 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_131 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_132 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_133 = tmp_130*tmp_5 + tmp_131*tmp_21 + tmp_132*tmp_27;
      real_t tmp_134 = tmp_133 - 1.0/4.0;
      real_t tmp_135 = tmp_130*tmp_35 + tmp_131*tmp_36 + tmp_132*tmp_37;
      real_t tmp_136 = tmp_135 - 1.0/4.0;
      real_t tmp_137 = tmp_130*tmp_40 + tmp_131*tmp_41 + tmp_132*tmp_42;
      real_t tmp_138 = tmp_137 - 1.0/4.0;
      real_t tmp_139 = tmp_0*tmp_134 + tmp_136*tmp_4 + tmp_138*tmp_2;
      real_t tmp_140 = tmp_10*tmp_136 + tmp_134*tmp_9 + tmp_138*tmp_7;
      real_t tmp_141 = tmp_1*tmp_136 + tmp_134*tmp_6 + tmp_138*tmp_3;
      real_t tmp_142 = -tmp_133 - tmp_135 - tmp_137 + 1;
      real_t tmp_143 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_72;
      real_t tmp_144 = 0.25*p_affine_13_0*tmp_15;
      real_t tmp_145 = tmp_144*tmp_5;
      real_t tmp_146 = tmp_144*tmp_21;
      real_t tmp_147 = 0.5*p_affine_13_0*tmp_63 + 0.5*p_affine_13_1*tmp_49 + 0.5*p_affine_13_2*tmp_56;
      real_t tmp_148 = tmp_144*tmp_35;
      real_t tmp_149 = tmp_144*tmp_36;
      real_t tmp_150 = 0.5*p_affine_13_0*tmp_62 + 0.5*p_affine_13_1*tmp_48 + 0.5*p_affine_13_2*tmp_55;
      real_t tmp_151 = tmp_144*tmp_40;
      real_t tmp_152 = tmp_144*tmp_41;
      real_t tmp_153 = 0.5*p_affine_13_0*tmp_61 + 0.5*p_affine_13_1*tmp_47 + 0.5*p_affine_13_2*tmp_54;
      real_t a_0_0 = tmp_101*(1.0*tmp_100*tmp_66*tmp_99 - tmp_100*tmp_71 - tmp_52*tmp_97 - tmp_58*tmp_98 - tmp_64*tmp_99) + tmp_115*(-tmp_111*tmp_52 - tmp_112*tmp_58 + 1.0*tmp_113*tmp_114*tmp_66 - tmp_113*tmp_64 - tmp_114*tmp_71) + tmp_129*(-tmp_125*tmp_52 - tmp_126*tmp_58 + 1.0*tmp_127*tmp_128*tmp_66 - tmp_127*tmp_64 - tmp_128*tmp_71) + tmp_143*(-tmp_139*tmp_52 - tmp_140*tmp_58 + 1.0*tmp_141*tmp_142*tmp_66 - tmp_141*tmp_64 - tmp_142*tmp_71) + tmp_73*(-tmp_45*tmp_52 - tmp_53*tmp_58 - tmp_59*tmp_64 + 1.0*tmp_59*tmp_66*tmp_67 - tmp_67*tmp_71) + tmp_87*(-tmp_52*tmp_83 - tmp_58*tmp_84 - tmp_64*tmp_85 + 1.0*tmp_66*tmp_85*tmp_86 - tmp_71*tmp_86);
      real_t a_1_0 = tmp_101*(-tmp_145*tmp_98 - tmp_146*tmp_97 - tmp_147*tmp_99 + 1.0*tmp_66*tmp_91*tmp_99 - tmp_71*tmp_91) + tmp_115*(1.0*tmp_105*tmp_113*tmp_66 - tmp_105*tmp_71 - tmp_111*tmp_146 - tmp_112*tmp_145 - tmp_113*tmp_147) + tmp_129*(1.0*tmp_119*tmp_127*tmp_66 - tmp_119*tmp_71 - tmp_125*tmp_146 - tmp_126*tmp_145 - tmp_127*tmp_147) + tmp_143*(1.0*tmp_133*tmp_141*tmp_66 - tmp_133*tmp_71 - tmp_139*tmp_146 - tmp_140*tmp_145 - tmp_141*tmp_147) + tmp_73*(-tmp_145*tmp_53 - tmp_146*tmp_45 - tmp_147*tmp_59 + 1.0*tmp_33*tmp_59*tmp_66 - tmp_33*tmp_71) + tmp_87*(-tmp_145*tmp_84 - tmp_146*tmp_83 - tmp_147*tmp_85 + 1.0*tmp_66*tmp_77*tmp_85 - tmp_71*tmp_77);
      real_t a_2_0 = tmp_101*(-tmp_148*tmp_98 - tmp_149*tmp_97 - tmp_150*tmp_99 + 1.0*tmp_66*tmp_93*tmp_99 - tmp_71*tmp_93) + tmp_115*(1.0*tmp_107*tmp_113*tmp_66 - tmp_107*tmp_71 - tmp_111*tmp_149 - tmp_112*tmp_148 - tmp_113*tmp_150) + tmp_129*(1.0*tmp_121*tmp_127*tmp_66 - tmp_121*tmp_71 - tmp_125*tmp_149 - tmp_126*tmp_148 - tmp_127*tmp_150) + tmp_143*(1.0*tmp_135*tmp_141*tmp_66 - tmp_135*tmp_71 - tmp_139*tmp_149 - tmp_140*tmp_148 - tmp_141*tmp_150) + tmp_73*(-tmp_148*tmp_53 - tmp_149*tmp_45 - tmp_150*tmp_59 + 1.0*tmp_38*tmp_59*tmp_66 - tmp_38*tmp_71) + tmp_87*(-tmp_148*tmp_84 - tmp_149*tmp_83 - tmp_150*tmp_85 + 1.0*tmp_66*tmp_79*tmp_85 - tmp_71*tmp_79);
      real_t a_3_0 = tmp_101*(-tmp_151*tmp_98 - tmp_152*tmp_97 - tmp_153*tmp_99 + 1.0*tmp_66*tmp_95*tmp_99 - tmp_71*tmp_95) + tmp_115*(1.0*tmp_109*tmp_113*tmp_66 - tmp_109*tmp_71 - tmp_111*tmp_152 - tmp_112*tmp_151 - tmp_113*tmp_153) + tmp_129*(1.0*tmp_123*tmp_127*tmp_66 - tmp_123*tmp_71 - tmp_125*tmp_152 - tmp_126*tmp_151 - tmp_127*tmp_153) + tmp_143*(1.0*tmp_137*tmp_141*tmp_66 - tmp_137*tmp_71 - tmp_139*tmp_152 - tmp_140*tmp_151 - tmp_141*tmp_153) + tmp_73*(-tmp_151*tmp_53 - tmp_152*tmp_45 - tmp_153*tmp_59 + 1.0*tmp_43*tmp_59*tmp_66 - tmp_43*tmp_71) + tmp_87*(-tmp_151*tmp_84 - tmp_152*tmp_83 - tmp_153*tmp_85 + 1.0*tmp_66*tmp_81*tmp_85 - tmp_71*tmp_81);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_1 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_2 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_5 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_9 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_10 = tmp_5*tmp_9;
      real_t tmp_11 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_12 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_13 = tmp_12*tmp_2;
      real_t tmp_14 = tmp_1*tmp_9;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_12*tmp_4 - tmp_0*tmp_14 + tmp_10*tmp_8 + tmp_11*tmp_3 - tmp_11*tmp_6 - tmp_13*tmp_8);
      real_t tmp_16 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_17 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_18 = -tmp_17;
      real_t tmp_19 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_20 = 0.091576213509770743*tmp_18 + 0.81684757298045851*tmp_19;
      real_t tmp_21 = tmp_15*(tmp_16 + tmp_20);
      real_t tmp_22 = tmp_12*tmp_4 - tmp_14;
      real_t tmp_23 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_24 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_25 = -tmp_24;
      real_t tmp_26 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_27 = 0.091576213509770743*tmp_25 + 0.81684757298045851*tmp_26;
      real_t tmp_28 = tmp_15*(tmp_23 + tmp_27);
      real_t tmp_29 = tmp_10 - tmp_13;
      real_t tmp_30 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_31 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_32 = -tmp_31;
      real_t tmp_33 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_34 = 0.091576213509770743*tmp_32 + 0.81684757298045851*tmp_33;
      real_t tmp_35 = tmp_15*(tmp_30 + tmp_34);
      real_t tmp_36 = tmp_21*tmp_7 + tmp_22*tmp_28 + tmp_29*tmp_35 - 1.0/4.0;
      real_t tmp_37 = tmp_0*tmp_4 - tmp_2*tmp_8;
      real_t tmp_38 = -tmp_11*tmp_4 + tmp_8*tmp_9;
      real_t tmp_39 = -tmp_0*tmp_9 + tmp_11*tmp_2;
      real_t tmp_40 = tmp_21*tmp_37 + tmp_28*tmp_38 + tmp_35*tmp_39 - 1.0/4.0;
      real_t tmp_41 = -tmp_0*tmp_1 + tmp_5*tmp_8;
      real_t tmp_42 = tmp_1*tmp_11 - tmp_12*tmp_8;
      real_t tmp_43 = tmp_0*tmp_12 - tmp_11*tmp_5;
      real_t tmp_44 = tmp_21*tmp_41 + tmp_28*tmp_42 + tmp_35*tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_36 + tmp_2*tmp_44 + tmp_40*tmp_5;
      real_t tmp_46 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_47 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_48 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_49 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_50 = -tmp_46*tmp_47 + tmp_48*tmp_49;
      real_t tmp_51 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_52 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_53 = tmp_51*tmp_52;
      real_t tmp_54 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_55 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_56 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_57 = tmp_47*tmp_56;
      real_t tmp_58 = tmp_47*tmp_54;
      real_t tmp_59 = tmp_52*tmp_56;
      real_t tmp_60 = tmp_49*tmp_51;
      real_t tmp_61 = 1.0 / (tmp_46*tmp_53 - tmp_46*tmp_58 + tmp_48*tmp_49*tmp_54 - tmp_48*tmp_59 + tmp_55*tmp_57 - tmp_55*tmp_60);
      real_t tmp_62 = 0.5*tmp_61;
      real_t tmp_63 = tmp_50*tmp_62;
      real_t tmp_64 = tmp_46*tmp_52 - tmp_49*tmp_55;
      real_t tmp_65 = tmp_62*tmp_64;
      real_t tmp_66 = tmp_47*tmp_55 - tmp_48*tmp_52;
      real_t tmp_67 = tmp_62*tmp_66;
      real_t tmp_68 = -tmp_63 - tmp_65 - tmp_67;
      real_t tmp_69 = 0.5*p_affine_13_0;
      real_t tmp_70 = tmp_68*tmp_69;
      real_t tmp_71 = tmp_11*tmp_36 + tmp_12*tmp_40 + tmp_44*tmp_9;
      real_t tmp_72 = tmp_46*tmp_51 - tmp_48*tmp_56;
      real_t tmp_73 = tmp_62*tmp_72;
      real_t tmp_74 = -tmp_46*tmp_54 + tmp_55*tmp_56;
      real_t tmp_75 = tmp_62*tmp_74;
      real_t tmp_76 = tmp_48*tmp_54 - tmp_51*tmp_55;
      real_t tmp_77 = tmp_62*tmp_76;
      real_t tmp_78 = -tmp_73 - tmp_75 - tmp_77;
      real_t tmp_79 = tmp_69*tmp_78;
      real_t tmp_80 = tmp_1*tmp_40 + tmp_36*tmp_8 + tmp_4*tmp_44;
      real_t tmp_81 = tmp_57 - tmp_60;
      real_t tmp_82 = 1.0*tmp_61;
      real_t tmp_83 = tmp_81*tmp_82;
      real_t tmp_84 = tmp_49*tmp_54 - tmp_59;
      real_t tmp_85 = tmp_82*tmp_84;
      real_t tmp_86 = tmp_53 - tmp_58;
      real_t tmp_87 = tmp_82*tmp_86;
      real_t tmp_88 = 0.5*p_affine_13_0*(-tmp_83 - tmp_85 - tmp_87) + 0.5*p_affine_13_1*tmp_68 + 0.5*p_affine_13_2*tmp_78;
      real_t tmp_89 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_90 = tmp_61*(tmp_20 + tmp_89);
      real_t tmp_91 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_92 = tmp_61*(tmp_27 + tmp_91);
      real_t tmp_93 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_94 = tmp_61*(tmp_34 + tmp_93);
      real_t tmp_95 = tmp_66*tmp_92 + tmp_76*tmp_90 + tmp_86*tmp_94;
      real_t tmp_96 = tmp_64*tmp_92 + tmp_74*tmp_90 + tmp_84*tmp_94;
      real_t tmp_97 = tmp_50*tmp_92 + tmp_72*tmp_90 + tmp_81*tmp_94;
      real_t tmp_98 = -tmp_95 - tmp_96 - tmp_97 + 1;
      real_t tmp_99 = (std::abs(tmp_17*tmp_26 - tmp_19*tmp_24)*std::abs(tmp_17*tmp_26 - tmp_19*tmp_24)) + (std::abs(tmp_17*tmp_33 - tmp_19*tmp_31)*std::abs(tmp_17*tmp_33 - tmp_19*tmp_31)) + (std::abs(tmp_24*tmp_33 - tmp_26*tmp_31)*std::abs(tmp_24*tmp_33 - tmp_26*tmp_31));
      real_t tmp_100 = 1.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_80;
      real_t tmp_102 = 1.0*tmp_15;
      real_t tmp_103 = 0.5*tmp_15;
      real_t tmp_104 = tmp_103*tmp_8;
      real_t tmp_105 = tmp_1*tmp_103;
      real_t tmp_106 = tmp_103*tmp_4;
      real_t tmp_107 = tmp_103*tmp_29;
      real_t tmp_108 = tmp_103*tmp_39;
      real_t tmp_109 = tmp_103*tmp_43;
      real_t tmp_110 = 0.5*p_affine_13_0*(tmp_1*tmp_102*tmp_39 + tmp_102*tmp_29*tmp_8 + tmp_102*tmp_4*tmp_43) + 0.5*p_affine_13_1*(tmp_0*tmp_107 + tmp_104*tmp_22 + tmp_105*tmp_38 + tmp_106*tmp_42 + tmp_108*tmp_5 + tmp_109*tmp_2) + 0.5*p_affine_13_2*(tmp_104*tmp_7 + tmp_105*tmp_37 + tmp_106*tmp_41 + tmp_107*tmp_11 + tmp_108*tmp_12 + tmp_109*tmp_9);
      real_t tmp_111 = 2.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_112 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_111;
      real_t tmp_113 = 0.44594849091596489*tmp_18 + 0.10810301816807022*tmp_19;
      real_t tmp_114 = tmp_15*(tmp_113 + tmp_16);
      real_t tmp_115 = 0.44594849091596489*tmp_25 + 0.10810301816807022*tmp_26;
      real_t tmp_116 = tmp_15*(tmp_115 + tmp_23);
      real_t tmp_117 = 0.44594849091596489*tmp_32 + 0.10810301816807022*tmp_33;
      real_t tmp_118 = tmp_15*(tmp_117 + tmp_30);
      real_t tmp_119 = tmp_114*tmp_7 + tmp_116*tmp_22 + tmp_118*tmp_29 - 1.0/4.0;
      real_t tmp_120 = tmp_114*tmp_37 + tmp_116*tmp_38 + tmp_118*tmp_39 - 1.0/4.0;
      real_t tmp_121 = tmp_114*tmp_41 + tmp_116*tmp_42 + tmp_118*tmp_43 - 1.0/4.0;
      real_t tmp_122 = tmp_0*tmp_119 + tmp_120*tmp_5 + tmp_121*tmp_2;
      real_t tmp_123 = tmp_11*tmp_119 + tmp_12*tmp_120 + tmp_121*tmp_9;
      real_t tmp_124 = tmp_1*tmp_120 + tmp_119*tmp_8 + tmp_121*tmp_4;
      real_t tmp_125 = tmp_61*(tmp_113 + tmp_89);
      real_t tmp_126 = tmp_61*(tmp_115 + tmp_91);
      real_t tmp_127 = tmp_61*(tmp_117 + tmp_93);
      real_t tmp_128 = tmp_125*tmp_76 + tmp_126*tmp_66 + tmp_127*tmp_86;
      real_t tmp_129 = tmp_125*tmp_74 + tmp_126*tmp_64 + tmp_127*tmp_84;
      real_t tmp_130 = tmp_125*tmp_72 + tmp_126*tmp_50 + tmp_127*tmp_81;
      real_t tmp_131 = -tmp_128 - tmp_129 - tmp_130 + 1;
      real_t tmp_132 = tmp_100*tmp_124;
      real_t tmp_133 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_111;
      real_t tmp_134 = 0.81684757298045851*tmp_18 + 0.091576213509770743*tmp_19;
      real_t tmp_135 = tmp_15*(tmp_134 + tmp_16);
      real_t tmp_136 = 0.81684757298045851*tmp_25 + 0.091576213509770743*tmp_26;
      real_t tmp_137 = tmp_15*(tmp_136 + tmp_23);
      real_t tmp_138 = 0.81684757298045851*tmp_32 + 0.091576213509770743*tmp_33;
      real_t tmp_139 = tmp_15*(tmp_138 + tmp_30);
      real_t tmp_140 = tmp_135*tmp_7 + tmp_137*tmp_22 + tmp_139*tmp_29 - 1.0/4.0;
      real_t tmp_141 = tmp_135*tmp_37 + tmp_137*tmp_38 + tmp_139*tmp_39 - 1.0/4.0;
      real_t tmp_142 = tmp_135*tmp_41 + tmp_137*tmp_42 + tmp_139*tmp_43 - 1.0/4.0;
      real_t tmp_143 = tmp_0*tmp_140 + tmp_141*tmp_5 + tmp_142*tmp_2;
      real_t tmp_144 = tmp_11*tmp_140 + tmp_12*tmp_141 + tmp_142*tmp_9;
      real_t tmp_145 = tmp_1*tmp_141 + tmp_140*tmp_8 + tmp_142*tmp_4;
      real_t tmp_146 = tmp_61*(tmp_134 + tmp_89);
      real_t tmp_147 = tmp_61*(tmp_136 + tmp_91);
      real_t tmp_148 = tmp_61*(tmp_138 + tmp_93);
      real_t tmp_149 = tmp_146*tmp_76 + tmp_147*tmp_66 + tmp_148*tmp_86;
      real_t tmp_150 = tmp_146*tmp_74 + tmp_147*tmp_64 + tmp_148*tmp_84;
      real_t tmp_151 = tmp_146*tmp_72 + tmp_147*tmp_50 + tmp_148*tmp_81;
      real_t tmp_152 = -tmp_149 - tmp_150 - tmp_151 + 1;
      real_t tmp_153 = tmp_100*tmp_145;
      real_t tmp_154 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_111;
      real_t tmp_155 = 0.10810301816807022*tmp_18 + 0.44594849091596489*tmp_19;
      real_t tmp_156 = tmp_15*(tmp_155 + tmp_16);
      real_t tmp_157 = 0.10810301816807022*tmp_25 + 0.44594849091596489*tmp_26;
      real_t tmp_158 = tmp_15*(tmp_157 + tmp_23);
      real_t tmp_159 = 0.10810301816807022*tmp_32 + 0.44594849091596489*tmp_33;
      real_t tmp_160 = tmp_15*(tmp_159 + tmp_30);
      real_t tmp_161 = tmp_156*tmp_7 + tmp_158*tmp_22 + tmp_160*tmp_29 - 1.0/4.0;
      real_t tmp_162 = tmp_156*tmp_37 + tmp_158*tmp_38 + tmp_160*tmp_39 - 1.0/4.0;
      real_t tmp_163 = tmp_156*tmp_41 + tmp_158*tmp_42 + tmp_160*tmp_43 - 1.0/4.0;
      real_t tmp_164 = tmp_0*tmp_161 + tmp_162*tmp_5 + tmp_163*tmp_2;
      real_t tmp_165 = tmp_11*tmp_161 + tmp_12*tmp_162 + tmp_163*tmp_9;
      real_t tmp_166 = tmp_1*tmp_162 + tmp_161*tmp_8 + tmp_163*tmp_4;
      real_t tmp_167 = tmp_61*(tmp_155 + tmp_89);
      real_t tmp_168 = tmp_61*(tmp_157 + tmp_91);
      real_t tmp_169 = tmp_61*(tmp_159 + tmp_93);
      real_t tmp_170 = tmp_167*tmp_76 + tmp_168*tmp_66 + tmp_169*tmp_86;
      real_t tmp_171 = tmp_167*tmp_74 + tmp_168*tmp_64 + tmp_169*tmp_84;
      real_t tmp_172 = tmp_167*tmp_72 + tmp_168*tmp_50 + tmp_169*tmp_81;
      real_t tmp_173 = -tmp_170 - tmp_171 - tmp_172 + 1;
      real_t tmp_174 = tmp_100*tmp_166;
      real_t tmp_175 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_111;
      real_t tmp_176 = 0.091576213509770743*tmp_18 + 0.091576213509770743*tmp_19;
      real_t tmp_177 = tmp_15*(tmp_16 + tmp_176);
      real_t tmp_178 = 0.091576213509770743*tmp_25 + 0.091576213509770743*tmp_26;
      real_t tmp_179 = tmp_15*(tmp_178 + tmp_23);
      real_t tmp_180 = 0.091576213509770743*tmp_32 + 0.091576213509770743*tmp_33;
      real_t tmp_181 = tmp_15*(tmp_180 + tmp_30);
      real_t tmp_182 = tmp_177*tmp_7 + tmp_179*tmp_22 + tmp_181*tmp_29 - 1.0/4.0;
      real_t tmp_183 = tmp_177*tmp_37 + tmp_179*tmp_38 + tmp_181*tmp_39 - 1.0/4.0;
      real_t tmp_184 = tmp_177*tmp_41 + tmp_179*tmp_42 + tmp_181*tmp_43 - 1.0/4.0;
      real_t tmp_185 = tmp_0*tmp_182 + tmp_183*tmp_5 + tmp_184*tmp_2;
      real_t tmp_186 = tmp_11*tmp_182 + tmp_12*tmp_183 + tmp_184*tmp_9;
      real_t tmp_187 = tmp_1*tmp_183 + tmp_182*tmp_8 + tmp_184*tmp_4;
      real_t tmp_188 = tmp_61*(tmp_176 + tmp_89);
      real_t tmp_189 = tmp_61*(tmp_178 + tmp_91);
      real_t tmp_190 = tmp_61*(tmp_180 + tmp_93);
      real_t tmp_191 = tmp_188*tmp_76 + tmp_189*tmp_66 + tmp_190*tmp_86;
      real_t tmp_192 = tmp_188*tmp_74 + tmp_189*tmp_64 + tmp_190*tmp_84;
      real_t tmp_193 = tmp_188*tmp_72 + tmp_189*tmp_50 + tmp_190*tmp_81;
      real_t tmp_194 = -tmp_191 - tmp_192 - tmp_193 + 1;
      real_t tmp_195 = tmp_100*tmp_187;
      real_t tmp_196 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_111;
      real_t tmp_197 = 0.44594849091596489*tmp_18 + 0.44594849091596489*tmp_19;
      real_t tmp_198 = tmp_15*(tmp_16 + tmp_197);
      real_t tmp_199 = 0.44594849091596489*tmp_25 + 0.44594849091596489*tmp_26;
      real_t tmp_200 = tmp_15*(tmp_199 + tmp_23);
      real_t tmp_201 = 0.44594849091596489*tmp_32 + 0.44594849091596489*tmp_33;
      real_t tmp_202 = tmp_15*(tmp_201 + tmp_30);
      real_t tmp_203 = tmp_198*tmp_7 + tmp_200*tmp_22 + tmp_202*tmp_29 - 1.0/4.0;
      real_t tmp_204 = tmp_198*tmp_37 + tmp_200*tmp_38 + tmp_202*tmp_39 - 1.0/4.0;
      real_t tmp_205 = tmp_198*tmp_41 + tmp_200*tmp_42 + tmp_202*tmp_43 - 1.0/4.0;
      real_t tmp_206 = tmp_0*tmp_203 + tmp_2*tmp_205 + tmp_204*tmp_5;
      real_t tmp_207 = tmp_11*tmp_203 + tmp_12*tmp_204 + tmp_205*tmp_9;
      real_t tmp_208 = tmp_1*tmp_204 + tmp_203*tmp_8 + tmp_205*tmp_4;
      real_t tmp_209 = tmp_61*(tmp_197 + tmp_89);
      real_t tmp_210 = tmp_61*(tmp_199 + tmp_91);
      real_t tmp_211 = tmp_61*(tmp_201 + tmp_93);
      real_t tmp_212 = tmp_209*tmp_76 + tmp_210*tmp_66 + tmp_211*tmp_86;
      real_t tmp_213 = tmp_209*tmp_74 + tmp_210*tmp_64 + tmp_211*tmp_84;
      real_t tmp_214 = tmp_209*tmp_72 + tmp_210*tmp_50 + tmp_211*tmp_81;
      real_t tmp_215 = -tmp_212 - tmp_213 - tmp_214 + 1;
      real_t tmp_216 = tmp_100*tmp_208;
      real_t tmp_217 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_111;
      real_t tmp_218 = 0.25*p_affine_13_0*tmp_61;
      real_t tmp_219 = tmp_218*tmp_76;
      real_t tmp_220 = tmp_218*tmp_66;
      real_t tmp_221 = 0.5*p_affine_13_0*tmp_87 + 0.5*p_affine_13_1*tmp_67 + 0.5*p_affine_13_2*tmp_77;
      real_t tmp_222 = tmp_218*tmp_74;
      real_t tmp_223 = tmp_218*tmp_64;
      real_t tmp_224 = 0.5*p_affine_13_0*tmp_85 + 0.5*p_affine_13_1*tmp_65 + 0.5*p_affine_13_2*tmp_75;
      real_t tmp_225 = tmp_218*tmp_72;
      real_t tmp_226 = tmp_218*tmp_50;
      real_t tmp_227 = 0.5*p_affine_13_0*tmp_83 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_73;
      real_t a_0_0 = tmp_112*(-tmp_101*tmp_98 - tmp_110*tmp_98 + tmp_45*tmp_70 + tmp_71*tmp_79 + tmp_80*tmp_88) + tmp_133*(-tmp_110*tmp_131 + tmp_122*tmp_70 + tmp_123*tmp_79 + tmp_124*tmp_88 - tmp_131*tmp_132) + tmp_154*(-tmp_110*tmp_152 + tmp_143*tmp_70 + tmp_144*tmp_79 + tmp_145*tmp_88 - tmp_152*tmp_153) + tmp_175*(-tmp_110*tmp_173 + tmp_164*tmp_70 + tmp_165*tmp_79 + tmp_166*tmp_88 - tmp_173*tmp_174) + tmp_196*(-tmp_110*tmp_194 + tmp_185*tmp_70 + tmp_186*tmp_79 + tmp_187*tmp_88 - tmp_194*tmp_195) + tmp_217*(-tmp_110*tmp_215 + tmp_206*tmp_70 + tmp_207*tmp_79 + tmp_208*tmp_88 - tmp_215*tmp_216);
      real_t a_1_0 = tmp_112*(-tmp_101*tmp_95 - tmp_110*tmp_95 + tmp_219*tmp_71 + tmp_220*tmp_45 + tmp_221*tmp_80) + tmp_133*(-tmp_110*tmp_128 + tmp_122*tmp_220 + tmp_123*tmp_219 + tmp_124*tmp_221 - tmp_128*tmp_132) + tmp_154*(-tmp_110*tmp_149 + tmp_143*tmp_220 + tmp_144*tmp_219 + tmp_145*tmp_221 - tmp_149*tmp_153) + tmp_175*(-tmp_110*tmp_170 + tmp_164*tmp_220 + tmp_165*tmp_219 + tmp_166*tmp_221 - tmp_170*tmp_174) + tmp_196*(-tmp_110*tmp_191 + tmp_185*tmp_220 + tmp_186*tmp_219 + tmp_187*tmp_221 - tmp_191*tmp_195) + tmp_217*(-tmp_110*tmp_212 + tmp_206*tmp_220 + tmp_207*tmp_219 + tmp_208*tmp_221 - tmp_212*tmp_216);
      real_t a_2_0 = tmp_112*(-tmp_101*tmp_96 - tmp_110*tmp_96 + tmp_222*tmp_71 + tmp_223*tmp_45 + tmp_224*tmp_80) + tmp_133*(-tmp_110*tmp_129 + tmp_122*tmp_223 + tmp_123*tmp_222 + tmp_124*tmp_224 - tmp_129*tmp_132) + tmp_154*(-tmp_110*tmp_150 + tmp_143*tmp_223 + tmp_144*tmp_222 + tmp_145*tmp_224 - tmp_150*tmp_153) + tmp_175*(-tmp_110*tmp_171 + tmp_164*tmp_223 + tmp_165*tmp_222 + tmp_166*tmp_224 - tmp_171*tmp_174) + tmp_196*(-tmp_110*tmp_192 + tmp_185*tmp_223 + tmp_186*tmp_222 + tmp_187*tmp_224 - tmp_192*tmp_195) + tmp_217*(-tmp_110*tmp_213 + tmp_206*tmp_223 + tmp_207*tmp_222 + tmp_208*tmp_224 - tmp_213*tmp_216);
      real_t a_3_0 = tmp_112*(-tmp_101*tmp_97 - tmp_110*tmp_97 + tmp_225*tmp_71 + tmp_226*tmp_45 + tmp_227*tmp_80) + tmp_133*(-tmp_110*tmp_130 + tmp_122*tmp_226 + tmp_123*tmp_225 + tmp_124*tmp_227 - tmp_130*tmp_132) + tmp_154*(-tmp_110*tmp_151 + tmp_143*tmp_226 + tmp_144*tmp_225 + tmp_145*tmp_227 - tmp_151*tmp_153) + tmp_175*(-tmp_110*tmp_172 + tmp_164*tmp_226 + tmp_165*tmp_225 + tmp_166*tmp_227 - tmp_172*tmp_174) + tmp_196*(-tmp_110*tmp_193 + tmp_185*tmp_226 + tmp_186*tmp_225 + tmp_187*tmp_227 - tmp_193*tmp_195) + tmp_217*(-tmp_110*tmp_214 + tmp_206*tmp_226 + tmp_207*tmp_225 + tmp_208*tmp_227 - tmp_214*tmp_216);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_0*tmp_10;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = tmp_4*tmp_9;
      real_t tmp_15 = 1.0 / (-tmp_1*tmp_13 + tmp_1*tmp_2*tmp_9 + tmp_11*tmp_3 - tmp_12*tmp_6 - tmp_14*tmp_3 + tmp_6*tmp_8);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_8;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = tmp_33 - 1.0/4.0;
      real_t tmp_35 = tmp_0*tmp_3 - tmp_2*tmp_6;
      real_t tmp_36 = -tmp_3*tmp_9 + tmp_6*tmp_7;
      real_t tmp_37 = -tmp_13 + tmp_2*tmp_9;
      real_t tmp_38 = tmp_20*tmp_35 + tmp_26*tmp_36 + tmp_32*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/4.0;
      real_t tmp_40 = -tmp_0*tmp_1 + tmp_4*tmp_6;
      real_t tmp_41 = tmp_1*tmp_9 - tmp_10*tmp_6;
      real_t tmp_42 = tmp_11 - tmp_14;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_26*tmp_41 + tmp_32*tmp_42;
      real_t tmp_44 = tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_34 + tmp_2*tmp_44 + tmp_39*tmp_4;
      real_t tmp_46 = 0.5*tmp_15;
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_36*tmp_46;
      real_t tmp_49 = tmp_21*tmp_46;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = p_affine_13_0*tmp_50;
      real_t tmp_52 = tmp_10*tmp_39 + tmp_34*tmp_9 + tmp_44*tmp_7;
      real_t tmp_53 = tmp_40*tmp_46;
      real_t tmp_54 = tmp_35*tmp_46;
      real_t tmp_55 = tmp_46*tmp_5;
      real_t tmp_56 = -tmp_53 - tmp_54 - tmp_55;
      real_t tmp_57 = p_affine_13_0*tmp_56;
      real_t tmp_58 = 1.0*tmp_15;
      real_t tmp_59 = tmp_42*tmp_58;
      real_t tmp_60 = tmp_37*tmp_58;
      real_t tmp_61 = tmp_27*tmp_58;
      real_t tmp_62 = p_affine_13_0*(-tmp_59 - tmp_60 - tmp_61) + p_affine_13_1*tmp_50 + p_affine_13_2*tmp_56;
      real_t tmp_63 = tmp_1*tmp_39 + tmp_3*tmp_44 + tmp_34*tmp_6;
      real_t tmp_64 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_65 = std::pow(tmp_64, -0.25);
      real_t tmp_66 = -tmp_33 - tmp_38 - tmp_43 + 1;
      real_t tmp_67 = tmp_27*tmp_46;
      real_t tmp_68 = tmp_37*tmp_46;
      real_t tmp_69 = tmp_42*tmp_46;
      real_t tmp_70 = p_affine_13_0*(tmp_1*tmp_60 + tmp_3*tmp_59 + tmp_6*tmp_61) + p_affine_13_1*(tmp_0*tmp_67 + tmp_1*tmp_48 + tmp_2*tmp_69 + tmp_3*tmp_47 + tmp_4*tmp_68 + tmp_49*tmp_6) + p_affine_13_2*(tmp_1*tmp_54 + tmp_10*tmp_68 + tmp_3*tmp_53 + tmp_55*tmp_6 + tmp_67*tmp_9 + tmp_69*tmp_7);
      real_t tmp_71 = 2.0*std::pow(tmp_64, 1.0/2.0);
      real_t tmp_72 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_71;
      real_t tmp_73 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_74 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_75 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_76 = tmp_21*tmp_74 + tmp_27*tmp_75 + tmp_5*tmp_73;
      real_t tmp_77 = tmp_76 - 1.0/4.0;
      real_t tmp_78 = tmp_35*tmp_73 + tmp_36*tmp_74 + tmp_37*tmp_75;
      real_t tmp_79 = tmp_78 - 1.0/4.0;
      real_t tmp_80 = tmp_40*tmp_73 + tmp_41*tmp_74 + tmp_42*tmp_75;
      real_t tmp_81 = tmp_80 - 1.0/4.0;
      real_t tmp_82 = tmp_0*tmp_77 + tmp_2*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_10*tmp_79 + tmp_7*tmp_81 + tmp_77*tmp_9;
      real_t tmp_84 = tmp_1*tmp_79 + tmp_3*tmp_81 + tmp_6*tmp_77;
      real_t tmp_85 = -tmp_76 - tmp_78 - tmp_80 + 1;
      real_t tmp_86 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_71;
      real_t tmp_87 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_88 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_89 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_90 = tmp_21*tmp_88 + tmp_27*tmp_89 + tmp_5*tmp_87;
      real_t tmp_91 = tmp_90 - 1.0/4.0;
      real_t tmp_92 = tmp_35*tmp_87 + tmp_36*tmp_88 + tmp_37*tmp_89;
      real_t tmp_93 = tmp_92 - 1.0/4.0;
      real_t tmp_94 = tmp_40*tmp_87 + tmp_41*tmp_88 + tmp_42*tmp_89;
      real_t tmp_95 = tmp_94 - 1.0/4.0;
      real_t tmp_96 = tmp_0*tmp_91 + tmp_2*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_10*tmp_93 + tmp_7*tmp_95 + tmp_9*tmp_91;
      real_t tmp_98 = tmp_1*tmp_93 + tmp_3*tmp_95 + tmp_6*tmp_91;
      real_t tmp_99 = -tmp_90 - tmp_92 - tmp_94 + 1;
      real_t tmp_100 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_71;
      real_t tmp_101 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_102 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_103 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_104 = tmp_101*tmp_5 + tmp_102*tmp_21 + tmp_103*tmp_27;
      real_t tmp_105 = tmp_104 - 1.0/4.0;
      real_t tmp_106 = tmp_101*tmp_35 + tmp_102*tmp_36 + tmp_103*tmp_37;
      real_t tmp_107 = tmp_106 - 1.0/4.0;
      real_t tmp_108 = tmp_101*tmp_40 + tmp_102*tmp_41 + tmp_103*tmp_42;
      real_t tmp_109 = tmp_108 - 1.0/4.0;
      real_t tmp_110 = tmp_0*tmp_105 + tmp_107*tmp_4 + tmp_109*tmp_2;
      real_t tmp_111 = tmp_10*tmp_107 + tmp_105*tmp_9 + tmp_109*tmp_7;
      real_t tmp_112 = tmp_1*tmp_107 + tmp_105*tmp_6 + tmp_109*tmp_3;
      real_t tmp_113 = -tmp_104 - tmp_106 - tmp_108 + 1;
      real_t tmp_114 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_71;
      real_t tmp_115 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_116 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_117 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_118 = tmp_115*tmp_5 + tmp_116*tmp_21 + tmp_117*tmp_27;
      real_t tmp_119 = tmp_118 - 1.0/4.0;
      real_t tmp_120 = tmp_115*tmp_35 + tmp_116*tmp_36 + tmp_117*tmp_37;
      real_t tmp_121 = tmp_120 - 1.0/4.0;
      real_t tmp_122 = tmp_115*tmp_40 + tmp_116*tmp_41 + tmp_117*tmp_42;
      real_t tmp_123 = tmp_122 - 1.0/4.0;
      real_t tmp_124 = tmp_0*tmp_119 + tmp_121*tmp_4 + tmp_123*tmp_2;
      real_t tmp_125 = tmp_10*tmp_121 + tmp_119*tmp_9 + tmp_123*tmp_7;
      real_t tmp_126 = tmp_1*tmp_121 + tmp_119*tmp_6 + tmp_123*tmp_3;
      real_t tmp_127 = -tmp_118 - tmp_120 - tmp_122 + 1;
      real_t tmp_128 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_71;
      real_t tmp_129 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_130 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_131 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_132 = tmp_129*tmp_5 + tmp_130*tmp_21 + tmp_131*tmp_27;
      real_t tmp_133 = tmp_132 - 1.0/4.0;
      real_t tmp_134 = tmp_129*tmp_35 + tmp_130*tmp_36 + tmp_131*tmp_37;
      real_t tmp_135 = tmp_134 - 1.0/4.0;
      real_t tmp_136 = tmp_129*tmp_40 + tmp_130*tmp_41 + tmp_131*tmp_42;
      real_t tmp_137 = tmp_136 - 1.0/4.0;
      real_t tmp_138 = tmp_0*tmp_133 + tmp_135*tmp_4 + tmp_137*tmp_2;
      real_t tmp_139 = tmp_10*tmp_135 + tmp_133*tmp_9 + tmp_137*tmp_7;
      real_t tmp_140 = tmp_1*tmp_135 + tmp_133*tmp_6 + tmp_137*tmp_3;
      real_t tmp_141 = -tmp_132 - tmp_134 - tmp_136 + 1;
      real_t tmp_142 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_71;
      real_t tmp_143 = p_affine_13_0*tmp_55;
      real_t tmp_144 = p_affine_13_0*tmp_49;
      real_t tmp_145 = p_affine_13_0*tmp_61 + p_affine_13_1*tmp_49 + p_affine_13_2*tmp_55;
      real_t tmp_146 = p_affine_13_0*tmp_54;
      real_t tmp_147 = p_affine_13_0*tmp_48;
      real_t tmp_148 = p_affine_13_0*tmp_60 + p_affine_13_1*tmp_48 + p_affine_13_2*tmp_54;
      real_t tmp_149 = p_affine_13_0*tmp_53;
      real_t tmp_150 = p_affine_13_0*tmp_47;
      real_t tmp_151 = p_affine_13_0*tmp_59 + p_affine_13_1*tmp_47 + p_affine_13_2*tmp_53;
      real_t a_0_0 = tmp_100*(-tmp_51*tmp_96 - tmp_57*tmp_97 - tmp_62*tmp_98 + 1.0*tmp_65*tmp_98*tmp_99 - tmp_70*tmp_99) + tmp_114*(-tmp_110*tmp_51 - tmp_111*tmp_57 + 1.0*tmp_112*tmp_113*tmp_65 - tmp_112*tmp_62 - tmp_113*tmp_70) + tmp_128*(-tmp_124*tmp_51 - tmp_125*tmp_57 + 1.0*tmp_126*tmp_127*tmp_65 - tmp_126*tmp_62 - tmp_127*tmp_70) + tmp_142*(-tmp_138*tmp_51 - tmp_139*tmp_57 + 1.0*tmp_140*tmp_141*tmp_65 - tmp_140*tmp_62 - tmp_141*tmp_70) + tmp_72*(-tmp_45*tmp_51 - tmp_52*tmp_57 - tmp_62*tmp_63 + 1.0*tmp_63*tmp_65*tmp_66 - tmp_66*tmp_70) + tmp_86*(-tmp_51*tmp_82 - tmp_57*tmp_83 - tmp_62*tmp_84 + 1.0*tmp_65*tmp_84*tmp_85 - tmp_70*tmp_85);
      real_t a_1_0 = tmp_100*(-tmp_143*tmp_97 - tmp_144*tmp_96 - tmp_145*tmp_98 + 1.0*tmp_65*tmp_90*tmp_98 - tmp_70*tmp_90) + tmp_114*(1.0*tmp_104*tmp_112*tmp_65 - tmp_104*tmp_70 - tmp_110*tmp_144 - tmp_111*tmp_143 - tmp_112*tmp_145) + tmp_128*(1.0*tmp_118*tmp_126*tmp_65 - tmp_118*tmp_70 - tmp_124*tmp_144 - tmp_125*tmp_143 - tmp_126*tmp_145) + tmp_142*(1.0*tmp_132*tmp_140*tmp_65 - tmp_132*tmp_70 - tmp_138*tmp_144 - tmp_139*tmp_143 - tmp_140*tmp_145) + tmp_72*(-tmp_143*tmp_52 - tmp_144*tmp_45 - tmp_145*tmp_63 + 1.0*tmp_33*tmp_63*tmp_65 - tmp_33*tmp_70) + tmp_86*(-tmp_143*tmp_83 - tmp_144*tmp_82 - tmp_145*tmp_84 + 1.0*tmp_65*tmp_76*tmp_84 - tmp_70*tmp_76);
      real_t a_2_0 = tmp_100*(-tmp_146*tmp_97 - tmp_147*tmp_96 - tmp_148*tmp_98 + 1.0*tmp_65*tmp_92*tmp_98 - tmp_70*tmp_92) + tmp_114*(1.0*tmp_106*tmp_112*tmp_65 - tmp_106*tmp_70 - tmp_110*tmp_147 - tmp_111*tmp_146 - tmp_112*tmp_148) + tmp_128*(1.0*tmp_120*tmp_126*tmp_65 - tmp_120*tmp_70 - tmp_124*tmp_147 - tmp_125*tmp_146 - tmp_126*tmp_148) + tmp_142*(1.0*tmp_134*tmp_140*tmp_65 - tmp_134*tmp_70 - tmp_138*tmp_147 - tmp_139*tmp_146 - tmp_140*tmp_148) + tmp_72*(-tmp_146*tmp_52 - tmp_147*tmp_45 - tmp_148*tmp_63 + 1.0*tmp_38*tmp_63*tmp_65 - tmp_38*tmp_70) + tmp_86*(-tmp_146*tmp_83 - tmp_147*tmp_82 - tmp_148*tmp_84 + 1.0*tmp_65*tmp_78*tmp_84 - tmp_70*tmp_78);
      real_t a_3_0 = tmp_100*(-tmp_149*tmp_97 - tmp_150*tmp_96 - tmp_151*tmp_98 + 1.0*tmp_65*tmp_94*tmp_98 - tmp_70*tmp_94) + tmp_114*(1.0*tmp_108*tmp_112*tmp_65 - tmp_108*tmp_70 - tmp_110*tmp_150 - tmp_111*tmp_149 - tmp_112*tmp_151) + tmp_128*(1.0*tmp_122*tmp_126*tmp_65 - tmp_122*tmp_70 - tmp_124*tmp_150 - tmp_125*tmp_149 - tmp_126*tmp_151) + tmp_142*(1.0*tmp_136*tmp_140*tmp_65 - tmp_136*tmp_70 - tmp_138*tmp_150 - tmp_139*tmp_149 - tmp_140*tmp_151) + tmp_72*(-tmp_149*tmp_52 - tmp_150*tmp_45 - tmp_151*tmp_63 + 1.0*tmp_43*tmp_63*tmp_65 - tmp_43*tmp_70) + tmp_86*(-tmp_149*tmp_83 - tmp_150*tmp_82 - tmp_151*tmp_84 + 1.0*tmp_65*tmp_80*tmp_84 - tmp_70*tmp_80);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }

};




class EGEpsilonFormEP1_1 : public hyteg::dg::DGForm
{

 public:
    EGEpsilonFormEP1_1(std::function< real_t ( const Point3D & ) > mu)
: callback_Scalar_Variable_Coefficient_3D_mu (mu)
, callback_Scalar_Variable_Coefficient_2D_mu (mu)
    {}

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_mu;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;

void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
}

      
 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = -tmp_3;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_2 - tmp_4*tmp_6);
      real_t tmp_8 = 1.0*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_2*tmp_8 + tmp_6*tmp_9;
      real_t tmp_11 = Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_10;
      real_t tmp_12 = -2*tmp_0*tmp_8 - 2*tmp_9;
      real_t tmp_13 = 0.5*tmp_7;
      real_t tmp_14 = tmp_0*tmp_13;
      real_t tmp_15 = tmp_1*tmp_13;
      real_t tmp_16 = tmp_13*tmp_5;
      real_t tmp_17 = tmp_1*tmp_16 + tmp_14*tmp_3 + tmp_14*tmp_4 + tmp_15*tmp_6;
      real_t tmp_18 = Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_17;
      real_t tmp_19 = -4*tmp_15 - 4*tmp_16;
      real_t tmp_20 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_21 = 0.054975871827660928*tmp_20;
      real_t tmp_22 = tmp_10*tmp_12;
      real_t tmp_23 = tmp_17*tmp_19;
      real_t tmp_24 = 0.11169079483900572*tmp_20;
      real_t tmp_25 = 0.054975871827660928*tmp_20;
      real_t tmp_26 = 0.11169079483900572*tmp_20;
      real_t tmp_27 = 0.054975871827660928*tmp_20;
      real_t tmp_28 = 0.11169079483900572*tmp_20;
      real_t tmp_29 = 2.0*tmp_7;
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_1*tmp_29;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_17*tmp_31;
      real_t tmp_34 = tmp_0*tmp_29;
      real_t tmp_35 = tmp_29*tmp_5;
      real_t tmp_36 = tmp_10*tmp_34;
      real_t tmp_37 = tmp_17*tmp_35;
      real_t a_0_0 = tmp_21*(tmp_11*tmp_12 + tmp_18*tmp_19) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_23) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_23) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_23) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_23) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_23);
      real_t a_0_1 = tmp_21*(tmp_11*tmp_30 + tmp_18*tmp_31) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_33) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_33) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_33) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_33) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_33);
      real_t a_0_2 = tmp_21*(tmp_11*tmp_34 + tmp_18*tmp_35) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_37) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_37) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_37) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_37) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_37);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = tmp_0*tmp_2;
      real_t tmp_4 = -tmp_1;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_3 - tmp_4*tmp_6);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_10 = tmp_7*(0.21132486540518713*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = tmp_7*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_1*tmp_10 + tmp_13*tmp_2;
      real_t tmp_15 = tmp_14 - 1.0/3.0;
      real_t tmp_16 = tmp_0*tmp_10 + tmp_13*tmp_5;
      real_t tmp_17 = tmp_16 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_1*(tmp_0*tmp_15 + tmp_17*tmp_4);
      real_t tmp_19 = 0.5*tmp_7;
      real_t tmp_20 = tmp_19*tmp_2;
      real_t tmp_21 = tmp_19*tmp_5;
      real_t tmp_22 = -tmp_20 - tmp_21;
      real_t tmp_23 = 0.5*tmp_22;
      real_t tmp_24 = tmp_15*tmp_6 + tmp_17*tmp_2;
      real_t tmp_25 = 1.0*tmp_7;
      real_t tmp_26 = tmp_0*tmp_25;
      real_t tmp_27 = tmp_1*tmp_25;
      real_t tmp_28 = 0.5*p_affine_10_0*tmp_22 + 0.5*p_affine_10_1*(-tmp_26 - tmp_27);
      real_t tmp_29 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_30 = 1.0 / (tmp_29);
      real_t tmp_31 = -tmp_14 - tmp_16 + 1;
      real_t tmp_32 = tmp_0*tmp_19;
      real_t tmp_33 = 0.5*p_affine_10_0*(tmp_1*tmp_32 + tmp_2*tmp_21 + tmp_20*tmp_6 + tmp_32*tmp_4) + 0.5*p_affine_10_1*(tmp_25*tmp_3 + tmp_27*tmp_6);
      real_t tmp_34 = 2*tmp_29;
      real_t tmp_35 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_34;
      real_t tmp_36 = tmp_7*(0.78867513459481287*tmp_8 + tmp_9);
      real_t tmp_37 = tmp_7*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_38 = tmp_1*tmp_36 + tmp_2*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/3.0;
      real_t tmp_40 = tmp_0*tmp_36 + tmp_37*tmp_5;
      real_t tmp_41 = tmp_40 - 1.0/3.0;
      real_t tmp_42 = p_affine_10_1*(tmp_0*tmp_39 + tmp_4*tmp_41);
      real_t tmp_43 = tmp_2*tmp_41 + tmp_39*tmp_6;
      real_t tmp_44 = -tmp_38 - tmp_40 + 1;
      real_t tmp_45 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_34;
      real_t tmp_46 = 0.25*tmp_7;
      real_t tmp_47 = tmp_2*tmp_46;
      real_t tmp_48 = 0.5*p_affine_10_0*tmp_20 + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_49 = tmp_46*tmp_5;
      real_t tmp_50 = 0.5*p_affine_10_0*tmp_21 + 0.5*p_affine_10_1*tmp_26;
      real_t a_0_0 = tmp_35*(-tmp_18*tmp_23 - tmp_24*tmp_28 + tmp_24*tmp_30*tmp_31 - tmp_31*tmp_33) + tmp_45*(-tmp_23*tmp_42 - tmp_28*tmp_43 + tmp_30*tmp_43*tmp_44 - tmp_33*tmp_44);
      real_t a_0_1 = tmp_35*(tmp_14*tmp_24*tmp_30 - tmp_14*tmp_33 - tmp_18*tmp_47 - tmp_24*tmp_48) + tmp_45*(tmp_30*tmp_38*tmp_43 - tmp_33*tmp_38 - tmp_42*tmp_47 - tmp_43*tmp_48);
      real_t a_0_2 = tmp_35*(tmp_16*tmp_24*tmp_30 - tmp_16*tmp_33 - tmp_18*tmp_49 - tmp_24*tmp_50) + tmp_45*(tmp_30*tmp_40*tmp_43 - tmp_33*tmp_40 - tmp_42*tmp_49 - tmp_43*tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElementInner[0]( 0 );
      const auto p_affine_0_1 = coordsElementInner[0]( 1 );

      const auto p_affine_1_0 = coordsElementInner[1]( 0 );
      const auto p_affine_1_1 = coordsElementInner[1]( 1 );

      const auto p_affine_2_0 = coordsElementInner[2]( 0 );
      const auto p_affine_2_1 = coordsElementInner[2]( 1 );

      const auto p_affine_3_0 = coordsElementOuter[0]( 0 );
      const auto p_affine_3_1 = coordsElementOuter[0]( 1 );

      const auto p_affine_4_0 = coordsElementOuter[1]( 0 );
      const auto p_affine_4_1 = coordsElementOuter[1]( 1 );

      const auto p_affine_5_0 = coordsElementOuter[2]( 0 );
      const auto p_affine_5_1 = coordsElementOuter[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = tmp_0*tmp_2;
      real_t tmp_4 = -tmp_1;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_3 - tmp_4*tmp_6);
      real_t tmp_8 = -p_affine_0_1;
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.21132486540518713*tmp_9;
      real_t tmp_11 = tmp_7*(tmp_10 + tmp_8);
      real_t tmp_12 = -p_affine_0_0;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + 0.21132486540518713*tmp_13;
      real_t tmp_15 = tmp_7*(tmp_12 + tmp_14);
      real_t tmp_16 = tmp_1*tmp_11 + tmp_15*tmp_2 - 1.0/3.0;
      real_t tmp_17 = tmp_0*tmp_11 + tmp_15*tmp_5 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_1*(tmp_0*tmp_16 + tmp_17*tmp_4);
      real_t tmp_19 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_20 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_21 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_22 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_23 = 1.0 / (tmp_19*tmp_20 - tmp_21*tmp_22);
      real_t tmp_24 = 0.5*tmp_23;
      real_t tmp_25 = tmp_19*tmp_24;
      real_t tmp_26 = tmp_22*tmp_24;
      real_t tmp_27 = -tmp_25 - tmp_26;
      real_t tmp_28 = 0.5*tmp_27;
      real_t tmp_29 = tmp_16*tmp_6 + tmp_17*tmp_2;
      real_t tmp_30 = 1.0*tmp_23;
      real_t tmp_31 = tmp_20*tmp_30;
      real_t tmp_32 = tmp_21*tmp_30;
      real_t tmp_33 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*(-tmp_31 - tmp_32);
      real_t tmp_34 = -p_affine_3_1;
      real_t tmp_35 = tmp_23*(tmp_10 + tmp_34);
      real_t tmp_36 = -p_affine_3_0;
      real_t tmp_37 = tmp_23*(tmp_14 + tmp_36);
      real_t tmp_38 = tmp_19*tmp_37 + tmp_21*tmp_35;
      real_t tmp_39 = tmp_20*tmp_35 + tmp_22*tmp_37;
      real_t tmp_40 = -tmp_38 - tmp_39 + 1;
      real_t tmp_41 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_42 = 1.0 / (tmp_41);
      real_t tmp_43 = tmp_29*tmp_42;
      real_t tmp_44 = 1.0*tmp_7;
      real_t tmp_45 = 0.5*tmp_7;
      real_t tmp_46 = tmp_0*tmp_45;
      real_t tmp_47 = tmp_2*tmp_45;
      real_t tmp_48 = p_affine_10_0*(tmp_1*tmp_46 + tmp_4*tmp_46 + tmp_47*tmp_5 + tmp_47*tmp_6) + p_affine_10_1*(tmp_1*tmp_44*tmp_6 + tmp_3*tmp_44);
      real_t tmp_49 = 2*tmp_41;
      real_t tmp_50 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_49;
      real_t tmp_51 = p_affine_6_1 + 0.78867513459481287*tmp_9;
      real_t tmp_52 = tmp_7*(tmp_51 + tmp_8);
      real_t tmp_53 = p_affine_6_0 + 0.78867513459481287*tmp_13;
      real_t tmp_54 = tmp_7*(tmp_12 + tmp_53);
      real_t tmp_55 = tmp_1*tmp_52 + tmp_2*tmp_54 - 1.0/3.0;
      real_t tmp_56 = tmp_0*tmp_52 + tmp_5*tmp_54 - 1.0/3.0;
      real_t tmp_57 = p_affine_10_1*(tmp_0*tmp_55 + tmp_4*tmp_56);
      real_t tmp_58 = tmp_2*tmp_56 + tmp_55*tmp_6;
      real_t tmp_59 = tmp_23*(tmp_34 + tmp_51);
      real_t tmp_60 = tmp_23*(tmp_36 + tmp_53);
      real_t tmp_61 = tmp_19*tmp_60 + tmp_21*tmp_59;
      real_t tmp_62 = tmp_20*tmp_59 + tmp_22*tmp_60;
      real_t tmp_63 = -tmp_61 - tmp_62 + 1;
      real_t tmp_64 = tmp_42*tmp_58;
      real_t tmp_65 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_49;
      real_t tmp_66 = 0.25*tmp_23;
      real_t tmp_67 = tmp_19*tmp_66;
      real_t tmp_68 = 0.5*p_affine_10_0*tmp_25 + 0.5*p_affine_10_1*tmp_32;
      real_t tmp_69 = tmp_22*tmp_66;
      real_t tmp_70 = 0.5*p_affine_10_0*tmp_26 + 0.5*p_affine_10_1*tmp_31;
      real_t a_0_0 = tmp_50*(-tmp_18*tmp_28 - tmp_29*tmp_33 - tmp_40*tmp_43 + 0.5*tmp_40*tmp_48) + tmp_65*(-tmp_28*tmp_57 - tmp_33*tmp_58 + 0.5*tmp_48*tmp_63 - tmp_63*tmp_64);
      real_t a_0_1 = tmp_50*(-tmp_18*tmp_67 - tmp_29*tmp_68 - tmp_38*tmp_43 + 0.5*tmp_38*tmp_48) + tmp_65*(0.5*tmp_48*tmp_61 - tmp_57*tmp_67 - tmp_58*tmp_68 - tmp_61*tmp_64);
      real_t a_0_2 = tmp_50*(-tmp_18*tmp_69 - tmp_29*tmp_70 - tmp_39*tmp_43 + 0.5*tmp_39*tmp_48) + tmp_65*(0.5*tmp_48*tmp_62 - tmp_57*tmp_69 - tmp_58*tmp_70 - tmp_62*tmp_64);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = -tmp_1;
      real_t tmp_4 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_5 = -tmp_4;
      real_t tmp_6 = 1.0 / (tmp_0*tmp_2 - tmp_3*tmp_5);
      real_t tmp_7 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_8 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_9 = tmp_6*(0.21132486540518713*tmp_7 + tmp_8);
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_12 = tmp_6*(0.21132486540518713*tmp_10 + tmp_11);
      real_t tmp_13 = tmp_1*tmp_9 + tmp_12*tmp_2 - 1.0/3.0;
      real_t tmp_14 = tmp_0*tmp_9 + tmp_12*tmp_4 - 1.0/3.0;
      real_t tmp_15 = tmp_0*tmp_13 + tmp_14*tmp_3;
      real_t tmp_16 = 0.5*tmp_6;
      real_t tmp_17 = tmp_16*tmp_2;
      real_t tmp_18 = tmp_16*tmp_4;
      real_t tmp_19 = -tmp_17 - tmp_18;
      real_t tmp_20 = p_affine_10_1*tmp_19;
      real_t tmp_21 = 1.0*tmp_6;
      real_t tmp_22 = tmp_0*tmp_21;
      real_t tmp_23 = tmp_1*tmp_21;
      real_t tmp_24 = p_affine_10_0*tmp_19 + p_affine_10_1*(-tmp_22 - tmp_23);
      real_t tmp_25 = tmp_13*tmp_5 + tmp_14*tmp_2;
      real_t tmp_26 = tmp_6*(0.78867513459481287*tmp_7 + tmp_8);
      real_t tmp_27 = tmp_6*(0.78867513459481287*tmp_10 + tmp_11);
      real_t tmp_28 = tmp_1*tmp_26 + tmp_2*tmp_27 - 1.0/3.0;
      real_t tmp_29 = tmp_0*tmp_26 + tmp_27*tmp_4 - 1.0/3.0;
      real_t tmp_30 = tmp_0*tmp_28 + tmp_29*tmp_3;
      real_t tmp_31 = tmp_2*tmp_29 + tmp_28*tmp_5;
      real_t tmp_32 = p_affine_10_1*tmp_17;
      real_t tmp_33 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_23;
      real_t tmp_34 = p_affine_10_1*tmp_18;
      real_t tmp_35 = p_affine_10_0*tmp_18 + p_affine_10_1*tmp_22;
      real_t a_0_0 = -0.5*tmp_15*tmp_20 - 0.5*tmp_20*tmp_30 - 0.5*tmp_24*tmp_25 - 0.5*tmp_24*tmp_31;
      real_t a_0_1 = -0.5*tmp_15*tmp_32 - 0.5*tmp_25*tmp_33 - 0.5*tmp_30*tmp_32 - 0.5*tmp_31*tmp_33;
      real_t a_0_2 = -0.5*tmp_15*tmp_34 - 0.5*tmp_25*tmp_35 - 0.5*tmp_30*tmp_34 - 0.5*tmp_31*tmp_35;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

  void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      elMat( 0, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

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

      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id6 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id7 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id8 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id9 );
      Scalar_Variable_Coefficient_3D_mu( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id10 );
      Scalar_Variable_Coefficient_3D_mu( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id11 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id12 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id13 );
      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_5 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_6 = -tmp_3 + tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_8 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_9 = tmp_2*tmp_8;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_12 = tmp_5*tmp_8;
      real_t tmp_13 = tmp_11*tmp_4;
      real_t tmp_14 = 1.0 / (-tmp_0*tmp_3 + tmp_0*tmp_4*tmp_5 + tmp_1*tmp_10*tmp_11 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_7*tmp_9);
      real_t tmp_15 = 1.0*tmp_14;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = -tmp_13 + tmp_9;
      real_t tmp_18 = tmp_15*tmp_17;
      real_t tmp_19 = tmp_1*tmp_11 - tmp_12;
      real_t tmp_20 = tmp_15*tmp_19;
      real_t tmp_21 = tmp_0*tmp_16 + tmp_10*tmp_20 + tmp_18*tmp_7;
      real_t tmp_22 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_21;
      real_t tmp_23 = -2*tmp_16 - 2*tmp_18 - 2*tmp_20;
      real_t tmp_24 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_25 = -tmp_0*tmp_1 + tmp_7*tmp_8;
      real_t tmp_26 = 0.5*tmp_14;
      real_t tmp_27 = tmp_25*tmp_26;
      real_t tmp_28 = tmp_0*tmp_4 - tmp_10*tmp_8;
      real_t tmp_29 = tmp_26*tmp_28;
      real_t tmp_30 = tmp_1*tmp_10 - tmp_4*tmp_7;
      real_t tmp_31 = tmp_26*tmp_30;
      real_t tmp_32 = tmp_26*tmp_6;
      real_t tmp_33 = tmp_17*tmp_26;
      real_t tmp_34 = tmp_19*tmp_26;
      real_t tmp_35 = tmp_0*tmp_31 + tmp_10*tmp_27 + tmp_11*tmp_32 + tmp_2*tmp_34 + tmp_29*tmp_7 + tmp_33*tmp_5;
      real_t tmp_36 = tmp_35*(-tmp_27 - tmp_29 - tmp_31);
      real_t tmp_37 = tmp_0*tmp_5 - tmp_11*tmp_7;
      real_t tmp_38 = tmp_26*tmp_37;
      real_t tmp_39 = -tmp_0*tmp_2 + tmp_10*tmp_11;
      real_t tmp_40 = tmp_26*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_5 + tmp_2*tmp_7;
      real_t tmp_42 = tmp_26*tmp_41;
      real_t tmp_43 = tmp_0*tmp_42 + tmp_1*tmp_33 + tmp_10*tmp_38 + tmp_32*tmp_8 + tmp_34*tmp_4 + tmp_40*tmp_7;
      real_t tmp_44 = tmp_43*(-tmp_38 - tmp_40 - tmp_42);
      real_t tmp_45 = p_affine_0_0*p_affine_1_1;
      real_t tmp_46 = p_affine_0_0*p_affine_1_2;
      real_t tmp_47 = p_affine_2_1*p_affine_3_2;
      real_t tmp_48 = p_affine_0_1*p_affine_1_0;
      real_t tmp_49 = p_affine_0_1*p_affine_1_2;
      real_t tmp_50 = p_affine_2_2*p_affine_3_0;
      real_t tmp_51 = p_affine_0_2*p_affine_1_0;
      real_t tmp_52 = p_affine_0_2*p_affine_1_1;
      real_t tmp_53 = p_affine_2_0*p_affine_3_1;
      real_t tmp_54 = p_affine_2_2*p_affine_3_1;
      real_t tmp_55 = p_affine_2_0*p_affine_3_2;
      real_t tmp_56 = p_affine_2_1*p_affine_3_0;
      real_t tmp_57 = std::abs(p_affine_0_0*tmp_47 - p_affine_0_0*tmp_54 + p_affine_0_1*tmp_50 - p_affine_0_1*tmp_55 + p_affine_0_2*tmp_53 - p_affine_0_2*tmp_56 - p_affine_1_0*tmp_47 + p_affine_1_0*tmp_54 - p_affine_1_1*tmp_50 + p_affine_1_1*tmp_55 - p_affine_1_2*tmp_53 + p_affine_1_2*tmp_56 + p_affine_2_0*tmp_49 - p_affine_2_0*tmp_52 - p_affine_2_1*tmp_46 + p_affine_2_1*tmp_51 + p_affine_2_2*tmp_45 - p_affine_2_2*tmp_48 - p_affine_3_0*tmp_49 + p_affine_3_0*tmp_52 + p_affine_3_1*tmp_46 - p_affine_3_1*tmp_51 - p_affine_3_2*tmp_45 + p_affine_3_2*tmp_48);
      real_t tmp_58 = 0.018781320953002646*tmp_57;
      real_t tmp_59 = tmp_21*tmp_23;
      real_t tmp_60 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id1;
      real_t tmp_61 = 0.012248840519393657*tmp_57;
      real_t tmp_62 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id2;
      real_t tmp_63 = 0.0070910034628469103*tmp_57;
      real_t tmp_64 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id3;
      real_t tmp_65 = 0.0070910034628469103*tmp_57;
      real_t tmp_66 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id4;
      real_t tmp_67 = 0.0070910034628469103*tmp_57;
      real_t tmp_68 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id5;
      real_t tmp_69 = 0.0070910034628469103*tmp_57;
      real_t tmp_70 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id6;
      real_t tmp_71 = 0.018781320953002646*tmp_57;
      real_t tmp_72 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id7;
      real_t tmp_73 = 0.012248840519393657*tmp_57;
      real_t tmp_74 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id8;
      real_t tmp_75 = 0.018781320953002646*tmp_57;
      real_t tmp_76 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id9;
      real_t tmp_77 = 0.012248840519393657*tmp_57;
      real_t tmp_78 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id10;
      real_t tmp_79 = 0.018781320953002646*tmp_57;
      real_t tmp_80 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id11;
      real_t tmp_81 = 0.012248840519393657*tmp_57;
      real_t tmp_82 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id12;
      real_t tmp_83 = 0.0070910034628469103*tmp_57;
      real_t tmp_84 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id13;
      real_t tmp_85 = 0.0070910034628469103*tmp_57;
      real_t tmp_86 = 2.0*tmp_14;
      real_t tmp_87 = tmp_6*tmp_86;
      real_t tmp_88 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_86;
      real_t tmp_89 = tmp_30*tmp_35;
      real_t tmp_90 = tmp_41*tmp_43;
      real_t tmp_91 = tmp_21*tmp_87;
      real_t tmp_92 = Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_86;
      real_t tmp_93 = Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_86;
      real_t tmp_94 = Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_86;
      real_t tmp_95 = Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_86;
      real_t tmp_96 = Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_86;
      real_t tmp_97 = Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_86;
      real_t tmp_98 = Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_86;
      real_t tmp_99 = Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_86;
      real_t tmp_100 = Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_86;
      real_t tmp_101 = Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_86;
      real_t tmp_102 = Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_86;
      real_t tmp_103 = Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_86;
      real_t tmp_104 = Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_86;
      real_t tmp_105 = tmp_22*tmp_86;
      real_t tmp_106 = tmp_28*tmp_35;
      real_t tmp_107 = tmp_39*tmp_43;
      real_t tmp_108 = tmp_17*tmp_21;
      real_t tmp_109 = tmp_25*tmp_35;
      real_t tmp_110 = tmp_37*tmp_43;
      real_t tmp_111 = tmp_19*tmp_21;
      real_t a_0_0 = tmp_58*(tmp_22*tmp_23 + tmp_24*tmp_36 + tmp_24*tmp_44) + tmp_61*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_59 + tmp_36*tmp_60 + tmp_44*tmp_60) + tmp_63*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_59 + tmp_36*tmp_62 + tmp_44*tmp_62) + tmp_65*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_59 + tmp_36*tmp_64 + tmp_44*tmp_64) + tmp_67*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_59 + tmp_36*tmp_66 + tmp_44*tmp_66) + tmp_69*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_59 + tmp_36*tmp_68 + tmp_44*tmp_68) + tmp_71*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_59 + tmp_36*tmp_70 + tmp_44*tmp_70) + tmp_73*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_59 + tmp_36*tmp_72 + tmp_44*tmp_72) + tmp_75*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_59 + tmp_36*tmp_74 + tmp_44*tmp_74) + tmp_77*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_59 + tmp_36*tmp_76 + tmp_44*tmp_76) + tmp_79*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_59 + tmp_36*tmp_78 + tmp_44*tmp_78) + tmp_81*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_59 + tmp_36*tmp_80 + tmp_44*tmp_80) + tmp_83*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_59 + tmp_36*tmp_82 + tmp_44*tmp_82) + tmp_85*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_59 + tmp_36*tmp_84 + tmp_44*tmp_84);
      real_t a_0_1 = tmp_58*(tmp_22*tmp_87 + tmp_88*tmp_89 + tmp_88*tmp_90) + tmp_61*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_91 + tmp_89*tmp_92 + tmp_90*tmp_92) + tmp_63*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_91 + tmp_89*tmp_93 + tmp_90*tmp_93) + tmp_65*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_91 + tmp_89*tmp_94 + tmp_90*tmp_94) + tmp_67*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_91 + tmp_89*tmp_95 + tmp_90*tmp_95) + tmp_69*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_91 + tmp_89*tmp_96 + tmp_90*tmp_96) + tmp_71*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_91 + tmp_89*tmp_97 + tmp_90*tmp_97) + tmp_73*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_91 + tmp_89*tmp_98 + tmp_90*tmp_98) + tmp_75*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_91 + tmp_89*tmp_99 + tmp_90*tmp_99) + tmp_77*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_91 + tmp_100*tmp_89 + tmp_100*tmp_90) + tmp_79*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_91 + tmp_101*tmp_89 + tmp_101*tmp_90) + tmp_81*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_91 + tmp_102*tmp_89 + tmp_102*tmp_90) + tmp_83*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_91 + tmp_103*tmp_89 + tmp_103*tmp_90) + tmp_85*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_91 + tmp_104*tmp_89 + tmp_104*tmp_90);
      real_t a_0_2 = tmp_58*(tmp_105*tmp_17 + tmp_106*tmp_88 + tmp_107*tmp_88) + tmp_61*(tmp_106*tmp_92 + tmp_107*tmp_92 + tmp_108*tmp_92) + tmp_63*(tmp_106*tmp_93 + tmp_107*tmp_93 + tmp_108*tmp_93) + tmp_65*(tmp_106*tmp_94 + tmp_107*tmp_94 + tmp_108*tmp_94) + tmp_67*(tmp_106*tmp_95 + tmp_107*tmp_95 + tmp_108*tmp_95) + tmp_69*(tmp_106*tmp_96 + tmp_107*tmp_96 + tmp_108*tmp_96) + tmp_71*(tmp_106*tmp_97 + tmp_107*tmp_97 + tmp_108*tmp_97) + tmp_73*(tmp_106*tmp_98 + tmp_107*tmp_98 + tmp_108*tmp_98) + tmp_75*(tmp_106*tmp_99 + tmp_107*tmp_99 + tmp_108*tmp_99) + tmp_77*(tmp_100*tmp_106 + tmp_100*tmp_107 + tmp_100*tmp_108) + tmp_79*(tmp_101*tmp_106 + tmp_101*tmp_107 + tmp_101*tmp_108) + tmp_81*(tmp_102*tmp_106 + tmp_102*tmp_107 + tmp_102*tmp_108) + tmp_83*(tmp_103*tmp_106 + tmp_103*tmp_107 + tmp_103*tmp_108) + tmp_85*(tmp_104*tmp_106 + tmp_104*tmp_107 + tmp_104*tmp_108);
      real_t a_0_3 = tmp_58*(tmp_105*tmp_19 + tmp_109*tmp_88 + tmp_110*tmp_88) + tmp_61*(tmp_109*tmp_92 + tmp_110*tmp_92 + tmp_111*tmp_92) + tmp_63*(tmp_109*tmp_93 + tmp_110*tmp_93 + tmp_111*tmp_93) + tmp_65*(tmp_109*tmp_94 + tmp_110*tmp_94 + tmp_111*tmp_94) + tmp_67*(tmp_109*tmp_95 + tmp_110*tmp_95 + tmp_111*tmp_95) + tmp_69*(tmp_109*tmp_96 + tmp_110*tmp_96 + tmp_111*tmp_96) + tmp_71*(tmp_109*tmp_97 + tmp_110*tmp_97 + tmp_111*tmp_97) + tmp_73*(tmp_109*tmp_98 + tmp_110*tmp_98 + tmp_111*tmp_98) + tmp_75*(tmp_109*tmp_99 + tmp_110*tmp_99 + tmp_111*tmp_99) + tmp_77*(tmp_100*tmp_109 + tmp_100*tmp_110 + tmp_100*tmp_111) + tmp_79*(tmp_101*tmp_109 + tmp_101*tmp_110 + tmp_101*tmp_111) + tmp_81*(tmp_102*tmp_109 + tmp_102*tmp_110 + tmp_102*tmp_111) + tmp_83*(tmp_103*tmp_109 + tmp_103*tmp_110 + tmp_103*tmp_111) + tmp_85*(tmp_104*tmp_109 + tmp_104*tmp_110 + tmp_104*tmp_111);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_6*tmp_9;
      real_t tmp_14 = tmp_4*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_12 + tmp_0*tmp_7 - tmp_1*tmp_13 + tmp_1*tmp_2*tmp_8 + tmp_11*tmp_3 - tmp_14*tmp_3);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_7;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = tmp_33 - 1.0/4.0;
      real_t tmp_35 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_36 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_37 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_38 = tmp_20*tmp_35 + tmp_26*tmp_36 + tmp_32*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/4.0;
      real_t tmp_40 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_41 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_42 = tmp_11 - tmp_14;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_26*tmp_41 + tmp_32*tmp_42;
      real_t tmp_44 = tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_34 + tmp_1*tmp_39 + tmp_3*tmp_44;
      real_t tmp_46 = 0.5*tmp_15;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = tmp_37*tmp_46;
      real_t tmp_49 = tmp_27*tmp_46;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = 0.5*p_affine_13_1;
      real_t tmp_52 = tmp_50*tmp_51;
      real_t tmp_53 = tmp_10*tmp_39 + tmp_34*tmp_8 + tmp_44*tmp_6;
      real_t tmp_54 = tmp_40*tmp_46;
      real_t tmp_55 = tmp_35*tmp_46;
      real_t tmp_56 = tmp_46*tmp_5;
      real_t tmp_57 = -tmp_54 - tmp_55 - tmp_56;
      real_t tmp_58 = tmp_51*tmp_57;
      real_t tmp_59 = tmp_2*tmp_44 + tmp_34*tmp_9 + tmp_39*tmp_4;
      real_t tmp_60 = 1.0*tmp_15;
      real_t tmp_61 = tmp_41*tmp_60;
      real_t tmp_62 = tmp_36*tmp_60;
      real_t tmp_63 = tmp_21*tmp_60;
      real_t tmp_64 = 0.5*p_affine_13_0*tmp_50 + 0.5*p_affine_13_1*(-tmp_61 - tmp_62 - tmp_63) + 0.5*p_affine_13_2*tmp_57;
      real_t tmp_65 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_66 = std::pow(tmp_65, -0.25);
      real_t tmp_67 = -tmp_33 - tmp_38 - tmp_43 + 1;
      real_t tmp_68 = tmp_21*tmp_46;
      real_t tmp_69 = tmp_36*tmp_46;
      real_t tmp_70 = tmp_41*tmp_46;
      real_t tmp_71 = 0.5*p_affine_13_0*(tmp_0*tmp_68 + tmp_1*tmp_69 + tmp_2*tmp_47 + tmp_3*tmp_70 + tmp_4*tmp_48 + tmp_49*tmp_9) + 0.5*p_affine_13_1*(tmp_2*tmp_61 + tmp_4*tmp_62 + tmp_63*tmp_9) + 0.5*p_affine_13_2*(tmp_10*tmp_69 + tmp_2*tmp_54 + tmp_4*tmp_55 + tmp_56*tmp_9 + tmp_6*tmp_70 + tmp_68*tmp_8);
      real_t tmp_72 = 2.0*std::pow(tmp_65, 1.0/2.0);
      real_t tmp_73 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_72;
      real_t tmp_74 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_75 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_76 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_77 = tmp_21*tmp_75 + tmp_27*tmp_76 + tmp_5*tmp_74;
      real_t tmp_78 = tmp_77 - 1.0/4.0;
      real_t tmp_79 = tmp_35*tmp_74 + tmp_36*tmp_75 + tmp_37*tmp_76;
      real_t tmp_80 = tmp_79 - 1.0/4.0;
      real_t tmp_81 = tmp_40*tmp_74 + tmp_41*tmp_75 + tmp_42*tmp_76;
      real_t tmp_82 = tmp_81 - 1.0/4.0;
      real_t tmp_83 = tmp_0*tmp_78 + tmp_1*tmp_80 + tmp_3*tmp_82;
      real_t tmp_84 = tmp_10*tmp_80 + tmp_6*tmp_82 + tmp_78*tmp_8;
      real_t tmp_85 = tmp_2*tmp_82 + tmp_4*tmp_80 + tmp_78*tmp_9;
      real_t tmp_86 = -tmp_77 - tmp_79 - tmp_81 + 1;
      real_t tmp_87 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_72;
      real_t tmp_88 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_89 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_90 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_91 = tmp_21*tmp_89 + tmp_27*tmp_90 + tmp_5*tmp_88;
      real_t tmp_92 = tmp_91 - 1.0/4.0;
      real_t tmp_93 = tmp_35*tmp_88 + tmp_36*tmp_89 + tmp_37*tmp_90;
      real_t tmp_94 = tmp_93 - 1.0/4.0;
      real_t tmp_95 = tmp_40*tmp_88 + tmp_41*tmp_89 + tmp_42*tmp_90;
      real_t tmp_96 = tmp_95 - 1.0/4.0;
      real_t tmp_97 = tmp_0*tmp_92 + tmp_1*tmp_94 + tmp_3*tmp_96;
      real_t tmp_98 = tmp_10*tmp_94 + tmp_6*tmp_96 + tmp_8*tmp_92;
      real_t tmp_99 = tmp_2*tmp_96 + tmp_4*tmp_94 + tmp_9*tmp_92;
      real_t tmp_100 = -tmp_91 - tmp_93 - tmp_95 + 1;
      real_t tmp_101 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_72;
      real_t tmp_102 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_103 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_104 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_105 = tmp_102*tmp_5 + tmp_103*tmp_21 + tmp_104*tmp_27;
      real_t tmp_106 = tmp_105 - 1.0/4.0;
      real_t tmp_107 = tmp_102*tmp_35 + tmp_103*tmp_36 + tmp_104*tmp_37;
      real_t tmp_108 = tmp_107 - 1.0/4.0;
      real_t tmp_109 = tmp_102*tmp_40 + tmp_103*tmp_41 + tmp_104*tmp_42;
      real_t tmp_110 = tmp_109 - 1.0/4.0;
      real_t tmp_111 = tmp_0*tmp_106 + tmp_1*tmp_108 + tmp_110*tmp_3;
      real_t tmp_112 = tmp_10*tmp_108 + tmp_106*tmp_8 + tmp_110*tmp_6;
      real_t tmp_113 = tmp_106*tmp_9 + tmp_108*tmp_4 + tmp_110*tmp_2;
      real_t tmp_114 = -tmp_105 - tmp_107 - tmp_109 + 1;
      real_t tmp_115 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_72;
      real_t tmp_116 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_117 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_118 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_119 = tmp_116*tmp_5 + tmp_117*tmp_21 + tmp_118*tmp_27;
      real_t tmp_120 = tmp_119 - 1.0/4.0;
      real_t tmp_121 = tmp_116*tmp_35 + tmp_117*tmp_36 + tmp_118*tmp_37;
      real_t tmp_122 = tmp_121 - 1.0/4.0;
      real_t tmp_123 = tmp_116*tmp_40 + tmp_117*tmp_41 + tmp_118*tmp_42;
      real_t tmp_124 = tmp_123 - 1.0/4.0;
      real_t tmp_125 = tmp_0*tmp_120 + tmp_1*tmp_122 + tmp_124*tmp_3;
      real_t tmp_126 = tmp_10*tmp_122 + tmp_120*tmp_8 + tmp_124*tmp_6;
      real_t tmp_127 = tmp_120*tmp_9 + tmp_122*tmp_4 + tmp_124*tmp_2;
      real_t tmp_128 = -tmp_119 - tmp_121 - tmp_123 + 1;
      real_t tmp_129 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_72;
      real_t tmp_130 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_131 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_132 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_133 = tmp_130*tmp_5 + tmp_131*tmp_21 + tmp_132*tmp_27;
      real_t tmp_134 = tmp_133 - 1.0/4.0;
      real_t tmp_135 = tmp_130*tmp_35 + tmp_131*tmp_36 + tmp_132*tmp_37;
      real_t tmp_136 = tmp_135 - 1.0/4.0;
      real_t tmp_137 = tmp_130*tmp_40 + tmp_131*tmp_41 + tmp_132*tmp_42;
      real_t tmp_138 = tmp_137 - 1.0/4.0;
      real_t tmp_139 = tmp_0*tmp_134 + tmp_1*tmp_136 + tmp_138*tmp_3;
      real_t tmp_140 = tmp_10*tmp_136 + tmp_134*tmp_8 + tmp_138*tmp_6;
      real_t tmp_141 = tmp_134*tmp_9 + tmp_136*tmp_4 + tmp_138*tmp_2;
      real_t tmp_142 = -tmp_133 - tmp_135 - tmp_137 + 1;
      real_t tmp_143 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_72;
      real_t tmp_144 = 0.25*p_affine_13_1*tmp_15;
      real_t tmp_145 = tmp_144*tmp_5;
      real_t tmp_146 = tmp_144*tmp_27;
      real_t tmp_147 = 0.5*p_affine_13_0*tmp_49 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_56;
      real_t tmp_148 = tmp_144*tmp_35;
      real_t tmp_149 = tmp_144*tmp_37;
      real_t tmp_150 = 0.5*p_affine_13_0*tmp_48 + 0.5*p_affine_13_1*tmp_62 + 0.5*p_affine_13_2*tmp_55;
      real_t tmp_151 = tmp_144*tmp_40;
      real_t tmp_152 = tmp_144*tmp_42;
      real_t tmp_153 = 0.5*p_affine_13_0*tmp_47 + 0.5*p_affine_13_1*tmp_61 + 0.5*p_affine_13_2*tmp_54;
      real_t a_0_0 = tmp_101*(1.0*tmp_100*tmp_66*tmp_99 - tmp_100*tmp_71 - tmp_52*tmp_97 - tmp_58*tmp_98 - tmp_64*tmp_99) + tmp_115*(-tmp_111*tmp_52 - tmp_112*tmp_58 + 1.0*tmp_113*tmp_114*tmp_66 - tmp_113*tmp_64 - tmp_114*tmp_71) + tmp_129*(-tmp_125*tmp_52 - tmp_126*tmp_58 + 1.0*tmp_127*tmp_128*tmp_66 - tmp_127*tmp_64 - tmp_128*tmp_71) + tmp_143*(-tmp_139*tmp_52 - tmp_140*tmp_58 + 1.0*tmp_141*tmp_142*tmp_66 - tmp_141*tmp_64 - tmp_142*tmp_71) + tmp_73*(-tmp_45*tmp_52 - tmp_53*tmp_58 - tmp_59*tmp_64 + 1.0*tmp_59*tmp_66*tmp_67 - tmp_67*tmp_71) + tmp_87*(-tmp_52*tmp_83 - tmp_58*tmp_84 - tmp_64*tmp_85 + 1.0*tmp_66*tmp_85*tmp_86 - tmp_71*tmp_86);
      real_t a_0_1 = tmp_101*(-tmp_145*tmp_98 - tmp_146*tmp_97 - tmp_147*tmp_99 + 1.0*tmp_66*tmp_91*tmp_99 - tmp_71*tmp_91) + tmp_115*(1.0*tmp_105*tmp_113*tmp_66 - tmp_105*tmp_71 - tmp_111*tmp_146 - tmp_112*tmp_145 - tmp_113*tmp_147) + tmp_129*(1.0*tmp_119*tmp_127*tmp_66 - tmp_119*tmp_71 - tmp_125*tmp_146 - tmp_126*tmp_145 - tmp_127*tmp_147) + tmp_143*(1.0*tmp_133*tmp_141*tmp_66 - tmp_133*tmp_71 - tmp_139*tmp_146 - tmp_140*tmp_145 - tmp_141*tmp_147) + tmp_73*(-tmp_145*tmp_53 - tmp_146*tmp_45 - tmp_147*tmp_59 + 1.0*tmp_33*tmp_59*tmp_66 - tmp_33*tmp_71) + tmp_87*(-tmp_145*tmp_84 - tmp_146*tmp_83 - tmp_147*tmp_85 + 1.0*tmp_66*tmp_77*tmp_85 - tmp_71*tmp_77);
      real_t a_0_2 = tmp_101*(-tmp_148*tmp_98 - tmp_149*tmp_97 - tmp_150*tmp_99 + 1.0*tmp_66*tmp_93*tmp_99 - tmp_71*tmp_93) + tmp_115*(1.0*tmp_107*tmp_113*tmp_66 - tmp_107*tmp_71 - tmp_111*tmp_149 - tmp_112*tmp_148 - tmp_113*tmp_150) + tmp_129*(1.0*tmp_121*tmp_127*tmp_66 - tmp_121*tmp_71 - tmp_125*tmp_149 - tmp_126*tmp_148 - tmp_127*tmp_150) + tmp_143*(1.0*tmp_135*tmp_141*tmp_66 - tmp_135*tmp_71 - tmp_139*tmp_149 - tmp_140*tmp_148 - tmp_141*tmp_150) + tmp_73*(-tmp_148*tmp_53 - tmp_149*tmp_45 - tmp_150*tmp_59 + 1.0*tmp_38*tmp_59*tmp_66 - tmp_38*tmp_71) + tmp_87*(-tmp_148*tmp_84 - tmp_149*tmp_83 - tmp_150*tmp_85 + 1.0*tmp_66*tmp_79*tmp_85 - tmp_71*tmp_79);
      real_t a_0_3 = tmp_101*(-tmp_151*tmp_98 - tmp_152*tmp_97 - tmp_153*tmp_99 + 1.0*tmp_66*tmp_95*tmp_99 - tmp_71*tmp_95) + tmp_115*(1.0*tmp_109*tmp_113*tmp_66 - tmp_109*tmp_71 - tmp_111*tmp_152 - tmp_112*tmp_151 - tmp_113*tmp_153) + tmp_129*(1.0*tmp_123*tmp_127*tmp_66 - tmp_123*tmp_71 - tmp_125*tmp_152 - tmp_126*tmp_151 - tmp_127*tmp_153) + tmp_143*(1.0*tmp_137*tmp_141*tmp_66 - tmp_137*tmp_71 - tmp_139*tmp_152 - tmp_140*tmp_151 - tmp_141*tmp_153) + tmp_73*(-tmp_151*tmp_53 - tmp_152*tmp_45 - tmp_153*tmp_59 + 1.0*tmp_43*tmp_59*tmp_66 - tmp_43*tmp_71) + tmp_87*(-tmp_151*tmp_84 - tmp_152*tmp_83 - tmp_153*tmp_85 + 1.0*tmp_66*tmp_81*tmp_85 - tmp_71*tmp_81);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_2;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_3*tmp_8;
      real_t tmp_14 = 1.0 / (-tmp_0*tmp_11 + tmp_0*tmp_7 + tmp_1*tmp_2*tmp_8 + tmp_10*tmp_3*tmp_9 - tmp_12*tmp_9 - tmp_13*tmp_4);
      real_t tmp_15 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_16 = -tmp_15;
      real_t tmp_17 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_18 = 0.091576213509770743*tmp_16 + 0.81684757298045851*tmp_17;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_14*(tmp_18 + tmp_19);
      real_t tmp_21 = tmp_10*tmp_3 - tmp_12;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = 0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24;
      real_t tmp_26 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_27 = tmp_14*(tmp_25 + tmp_26);
      real_t tmp_28 = -tmp_11 + tmp_7;
      real_t tmp_29 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_30 = -tmp_29;
      real_t tmp_31 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_32 = 0.091576213509770743*tmp_30 + 0.81684757298045851*tmp_31;
      real_t tmp_33 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_34 = tmp_14*(tmp_32 + tmp_33);
      real_t tmp_35 = tmp_20*tmp_5 + tmp_21*tmp_27 + tmp_28*tmp_34 - 1.0/4.0;
      real_t tmp_36 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_37 = tmp_0*tmp_6 - tmp_13;
      real_t tmp_38 = tmp_2*tmp_8 - tmp_6*tmp_9;
      real_t tmp_39 = tmp_20*tmp_36 + tmp_27*tmp_37 + tmp_34*tmp_38 - 1.0/4.0;
      real_t tmp_40 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_41 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_42 = tmp_10*tmp_9 - tmp_4*tmp_8;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_27*tmp_41 + tmp_34*tmp_42 - 1.0/4.0;
      real_t tmp_44 = tmp_0*tmp_35 + tmp_1*tmp_39 + tmp_3*tmp_43;
      real_t tmp_45 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_46 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_47 = tmp_45*tmp_46;
      real_t tmp_48 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_49 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_50 = tmp_47 - tmp_48*tmp_49;
      real_t tmp_51 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_52 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_53 = tmp_48*tmp_52;
      real_t tmp_54 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_55 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_56 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_57 = tmp_52*tmp_54;
      real_t tmp_58 = tmp_45*tmp_55;
      real_t tmp_59 = tmp_49*tmp_56;
      real_t tmp_60 = 1.0 / (-tmp_46*tmp_57 + tmp_47*tmp_56 - tmp_48*tmp_59 + tmp_49*tmp_54*tmp_55 + tmp_51*tmp_53 - tmp_51*tmp_58);
      real_t tmp_61 = 0.5*tmp_60;
      real_t tmp_62 = tmp_50*tmp_61;
      real_t tmp_63 = -tmp_45*tmp_51 + tmp_49*tmp_54;
      real_t tmp_64 = tmp_61*tmp_63;
      real_t tmp_65 = -tmp_46*tmp_54 + tmp_48*tmp_51;
      real_t tmp_66 = tmp_61*tmp_65;
      real_t tmp_67 = -tmp_62 - tmp_64 - tmp_66;
      real_t tmp_68 = 0.5*p_affine_13_1;
      real_t tmp_69 = tmp_67*tmp_68;
      real_t tmp_70 = tmp_10*tmp_39 + tmp_35*tmp_8 + tmp_43*tmp_6;
      real_t tmp_71 = tmp_53 - tmp_58;
      real_t tmp_72 = tmp_61*tmp_71;
      real_t tmp_73 = tmp_45*tmp_56 - tmp_57;
      real_t tmp_74 = tmp_61*tmp_73;
      real_t tmp_75 = -tmp_48*tmp_56 + tmp_54*tmp_55;
      real_t tmp_76 = tmp_61*tmp_75;
      real_t tmp_77 = -tmp_72 - tmp_74 - tmp_76;
      real_t tmp_78 = tmp_68*tmp_77;
      real_t tmp_79 = tmp_2*tmp_43 + tmp_35*tmp_9 + tmp_39*tmp_4;
      real_t tmp_80 = -tmp_46*tmp_52 + tmp_49*tmp_55;
      real_t tmp_81 = 1.0*tmp_60;
      real_t tmp_82 = tmp_80*tmp_81;
      real_t tmp_83 = tmp_51*tmp_52 - tmp_59;
      real_t tmp_84 = tmp_81*tmp_83;
      real_t tmp_85 = tmp_46*tmp_56 - tmp_51*tmp_55;
      real_t tmp_86 = tmp_81*tmp_85;
      real_t tmp_87 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*(-tmp_82 - tmp_84 - tmp_86) + 0.5*p_affine_13_2*tmp_77;
      real_t tmp_88 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_89 = tmp_60*(tmp_18 + tmp_88);
      real_t tmp_90 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_91 = tmp_60*(tmp_25 + tmp_90);
      real_t tmp_92 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_93 = tmp_60*(tmp_32 + tmp_92);
      real_t tmp_94 = tmp_65*tmp_93 + tmp_75*tmp_89 + tmp_85*tmp_91;
      real_t tmp_95 = tmp_63*tmp_93 + tmp_73*tmp_89 + tmp_83*tmp_91;
      real_t tmp_96 = tmp_50*tmp_93 + tmp_71*tmp_89 + tmp_80*tmp_91;
      real_t tmp_97 = -tmp_94 - tmp_95 - tmp_96 + 1;
      real_t tmp_98 = (std::abs(tmp_15*tmp_24 - tmp_17*tmp_22)*std::abs(tmp_15*tmp_24 - tmp_17*tmp_22)) + (std::abs(tmp_15*tmp_31 - tmp_17*tmp_29)*std::abs(tmp_15*tmp_31 - tmp_17*tmp_29)) + (std::abs(tmp_22*tmp_31 - tmp_24*tmp_29)*std::abs(tmp_22*tmp_31 - tmp_24*tmp_29));
      real_t tmp_99 = 1.0*std::pow(tmp_98, -0.25);
      real_t tmp_100 = tmp_79*tmp_99;
      real_t tmp_101 = 1.0*tmp_14;
      real_t tmp_102 = 0.5*tmp_14;
      real_t tmp_103 = tmp_102*tmp_21;
      real_t tmp_104 = tmp_102*tmp_37;
      real_t tmp_105 = tmp_102*tmp_41;
      real_t tmp_106 = tmp_102*tmp_9;
      real_t tmp_107 = tmp_102*tmp_4;
      real_t tmp_108 = tmp_102*tmp_2;
      real_t tmp_109 = p_affine_13_0*(tmp_0*tmp_103 + tmp_1*tmp_104 + tmp_105*tmp_3 + tmp_106*tmp_28 + tmp_107*tmp_38 + tmp_108*tmp_42) + p_affine_13_1*(tmp_101*tmp_2*tmp_41 + tmp_101*tmp_21*tmp_9 + tmp_101*tmp_37*tmp_4) + p_affine_13_2*(tmp_10*tmp_104 + tmp_103*tmp_8 + tmp_105*tmp_6 + tmp_106*tmp_5 + tmp_107*tmp_36 + tmp_108*tmp_40);
      real_t tmp_110 = 2.0*std::pow(tmp_98, 1.0/2.0);
      real_t tmp_111 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_110;
      real_t tmp_112 = 0.44594849091596489*tmp_16 + 0.10810301816807022*tmp_17;
      real_t tmp_113 = tmp_14*(tmp_112 + tmp_19);
      real_t tmp_114 = 0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24;
      real_t tmp_115 = tmp_14*(tmp_114 + tmp_26);
      real_t tmp_116 = 0.44594849091596489*tmp_30 + 0.10810301816807022*tmp_31;
      real_t tmp_117 = tmp_14*(tmp_116 + tmp_33);
      real_t tmp_118 = tmp_113*tmp_5 + tmp_115*tmp_21 + tmp_117*tmp_28 - 1.0/4.0;
      real_t tmp_119 = tmp_113*tmp_36 + tmp_115*tmp_37 + tmp_117*tmp_38 - 1.0/4.0;
      real_t tmp_120 = tmp_113*tmp_40 + tmp_115*tmp_41 + tmp_117*tmp_42 - 1.0/4.0;
      real_t tmp_121 = tmp_0*tmp_118 + tmp_1*tmp_119 + tmp_120*tmp_3;
      real_t tmp_122 = tmp_10*tmp_119 + tmp_118*tmp_8 + tmp_120*tmp_6;
      real_t tmp_123 = tmp_118*tmp_9 + tmp_119*tmp_4 + tmp_120*tmp_2;
      real_t tmp_124 = tmp_60*(tmp_112 + tmp_88);
      real_t tmp_125 = tmp_60*(tmp_114 + tmp_90);
      real_t tmp_126 = tmp_60*(tmp_116 + tmp_92);
      real_t tmp_127 = tmp_124*tmp_75 + tmp_125*tmp_85 + tmp_126*tmp_65;
      real_t tmp_128 = tmp_124*tmp_73 + tmp_125*tmp_83 + tmp_126*tmp_63;
      real_t tmp_129 = tmp_124*tmp_71 + tmp_125*tmp_80 + tmp_126*tmp_50;
      real_t tmp_130 = -tmp_127 - tmp_128 - tmp_129 + 1;
      real_t tmp_131 = tmp_123*tmp_99;
      real_t tmp_132 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_110;
      real_t tmp_133 = 0.81684757298045851*tmp_16 + 0.091576213509770743*tmp_17;
      real_t tmp_134 = tmp_14*(tmp_133 + tmp_19);
      real_t tmp_135 = 0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24;
      real_t tmp_136 = tmp_14*(tmp_135 + tmp_26);
      real_t tmp_137 = 0.81684757298045851*tmp_30 + 0.091576213509770743*tmp_31;
      real_t tmp_138 = tmp_14*(tmp_137 + tmp_33);
      real_t tmp_139 = tmp_134*tmp_5 + tmp_136*tmp_21 + tmp_138*tmp_28 - 1.0/4.0;
      real_t tmp_140 = tmp_134*tmp_36 + tmp_136*tmp_37 + tmp_138*tmp_38 - 1.0/4.0;
      real_t tmp_141 = tmp_134*tmp_40 + tmp_136*tmp_41 + tmp_138*tmp_42 - 1.0/4.0;
      real_t tmp_142 = tmp_0*tmp_139 + tmp_1*tmp_140 + tmp_141*tmp_3;
      real_t tmp_143 = tmp_10*tmp_140 + tmp_139*tmp_8 + tmp_141*tmp_6;
      real_t tmp_144 = tmp_139*tmp_9 + tmp_140*tmp_4 + tmp_141*tmp_2;
      real_t tmp_145 = tmp_60*(tmp_133 + tmp_88);
      real_t tmp_146 = tmp_60*(tmp_135 + tmp_90);
      real_t tmp_147 = tmp_60*(tmp_137 + tmp_92);
      real_t tmp_148 = tmp_145*tmp_75 + tmp_146*tmp_85 + tmp_147*tmp_65;
      real_t tmp_149 = tmp_145*tmp_73 + tmp_146*tmp_83 + tmp_147*tmp_63;
      real_t tmp_150 = tmp_145*tmp_71 + tmp_146*tmp_80 + tmp_147*tmp_50;
      real_t tmp_151 = -tmp_148 - tmp_149 - tmp_150 + 1;
      real_t tmp_152 = tmp_144*tmp_99;
      real_t tmp_153 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_110;
      real_t tmp_154 = 0.10810301816807022*tmp_16 + 0.44594849091596489*tmp_17;
      real_t tmp_155 = tmp_14*(tmp_154 + tmp_19);
      real_t tmp_156 = 0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24;
      real_t tmp_157 = tmp_14*(tmp_156 + tmp_26);
      real_t tmp_158 = 0.10810301816807022*tmp_30 + 0.44594849091596489*tmp_31;
      real_t tmp_159 = tmp_14*(tmp_158 + tmp_33);
      real_t tmp_160 = tmp_155*tmp_5 + tmp_157*tmp_21 + tmp_159*tmp_28 - 1.0/4.0;
      real_t tmp_161 = tmp_155*tmp_36 + tmp_157*tmp_37 + tmp_159*tmp_38 - 1.0/4.0;
      real_t tmp_162 = tmp_155*tmp_40 + tmp_157*tmp_41 + tmp_159*tmp_42 - 1.0/4.0;
      real_t tmp_163 = tmp_0*tmp_160 + tmp_1*tmp_161 + tmp_162*tmp_3;
      real_t tmp_164 = tmp_10*tmp_161 + tmp_160*tmp_8 + tmp_162*tmp_6;
      real_t tmp_165 = tmp_160*tmp_9 + tmp_161*tmp_4 + tmp_162*tmp_2;
      real_t tmp_166 = tmp_60*(tmp_154 + tmp_88);
      real_t tmp_167 = tmp_60*(tmp_156 + tmp_90);
      real_t tmp_168 = tmp_60*(tmp_158 + tmp_92);
      real_t tmp_169 = tmp_166*tmp_75 + tmp_167*tmp_85 + tmp_168*tmp_65;
      real_t tmp_170 = tmp_166*tmp_73 + tmp_167*tmp_83 + tmp_168*tmp_63;
      real_t tmp_171 = tmp_166*tmp_71 + tmp_167*tmp_80 + tmp_168*tmp_50;
      real_t tmp_172 = -tmp_169 - tmp_170 - tmp_171 + 1;
      real_t tmp_173 = tmp_165*tmp_99;
      real_t tmp_174 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_110;
      real_t tmp_175 = 0.091576213509770743*tmp_16 + 0.091576213509770743*tmp_17;
      real_t tmp_176 = tmp_14*(tmp_175 + tmp_19);
      real_t tmp_177 = 0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24;
      real_t tmp_178 = tmp_14*(tmp_177 + tmp_26);
      real_t tmp_179 = 0.091576213509770743*tmp_30 + 0.091576213509770743*tmp_31;
      real_t tmp_180 = tmp_14*(tmp_179 + tmp_33);
      real_t tmp_181 = tmp_176*tmp_5 + tmp_178*tmp_21 + tmp_180*tmp_28 - 1.0/4.0;
      real_t tmp_182 = tmp_176*tmp_36 + tmp_178*tmp_37 + tmp_180*tmp_38 - 1.0/4.0;
      real_t tmp_183 = tmp_176*tmp_40 + tmp_178*tmp_41 + tmp_180*tmp_42 - 1.0/4.0;
      real_t tmp_184 = tmp_0*tmp_181 + tmp_1*tmp_182 + tmp_183*tmp_3;
      real_t tmp_185 = tmp_10*tmp_182 + tmp_181*tmp_8 + tmp_183*tmp_6;
      real_t tmp_186 = tmp_181*tmp_9 + tmp_182*tmp_4 + tmp_183*tmp_2;
      real_t tmp_187 = tmp_60*(tmp_175 + tmp_88);
      real_t tmp_188 = tmp_60*(tmp_177 + tmp_90);
      real_t tmp_189 = tmp_60*(tmp_179 + tmp_92);
      real_t tmp_190 = tmp_187*tmp_75 + tmp_188*tmp_85 + tmp_189*tmp_65;
      real_t tmp_191 = tmp_187*tmp_73 + tmp_188*tmp_83 + tmp_189*tmp_63;
      real_t tmp_192 = tmp_187*tmp_71 + tmp_188*tmp_80 + tmp_189*tmp_50;
      real_t tmp_193 = -tmp_190 - tmp_191 - tmp_192 + 1;
      real_t tmp_194 = tmp_186*tmp_99;
      real_t tmp_195 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_110;
      real_t tmp_196 = 0.44594849091596489*tmp_16 + 0.44594849091596489*tmp_17;
      real_t tmp_197 = tmp_14*(tmp_19 + tmp_196);
      real_t tmp_198 = 0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24;
      real_t tmp_199 = tmp_14*(tmp_198 + tmp_26);
      real_t tmp_200 = 0.44594849091596489*tmp_30 + 0.44594849091596489*tmp_31;
      real_t tmp_201 = tmp_14*(tmp_200 + tmp_33);
      real_t tmp_202 = tmp_197*tmp_5 + tmp_199*tmp_21 + tmp_201*tmp_28 - 1.0/4.0;
      real_t tmp_203 = tmp_197*tmp_36 + tmp_199*tmp_37 + tmp_201*tmp_38 - 1.0/4.0;
      real_t tmp_204 = tmp_197*tmp_40 + tmp_199*tmp_41 + tmp_201*tmp_42 - 1.0/4.0;
      real_t tmp_205 = tmp_0*tmp_202 + tmp_1*tmp_203 + tmp_204*tmp_3;
      real_t tmp_206 = tmp_10*tmp_203 + tmp_202*tmp_8 + tmp_204*tmp_6;
      real_t tmp_207 = tmp_2*tmp_204 + tmp_202*tmp_9 + tmp_203*tmp_4;
      real_t tmp_208 = tmp_60*(tmp_196 + tmp_88);
      real_t tmp_209 = tmp_60*(tmp_198 + tmp_90);
      real_t tmp_210 = tmp_60*(tmp_200 + tmp_92);
      real_t tmp_211 = tmp_208*tmp_75 + tmp_209*tmp_85 + tmp_210*tmp_65;
      real_t tmp_212 = tmp_208*tmp_73 + tmp_209*tmp_83 + tmp_210*tmp_63;
      real_t tmp_213 = tmp_208*tmp_71 + tmp_209*tmp_80 + tmp_210*tmp_50;
      real_t tmp_214 = -tmp_211 - tmp_212 - tmp_213 + 1;
      real_t tmp_215 = tmp_207*tmp_99;
      real_t tmp_216 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_110;
      real_t tmp_217 = 0.25*p_affine_13_1*tmp_60;
      real_t tmp_218 = tmp_217*tmp_75;
      real_t tmp_219 = tmp_217*tmp_65;
      real_t tmp_220 = 0.5*p_affine_13_0*tmp_66 + 0.5*p_affine_13_1*tmp_86 + 0.5*p_affine_13_2*tmp_76;
      real_t tmp_221 = tmp_217*tmp_73;
      real_t tmp_222 = tmp_217*tmp_63;
      real_t tmp_223 = 0.5*p_affine_13_0*tmp_64 + 0.5*p_affine_13_1*tmp_84 + 0.5*p_affine_13_2*tmp_74;
      real_t tmp_224 = tmp_217*tmp_71;
      real_t tmp_225 = tmp_217*tmp_50;
      real_t tmp_226 = 0.5*p_affine_13_0*tmp_62 + 0.5*p_affine_13_1*tmp_82 + 0.5*p_affine_13_2*tmp_72;
      real_t a_0_0 = tmp_111*(-tmp_100*tmp_97 + 0.5*tmp_109*tmp_97 - tmp_44*tmp_69 - tmp_70*tmp_78 - tmp_79*tmp_87) + tmp_132*(0.5*tmp_109*tmp_130 - tmp_121*tmp_69 - tmp_122*tmp_78 - tmp_123*tmp_87 - tmp_130*tmp_131) + tmp_153*(0.5*tmp_109*tmp_151 - tmp_142*tmp_69 - tmp_143*tmp_78 - tmp_144*tmp_87 - tmp_151*tmp_152) + tmp_174*(0.5*tmp_109*tmp_172 - tmp_163*tmp_69 - tmp_164*tmp_78 - tmp_165*tmp_87 - tmp_172*tmp_173) + tmp_195*(0.5*tmp_109*tmp_193 - tmp_184*tmp_69 - tmp_185*tmp_78 - tmp_186*tmp_87 - tmp_193*tmp_194) + tmp_216*(0.5*tmp_109*tmp_214 - tmp_205*tmp_69 - tmp_206*tmp_78 - tmp_207*tmp_87 - tmp_214*tmp_215);
      real_t a_0_1 = tmp_111*(-tmp_100*tmp_94 + 0.5*tmp_109*tmp_94 - tmp_218*tmp_70 - tmp_219*tmp_44 - tmp_220*tmp_79) + tmp_132*(0.5*tmp_109*tmp_127 - tmp_121*tmp_219 - tmp_122*tmp_218 - tmp_123*tmp_220 - tmp_127*tmp_131) + tmp_153*(0.5*tmp_109*tmp_148 - tmp_142*tmp_219 - tmp_143*tmp_218 - tmp_144*tmp_220 - tmp_148*tmp_152) + tmp_174*(0.5*tmp_109*tmp_169 - tmp_163*tmp_219 - tmp_164*tmp_218 - tmp_165*tmp_220 - tmp_169*tmp_173) + tmp_195*(0.5*tmp_109*tmp_190 - tmp_184*tmp_219 - tmp_185*tmp_218 - tmp_186*tmp_220 - tmp_190*tmp_194) + tmp_216*(0.5*tmp_109*tmp_211 - tmp_205*tmp_219 - tmp_206*tmp_218 - tmp_207*tmp_220 - tmp_211*tmp_215);
      real_t a_0_2 = tmp_111*(-tmp_100*tmp_95 + 0.5*tmp_109*tmp_95 - tmp_221*tmp_70 - tmp_222*tmp_44 - tmp_223*tmp_79) + tmp_132*(0.5*tmp_109*tmp_128 - tmp_121*tmp_222 - tmp_122*tmp_221 - tmp_123*tmp_223 - tmp_128*tmp_131) + tmp_153*(0.5*tmp_109*tmp_149 - tmp_142*tmp_222 - tmp_143*tmp_221 - tmp_144*tmp_223 - tmp_149*tmp_152) + tmp_174*(0.5*tmp_109*tmp_170 - tmp_163*tmp_222 - tmp_164*tmp_221 - tmp_165*tmp_223 - tmp_170*tmp_173) + tmp_195*(0.5*tmp_109*tmp_191 - tmp_184*tmp_222 - tmp_185*tmp_221 - tmp_186*tmp_223 - tmp_191*tmp_194) + tmp_216*(0.5*tmp_109*tmp_212 - tmp_205*tmp_222 - tmp_206*tmp_221 - tmp_207*tmp_223 - tmp_212*tmp_215);
      real_t a_0_3 = tmp_111*(-tmp_100*tmp_96 + 0.5*tmp_109*tmp_96 - tmp_224*tmp_70 - tmp_225*tmp_44 - tmp_226*tmp_79) + tmp_132*(0.5*tmp_109*tmp_129 - tmp_121*tmp_225 - tmp_122*tmp_224 - tmp_123*tmp_226 - tmp_129*tmp_131) + tmp_153*(0.5*tmp_109*tmp_150 - tmp_142*tmp_225 - tmp_143*tmp_224 - tmp_144*tmp_226 - tmp_150*tmp_152) + tmp_174*(0.5*tmp_109*tmp_171 - tmp_163*tmp_225 - tmp_164*tmp_224 - tmp_165*tmp_226 - tmp_171*tmp_173) + tmp_195*(0.5*tmp_109*tmp_192 - tmp_184*tmp_225 - tmp_185*tmp_224 - tmp_186*tmp_226 - tmp_192*tmp_194) + tmp_216*(0.5*tmp_109*tmp_213 - tmp_205*tmp_225 - tmp_206*tmp_224 - tmp_207*tmp_226 - tmp_213*tmp_215);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_6*tmp_9;
      real_t tmp_14 = tmp_4*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_12 + tmp_0*tmp_7 - tmp_1*tmp_13 + tmp_1*tmp_2*tmp_8 + tmp_11*tmp_3 - tmp_14*tmp_3);
      real_t tmp_16 = -p_affine_8_2 + p_affine_9_2;
      real_t tmp_17 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_18 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_19 = tmp_15*(0.091576213509770743*tmp_16 + 0.81684757298045851*tmp_17 + tmp_18);
      real_t tmp_20 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_21 = -p_affine_8_1 + p_affine_9_1;
      real_t tmp_22 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_23 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_24 = tmp_15*(0.091576213509770743*tmp_21 + 0.81684757298045851*tmp_22 + tmp_23);
      real_t tmp_25 = -tmp_12 + tmp_7;
      real_t tmp_26 = -p_affine_8_0 + p_affine_9_0;
      real_t tmp_27 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_28 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_29 = tmp_15*(0.091576213509770743*tmp_26 + 0.81684757298045851*tmp_27 + tmp_28);
      real_t tmp_30 = tmp_19*tmp_5 + tmp_20*tmp_24 + tmp_25*tmp_29 - 1.0/4.0;
      real_t tmp_31 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_32 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_33 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_34 = tmp_19*tmp_31 + tmp_24*tmp_32 + tmp_29*tmp_33 - 1.0/4.0;
      real_t tmp_35 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_36 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_37 = tmp_11 - tmp_14;
      real_t tmp_38 = tmp_19*tmp_35 + tmp_24*tmp_36 + tmp_29*tmp_37 - 1.0/4.0;
      real_t tmp_39 = tmp_0*tmp_30 + tmp_1*tmp_34 + tmp_3*tmp_38;
      real_t tmp_40 = 0.5*tmp_15;
      real_t tmp_41 = tmp_37*tmp_40;
      real_t tmp_42 = tmp_33*tmp_40;
      real_t tmp_43 = tmp_25*tmp_40;
      real_t tmp_44 = -tmp_41 - tmp_42 - tmp_43;
      real_t tmp_45 = p_affine_13_1*tmp_44;
      real_t tmp_46 = tmp_10*tmp_34 + tmp_30*tmp_8 + tmp_38*tmp_6;
      real_t tmp_47 = tmp_35*tmp_40;
      real_t tmp_48 = tmp_31*tmp_40;
      real_t tmp_49 = tmp_40*tmp_5;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = p_affine_13_1*tmp_50;
      real_t tmp_52 = 1.0*tmp_15;
      real_t tmp_53 = tmp_36*tmp_52;
      real_t tmp_54 = tmp_32*tmp_52;
      real_t tmp_55 = tmp_20*tmp_52;
      real_t tmp_56 = p_affine_13_0*tmp_44 + p_affine_13_1*(-tmp_53 - tmp_54 - tmp_55) + p_affine_13_2*tmp_50;
      real_t tmp_57 = tmp_2*tmp_38 + tmp_30*tmp_9 + tmp_34*tmp_4;
      real_t tmp_58 = tmp_15*(0.44594849091596489*tmp_16 + 0.10810301816807022*tmp_17 + tmp_18);
      real_t tmp_59 = tmp_15*(0.44594849091596489*tmp_21 + 0.10810301816807022*tmp_22 + tmp_23);
      real_t tmp_60 = tmp_15*(0.44594849091596489*tmp_26 + 0.10810301816807022*tmp_27 + tmp_28);
      real_t tmp_61 = tmp_20*tmp_59 + tmp_25*tmp_60 + tmp_5*tmp_58 - 1.0/4.0;
      real_t tmp_62 = tmp_31*tmp_58 + tmp_32*tmp_59 + tmp_33*tmp_60 - 1.0/4.0;
      real_t tmp_63 = tmp_35*tmp_58 + tmp_36*tmp_59 + tmp_37*tmp_60 - 1.0/4.0;
      real_t tmp_64 = tmp_0*tmp_61 + tmp_1*tmp_62 + tmp_3*tmp_63;
      real_t tmp_65 = tmp_10*tmp_62 + tmp_6*tmp_63 + tmp_61*tmp_8;
      real_t tmp_66 = tmp_2*tmp_63 + tmp_4*tmp_62 + tmp_61*tmp_9;
      real_t tmp_67 = tmp_15*(0.81684757298045851*tmp_16 + 0.091576213509770743*tmp_17 + tmp_18);
      real_t tmp_68 = tmp_15*(0.81684757298045851*tmp_21 + 0.091576213509770743*tmp_22 + tmp_23);
      real_t tmp_69 = tmp_15*(0.81684757298045851*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28);
      real_t tmp_70 = tmp_20*tmp_68 + tmp_25*tmp_69 + tmp_5*tmp_67 - 1.0/4.0;
      real_t tmp_71 = tmp_31*tmp_67 + tmp_32*tmp_68 + tmp_33*tmp_69 - 1.0/4.0;
      real_t tmp_72 = tmp_35*tmp_67 + tmp_36*tmp_68 + tmp_37*tmp_69 - 1.0/4.0;
      real_t tmp_73 = tmp_0*tmp_70 + tmp_1*tmp_71 + tmp_3*tmp_72;
      real_t tmp_74 = tmp_10*tmp_71 + tmp_6*tmp_72 + tmp_70*tmp_8;
      real_t tmp_75 = tmp_2*tmp_72 + tmp_4*tmp_71 + tmp_70*tmp_9;
      real_t tmp_76 = tmp_15*(0.10810301816807022*tmp_16 + 0.44594849091596489*tmp_17 + tmp_18);
      real_t tmp_77 = tmp_15*(0.10810301816807022*tmp_21 + 0.44594849091596489*tmp_22 + tmp_23);
      real_t tmp_78 = tmp_15*(0.10810301816807022*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28);
      real_t tmp_79 = tmp_20*tmp_77 + tmp_25*tmp_78 + tmp_5*tmp_76 - 1.0/4.0;
      real_t tmp_80 = tmp_31*tmp_76 + tmp_32*tmp_77 + tmp_33*tmp_78 - 1.0/4.0;
      real_t tmp_81 = tmp_35*tmp_76 + tmp_36*tmp_77 + tmp_37*tmp_78 - 1.0/4.0;
      real_t tmp_82 = tmp_0*tmp_79 + tmp_1*tmp_80 + tmp_3*tmp_81;
      real_t tmp_83 = tmp_10*tmp_80 + tmp_6*tmp_81 + tmp_79*tmp_8;
      real_t tmp_84 = tmp_2*tmp_81 + tmp_4*tmp_80 + tmp_79*tmp_9;
      real_t tmp_85 = tmp_15*(0.091576213509770743*tmp_16 + 0.091576213509770743*tmp_17 + tmp_18);
      real_t tmp_86 = tmp_15*(0.091576213509770743*tmp_21 + 0.091576213509770743*tmp_22 + tmp_23);
      real_t tmp_87 = tmp_15*(0.091576213509770743*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28);
      real_t tmp_88 = tmp_20*tmp_86 + tmp_25*tmp_87 + tmp_5*tmp_85 - 1.0/4.0;
      real_t tmp_89 = tmp_31*tmp_85 + tmp_32*tmp_86 + tmp_33*tmp_87 - 1.0/4.0;
      real_t tmp_90 = tmp_35*tmp_85 + tmp_36*tmp_86 + tmp_37*tmp_87 - 1.0/4.0;
      real_t tmp_91 = tmp_0*tmp_88 + tmp_1*tmp_89 + tmp_3*tmp_90;
      real_t tmp_92 = tmp_10*tmp_89 + tmp_6*tmp_90 + tmp_8*tmp_88;
      real_t tmp_93 = tmp_2*tmp_90 + tmp_4*tmp_89 + tmp_88*tmp_9;
      real_t tmp_94 = tmp_15*(0.44594849091596489*tmp_16 + 0.44594849091596489*tmp_17 + tmp_18);
      real_t tmp_95 = tmp_15*(0.44594849091596489*tmp_21 + 0.44594849091596489*tmp_22 + tmp_23);
      real_t tmp_96 = tmp_15*(0.44594849091596489*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28);
      real_t tmp_97 = tmp_20*tmp_95 + tmp_25*tmp_96 + tmp_5*tmp_94 - 1.0/4.0;
      real_t tmp_98 = tmp_31*tmp_94 + tmp_32*tmp_95 + tmp_33*tmp_96 - 1.0/4.0;
      real_t tmp_99 = tmp_35*tmp_94 + tmp_36*tmp_95 + tmp_37*tmp_96 - 1.0/4.0;
      real_t tmp_100 = tmp_0*tmp_97 + tmp_1*tmp_98 + tmp_3*tmp_99;
      real_t tmp_101 = tmp_10*tmp_98 + tmp_6*tmp_99 + tmp_8*tmp_97;
      real_t tmp_102 = tmp_2*tmp_99 + tmp_4*tmp_98 + tmp_9*tmp_97;
      real_t tmp_103 = p_affine_13_1*tmp_49;
      real_t tmp_104 = p_affine_13_1*tmp_43;
      real_t tmp_105 = p_affine_13_0*tmp_43 + p_affine_13_1*tmp_55 + p_affine_13_2*tmp_49;
      real_t tmp_106 = p_affine_13_1*tmp_48;
      real_t tmp_107 = p_affine_13_1*tmp_42;
      real_t tmp_108 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_54 + p_affine_13_2*tmp_48;
      real_t tmp_109 = p_affine_13_1*tmp_47;
      real_t tmp_110 = p_affine_13_1*tmp_41;
      real_t tmp_111 = p_affine_13_0*tmp_41 + p_affine_13_1*tmp_53 + p_affine_13_2*tmp_47;
      real_t a_0_0 = -0.11169079483900572*tmp_100*tmp_45 - 0.11169079483900572*tmp_101*tmp_51 - 0.11169079483900572*tmp_102*tmp_56 - 0.054975871827660928*tmp_39*tmp_45 - 0.11169079483900572*tmp_45*tmp_64 - 0.054975871827660928*tmp_45*tmp_73 - 0.11169079483900572*tmp_45*tmp_82 - 0.054975871827660928*tmp_45*tmp_91 - 0.054975871827660928*tmp_46*tmp_51 - 0.11169079483900572*tmp_51*tmp_65 - 0.054975871827660928*tmp_51*tmp_74 - 0.11169079483900572*tmp_51*tmp_83 - 0.054975871827660928*tmp_51*tmp_92 - 0.054975871827660928*tmp_56*tmp_57 - 0.11169079483900572*tmp_56*tmp_66 - 0.054975871827660928*tmp_56*tmp_75 - 0.11169079483900572*tmp_56*tmp_84 - 0.054975871827660928*tmp_56*tmp_93;
      real_t a_0_1 = -0.11169079483900572*tmp_100*tmp_104 - 0.11169079483900572*tmp_101*tmp_103 - 0.11169079483900572*tmp_102*tmp_105 - 0.054975871827660928*tmp_103*tmp_46 - 0.11169079483900572*tmp_103*tmp_65 - 0.054975871827660928*tmp_103*tmp_74 - 0.11169079483900572*tmp_103*tmp_83 - 0.054975871827660928*tmp_103*tmp_92 - 0.054975871827660928*tmp_104*tmp_39 - 0.11169079483900572*tmp_104*tmp_64 - 0.054975871827660928*tmp_104*tmp_73 - 0.11169079483900572*tmp_104*tmp_82 - 0.054975871827660928*tmp_104*tmp_91 - 0.054975871827660928*tmp_105*tmp_57 - 0.11169079483900572*tmp_105*tmp_66 - 0.054975871827660928*tmp_105*tmp_75 - 0.11169079483900572*tmp_105*tmp_84 - 0.054975871827660928*tmp_105*tmp_93;
      real_t a_0_2 = -0.11169079483900572*tmp_100*tmp_107 - 0.11169079483900572*tmp_101*tmp_106 - 0.11169079483900572*tmp_102*tmp_108 - 0.054975871827660928*tmp_106*tmp_46 - 0.11169079483900572*tmp_106*tmp_65 - 0.054975871827660928*tmp_106*tmp_74 - 0.11169079483900572*tmp_106*tmp_83 - 0.054975871827660928*tmp_106*tmp_92 - 0.054975871827660928*tmp_107*tmp_39 - 0.11169079483900572*tmp_107*tmp_64 - 0.054975871827660928*tmp_107*tmp_73 - 0.11169079483900572*tmp_107*tmp_82 - 0.054975871827660928*tmp_107*tmp_91 - 0.054975871827660928*tmp_108*tmp_57 - 0.11169079483900572*tmp_108*tmp_66 - 0.054975871827660928*tmp_108*tmp_75 - 0.11169079483900572*tmp_108*tmp_84 - 0.054975871827660928*tmp_108*tmp_93;
      real_t a_0_3 = -0.11169079483900572*tmp_100*tmp_110 - 0.11169079483900572*tmp_101*tmp_109 - 0.11169079483900572*tmp_102*tmp_111 - 0.054975871827660928*tmp_109*tmp_46 - 0.11169079483900572*tmp_109*tmp_65 - 0.054975871827660928*tmp_109*tmp_74 - 0.11169079483900572*tmp_109*tmp_83 - 0.054975871827660928*tmp_109*tmp_92 - 0.054975871827660928*tmp_110*tmp_39 - 0.11169079483900572*tmp_110*tmp_64 - 0.054975871827660928*tmp_110*tmp_73 - 0.11169079483900572*tmp_110*tmp_82 - 0.054975871827660928*tmp_110*tmp_91 - 0.054975871827660928*tmp_111*tmp_57 - 0.11169079483900572*tmp_111*tmp_66 - 0.054975871827660928*tmp_111*tmp_75 - 0.11169079483900572*tmp_111*tmp_84 - 0.054975871827660928*tmp_111*tmp_93;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }

};




class EGEpsilonFormP1E_1 : public hyteg::dg::DGForm
{

 public:
    EGEpsilonFormP1E_1(std::function< real_t ( const Point3D & ) > mu)
: callback_Scalar_Variable_Coefficient_3D_mu (mu)
, callback_Scalar_Variable_Coefficient_2D_mu (mu)
    {}

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_mu;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;

void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
}

      
 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = -tmp_3;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_2 - tmp_4*tmp_6);
      real_t tmp_8 = 1.0*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_2*tmp_8 + tmp_6*tmp_9;
      real_t tmp_11 = Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_10;
      real_t tmp_12 = -2*tmp_0*tmp_8 - 2*tmp_9;
      real_t tmp_13 = 0.5*tmp_7;
      real_t tmp_14 = tmp_0*tmp_13;
      real_t tmp_15 = tmp_1*tmp_13;
      real_t tmp_16 = tmp_13*tmp_5;
      real_t tmp_17 = tmp_1*tmp_16 + tmp_14*tmp_3 + tmp_14*tmp_4 + tmp_15*tmp_6;
      real_t tmp_18 = Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_17;
      real_t tmp_19 = -4*tmp_15 - 4*tmp_16;
      real_t tmp_20 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_21 = 0.054975871827660928*tmp_20;
      real_t tmp_22 = tmp_10*tmp_12;
      real_t tmp_23 = tmp_17*tmp_19;
      real_t tmp_24 = 0.11169079483900572*tmp_20;
      real_t tmp_25 = 0.054975871827660928*tmp_20;
      real_t tmp_26 = 0.11169079483900572*tmp_20;
      real_t tmp_27 = 0.054975871827660928*tmp_20;
      real_t tmp_28 = 0.11169079483900572*tmp_20;
      real_t tmp_29 = 2.0*tmp_7;
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_1*tmp_29;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_17*tmp_31;
      real_t tmp_34 = tmp_0*tmp_29;
      real_t tmp_35 = tmp_29*tmp_5;
      real_t tmp_36 = tmp_10*tmp_34;
      real_t tmp_37 = tmp_17*tmp_35;
      real_t a_0_0 = tmp_21*(tmp_11*tmp_12 + tmp_18*tmp_19) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_23) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_23) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_23) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_23) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_22 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_23);
      real_t a_1_0 = tmp_21*(tmp_11*tmp_30 + tmp_18*tmp_31) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_33) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_33) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_33) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_33) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_32 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_33);
      real_t a_2_0 = tmp_21*(tmp_11*tmp_34 + tmp_18*tmp_35) + tmp_24*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_37) + tmp_25*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_37) + tmp_26*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_37) + tmp_27*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_37) + tmp_28*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_36 + Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_37);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = tmp_0*tmp_2;
      real_t tmp_4 = -tmp_1;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_3 - tmp_4*tmp_6);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_10 = tmp_7*(0.21132486540518713*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = tmp_7*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_1*tmp_10 + tmp_13*tmp_2;
      real_t tmp_15 = tmp_14 - 1.0/3.0;
      real_t tmp_16 = tmp_0*tmp_10 + tmp_13*tmp_5;
      real_t tmp_17 = tmp_16 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_1*(tmp_0*tmp_15 + tmp_17*tmp_4);
      real_t tmp_19 = 0.5*tmp_7;
      real_t tmp_20 = tmp_19*tmp_2;
      real_t tmp_21 = tmp_19*tmp_5;
      real_t tmp_22 = -tmp_20 - tmp_21;
      real_t tmp_23 = 0.5*tmp_22;
      real_t tmp_24 = tmp_15*tmp_6 + tmp_17*tmp_2;
      real_t tmp_25 = 1.0*tmp_7;
      real_t tmp_26 = tmp_0*tmp_25;
      real_t tmp_27 = tmp_1*tmp_25;
      real_t tmp_28 = 0.5*p_affine_10_0*tmp_22 + 0.5*p_affine_10_1*(-tmp_26 - tmp_27);
      real_t tmp_29 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_30 = 1.0 / (tmp_29);
      real_t tmp_31 = -tmp_14 - tmp_16 + 1;
      real_t tmp_32 = tmp_0*tmp_19;
      real_t tmp_33 = 0.5*p_affine_10_0*(tmp_1*tmp_32 + tmp_2*tmp_21 + tmp_20*tmp_6 + tmp_32*tmp_4) + 0.5*p_affine_10_1*(tmp_25*tmp_3 + tmp_27*tmp_6);
      real_t tmp_34 = 2*tmp_29;
      real_t tmp_35 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_34;
      real_t tmp_36 = tmp_7*(0.78867513459481287*tmp_8 + tmp_9);
      real_t tmp_37 = tmp_7*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_38 = tmp_1*tmp_36 + tmp_2*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/3.0;
      real_t tmp_40 = tmp_0*tmp_36 + tmp_37*tmp_5;
      real_t tmp_41 = tmp_40 - 1.0/3.0;
      real_t tmp_42 = p_affine_10_1*(tmp_0*tmp_39 + tmp_4*tmp_41);
      real_t tmp_43 = tmp_2*tmp_41 + tmp_39*tmp_6;
      real_t tmp_44 = -tmp_38 - tmp_40 + 1;
      real_t tmp_45 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_34;
      real_t tmp_46 = 0.25*tmp_7;
      real_t tmp_47 = tmp_2*tmp_46;
      real_t tmp_48 = 0.5*p_affine_10_0*tmp_20 + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_49 = tmp_46*tmp_5;
      real_t tmp_50 = 0.5*p_affine_10_0*tmp_21 + 0.5*p_affine_10_1*tmp_26;
      real_t a_0_0 = tmp_35*(-tmp_18*tmp_23 - tmp_24*tmp_28 + tmp_24*tmp_30*tmp_31 - tmp_31*tmp_33) + tmp_45*(-tmp_23*tmp_42 - tmp_28*tmp_43 + tmp_30*tmp_43*tmp_44 - tmp_33*tmp_44);
      real_t a_1_0 = tmp_35*(tmp_14*tmp_24*tmp_30 - tmp_14*tmp_33 - tmp_18*tmp_47 - tmp_24*tmp_48) + tmp_45*(tmp_30*tmp_38*tmp_43 - tmp_33*tmp_38 - tmp_42*tmp_47 - tmp_43*tmp_48);
      real_t a_2_0 = tmp_35*(tmp_16*tmp_24*tmp_30 - tmp_16*tmp_33 - tmp_18*tmp_49 - tmp_24*tmp_50) + tmp_45*(tmp_30*tmp_40*tmp_43 - tmp_33*tmp_40 - tmp_42*tmp_49 - tmp_43*tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElementInner[0]( 0 );
      const auto p_affine_0_1 = coordsElementInner[0]( 1 );

      const auto p_affine_1_0 = coordsElementInner[1]( 0 );
      const auto p_affine_1_1 = coordsElementInner[1]( 1 );

      const auto p_affine_2_0 = coordsElementInner[2]( 0 );
      const auto p_affine_2_1 = coordsElementInner[2]( 1 );

      const auto p_affine_3_0 = coordsElementOuter[0]( 0 );
      const auto p_affine_3_1 = coordsElementOuter[0]( 1 );

      const auto p_affine_4_0 = coordsElementOuter[1]( 0 );
      const auto p_affine_4_1 = coordsElementOuter[1]( 1 );

      const auto p_affine_5_0 = coordsElementOuter[2]( 0 );
      const auto p_affine_5_1 = coordsElementOuter[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_1 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_2 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_3 = tmp_0*tmp_2;
      real_t tmp_4 = -tmp_1;
      real_t tmp_5 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_3 - tmp_4*tmp_6);
      real_t tmp_8 = -p_affine_3_1;
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.21132486540518713*tmp_9;
      real_t tmp_11 = tmp_7*(tmp_10 + tmp_8);
      real_t tmp_12 = -p_affine_3_0;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + 0.21132486540518713*tmp_13;
      real_t tmp_15 = tmp_7*(tmp_12 + tmp_14);
      real_t tmp_16 = tmp_1*tmp_11 + tmp_15*tmp_2 - 1.0/3.0;
      real_t tmp_17 = tmp_0*tmp_11 + tmp_15*tmp_5 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_1*(tmp_0*tmp_16 + tmp_17*tmp_4);
      real_t tmp_19 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_20 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_21 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_22 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_23 = 1.0 / (tmp_19*tmp_20 - tmp_21*tmp_22);
      real_t tmp_24 = 0.5*tmp_23;
      real_t tmp_25 = tmp_19*tmp_24;
      real_t tmp_26 = tmp_22*tmp_24;
      real_t tmp_27 = -tmp_25 - tmp_26;
      real_t tmp_28 = 0.5*tmp_27;
      real_t tmp_29 = tmp_16*tmp_6 + tmp_17*tmp_2;
      real_t tmp_30 = 1.0*tmp_23;
      real_t tmp_31 = tmp_20*tmp_30;
      real_t tmp_32 = tmp_21*tmp_30;
      real_t tmp_33 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*(-tmp_31 - tmp_32);
      real_t tmp_34 = -p_affine_0_1;
      real_t tmp_35 = tmp_23*(tmp_10 + tmp_34);
      real_t tmp_36 = -p_affine_0_0;
      real_t tmp_37 = tmp_23*(tmp_14 + tmp_36);
      real_t tmp_38 = tmp_19*tmp_37 + tmp_21*tmp_35;
      real_t tmp_39 = tmp_20*tmp_35 + tmp_22*tmp_37;
      real_t tmp_40 = -tmp_38 - tmp_39 + 1;
      real_t tmp_41 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_42 = 1.0 / (tmp_41);
      real_t tmp_43 = tmp_29*tmp_42;
      real_t tmp_44 = 1.0*tmp_7;
      real_t tmp_45 = 0.5*tmp_7;
      real_t tmp_46 = tmp_0*tmp_45;
      real_t tmp_47 = tmp_2*tmp_45;
      real_t tmp_48 = 0.5*p_affine_10_0*(tmp_1*tmp_46 + tmp_4*tmp_46 + tmp_47*tmp_5 + tmp_47*tmp_6) + 0.5*p_affine_10_1*(tmp_1*tmp_44*tmp_6 + tmp_3*tmp_44);
      real_t tmp_49 = 2*tmp_41;
      real_t tmp_50 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_49;
      real_t tmp_51 = p_affine_6_1 + 0.78867513459481287*tmp_9;
      real_t tmp_52 = tmp_7*(tmp_51 + tmp_8);
      real_t tmp_53 = p_affine_6_0 + 0.78867513459481287*tmp_13;
      real_t tmp_54 = tmp_7*(tmp_12 + tmp_53);
      real_t tmp_55 = tmp_1*tmp_52 + tmp_2*tmp_54 - 1.0/3.0;
      real_t tmp_56 = tmp_0*tmp_52 + tmp_5*tmp_54 - 1.0/3.0;
      real_t tmp_57 = p_affine_10_1*(tmp_0*tmp_55 + tmp_4*tmp_56);
      real_t tmp_58 = tmp_2*tmp_56 + tmp_55*tmp_6;
      real_t tmp_59 = tmp_23*(tmp_34 + tmp_51);
      real_t tmp_60 = tmp_23*(tmp_36 + tmp_53);
      real_t tmp_61 = tmp_19*tmp_60 + tmp_21*tmp_59;
      real_t tmp_62 = tmp_20*tmp_59 + tmp_22*tmp_60;
      real_t tmp_63 = -tmp_61 - tmp_62 + 1;
      real_t tmp_64 = tmp_42*tmp_58;
      real_t tmp_65 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_49;
      real_t tmp_66 = 0.25*tmp_23;
      real_t tmp_67 = tmp_19*tmp_66;
      real_t tmp_68 = 0.5*p_affine_10_0*tmp_25 + 0.5*p_affine_10_1*tmp_32;
      real_t tmp_69 = tmp_22*tmp_66;
      real_t tmp_70 = 0.5*p_affine_10_0*tmp_26 + 0.5*p_affine_10_1*tmp_31;
      real_t a_0_0 = tmp_50*(tmp_18*tmp_28 + tmp_29*tmp_33 - tmp_40*tmp_43 - tmp_40*tmp_48) + tmp_65*(tmp_28*tmp_57 + tmp_33*tmp_58 - tmp_48*tmp_63 - tmp_63*tmp_64);
      real_t a_1_0 = tmp_50*(tmp_18*tmp_67 + tmp_29*tmp_68 - tmp_38*tmp_43 - tmp_38*tmp_48) + tmp_65*(-tmp_48*tmp_61 + tmp_57*tmp_67 + tmp_58*tmp_68 - tmp_61*tmp_64);
      real_t a_2_0 = tmp_50*(tmp_18*tmp_69 + tmp_29*tmp_70 - tmp_39*tmp_43 - tmp_39*tmp_48) + tmp_65*(-tmp_48*tmp_62 + tmp_57*tmp_69 + tmp_58*tmp_70 - tmp_62*tmp_64);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = tmp_0*tmp_2;
      real_t tmp_4 = -tmp_1;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_3 - tmp_4*tmp_6);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_10 = tmp_7*(0.21132486540518713*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = tmp_7*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_1*tmp_10 + tmp_13*tmp_2;
      real_t tmp_15 = tmp_14 - 1.0/3.0;
      real_t tmp_16 = tmp_0*tmp_10 + tmp_13*tmp_5;
      real_t tmp_17 = tmp_16 - 1.0/3.0;
      real_t tmp_18 = tmp_0*tmp_15 + tmp_17*tmp_4;
      real_t tmp_19 = 0.5*tmp_7;
      real_t tmp_20 = tmp_19*tmp_2;
      real_t tmp_21 = tmp_19*tmp_5;
      real_t tmp_22 = -tmp_20 - tmp_21;
      real_t tmp_23 = p_affine_10_1*tmp_22;
      real_t tmp_24 = 1.0*tmp_7;
      real_t tmp_25 = tmp_0*tmp_24;
      real_t tmp_26 = tmp_1*tmp_24;
      real_t tmp_27 = p_affine_10_0*tmp_22 + p_affine_10_1*(-tmp_25 - tmp_26);
      real_t tmp_28 = tmp_15*tmp_6 + tmp_17*tmp_2;
      real_t tmp_29 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_30 = 1.0 / (tmp_29);
      real_t tmp_31 = -tmp_14 - tmp_16 + 1;
      real_t tmp_32 = tmp_0*tmp_19;
      real_t tmp_33 = p_affine_10_0*(tmp_1*tmp_32 + tmp_2*tmp_21 + tmp_20*tmp_6 + tmp_32*tmp_4) + p_affine_10_1*(tmp_24*tmp_3 + tmp_26*tmp_6);
      real_t tmp_34 = 2*tmp_29;
      real_t tmp_35 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_34;
      real_t tmp_36 = tmp_7*(0.78867513459481287*tmp_8 + tmp_9);
      real_t tmp_37 = tmp_7*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_38 = tmp_1*tmp_36 + tmp_2*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/3.0;
      real_t tmp_40 = tmp_0*tmp_36 + tmp_37*tmp_5;
      real_t tmp_41 = tmp_40 - 1.0/3.0;
      real_t tmp_42 = tmp_0*tmp_39 + tmp_4*tmp_41;
      real_t tmp_43 = tmp_2*tmp_41 + tmp_39*tmp_6;
      real_t tmp_44 = -tmp_38 - tmp_40 + 1;
      real_t tmp_45 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_34;
      real_t tmp_46 = p_affine_10_1*tmp_20;
      real_t tmp_47 = p_affine_10_0*tmp_20 + p_affine_10_1*tmp_26;
      real_t tmp_48 = p_affine_10_1*tmp_21;
      real_t tmp_49 = p_affine_10_0*tmp_21 + p_affine_10_1*tmp_25;
      real_t a_0_0 = tmp_35*(-tmp_18*tmp_23 - tmp_27*tmp_28 + tmp_28*tmp_30*tmp_31 - tmp_31*tmp_33) + tmp_45*(-tmp_23*tmp_42 - tmp_27*tmp_43 + tmp_30*tmp_43*tmp_44 - tmp_33*tmp_44);
      real_t a_1_0 = tmp_35*(tmp_14*tmp_28*tmp_30 - tmp_14*tmp_33 - tmp_18*tmp_46 - tmp_28*tmp_47) + tmp_45*(tmp_30*tmp_38*tmp_43 - tmp_33*tmp_38 - tmp_42*tmp_46 - tmp_43*tmp_47);
      real_t a_2_0 = tmp_35*(tmp_16*tmp_28*tmp_30 - tmp_16*tmp_33 - tmp_18*tmp_48 - tmp_28*tmp_49) + tmp_45*(tmp_30*tmp_40*tmp_43 - tmp_33*tmp_40 - tmp_42*tmp_48 - tmp_43*tmp_49);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

  void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

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

      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id6 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id7 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id8 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id9 );
      Scalar_Variable_Coefficient_3D_mu( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id10 );
      Scalar_Variable_Coefficient_3D_mu( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id11 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id12 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id13 );
      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_5 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_6 = -tmp_3 + tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_8 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_9 = tmp_2*tmp_8;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_12 = tmp_5*tmp_8;
      real_t tmp_13 = tmp_11*tmp_4;
      real_t tmp_14 = 1.0 / (-tmp_0*tmp_3 + tmp_0*tmp_4*tmp_5 + tmp_1*tmp_10*tmp_11 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_7*tmp_9);
      real_t tmp_15 = 1.0*tmp_14;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = -tmp_13 + tmp_9;
      real_t tmp_18 = tmp_15*tmp_17;
      real_t tmp_19 = tmp_1*tmp_11 - tmp_12;
      real_t tmp_20 = tmp_15*tmp_19;
      real_t tmp_21 = tmp_0*tmp_16 + tmp_10*tmp_20 + tmp_18*tmp_7;
      real_t tmp_22 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_21;
      real_t tmp_23 = -2*tmp_16 - 2*tmp_18 - 2*tmp_20;
      real_t tmp_24 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_25 = -tmp_0*tmp_1 + tmp_7*tmp_8;
      real_t tmp_26 = 0.5*tmp_14;
      real_t tmp_27 = tmp_25*tmp_26;
      real_t tmp_28 = tmp_0*tmp_4 - tmp_10*tmp_8;
      real_t tmp_29 = tmp_26*tmp_28;
      real_t tmp_30 = tmp_1*tmp_10 - tmp_4*tmp_7;
      real_t tmp_31 = tmp_26*tmp_30;
      real_t tmp_32 = tmp_26*tmp_6;
      real_t tmp_33 = tmp_17*tmp_26;
      real_t tmp_34 = tmp_19*tmp_26;
      real_t tmp_35 = tmp_0*tmp_31 + tmp_10*tmp_27 + tmp_11*tmp_32 + tmp_2*tmp_34 + tmp_29*tmp_7 + tmp_33*tmp_5;
      real_t tmp_36 = tmp_35*(-tmp_27 - tmp_29 - tmp_31);
      real_t tmp_37 = tmp_0*tmp_5 - tmp_11*tmp_7;
      real_t tmp_38 = tmp_26*tmp_37;
      real_t tmp_39 = -tmp_0*tmp_2 + tmp_10*tmp_11;
      real_t tmp_40 = tmp_26*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_5 + tmp_2*tmp_7;
      real_t tmp_42 = tmp_26*tmp_41;
      real_t tmp_43 = tmp_0*tmp_42 + tmp_1*tmp_33 + tmp_10*tmp_38 + tmp_32*tmp_8 + tmp_34*tmp_4 + tmp_40*tmp_7;
      real_t tmp_44 = tmp_43*(-tmp_38 - tmp_40 - tmp_42);
      real_t tmp_45 = p_affine_0_0*p_affine_1_1;
      real_t tmp_46 = p_affine_0_0*p_affine_1_2;
      real_t tmp_47 = p_affine_2_1*p_affine_3_2;
      real_t tmp_48 = p_affine_0_1*p_affine_1_0;
      real_t tmp_49 = p_affine_0_1*p_affine_1_2;
      real_t tmp_50 = p_affine_2_2*p_affine_3_0;
      real_t tmp_51 = p_affine_0_2*p_affine_1_0;
      real_t tmp_52 = p_affine_0_2*p_affine_1_1;
      real_t tmp_53 = p_affine_2_0*p_affine_3_1;
      real_t tmp_54 = p_affine_2_2*p_affine_3_1;
      real_t tmp_55 = p_affine_2_0*p_affine_3_2;
      real_t tmp_56 = p_affine_2_1*p_affine_3_0;
      real_t tmp_57 = std::abs(p_affine_0_0*tmp_47 - p_affine_0_0*tmp_54 + p_affine_0_1*tmp_50 - p_affine_0_1*tmp_55 + p_affine_0_2*tmp_53 - p_affine_0_2*tmp_56 - p_affine_1_0*tmp_47 + p_affine_1_0*tmp_54 - p_affine_1_1*tmp_50 + p_affine_1_1*tmp_55 - p_affine_1_2*tmp_53 + p_affine_1_2*tmp_56 + p_affine_2_0*tmp_49 - p_affine_2_0*tmp_52 - p_affine_2_1*tmp_46 + p_affine_2_1*tmp_51 + p_affine_2_2*tmp_45 - p_affine_2_2*tmp_48 - p_affine_3_0*tmp_49 + p_affine_3_0*tmp_52 + p_affine_3_1*tmp_46 - p_affine_3_1*tmp_51 - p_affine_3_2*tmp_45 + p_affine_3_2*tmp_48);
      real_t tmp_58 = 0.018781320953002646*tmp_57;
      real_t tmp_59 = tmp_21*tmp_23;
      real_t tmp_60 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id1;
      real_t tmp_61 = 0.012248840519393657*tmp_57;
      real_t tmp_62 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id2;
      real_t tmp_63 = 0.0070910034628469103*tmp_57;
      real_t tmp_64 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id3;
      real_t tmp_65 = 0.0070910034628469103*tmp_57;
      real_t tmp_66 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id4;
      real_t tmp_67 = 0.0070910034628469103*tmp_57;
      real_t tmp_68 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id5;
      real_t tmp_69 = 0.0070910034628469103*tmp_57;
      real_t tmp_70 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id6;
      real_t tmp_71 = 0.018781320953002646*tmp_57;
      real_t tmp_72 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id7;
      real_t tmp_73 = 0.012248840519393657*tmp_57;
      real_t tmp_74 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id8;
      real_t tmp_75 = 0.018781320953002646*tmp_57;
      real_t tmp_76 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id9;
      real_t tmp_77 = 0.012248840519393657*tmp_57;
      real_t tmp_78 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id10;
      real_t tmp_79 = 0.018781320953002646*tmp_57;
      real_t tmp_80 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id11;
      real_t tmp_81 = 0.012248840519393657*tmp_57;
      real_t tmp_82 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id12;
      real_t tmp_83 = 0.0070910034628469103*tmp_57;
      real_t tmp_84 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id13;
      real_t tmp_85 = 0.0070910034628469103*tmp_57;
      real_t tmp_86 = 2.0*tmp_14;
      real_t tmp_87 = tmp_6*tmp_86;
      real_t tmp_88 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_86;
      real_t tmp_89 = tmp_30*tmp_35;
      real_t tmp_90 = tmp_41*tmp_43;
      real_t tmp_91 = tmp_21*tmp_87;
      real_t tmp_92 = Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_86;
      real_t tmp_93 = Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_86;
      real_t tmp_94 = Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_86;
      real_t tmp_95 = Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_86;
      real_t tmp_96 = Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_86;
      real_t tmp_97 = Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_86;
      real_t tmp_98 = Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_86;
      real_t tmp_99 = Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_86;
      real_t tmp_100 = Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_86;
      real_t tmp_101 = Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_86;
      real_t tmp_102 = Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_86;
      real_t tmp_103 = Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_86;
      real_t tmp_104 = Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_86;
      real_t tmp_105 = tmp_22*tmp_86;
      real_t tmp_106 = tmp_28*tmp_35;
      real_t tmp_107 = tmp_39*tmp_43;
      real_t tmp_108 = tmp_17*tmp_21;
      real_t tmp_109 = tmp_25*tmp_35;
      real_t tmp_110 = tmp_37*tmp_43;
      real_t tmp_111 = tmp_19*tmp_21;
      real_t a_0_0 = tmp_58*(tmp_22*tmp_23 + tmp_24*tmp_36 + tmp_24*tmp_44) + tmp_61*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_59 + tmp_36*tmp_60 + tmp_44*tmp_60) + tmp_63*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_59 + tmp_36*tmp_62 + tmp_44*tmp_62) + tmp_65*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_59 + tmp_36*tmp_64 + tmp_44*tmp_64) + tmp_67*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_59 + tmp_36*tmp_66 + tmp_44*tmp_66) + tmp_69*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_59 + tmp_36*tmp_68 + tmp_44*tmp_68) + tmp_71*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_59 + tmp_36*tmp_70 + tmp_44*tmp_70) + tmp_73*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_59 + tmp_36*tmp_72 + tmp_44*tmp_72) + tmp_75*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_59 + tmp_36*tmp_74 + tmp_44*tmp_74) + tmp_77*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_59 + tmp_36*tmp_76 + tmp_44*tmp_76) + tmp_79*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_59 + tmp_36*tmp_78 + tmp_44*tmp_78) + tmp_81*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_59 + tmp_36*tmp_80 + tmp_44*tmp_80) + tmp_83*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_59 + tmp_36*tmp_82 + tmp_44*tmp_82) + tmp_85*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_59 + tmp_36*tmp_84 + tmp_44*tmp_84);
      real_t a_1_0 = tmp_58*(tmp_22*tmp_87 + tmp_88*tmp_89 + tmp_88*tmp_90) + tmp_61*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_91 + tmp_89*tmp_92 + tmp_90*tmp_92) + tmp_63*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_91 + tmp_89*tmp_93 + tmp_90*tmp_93) + tmp_65*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_91 + tmp_89*tmp_94 + tmp_90*tmp_94) + tmp_67*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_91 + tmp_89*tmp_95 + tmp_90*tmp_95) + tmp_69*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_91 + tmp_89*tmp_96 + tmp_90*tmp_96) + tmp_71*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_91 + tmp_89*tmp_97 + tmp_90*tmp_97) + tmp_73*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_91 + tmp_89*tmp_98 + tmp_90*tmp_98) + tmp_75*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_91 + tmp_89*tmp_99 + tmp_90*tmp_99) + tmp_77*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_91 + tmp_100*tmp_89 + tmp_100*tmp_90) + tmp_79*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_91 + tmp_101*tmp_89 + tmp_101*tmp_90) + tmp_81*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_91 + tmp_102*tmp_89 + tmp_102*tmp_90) + tmp_83*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_91 + tmp_103*tmp_89 + tmp_103*tmp_90) + tmp_85*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_91 + tmp_104*tmp_89 + tmp_104*tmp_90);
      real_t a_2_0 = tmp_58*(tmp_105*tmp_17 + tmp_106*tmp_88 + tmp_107*tmp_88) + tmp_61*(tmp_106*tmp_92 + tmp_107*tmp_92 + tmp_108*tmp_92) + tmp_63*(tmp_106*tmp_93 + tmp_107*tmp_93 + tmp_108*tmp_93) + tmp_65*(tmp_106*tmp_94 + tmp_107*tmp_94 + tmp_108*tmp_94) + tmp_67*(tmp_106*tmp_95 + tmp_107*tmp_95 + tmp_108*tmp_95) + tmp_69*(tmp_106*tmp_96 + tmp_107*tmp_96 + tmp_108*tmp_96) + tmp_71*(tmp_106*tmp_97 + tmp_107*tmp_97 + tmp_108*tmp_97) + tmp_73*(tmp_106*tmp_98 + tmp_107*tmp_98 + tmp_108*tmp_98) + tmp_75*(tmp_106*tmp_99 + tmp_107*tmp_99 + tmp_108*tmp_99) + tmp_77*(tmp_100*tmp_106 + tmp_100*tmp_107 + tmp_100*tmp_108) + tmp_79*(tmp_101*tmp_106 + tmp_101*tmp_107 + tmp_101*tmp_108) + tmp_81*(tmp_102*tmp_106 + tmp_102*tmp_107 + tmp_102*tmp_108) + tmp_83*(tmp_103*tmp_106 + tmp_103*tmp_107 + tmp_103*tmp_108) + tmp_85*(tmp_104*tmp_106 + tmp_104*tmp_107 + tmp_104*tmp_108);
      real_t a_3_0 = tmp_58*(tmp_105*tmp_19 + tmp_109*tmp_88 + tmp_110*tmp_88) + tmp_61*(tmp_109*tmp_92 + tmp_110*tmp_92 + tmp_111*tmp_92) + tmp_63*(tmp_109*tmp_93 + tmp_110*tmp_93 + tmp_111*tmp_93) + tmp_65*(tmp_109*tmp_94 + tmp_110*tmp_94 + tmp_111*tmp_94) + tmp_67*(tmp_109*tmp_95 + tmp_110*tmp_95 + tmp_111*tmp_95) + tmp_69*(tmp_109*tmp_96 + tmp_110*tmp_96 + tmp_111*tmp_96) + tmp_71*(tmp_109*tmp_97 + tmp_110*tmp_97 + tmp_111*tmp_97) + tmp_73*(tmp_109*tmp_98 + tmp_110*tmp_98 + tmp_111*tmp_98) + tmp_75*(tmp_109*tmp_99 + tmp_110*tmp_99 + tmp_111*tmp_99) + tmp_77*(tmp_100*tmp_109 + tmp_100*tmp_110 + tmp_100*tmp_111) + tmp_79*(tmp_101*tmp_109 + tmp_101*tmp_110 + tmp_101*tmp_111) + tmp_81*(tmp_102*tmp_109 + tmp_102*tmp_110 + tmp_102*tmp_111) + tmp_83*(tmp_103*tmp_109 + tmp_103*tmp_110 + tmp_103*tmp_111) + tmp_85*(tmp_104*tmp_109 + tmp_104*tmp_110 + tmp_104*tmp_111);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_6*tmp_9;
      real_t tmp_14 = tmp_4*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_12 + tmp_0*tmp_7 - tmp_1*tmp_13 + tmp_1*tmp_2*tmp_8 + tmp_11*tmp_3 - tmp_14*tmp_3);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_7;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = tmp_33 - 1.0/4.0;
      real_t tmp_35 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_36 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_37 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_38 = tmp_20*tmp_35 + tmp_26*tmp_36 + tmp_32*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/4.0;
      real_t tmp_40 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_41 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_42 = tmp_11 - tmp_14;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_26*tmp_41 + tmp_32*tmp_42;
      real_t tmp_44 = tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_34 + tmp_1*tmp_39 + tmp_3*tmp_44;
      real_t tmp_46 = 0.5*tmp_15;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = tmp_37*tmp_46;
      real_t tmp_49 = tmp_27*tmp_46;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = 0.5*p_affine_13_1;
      real_t tmp_52 = tmp_50*tmp_51;
      real_t tmp_53 = tmp_10*tmp_39 + tmp_34*tmp_8 + tmp_44*tmp_6;
      real_t tmp_54 = tmp_40*tmp_46;
      real_t tmp_55 = tmp_35*tmp_46;
      real_t tmp_56 = tmp_46*tmp_5;
      real_t tmp_57 = -tmp_54 - tmp_55 - tmp_56;
      real_t tmp_58 = tmp_51*tmp_57;
      real_t tmp_59 = tmp_2*tmp_44 + tmp_34*tmp_9 + tmp_39*tmp_4;
      real_t tmp_60 = 1.0*tmp_15;
      real_t tmp_61 = tmp_41*tmp_60;
      real_t tmp_62 = tmp_36*tmp_60;
      real_t tmp_63 = tmp_21*tmp_60;
      real_t tmp_64 = 0.5*p_affine_13_0*tmp_50 + 0.5*p_affine_13_1*(-tmp_61 - tmp_62 - tmp_63) + 0.5*p_affine_13_2*tmp_57;
      real_t tmp_65 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_66 = std::pow(tmp_65, -0.25);
      real_t tmp_67 = -tmp_33 - tmp_38 - tmp_43 + 1;
      real_t tmp_68 = tmp_21*tmp_46;
      real_t tmp_69 = tmp_36*tmp_46;
      real_t tmp_70 = tmp_41*tmp_46;
      real_t tmp_71 = 0.5*p_affine_13_0*(tmp_0*tmp_68 + tmp_1*tmp_69 + tmp_2*tmp_47 + tmp_3*tmp_70 + tmp_4*tmp_48 + tmp_49*tmp_9) + 0.5*p_affine_13_1*(tmp_2*tmp_61 + tmp_4*tmp_62 + tmp_63*tmp_9) + 0.5*p_affine_13_2*(tmp_10*tmp_69 + tmp_2*tmp_54 + tmp_4*tmp_55 + tmp_56*tmp_9 + tmp_6*tmp_70 + tmp_68*tmp_8);
      real_t tmp_72 = 2.0*std::pow(tmp_65, 1.0/2.0);
      real_t tmp_73 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_72;
      real_t tmp_74 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_75 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_76 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_77 = tmp_21*tmp_75 + tmp_27*tmp_76 + tmp_5*tmp_74;
      real_t tmp_78 = tmp_77 - 1.0/4.0;
      real_t tmp_79 = tmp_35*tmp_74 + tmp_36*tmp_75 + tmp_37*tmp_76;
      real_t tmp_80 = tmp_79 - 1.0/4.0;
      real_t tmp_81 = tmp_40*tmp_74 + tmp_41*tmp_75 + tmp_42*tmp_76;
      real_t tmp_82 = tmp_81 - 1.0/4.0;
      real_t tmp_83 = tmp_0*tmp_78 + tmp_1*tmp_80 + tmp_3*tmp_82;
      real_t tmp_84 = tmp_10*tmp_80 + tmp_6*tmp_82 + tmp_78*tmp_8;
      real_t tmp_85 = tmp_2*tmp_82 + tmp_4*tmp_80 + tmp_78*tmp_9;
      real_t tmp_86 = -tmp_77 - tmp_79 - tmp_81 + 1;
      real_t tmp_87 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_72;
      real_t tmp_88 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_89 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_90 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_91 = tmp_21*tmp_89 + tmp_27*tmp_90 + tmp_5*tmp_88;
      real_t tmp_92 = tmp_91 - 1.0/4.0;
      real_t tmp_93 = tmp_35*tmp_88 + tmp_36*tmp_89 + tmp_37*tmp_90;
      real_t tmp_94 = tmp_93 - 1.0/4.0;
      real_t tmp_95 = tmp_40*tmp_88 + tmp_41*tmp_89 + tmp_42*tmp_90;
      real_t tmp_96 = tmp_95 - 1.0/4.0;
      real_t tmp_97 = tmp_0*tmp_92 + tmp_1*tmp_94 + tmp_3*tmp_96;
      real_t tmp_98 = tmp_10*tmp_94 + tmp_6*tmp_96 + tmp_8*tmp_92;
      real_t tmp_99 = tmp_2*tmp_96 + tmp_4*tmp_94 + tmp_9*tmp_92;
      real_t tmp_100 = -tmp_91 - tmp_93 - tmp_95 + 1;
      real_t tmp_101 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_72;
      real_t tmp_102 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_103 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_104 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_105 = tmp_102*tmp_5 + tmp_103*tmp_21 + tmp_104*tmp_27;
      real_t tmp_106 = tmp_105 - 1.0/4.0;
      real_t tmp_107 = tmp_102*tmp_35 + tmp_103*tmp_36 + tmp_104*tmp_37;
      real_t tmp_108 = tmp_107 - 1.0/4.0;
      real_t tmp_109 = tmp_102*tmp_40 + tmp_103*tmp_41 + tmp_104*tmp_42;
      real_t tmp_110 = tmp_109 - 1.0/4.0;
      real_t tmp_111 = tmp_0*tmp_106 + tmp_1*tmp_108 + tmp_110*tmp_3;
      real_t tmp_112 = tmp_10*tmp_108 + tmp_106*tmp_8 + tmp_110*tmp_6;
      real_t tmp_113 = tmp_106*tmp_9 + tmp_108*tmp_4 + tmp_110*tmp_2;
      real_t tmp_114 = -tmp_105 - tmp_107 - tmp_109 + 1;
      real_t tmp_115 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_72;
      real_t tmp_116 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_117 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_118 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_119 = tmp_116*tmp_5 + tmp_117*tmp_21 + tmp_118*tmp_27;
      real_t tmp_120 = tmp_119 - 1.0/4.0;
      real_t tmp_121 = tmp_116*tmp_35 + tmp_117*tmp_36 + tmp_118*tmp_37;
      real_t tmp_122 = tmp_121 - 1.0/4.0;
      real_t tmp_123 = tmp_116*tmp_40 + tmp_117*tmp_41 + tmp_118*tmp_42;
      real_t tmp_124 = tmp_123 - 1.0/4.0;
      real_t tmp_125 = tmp_0*tmp_120 + tmp_1*tmp_122 + tmp_124*tmp_3;
      real_t tmp_126 = tmp_10*tmp_122 + tmp_120*tmp_8 + tmp_124*tmp_6;
      real_t tmp_127 = tmp_120*tmp_9 + tmp_122*tmp_4 + tmp_124*tmp_2;
      real_t tmp_128 = -tmp_119 - tmp_121 - tmp_123 + 1;
      real_t tmp_129 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_72;
      real_t tmp_130 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_131 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_132 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_133 = tmp_130*tmp_5 + tmp_131*tmp_21 + tmp_132*tmp_27;
      real_t tmp_134 = tmp_133 - 1.0/4.0;
      real_t tmp_135 = tmp_130*tmp_35 + tmp_131*tmp_36 + tmp_132*tmp_37;
      real_t tmp_136 = tmp_135 - 1.0/4.0;
      real_t tmp_137 = tmp_130*tmp_40 + tmp_131*tmp_41 + tmp_132*tmp_42;
      real_t tmp_138 = tmp_137 - 1.0/4.0;
      real_t tmp_139 = tmp_0*tmp_134 + tmp_1*tmp_136 + tmp_138*tmp_3;
      real_t tmp_140 = tmp_10*tmp_136 + tmp_134*tmp_8 + tmp_138*tmp_6;
      real_t tmp_141 = tmp_134*tmp_9 + tmp_136*tmp_4 + tmp_138*tmp_2;
      real_t tmp_142 = -tmp_133 - tmp_135 - tmp_137 + 1;
      real_t tmp_143 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_72;
      real_t tmp_144 = 0.25*p_affine_13_1*tmp_15;
      real_t tmp_145 = tmp_144*tmp_5;
      real_t tmp_146 = tmp_144*tmp_27;
      real_t tmp_147 = 0.5*p_affine_13_0*tmp_49 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_56;
      real_t tmp_148 = tmp_144*tmp_35;
      real_t tmp_149 = tmp_144*tmp_37;
      real_t tmp_150 = 0.5*p_affine_13_0*tmp_48 + 0.5*p_affine_13_1*tmp_62 + 0.5*p_affine_13_2*tmp_55;
      real_t tmp_151 = tmp_144*tmp_40;
      real_t tmp_152 = tmp_144*tmp_42;
      real_t tmp_153 = 0.5*p_affine_13_0*tmp_47 + 0.5*p_affine_13_1*tmp_61 + 0.5*p_affine_13_2*tmp_54;
      real_t a_0_0 = tmp_101*(1.0*tmp_100*tmp_66*tmp_99 - tmp_100*tmp_71 - tmp_52*tmp_97 - tmp_58*tmp_98 - tmp_64*tmp_99) + tmp_115*(-tmp_111*tmp_52 - tmp_112*tmp_58 + 1.0*tmp_113*tmp_114*tmp_66 - tmp_113*tmp_64 - tmp_114*tmp_71) + tmp_129*(-tmp_125*tmp_52 - tmp_126*tmp_58 + 1.0*tmp_127*tmp_128*tmp_66 - tmp_127*tmp_64 - tmp_128*tmp_71) + tmp_143*(-tmp_139*tmp_52 - tmp_140*tmp_58 + 1.0*tmp_141*tmp_142*tmp_66 - tmp_141*tmp_64 - tmp_142*tmp_71) + tmp_73*(-tmp_45*tmp_52 - tmp_53*tmp_58 - tmp_59*tmp_64 + 1.0*tmp_59*tmp_66*tmp_67 - tmp_67*tmp_71) + tmp_87*(-tmp_52*tmp_83 - tmp_58*tmp_84 - tmp_64*tmp_85 + 1.0*tmp_66*tmp_85*tmp_86 - tmp_71*tmp_86);
      real_t a_1_0 = tmp_101*(-tmp_145*tmp_98 - tmp_146*tmp_97 - tmp_147*tmp_99 + 1.0*tmp_66*tmp_91*tmp_99 - tmp_71*tmp_91) + tmp_115*(1.0*tmp_105*tmp_113*tmp_66 - tmp_105*tmp_71 - tmp_111*tmp_146 - tmp_112*tmp_145 - tmp_113*tmp_147) + tmp_129*(1.0*tmp_119*tmp_127*tmp_66 - tmp_119*tmp_71 - tmp_125*tmp_146 - tmp_126*tmp_145 - tmp_127*tmp_147) + tmp_143*(1.0*tmp_133*tmp_141*tmp_66 - tmp_133*tmp_71 - tmp_139*tmp_146 - tmp_140*tmp_145 - tmp_141*tmp_147) + tmp_73*(-tmp_145*tmp_53 - tmp_146*tmp_45 - tmp_147*tmp_59 + 1.0*tmp_33*tmp_59*tmp_66 - tmp_33*tmp_71) + tmp_87*(-tmp_145*tmp_84 - tmp_146*tmp_83 - tmp_147*tmp_85 + 1.0*tmp_66*tmp_77*tmp_85 - tmp_71*tmp_77);
      real_t a_2_0 = tmp_101*(-tmp_148*tmp_98 - tmp_149*tmp_97 - tmp_150*tmp_99 + 1.0*tmp_66*tmp_93*tmp_99 - tmp_71*tmp_93) + tmp_115*(1.0*tmp_107*tmp_113*tmp_66 - tmp_107*tmp_71 - tmp_111*tmp_149 - tmp_112*tmp_148 - tmp_113*tmp_150) + tmp_129*(1.0*tmp_121*tmp_127*tmp_66 - tmp_121*tmp_71 - tmp_125*tmp_149 - tmp_126*tmp_148 - tmp_127*tmp_150) + tmp_143*(1.0*tmp_135*tmp_141*tmp_66 - tmp_135*tmp_71 - tmp_139*tmp_149 - tmp_140*tmp_148 - tmp_141*tmp_150) + tmp_73*(-tmp_148*tmp_53 - tmp_149*tmp_45 - tmp_150*tmp_59 + 1.0*tmp_38*tmp_59*tmp_66 - tmp_38*tmp_71) + tmp_87*(-tmp_148*tmp_84 - tmp_149*tmp_83 - tmp_150*tmp_85 + 1.0*tmp_66*tmp_79*tmp_85 - tmp_71*tmp_79);
      real_t a_3_0 = tmp_101*(-tmp_151*tmp_98 - tmp_152*tmp_97 - tmp_153*tmp_99 + 1.0*tmp_66*tmp_95*tmp_99 - tmp_71*tmp_95) + tmp_115*(1.0*tmp_109*tmp_113*tmp_66 - tmp_109*tmp_71 - tmp_111*tmp_152 - tmp_112*tmp_151 - tmp_113*tmp_153) + tmp_129*(1.0*tmp_123*tmp_127*tmp_66 - tmp_123*tmp_71 - tmp_125*tmp_152 - tmp_126*tmp_151 - tmp_127*tmp_153) + tmp_143*(1.0*tmp_137*tmp_141*tmp_66 - tmp_137*tmp_71 - tmp_139*tmp_152 - tmp_140*tmp_151 - tmp_141*tmp_153) + tmp_73*(-tmp_151*tmp_53 - tmp_152*tmp_45 - tmp_153*tmp_59 + 1.0*tmp_43*tmp_59*tmp_66 - tmp_43*tmp_71) + tmp_87*(-tmp_151*tmp_84 - tmp_152*tmp_83 - tmp_153*tmp_85 + 1.0*tmp_66*tmp_81*tmp_85 - tmp_71*tmp_81);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_1 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_2 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_5 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_11 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_12 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_13 = tmp_12*tmp_2;
      real_t tmp_14 = tmp_1*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_13 + tmp_0*tmp_9 + tmp_10*tmp_3 - tmp_10*tmp_6 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_14);
      real_t tmp_16 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_17 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_18 = -tmp_17;
      real_t tmp_19 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_20 = 0.091576213509770743*tmp_18 + 0.81684757298045851*tmp_19;
      real_t tmp_21 = tmp_15*(tmp_16 + tmp_20);
      real_t tmp_22 = tmp_12*tmp_4 - tmp_14;
      real_t tmp_23 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_24 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_25 = -tmp_24;
      real_t tmp_26 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_27 = 0.091576213509770743*tmp_25 + 0.81684757298045851*tmp_26;
      real_t tmp_28 = tmp_15*(tmp_23 + tmp_27);
      real_t tmp_29 = -tmp_13 + tmp_9;
      real_t tmp_30 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_31 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_32 = -tmp_31;
      real_t tmp_33 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_34 = 0.091576213509770743*tmp_32 + 0.81684757298045851*tmp_33;
      real_t tmp_35 = tmp_15*(tmp_30 + tmp_34);
      real_t tmp_36 = tmp_21*tmp_7 + tmp_22*tmp_28 + tmp_29*tmp_35 - 1.0/4.0;
      real_t tmp_37 = -tmp_0*tmp_2 + tmp_11*tmp_4;
      real_t tmp_38 = tmp_0*tmp_8 - tmp_10*tmp_4;
      real_t tmp_39 = tmp_10*tmp_2 - tmp_11*tmp_8;
      real_t tmp_40 = tmp_21*tmp_37 + tmp_28*tmp_38 + tmp_35*tmp_39 - 1.0/4.0;
      real_t tmp_41 = tmp_0*tmp_5 - tmp_1*tmp_11;
      real_t tmp_42 = -tmp_0*tmp_12 + tmp_1*tmp_10;
      real_t tmp_43 = -tmp_10*tmp_5 + tmp_11*tmp_12;
      real_t tmp_44 = tmp_21*tmp_41 + tmp_28*tmp_42 + tmp_35*tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_36 + tmp_1*tmp_40 + tmp_4*tmp_44;
      real_t tmp_46 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_47 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_48 = tmp_46*tmp_47;
      real_t tmp_49 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_50 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_51 = tmp_49*tmp_50;
      real_t tmp_52 = tmp_48 - tmp_51;
      real_t tmp_53 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_54 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_55 = tmp_49*tmp_54;
      real_t tmp_56 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_57 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_58 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_59 = tmp_47*tmp_57;
      real_t tmp_60 = tmp_46*tmp_54;
      real_t tmp_61 = 1.0 / (tmp_48*tmp_58 + tmp_50*tmp_56*tmp_57 - tmp_51*tmp_58 + tmp_53*tmp_55 - tmp_53*tmp_59 - tmp_56*tmp_60);
      real_t tmp_62 = 0.5*tmp_61;
      real_t tmp_63 = tmp_52*tmp_62;
      real_t tmp_64 = tmp_50*tmp_57 - tmp_60;
      real_t tmp_65 = tmp_62*tmp_64;
      real_t tmp_66 = tmp_55 - tmp_59;
      real_t tmp_67 = tmp_62*tmp_66;
      real_t tmp_68 = -tmp_63 - tmp_65 - tmp_67;
      real_t tmp_69 = 0.5*p_affine_13_1;
      real_t tmp_70 = tmp_68*tmp_69;
      real_t tmp_71 = tmp_10*tmp_36 + tmp_12*tmp_40 + tmp_44*tmp_8;
      real_t tmp_72 = -tmp_46*tmp_56 + tmp_49*tmp_53;
      real_t tmp_73 = tmp_62*tmp_72;
      real_t tmp_74 = tmp_46*tmp_58 - tmp_53*tmp_57;
      real_t tmp_75 = tmp_62*tmp_74;
      real_t tmp_76 = -tmp_49*tmp_58 + tmp_56*tmp_57;
      real_t tmp_77 = tmp_62*tmp_76;
      real_t tmp_78 = -tmp_73 - tmp_75 - tmp_77;
      real_t tmp_79 = tmp_69*tmp_78;
      real_t tmp_80 = tmp_11*tmp_36 + tmp_2*tmp_44 + tmp_40*tmp_5;
      real_t tmp_81 = -tmp_47*tmp_53 + tmp_50*tmp_56;
      real_t tmp_82 = 1.0*tmp_61;
      real_t tmp_83 = tmp_81*tmp_82;
      real_t tmp_84 = -tmp_50*tmp_58 + tmp_53*tmp_54;
      real_t tmp_85 = tmp_82*tmp_84;
      real_t tmp_86 = tmp_47*tmp_58 - tmp_54*tmp_56;
      real_t tmp_87 = tmp_82*tmp_86;
      real_t tmp_88 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*(-tmp_83 - tmp_85 - tmp_87) + 0.5*p_affine_13_2*tmp_78;
      real_t tmp_89 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_90 = tmp_61*(tmp_20 + tmp_89);
      real_t tmp_91 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_92 = tmp_61*(tmp_27 + tmp_91);
      real_t tmp_93 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_94 = tmp_61*(tmp_34 + tmp_93);
      real_t tmp_95 = tmp_66*tmp_94 + tmp_76*tmp_90 + tmp_86*tmp_92;
      real_t tmp_96 = tmp_64*tmp_94 + tmp_74*tmp_90 + tmp_84*tmp_92;
      real_t tmp_97 = tmp_52*tmp_94 + tmp_72*tmp_90 + tmp_81*tmp_92;
      real_t tmp_98 = -tmp_95 - tmp_96 - tmp_97 + 1;
      real_t tmp_99 = (std::abs(tmp_17*tmp_26 - tmp_19*tmp_24)*std::abs(tmp_17*tmp_26 - tmp_19*tmp_24)) + (std::abs(tmp_17*tmp_33 - tmp_19*tmp_31)*std::abs(tmp_17*tmp_33 - tmp_19*tmp_31)) + (std::abs(tmp_24*tmp_33 - tmp_26*tmp_31)*std::abs(tmp_24*tmp_33 - tmp_26*tmp_31));
      real_t tmp_100 = 1.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_80;
      real_t tmp_102 = 1.0*tmp_15;
      real_t tmp_103 = 0.5*tmp_15;
      real_t tmp_104 = tmp_103*tmp_22;
      real_t tmp_105 = tmp_103*tmp_38;
      real_t tmp_106 = tmp_103*tmp_42;
      real_t tmp_107 = tmp_103*tmp_11;
      real_t tmp_108 = tmp_103*tmp_5;
      real_t tmp_109 = tmp_103*tmp_2;
      real_t tmp_110 = 0.5*p_affine_13_0*(tmp_0*tmp_104 + tmp_1*tmp_105 + tmp_106*tmp_4 + tmp_107*tmp_29 + tmp_108*tmp_39 + tmp_109*tmp_43) + 0.5*p_affine_13_1*(tmp_102*tmp_11*tmp_22 + tmp_102*tmp_2*tmp_42 + tmp_102*tmp_38*tmp_5) + 0.5*p_affine_13_2*(tmp_10*tmp_104 + tmp_105*tmp_12 + tmp_106*tmp_8 + tmp_107*tmp_7 + tmp_108*tmp_37 + tmp_109*tmp_41);
      real_t tmp_111 = 2.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_112 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_111;
      real_t tmp_113 = 0.44594849091596489*tmp_18 + 0.10810301816807022*tmp_19;
      real_t tmp_114 = tmp_15*(tmp_113 + tmp_16);
      real_t tmp_115 = 0.44594849091596489*tmp_25 + 0.10810301816807022*tmp_26;
      real_t tmp_116 = tmp_15*(tmp_115 + tmp_23);
      real_t tmp_117 = 0.44594849091596489*tmp_32 + 0.10810301816807022*tmp_33;
      real_t tmp_118 = tmp_15*(tmp_117 + tmp_30);
      real_t tmp_119 = tmp_114*tmp_7 + tmp_116*tmp_22 + tmp_118*tmp_29 - 1.0/4.0;
      real_t tmp_120 = tmp_114*tmp_37 + tmp_116*tmp_38 + tmp_118*tmp_39 - 1.0/4.0;
      real_t tmp_121 = tmp_114*tmp_41 + tmp_116*tmp_42 + tmp_118*tmp_43 - 1.0/4.0;
      real_t tmp_122 = tmp_0*tmp_119 + tmp_1*tmp_120 + tmp_121*tmp_4;
      real_t tmp_123 = tmp_10*tmp_119 + tmp_12*tmp_120 + tmp_121*tmp_8;
      real_t tmp_124 = tmp_11*tmp_119 + tmp_120*tmp_5 + tmp_121*tmp_2;
      real_t tmp_125 = tmp_61*(tmp_113 + tmp_89);
      real_t tmp_126 = tmp_61*(tmp_115 + tmp_91);
      real_t tmp_127 = tmp_61*(tmp_117 + tmp_93);
      real_t tmp_128 = tmp_125*tmp_76 + tmp_126*tmp_86 + tmp_127*tmp_66;
      real_t tmp_129 = tmp_125*tmp_74 + tmp_126*tmp_84 + tmp_127*tmp_64;
      real_t tmp_130 = tmp_125*tmp_72 + tmp_126*tmp_81 + tmp_127*tmp_52;
      real_t tmp_131 = -tmp_128 - tmp_129 - tmp_130 + 1;
      real_t tmp_132 = tmp_100*tmp_124;
      real_t tmp_133 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_111;
      real_t tmp_134 = 0.81684757298045851*tmp_18 + 0.091576213509770743*tmp_19;
      real_t tmp_135 = tmp_15*(tmp_134 + tmp_16);
      real_t tmp_136 = 0.81684757298045851*tmp_25 + 0.091576213509770743*tmp_26;
      real_t tmp_137 = tmp_15*(tmp_136 + tmp_23);
      real_t tmp_138 = 0.81684757298045851*tmp_32 + 0.091576213509770743*tmp_33;
      real_t tmp_139 = tmp_15*(tmp_138 + tmp_30);
      real_t tmp_140 = tmp_135*tmp_7 + tmp_137*tmp_22 + tmp_139*tmp_29 - 1.0/4.0;
      real_t tmp_141 = tmp_135*tmp_37 + tmp_137*tmp_38 + tmp_139*tmp_39 - 1.0/4.0;
      real_t tmp_142 = tmp_135*tmp_41 + tmp_137*tmp_42 + tmp_139*tmp_43 - 1.0/4.0;
      real_t tmp_143 = tmp_0*tmp_140 + tmp_1*tmp_141 + tmp_142*tmp_4;
      real_t tmp_144 = tmp_10*tmp_140 + tmp_12*tmp_141 + tmp_142*tmp_8;
      real_t tmp_145 = tmp_11*tmp_140 + tmp_141*tmp_5 + tmp_142*tmp_2;
      real_t tmp_146 = tmp_61*(tmp_134 + tmp_89);
      real_t tmp_147 = tmp_61*(tmp_136 + tmp_91);
      real_t tmp_148 = tmp_61*(tmp_138 + tmp_93);
      real_t tmp_149 = tmp_146*tmp_76 + tmp_147*tmp_86 + tmp_148*tmp_66;
      real_t tmp_150 = tmp_146*tmp_74 + tmp_147*tmp_84 + tmp_148*tmp_64;
      real_t tmp_151 = tmp_146*tmp_72 + tmp_147*tmp_81 + tmp_148*tmp_52;
      real_t tmp_152 = -tmp_149 - tmp_150 - tmp_151 + 1;
      real_t tmp_153 = tmp_100*tmp_145;
      real_t tmp_154 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_111;
      real_t tmp_155 = 0.10810301816807022*tmp_18 + 0.44594849091596489*tmp_19;
      real_t tmp_156 = tmp_15*(tmp_155 + tmp_16);
      real_t tmp_157 = 0.10810301816807022*tmp_25 + 0.44594849091596489*tmp_26;
      real_t tmp_158 = tmp_15*(tmp_157 + tmp_23);
      real_t tmp_159 = 0.10810301816807022*tmp_32 + 0.44594849091596489*tmp_33;
      real_t tmp_160 = tmp_15*(tmp_159 + tmp_30);
      real_t tmp_161 = tmp_156*tmp_7 + tmp_158*tmp_22 + tmp_160*tmp_29 - 1.0/4.0;
      real_t tmp_162 = tmp_156*tmp_37 + tmp_158*tmp_38 + tmp_160*tmp_39 - 1.0/4.0;
      real_t tmp_163 = tmp_156*tmp_41 + tmp_158*tmp_42 + tmp_160*tmp_43 - 1.0/4.0;
      real_t tmp_164 = tmp_0*tmp_161 + tmp_1*tmp_162 + tmp_163*tmp_4;
      real_t tmp_165 = tmp_10*tmp_161 + tmp_12*tmp_162 + tmp_163*tmp_8;
      real_t tmp_166 = tmp_11*tmp_161 + tmp_162*tmp_5 + tmp_163*tmp_2;
      real_t tmp_167 = tmp_61*(tmp_155 + tmp_89);
      real_t tmp_168 = tmp_61*(tmp_157 + tmp_91);
      real_t tmp_169 = tmp_61*(tmp_159 + tmp_93);
      real_t tmp_170 = tmp_167*tmp_76 + tmp_168*tmp_86 + tmp_169*tmp_66;
      real_t tmp_171 = tmp_167*tmp_74 + tmp_168*tmp_84 + tmp_169*tmp_64;
      real_t tmp_172 = tmp_167*tmp_72 + tmp_168*tmp_81 + tmp_169*tmp_52;
      real_t tmp_173 = -tmp_170 - tmp_171 - tmp_172 + 1;
      real_t tmp_174 = tmp_100*tmp_166;
      real_t tmp_175 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_111;
      real_t tmp_176 = 0.091576213509770743*tmp_18 + 0.091576213509770743*tmp_19;
      real_t tmp_177 = tmp_15*(tmp_16 + tmp_176);
      real_t tmp_178 = 0.091576213509770743*tmp_25 + 0.091576213509770743*tmp_26;
      real_t tmp_179 = tmp_15*(tmp_178 + tmp_23);
      real_t tmp_180 = 0.091576213509770743*tmp_32 + 0.091576213509770743*tmp_33;
      real_t tmp_181 = tmp_15*(tmp_180 + tmp_30);
      real_t tmp_182 = tmp_177*tmp_7 + tmp_179*tmp_22 + tmp_181*tmp_29 - 1.0/4.0;
      real_t tmp_183 = tmp_177*tmp_37 + tmp_179*tmp_38 + tmp_181*tmp_39 - 1.0/4.0;
      real_t tmp_184 = tmp_177*tmp_41 + tmp_179*tmp_42 + tmp_181*tmp_43 - 1.0/4.0;
      real_t tmp_185 = tmp_0*tmp_182 + tmp_1*tmp_183 + tmp_184*tmp_4;
      real_t tmp_186 = tmp_10*tmp_182 + tmp_12*tmp_183 + tmp_184*tmp_8;
      real_t tmp_187 = tmp_11*tmp_182 + tmp_183*tmp_5 + tmp_184*tmp_2;
      real_t tmp_188 = tmp_61*(tmp_176 + tmp_89);
      real_t tmp_189 = tmp_61*(tmp_178 + tmp_91);
      real_t tmp_190 = tmp_61*(tmp_180 + tmp_93);
      real_t tmp_191 = tmp_188*tmp_76 + tmp_189*tmp_86 + tmp_190*tmp_66;
      real_t tmp_192 = tmp_188*tmp_74 + tmp_189*tmp_84 + tmp_190*tmp_64;
      real_t tmp_193 = tmp_188*tmp_72 + tmp_189*tmp_81 + tmp_190*tmp_52;
      real_t tmp_194 = -tmp_191 - tmp_192 - tmp_193 + 1;
      real_t tmp_195 = tmp_100*tmp_187;
      real_t tmp_196 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_111;
      real_t tmp_197 = 0.44594849091596489*tmp_18 + 0.44594849091596489*tmp_19;
      real_t tmp_198 = tmp_15*(tmp_16 + tmp_197);
      real_t tmp_199 = 0.44594849091596489*tmp_25 + 0.44594849091596489*tmp_26;
      real_t tmp_200 = tmp_15*(tmp_199 + tmp_23);
      real_t tmp_201 = 0.44594849091596489*tmp_32 + 0.44594849091596489*tmp_33;
      real_t tmp_202 = tmp_15*(tmp_201 + tmp_30);
      real_t tmp_203 = tmp_198*tmp_7 + tmp_200*tmp_22 + tmp_202*tmp_29 - 1.0/4.0;
      real_t tmp_204 = tmp_198*tmp_37 + tmp_200*tmp_38 + tmp_202*tmp_39 - 1.0/4.0;
      real_t tmp_205 = tmp_198*tmp_41 + tmp_200*tmp_42 + tmp_202*tmp_43 - 1.0/4.0;
      real_t tmp_206 = tmp_0*tmp_203 + tmp_1*tmp_204 + tmp_205*tmp_4;
      real_t tmp_207 = tmp_10*tmp_203 + tmp_12*tmp_204 + tmp_205*tmp_8;
      real_t tmp_208 = tmp_11*tmp_203 + tmp_2*tmp_205 + tmp_204*tmp_5;
      real_t tmp_209 = tmp_61*(tmp_197 + tmp_89);
      real_t tmp_210 = tmp_61*(tmp_199 + tmp_91);
      real_t tmp_211 = tmp_61*(tmp_201 + tmp_93);
      real_t tmp_212 = tmp_209*tmp_76 + tmp_210*tmp_86 + tmp_211*tmp_66;
      real_t tmp_213 = tmp_209*tmp_74 + tmp_210*tmp_84 + tmp_211*tmp_64;
      real_t tmp_214 = tmp_209*tmp_72 + tmp_210*tmp_81 + tmp_211*tmp_52;
      real_t tmp_215 = -tmp_212 - tmp_213 - tmp_214 + 1;
      real_t tmp_216 = tmp_100*tmp_208;
      real_t tmp_217 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_111;
      real_t tmp_218 = 0.25*p_affine_13_1*tmp_61;
      real_t tmp_219 = tmp_218*tmp_76;
      real_t tmp_220 = tmp_218*tmp_66;
      real_t tmp_221 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*tmp_87 + 0.5*p_affine_13_2*tmp_77;
      real_t tmp_222 = tmp_218*tmp_74;
      real_t tmp_223 = tmp_218*tmp_64;
      real_t tmp_224 = 0.5*p_affine_13_0*tmp_65 + 0.5*p_affine_13_1*tmp_85 + 0.5*p_affine_13_2*tmp_75;
      real_t tmp_225 = tmp_218*tmp_72;
      real_t tmp_226 = tmp_218*tmp_52;
      real_t tmp_227 = 0.5*p_affine_13_0*tmp_63 + 0.5*p_affine_13_1*tmp_83 + 0.5*p_affine_13_2*tmp_73;
      real_t a_0_0 = tmp_112*(-tmp_101*tmp_98 - tmp_110*tmp_98 + tmp_45*tmp_70 + tmp_71*tmp_79 + tmp_80*tmp_88) + tmp_133*(-tmp_110*tmp_131 + tmp_122*tmp_70 + tmp_123*tmp_79 + tmp_124*tmp_88 - tmp_131*tmp_132) + tmp_154*(-tmp_110*tmp_152 + tmp_143*tmp_70 + tmp_144*tmp_79 + tmp_145*tmp_88 - tmp_152*tmp_153) + tmp_175*(-tmp_110*tmp_173 + tmp_164*tmp_70 + tmp_165*tmp_79 + tmp_166*tmp_88 - tmp_173*tmp_174) + tmp_196*(-tmp_110*tmp_194 + tmp_185*tmp_70 + tmp_186*tmp_79 + tmp_187*tmp_88 - tmp_194*tmp_195) + tmp_217*(-tmp_110*tmp_215 + tmp_206*tmp_70 + tmp_207*tmp_79 + tmp_208*tmp_88 - tmp_215*tmp_216);
      real_t a_1_0 = tmp_112*(-tmp_101*tmp_95 - tmp_110*tmp_95 + tmp_219*tmp_71 + tmp_220*tmp_45 + tmp_221*tmp_80) + tmp_133*(-tmp_110*tmp_128 + tmp_122*tmp_220 + tmp_123*tmp_219 + tmp_124*tmp_221 - tmp_128*tmp_132) + tmp_154*(-tmp_110*tmp_149 + tmp_143*tmp_220 + tmp_144*tmp_219 + tmp_145*tmp_221 - tmp_149*tmp_153) + tmp_175*(-tmp_110*tmp_170 + tmp_164*tmp_220 + tmp_165*tmp_219 + tmp_166*tmp_221 - tmp_170*tmp_174) + tmp_196*(-tmp_110*tmp_191 + tmp_185*tmp_220 + tmp_186*tmp_219 + tmp_187*tmp_221 - tmp_191*tmp_195) + tmp_217*(-tmp_110*tmp_212 + tmp_206*tmp_220 + tmp_207*tmp_219 + tmp_208*tmp_221 - tmp_212*tmp_216);
      real_t a_2_0 = tmp_112*(-tmp_101*tmp_96 - tmp_110*tmp_96 + tmp_222*tmp_71 + tmp_223*tmp_45 + tmp_224*tmp_80) + tmp_133*(-tmp_110*tmp_129 + tmp_122*tmp_223 + tmp_123*tmp_222 + tmp_124*tmp_224 - tmp_129*tmp_132) + tmp_154*(-tmp_110*tmp_150 + tmp_143*tmp_223 + tmp_144*tmp_222 + tmp_145*tmp_224 - tmp_150*tmp_153) + tmp_175*(-tmp_110*tmp_171 + tmp_164*tmp_223 + tmp_165*tmp_222 + tmp_166*tmp_224 - tmp_171*tmp_174) + tmp_196*(-tmp_110*tmp_192 + tmp_185*tmp_223 + tmp_186*tmp_222 + tmp_187*tmp_224 - tmp_192*tmp_195) + tmp_217*(-tmp_110*tmp_213 + tmp_206*tmp_223 + tmp_207*tmp_222 + tmp_208*tmp_224 - tmp_213*tmp_216);
      real_t a_3_0 = tmp_112*(-tmp_101*tmp_97 - tmp_110*tmp_97 + tmp_225*tmp_71 + tmp_226*tmp_45 + tmp_227*tmp_80) + tmp_133*(-tmp_110*tmp_130 + tmp_122*tmp_226 + tmp_123*tmp_225 + tmp_124*tmp_227 - tmp_130*tmp_132) + tmp_154*(-tmp_110*tmp_151 + tmp_143*tmp_226 + tmp_144*tmp_225 + tmp_145*tmp_227 - tmp_151*tmp_153) + tmp_175*(-tmp_110*tmp_172 + tmp_164*tmp_226 + tmp_165*tmp_225 + tmp_166*tmp_227 - tmp_172*tmp_174) + tmp_196*(-tmp_110*tmp_193 + tmp_185*tmp_226 + tmp_186*tmp_225 + tmp_187*tmp_227 - tmp_193*tmp_195) + tmp_217*(-tmp_110*tmp_214 + tmp_206*tmp_226 + tmp_207*tmp_225 + tmp_208*tmp_227 - tmp_214*tmp_216);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_6*tmp_9;
      real_t tmp_14 = tmp_4*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_12 + tmp_0*tmp_7 - tmp_1*tmp_13 + tmp_1*tmp_2*tmp_8 + tmp_11*tmp_3 - tmp_14*tmp_3);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_7;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = tmp_33 - 1.0/4.0;
      real_t tmp_35 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_36 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_37 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_38 = tmp_20*tmp_35 + tmp_26*tmp_36 + tmp_32*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/4.0;
      real_t tmp_40 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_41 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_42 = tmp_11 - tmp_14;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_26*tmp_41 + tmp_32*tmp_42;
      real_t tmp_44 = tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_34 + tmp_1*tmp_39 + tmp_3*tmp_44;
      real_t tmp_46 = 0.5*tmp_15;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = tmp_37*tmp_46;
      real_t tmp_49 = tmp_27*tmp_46;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = p_affine_13_1*tmp_50;
      real_t tmp_52 = tmp_10*tmp_39 + tmp_34*tmp_8 + tmp_44*tmp_6;
      real_t tmp_53 = tmp_40*tmp_46;
      real_t tmp_54 = tmp_35*tmp_46;
      real_t tmp_55 = tmp_46*tmp_5;
      real_t tmp_56 = -tmp_53 - tmp_54 - tmp_55;
      real_t tmp_57 = p_affine_13_1*tmp_56;
      real_t tmp_58 = 1.0*tmp_15;
      real_t tmp_59 = tmp_41*tmp_58;
      real_t tmp_60 = tmp_36*tmp_58;
      real_t tmp_61 = tmp_21*tmp_58;
      real_t tmp_62 = p_affine_13_0*tmp_50 + p_affine_13_1*(-tmp_59 - tmp_60 - tmp_61) + p_affine_13_2*tmp_56;
      real_t tmp_63 = tmp_2*tmp_44 + tmp_34*tmp_9 + tmp_39*tmp_4;
      real_t tmp_64 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_65 = std::pow(tmp_64, -0.25);
      real_t tmp_66 = -tmp_33 - tmp_38 - tmp_43 + 1;
      real_t tmp_67 = tmp_21*tmp_46;
      real_t tmp_68 = tmp_36*tmp_46;
      real_t tmp_69 = tmp_41*tmp_46;
      real_t tmp_70 = p_affine_13_0*(tmp_0*tmp_67 + tmp_1*tmp_68 + tmp_2*tmp_47 + tmp_3*tmp_69 + tmp_4*tmp_48 + tmp_49*tmp_9) + p_affine_13_1*(tmp_2*tmp_59 + tmp_4*tmp_60 + tmp_61*tmp_9) + p_affine_13_2*(tmp_10*tmp_68 + tmp_2*tmp_53 + tmp_4*tmp_54 + tmp_55*tmp_9 + tmp_6*tmp_69 + tmp_67*tmp_8);
      real_t tmp_71 = 2.0*std::pow(tmp_64, 1.0/2.0);
      real_t tmp_72 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_71;
      real_t tmp_73 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_74 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_75 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_76 = tmp_21*tmp_74 + tmp_27*tmp_75 + tmp_5*tmp_73;
      real_t tmp_77 = tmp_76 - 1.0/4.0;
      real_t tmp_78 = tmp_35*tmp_73 + tmp_36*tmp_74 + tmp_37*tmp_75;
      real_t tmp_79 = tmp_78 - 1.0/4.0;
      real_t tmp_80 = tmp_40*tmp_73 + tmp_41*tmp_74 + tmp_42*tmp_75;
      real_t tmp_81 = tmp_80 - 1.0/4.0;
      real_t tmp_82 = tmp_0*tmp_77 + tmp_1*tmp_79 + tmp_3*tmp_81;
      real_t tmp_83 = tmp_10*tmp_79 + tmp_6*tmp_81 + tmp_77*tmp_8;
      real_t tmp_84 = tmp_2*tmp_81 + tmp_4*tmp_79 + tmp_77*tmp_9;
      real_t tmp_85 = -tmp_76 - tmp_78 - tmp_80 + 1;
      real_t tmp_86 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_71;
      real_t tmp_87 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_88 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_89 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_90 = tmp_21*tmp_88 + tmp_27*tmp_89 + tmp_5*tmp_87;
      real_t tmp_91 = tmp_90 - 1.0/4.0;
      real_t tmp_92 = tmp_35*tmp_87 + tmp_36*tmp_88 + tmp_37*tmp_89;
      real_t tmp_93 = tmp_92 - 1.0/4.0;
      real_t tmp_94 = tmp_40*tmp_87 + tmp_41*tmp_88 + tmp_42*tmp_89;
      real_t tmp_95 = tmp_94 - 1.0/4.0;
      real_t tmp_96 = tmp_0*tmp_91 + tmp_1*tmp_93 + tmp_3*tmp_95;
      real_t tmp_97 = tmp_10*tmp_93 + tmp_6*tmp_95 + tmp_8*tmp_91;
      real_t tmp_98 = tmp_2*tmp_95 + tmp_4*tmp_93 + tmp_9*tmp_91;
      real_t tmp_99 = -tmp_90 - tmp_92 - tmp_94 + 1;
      real_t tmp_100 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_71;
      real_t tmp_101 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_102 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_103 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_104 = tmp_101*tmp_5 + tmp_102*tmp_21 + tmp_103*tmp_27;
      real_t tmp_105 = tmp_104 - 1.0/4.0;
      real_t tmp_106 = tmp_101*tmp_35 + tmp_102*tmp_36 + tmp_103*tmp_37;
      real_t tmp_107 = tmp_106 - 1.0/4.0;
      real_t tmp_108 = tmp_101*tmp_40 + tmp_102*tmp_41 + tmp_103*tmp_42;
      real_t tmp_109 = tmp_108 - 1.0/4.0;
      real_t tmp_110 = tmp_0*tmp_105 + tmp_1*tmp_107 + tmp_109*tmp_3;
      real_t tmp_111 = tmp_10*tmp_107 + tmp_105*tmp_8 + tmp_109*tmp_6;
      real_t tmp_112 = tmp_105*tmp_9 + tmp_107*tmp_4 + tmp_109*tmp_2;
      real_t tmp_113 = -tmp_104 - tmp_106 - tmp_108 + 1;
      real_t tmp_114 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_71;
      real_t tmp_115 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_116 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_117 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_118 = tmp_115*tmp_5 + tmp_116*tmp_21 + tmp_117*tmp_27;
      real_t tmp_119 = tmp_118 - 1.0/4.0;
      real_t tmp_120 = tmp_115*tmp_35 + tmp_116*tmp_36 + tmp_117*tmp_37;
      real_t tmp_121 = tmp_120 - 1.0/4.0;
      real_t tmp_122 = tmp_115*tmp_40 + tmp_116*tmp_41 + tmp_117*tmp_42;
      real_t tmp_123 = tmp_122 - 1.0/4.0;
      real_t tmp_124 = tmp_0*tmp_119 + tmp_1*tmp_121 + tmp_123*tmp_3;
      real_t tmp_125 = tmp_10*tmp_121 + tmp_119*tmp_8 + tmp_123*tmp_6;
      real_t tmp_126 = tmp_119*tmp_9 + tmp_121*tmp_4 + tmp_123*tmp_2;
      real_t tmp_127 = -tmp_118 - tmp_120 - tmp_122 + 1;
      real_t tmp_128 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_71;
      real_t tmp_129 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_130 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_131 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_132 = tmp_129*tmp_5 + tmp_130*tmp_21 + tmp_131*tmp_27;
      real_t tmp_133 = tmp_132 - 1.0/4.0;
      real_t tmp_134 = tmp_129*tmp_35 + tmp_130*tmp_36 + tmp_131*tmp_37;
      real_t tmp_135 = tmp_134 - 1.0/4.0;
      real_t tmp_136 = tmp_129*tmp_40 + tmp_130*tmp_41 + tmp_131*tmp_42;
      real_t tmp_137 = tmp_136 - 1.0/4.0;
      real_t tmp_138 = tmp_0*tmp_133 + tmp_1*tmp_135 + tmp_137*tmp_3;
      real_t tmp_139 = tmp_10*tmp_135 + tmp_133*tmp_8 + tmp_137*tmp_6;
      real_t tmp_140 = tmp_133*tmp_9 + tmp_135*tmp_4 + tmp_137*tmp_2;
      real_t tmp_141 = -tmp_132 - tmp_134 - tmp_136 + 1;
      real_t tmp_142 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_71;
      real_t tmp_143 = p_affine_13_1*tmp_55;
      real_t tmp_144 = p_affine_13_1*tmp_49;
      real_t tmp_145 = p_affine_13_0*tmp_49 + p_affine_13_1*tmp_61 + p_affine_13_2*tmp_55;
      real_t tmp_146 = p_affine_13_1*tmp_54;
      real_t tmp_147 = p_affine_13_1*tmp_48;
      real_t tmp_148 = p_affine_13_0*tmp_48 + p_affine_13_1*tmp_60 + p_affine_13_2*tmp_54;
      real_t tmp_149 = p_affine_13_1*tmp_53;
      real_t tmp_150 = p_affine_13_1*tmp_47;
      real_t tmp_151 = p_affine_13_0*tmp_47 + p_affine_13_1*tmp_59 + p_affine_13_2*tmp_53;
      real_t a_0_0 = tmp_100*(-tmp_51*tmp_96 - tmp_57*tmp_97 - tmp_62*tmp_98 + 1.0*tmp_65*tmp_98*tmp_99 - tmp_70*tmp_99) + tmp_114*(-tmp_110*tmp_51 - tmp_111*tmp_57 + 1.0*tmp_112*tmp_113*tmp_65 - tmp_112*tmp_62 - tmp_113*tmp_70) + tmp_128*(-tmp_124*tmp_51 - tmp_125*tmp_57 + 1.0*tmp_126*tmp_127*tmp_65 - tmp_126*tmp_62 - tmp_127*tmp_70) + tmp_142*(-tmp_138*tmp_51 - tmp_139*tmp_57 + 1.0*tmp_140*tmp_141*tmp_65 - tmp_140*tmp_62 - tmp_141*tmp_70) + tmp_72*(-tmp_45*tmp_51 - tmp_52*tmp_57 - tmp_62*tmp_63 + 1.0*tmp_63*tmp_65*tmp_66 - tmp_66*tmp_70) + tmp_86*(-tmp_51*tmp_82 - tmp_57*tmp_83 - tmp_62*tmp_84 + 1.0*tmp_65*tmp_84*tmp_85 - tmp_70*tmp_85);
      real_t a_1_0 = tmp_100*(-tmp_143*tmp_97 - tmp_144*tmp_96 - tmp_145*tmp_98 + 1.0*tmp_65*tmp_90*tmp_98 - tmp_70*tmp_90) + tmp_114*(1.0*tmp_104*tmp_112*tmp_65 - tmp_104*tmp_70 - tmp_110*tmp_144 - tmp_111*tmp_143 - tmp_112*tmp_145) + tmp_128*(1.0*tmp_118*tmp_126*tmp_65 - tmp_118*tmp_70 - tmp_124*tmp_144 - tmp_125*tmp_143 - tmp_126*tmp_145) + tmp_142*(1.0*tmp_132*tmp_140*tmp_65 - tmp_132*tmp_70 - tmp_138*tmp_144 - tmp_139*tmp_143 - tmp_140*tmp_145) + tmp_72*(-tmp_143*tmp_52 - tmp_144*tmp_45 - tmp_145*tmp_63 + 1.0*tmp_33*tmp_63*tmp_65 - tmp_33*tmp_70) + tmp_86*(-tmp_143*tmp_83 - tmp_144*tmp_82 - tmp_145*tmp_84 + 1.0*tmp_65*tmp_76*tmp_84 - tmp_70*tmp_76);
      real_t a_2_0 = tmp_100*(-tmp_146*tmp_97 - tmp_147*tmp_96 - tmp_148*tmp_98 + 1.0*tmp_65*tmp_92*tmp_98 - tmp_70*tmp_92) + tmp_114*(1.0*tmp_106*tmp_112*tmp_65 - tmp_106*tmp_70 - tmp_110*tmp_147 - tmp_111*tmp_146 - tmp_112*tmp_148) + tmp_128*(1.0*tmp_120*tmp_126*tmp_65 - tmp_120*tmp_70 - tmp_124*tmp_147 - tmp_125*tmp_146 - tmp_126*tmp_148) + tmp_142*(1.0*tmp_134*tmp_140*tmp_65 - tmp_134*tmp_70 - tmp_138*tmp_147 - tmp_139*tmp_146 - tmp_140*tmp_148) + tmp_72*(-tmp_146*tmp_52 - tmp_147*tmp_45 - tmp_148*tmp_63 + 1.0*tmp_38*tmp_63*tmp_65 - tmp_38*tmp_70) + tmp_86*(-tmp_146*tmp_83 - tmp_147*tmp_82 - tmp_148*tmp_84 + 1.0*tmp_65*tmp_78*tmp_84 - tmp_70*tmp_78);
      real_t a_3_0 = tmp_100*(-tmp_149*tmp_97 - tmp_150*tmp_96 - tmp_151*tmp_98 + 1.0*tmp_65*tmp_94*tmp_98 - tmp_70*tmp_94) + tmp_114*(1.0*tmp_108*tmp_112*tmp_65 - tmp_108*tmp_70 - tmp_110*tmp_150 - tmp_111*tmp_149 - tmp_112*tmp_151) + tmp_128*(1.0*tmp_122*tmp_126*tmp_65 - tmp_122*tmp_70 - tmp_124*tmp_150 - tmp_125*tmp_149 - tmp_126*tmp_151) + tmp_142*(1.0*tmp_136*tmp_140*tmp_65 - tmp_136*tmp_70 - tmp_138*tmp_150 - tmp_139*tmp_149 - tmp_140*tmp_151) + tmp_72*(-tmp_149*tmp_52 - tmp_150*tmp_45 - tmp_151*tmp_63 + 1.0*tmp_43*tmp_63*tmp_65 - tmp_43*tmp_70) + tmp_86*(-tmp_149*tmp_83 - tmp_150*tmp_82 - tmp_151*tmp_84 + 1.0*tmp_65*tmp_80*tmp_84 - tmp_70*tmp_80);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }

};




class EGEpsilonFormEP1_2 : public hyteg::dg::DGForm
{

 public:
    EGEpsilonFormEP1_2(std::function< real_t ( const Point3D & ) > mu)
: callback_Scalar_Variable_Coefficient_3D_mu (mu)
    {}

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;



void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
}

      
 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElementInner[0]( 0 );
      const auto p_affine_0_1 = coordsElementInner[0]( 1 );

      const auto p_affine_1_0 = coordsElementInner[1]( 0 );
      const auto p_affine_1_1 = coordsElementInner[1]( 1 );

      const auto p_affine_2_0 = coordsElementInner[2]( 0 );
      const auto p_affine_2_1 = coordsElementInner[2]( 1 );

      const auto p_affine_3_0 = coordsElementOuter[0]( 0 );
      const auto p_affine_3_1 = coordsElementOuter[0]( 1 );

      const auto p_affine_4_0 = coordsElementOuter[1]( 0 );
      const auto p_affine_4_1 = coordsElementOuter[1]( 1 );

      const auto p_affine_5_0 = coordsElementOuter[2]( 0 );
      const auto p_affine_5_1 = coordsElementOuter[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

  void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      elMat( 0, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

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

      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id6 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id7 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id8 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id9 );
      Scalar_Variable_Coefficient_3D_mu( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id10 );
      Scalar_Variable_Coefficient_3D_mu( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id11 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id12 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id13 );
      real_t tmp_0 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_9 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_10 = tmp_5*tmp_9;
      real_t tmp_11 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_12 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_13 = tmp_2*tmp_9;
      real_t tmp_14 = tmp_1*tmp_12;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_3 - tmp_0*tmp_6 + tmp_10*tmp_8 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_13 - tmp_14*tmp_8);
      real_t tmp_16 = 1.0*tmp_15;
      real_t tmp_17 = tmp_16*tmp_7;
      real_t tmp_18 = tmp_12*tmp_4 - tmp_13;
      real_t tmp_19 = tmp_16*tmp_18;
      real_t tmp_20 = tmp_10 - tmp_14;
      real_t tmp_21 = tmp_16*tmp_20;
      real_t tmp_22 = tmp_0*tmp_17 + tmp_11*tmp_19 + tmp_21*tmp_8;
      real_t tmp_23 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_22;
      real_t tmp_24 = -2*tmp_17 - 2*tmp_19 - 2*tmp_21;
      real_t tmp_25 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_26 = tmp_0*tmp_1 - tmp_11*tmp_9;
      real_t tmp_27 = 0.5*tmp_15;
      real_t tmp_28 = tmp_26*tmp_27;
      real_t tmp_29 = -tmp_0*tmp_4 + tmp_8*tmp_9;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = -tmp_1*tmp_8 + tmp_11*tmp_4;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_27*tmp_7;
      real_t tmp_34 = tmp_18*tmp_27;
      real_t tmp_35 = tmp_20*tmp_27;
      real_t tmp_36 = tmp_0*tmp_32 + tmp_11*tmp_30 + tmp_12*tmp_33 + tmp_2*tmp_35 + tmp_28*tmp_8 + tmp_34*tmp_5;
      real_t tmp_37 = tmp_36*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_38 = -tmp_0*tmp_5 + tmp_11*tmp_12;
      real_t tmp_39 = tmp_27*tmp_38;
      real_t tmp_40 = tmp_0*tmp_2 - tmp_12*tmp_8;
      real_t tmp_41 = tmp_27*tmp_40;
      real_t tmp_42 = -tmp_11*tmp_2 + tmp_5*tmp_8;
      real_t tmp_43 = tmp_27*tmp_42;
      real_t tmp_44 = tmp_0*tmp_43 + tmp_1*tmp_34 + tmp_11*tmp_41 + tmp_33*tmp_9 + tmp_35*tmp_4 + tmp_39*tmp_8;
      real_t tmp_45 = tmp_44*(-tmp_39 - tmp_41 - tmp_43);
      real_t tmp_46 = p_affine_0_0*p_affine_1_1;
      real_t tmp_47 = p_affine_0_0*p_affine_1_2;
      real_t tmp_48 = p_affine_2_1*p_affine_3_2;
      real_t tmp_49 = p_affine_0_1*p_affine_1_0;
      real_t tmp_50 = p_affine_0_1*p_affine_1_2;
      real_t tmp_51 = p_affine_2_2*p_affine_3_0;
      real_t tmp_52 = p_affine_0_2*p_affine_1_0;
      real_t tmp_53 = p_affine_0_2*p_affine_1_1;
      real_t tmp_54 = p_affine_2_0*p_affine_3_1;
      real_t tmp_55 = p_affine_2_2*p_affine_3_1;
      real_t tmp_56 = p_affine_2_0*p_affine_3_2;
      real_t tmp_57 = p_affine_2_1*p_affine_3_0;
      real_t tmp_58 = std::abs(p_affine_0_0*tmp_48 - p_affine_0_0*tmp_55 + p_affine_0_1*tmp_51 - p_affine_0_1*tmp_56 + p_affine_0_2*tmp_54 - p_affine_0_2*tmp_57 - p_affine_1_0*tmp_48 + p_affine_1_0*tmp_55 - p_affine_1_1*tmp_51 + p_affine_1_1*tmp_56 - p_affine_1_2*tmp_54 + p_affine_1_2*tmp_57 + p_affine_2_0*tmp_50 - p_affine_2_0*tmp_53 - p_affine_2_1*tmp_47 + p_affine_2_1*tmp_52 + p_affine_2_2*tmp_46 - p_affine_2_2*tmp_49 - p_affine_3_0*tmp_50 + p_affine_3_0*tmp_53 + p_affine_3_1*tmp_47 - p_affine_3_1*tmp_52 - p_affine_3_2*tmp_46 + p_affine_3_2*tmp_49);
      real_t tmp_59 = 0.018781320953002646*tmp_58;
      real_t tmp_60 = tmp_22*tmp_24;
      real_t tmp_61 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id1;
      real_t tmp_62 = 0.012248840519393657*tmp_58;
      real_t tmp_63 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id2;
      real_t tmp_64 = 0.0070910034628469103*tmp_58;
      real_t tmp_65 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id3;
      real_t tmp_66 = 0.0070910034628469103*tmp_58;
      real_t tmp_67 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id4;
      real_t tmp_68 = 0.0070910034628469103*tmp_58;
      real_t tmp_69 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id5;
      real_t tmp_70 = 0.0070910034628469103*tmp_58;
      real_t tmp_71 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id6;
      real_t tmp_72 = 0.018781320953002646*tmp_58;
      real_t tmp_73 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id7;
      real_t tmp_74 = 0.012248840519393657*tmp_58;
      real_t tmp_75 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id8;
      real_t tmp_76 = 0.018781320953002646*tmp_58;
      real_t tmp_77 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id9;
      real_t tmp_78 = 0.012248840519393657*tmp_58;
      real_t tmp_79 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id10;
      real_t tmp_80 = 0.018781320953002646*tmp_58;
      real_t tmp_81 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id11;
      real_t tmp_82 = 0.012248840519393657*tmp_58;
      real_t tmp_83 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id12;
      real_t tmp_84 = 0.0070910034628469103*tmp_58;
      real_t tmp_85 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id13;
      real_t tmp_86 = 0.0070910034628469103*tmp_58;
      real_t tmp_87 = 2.0*tmp_15;
      real_t tmp_88 = tmp_7*tmp_87;
      real_t tmp_89 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_87;
      real_t tmp_90 = tmp_31*tmp_36;
      real_t tmp_91 = tmp_42*tmp_44;
      real_t tmp_92 = tmp_22*tmp_88;
      real_t tmp_93 = Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_87;
      real_t tmp_94 = Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_87;
      real_t tmp_95 = Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_87;
      real_t tmp_96 = Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_87;
      real_t tmp_97 = Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_87;
      real_t tmp_98 = Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_87;
      real_t tmp_99 = Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_87;
      real_t tmp_100 = Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_87;
      real_t tmp_101 = Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_87;
      real_t tmp_102 = Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_87;
      real_t tmp_103 = Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_87;
      real_t tmp_104 = Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_87;
      real_t tmp_105 = Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_87;
      real_t tmp_106 = tmp_23*tmp_87;
      real_t tmp_107 = tmp_29*tmp_36;
      real_t tmp_108 = tmp_40*tmp_44;
      real_t tmp_109 = tmp_18*tmp_22;
      real_t tmp_110 = tmp_26*tmp_36;
      real_t tmp_111 = tmp_38*tmp_44;
      real_t tmp_112 = tmp_20*tmp_22;
      real_t a_0_0 = tmp_59*(tmp_23*tmp_24 + tmp_25*tmp_37 + tmp_25*tmp_45) + tmp_62*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_60 + tmp_37*tmp_61 + tmp_45*tmp_61) + tmp_64*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_60 + tmp_37*tmp_63 + tmp_45*tmp_63) + tmp_66*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_60 + tmp_37*tmp_65 + tmp_45*tmp_65) + tmp_68*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_60 + tmp_37*tmp_67 + tmp_45*tmp_67) + tmp_70*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_60 + tmp_37*tmp_69 + tmp_45*tmp_69) + tmp_72*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_60 + tmp_37*tmp_71 + tmp_45*tmp_71) + tmp_74*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_60 + tmp_37*tmp_73 + tmp_45*tmp_73) + tmp_76*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_60 + tmp_37*tmp_75 + tmp_45*tmp_75) + tmp_78*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_60 + tmp_37*tmp_77 + tmp_45*tmp_77) + tmp_80*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_60 + tmp_37*tmp_79 + tmp_45*tmp_79) + tmp_82*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_60 + tmp_37*tmp_81 + tmp_45*tmp_81) + tmp_84*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_60 + tmp_37*tmp_83 + tmp_45*tmp_83) + tmp_86*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_60 + tmp_37*tmp_85 + tmp_45*tmp_85);
      real_t a_0_1 = tmp_59*(tmp_23*tmp_88 + tmp_89*tmp_90 + tmp_89*tmp_91) + tmp_62*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_92 + tmp_90*tmp_93 + tmp_91*tmp_93) + tmp_64*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_92 + tmp_90*tmp_94 + tmp_91*tmp_94) + tmp_66*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_92 + tmp_90*tmp_95 + tmp_91*tmp_95) + tmp_68*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_92 + tmp_90*tmp_96 + tmp_91*tmp_96) + tmp_70*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_92 + tmp_90*tmp_97 + tmp_91*tmp_97) + tmp_72*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_92 + tmp_90*tmp_98 + tmp_91*tmp_98) + tmp_74*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_92 + tmp_90*tmp_99 + tmp_91*tmp_99) + tmp_76*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_92 + tmp_100*tmp_90 + tmp_100*tmp_91) + tmp_78*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_92 + tmp_101*tmp_90 + tmp_101*tmp_91) + tmp_80*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_92 + tmp_102*tmp_90 + tmp_102*tmp_91) + tmp_82*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_92 + tmp_103*tmp_90 + tmp_103*tmp_91) + tmp_84*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_92 + tmp_104*tmp_90 + tmp_104*tmp_91) + tmp_86*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_92 + tmp_105*tmp_90 + tmp_105*tmp_91);
      real_t a_0_2 = tmp_59*(tmp_106*tmp_18 + tmp_107*tmp_89 + tmp_108*tmp_89) + tmp_62*(tmp_107*tmp_93 + tmp_108*tmp_93 + tmp_109*tmp_93) + tmp_64*(tmp_107*tmp_94 + tmp_108*tmp_94 + tmp_109*tmp_94) + tmp_66*(tmp_107*tmp_95 + tmp_108*tmp_95 + tmp_109*tmp_95) + tmp_68*(tmp_107*tmp_96 + tmp_108*tmp_96 + tmp_109*tmp_96) + tmp_70*(tmp_107*tmp_97 + tmp_108*tmp_97 + tmp_109*tmp_97) + tmp_72*(tmp_107*tmp_98 + tmp_108*tmp_98 + tmp_109*tmp_98) + tmp_74*(tmp_107*tmp_99 + tmp_108*tmp_99 + tmp_109*tmp_99) + tmp_76*(tmp_100*tmp_107 + tmp_100*tmp_108 + tmp_100*tmp_109) + tmp_78*(tmp_101*tmp_107 + tmp_101*tmp_108 + tmp_101*tmp_109) + tmp_80*(tmp_102*tmp_107 + tmp_102*tmp_108 + tmp_102*tmp_109) + tmp_82*(tmp_103*tmp_107 + tmp_103*tmp_108 + tmp_103*tmp_109) + tmp_84*(tmp_104*tmp_107 + tmp_104*tmp_108 + tmp_104*tmp_109) + tmp_86*(tmp_105*tmp_107 + tmp_105*tmp_108 + tmp_105*tmp_109);
      real_t a_0_3 = tmp_59*(tmp_106*tmp_20 + tmp_110*tmp_89 + tmp_111*tmp_89) + tmp_62*(tmp_110*tmp_93 + tmp_111*tmp_93 + tmp_112*tmp_93) + tmp_64*(tmp_110*tmp_94 + tmp_111*tmp_94 + tmp_112*tmp_94) + tmp_66*(tmp_110*tmp_95 + tmp_111*tmp_95 + tmp_112*tmp_95) + tmp_68*(tmp_110*tmp_96 + tmp_111*tmp_96 + tmp_112*tmp_96) + tmp_70*(tmp_110*tmp_97 + tmp_111*tmp_97 + tmp_112*tmp_97) + tmp_72*(tmp_110*tmp_98 + tmp_111*tmp_98 + tmp_112*tmp_98) + tmp_74*(tmp_110*tmp_99 + tmp_111*tmp_99 + tmp_112*tmp_99) + tmp_76*(tmp_100*tmp_110 + tmp_100*tmp_111 + tmp_100*tmp_112) + tmp_78*(tmp_101*tmp_110 + tmp_101*tmp_111 + tmp_101*tmp_112) + tmp_80*(tmp_102*tmp_110 + tmp_102*tmp_111 + tmp_102*tmp_112) + tmp_82*(tmp_103*tmp_110 + tmp_103*tmp_111 + tmp_103*tmp_112) + tmp_84*(tmp_104*tmp_110 + tmp_104*tmp_111 + tmp_104*tmp_112) + tmp_86*(tmp_105*tmp_110 + tmp_105*tmp_111 + tmp_105*tmp_112);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_6*tmp_9;
      real_t tmp_14 = tmp_4*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_12 + tmp_0*tmp_7 - tmp_1*tmp_13 + tmp_1*tmp_2*tmp_8 + tmp_11*tmp_3 - tmp_14*tmp_3);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_7;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = tmp_33 - 1.0/4.0;
      real_t tmp_35 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_36 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_37 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_38 = tmp_20*tmp_35 + tmp_26*tmp_36 + tmp_32*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/4.0;
      real_t tmp_40 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_41 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_42 = tmp_11 - tmp_14;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_26*tmp_41 + tmp_32*tmp_42;
      real_t tmp_44 = tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_34 + tmp_1*tmp_39 + tmp_3*tmp_44;
      real_t tmp_46 = 0.5*tmp_15;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = tmp_37*tmp_46;
      real_t tmp_49 = tmp_27*tmp_46;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = 0.5*p_affine_13_2;
      real_t tmp_52 = tmp_50*tmp_51;
      real_t tmp_53 = tmp_2*tmp_44 + tmp_34*tmp_9 + tmp_39*tmp_4;
      real_t tmp_54 = tmp_41*tmp_46;
      real_t tmp_55 = tmp_36*tmp_46;
      real_t tmp_56 = tmp_21*tmp_46;
      real_t tmp_57 = -tmp_54 - tmp_55 - tmp_56;
      real_t tmp_58 = tmp_51*tmp_57;
      real_t tmp_59 = tmp_10*tmp_39 + tmp_34*tmp_8 + tmp_44*tmp_6;
      real_t tmp_60 = 1.0*tmp_15;
      real_t tmp_61 = tmp_40*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_5*tmp_60;
      real_t tmp_64 = 0.5*p_affine_13_0*tmp_50 + 0.5*p_affine_13_1*tmp_57 + 0.5*p_affine_13_2*(-tmp_61 - tmp_62 - tmp_63);
      real_t tmp_65 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_66 = std::pow(tmp_65, -0.25);
      real_t tmp_67 = -tmp_33 - tmp_38 - tmp_43 + 1;
      real_t tmp_68 = tmp_46*tmp_5;
      real_t tmp_69 = tmp_35*tmp_46;
      real_t tmp_70 = tmp_40*tmp_46;
      real_t tmp_71 = 0.5*p_affine_13_0*(tmp_0*tmp_68 + tmp_1*tmp_69 + tmp_10*tmp_48 + tmp_3*tmp_70 + tmp_47*tmp_6 + tmp_49*tmp_8) + 0.5*p_affine_13_1*(tmp_10*tmp_55 + tmp_2*tmp_70 + tmp_4*tmp_69 + tmp_54*tmp_6 + tmp_56*tmp_8 + tmp_68*tmp_9) + 0.5*p_affine_13_2*(tmp_10*tmp_62 + tmp_6*tmp_61 + tmp_63*tmp_8);
      real_t tmp_72 = 2.0*std::pow(tmp_65, 1.0/2.0);
      real_t tmp_73 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_72;
      real_t tmp_74 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_75 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_76 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_77 = tmp_21*tmp_75 + tmp_27*tmp_76 + tmp_5*tmp_74;
      real_t tmp_78 = tmp_77 - 1.0/4.0;
      real_t tmp_79 = tmp_35*tmp_74 + tmp_36*tmp_75 + tmp_37*tmp_76;
      real_t tmp_80 = tmp_79 - 1.0/4.0;
      real_t tmp_81 = tmp_40*tmp_74 + tmp_41*tmp_75 + tmp_42*tmp_76;
      real_t tmp_82 = tmp_81 - 1.0/4.0;
      real_t tmp_83 = tmp_0*tmp_78 + tmp_1*tmp_80 + tmp_3*tmp_82;
      real_t tmp_84 = tmp_2*tmp_82 + tmp_4*tmp_80 + tmp_78*tmp_9;
      real_t tmp_85 = tmp_10*tmp_80 + tmp_6*tmp_82 + tmp_78*tmp_8;
      real_t tmp_86 = -tmp_77 - tmp_79 - tmp_81 + 1;
      real_t tmp_87 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_72;
      real_t tmp_88 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_89 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_90 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_91 = tmp_21*tmp_89 + tmp_27*tmp_90 + tmp_5*tmp_88;
      real_t tmp_92 = tmp_91 - 1.0/4.0;
      real_t tmp_93 = tmp_35*tmp_88 + tmp_36*tmp_89 + tmp_37*tmp_90;
      real_t tmp_94 = tmp_93 - 1.0/4.0;
      real_t tmp_95 = tmp_40*tmp_88 + tmp_41*tmp_89 + tmp_42*tmp_90;
      real_t tmp_96 = tmp_95 - 1.0/4.0;
      real_t tmp_97 = tmp_0*tmp_92 + tmp_1*tmp_94 + tmp_3*tmp_96;
      real_t tmp_98 = tmp_2*tmp_96 + tmp_4*tmp_94 + tmp_9*tmp_92;
      real_t tmp_99 = tmp_10*tmp_94 + tmp_6*tmp_96 + tmp_8*tmp_92;
      real_t tmp_100 = -tmp_91 - tmp_93 - tmp_95 + 1;
      real_t tmp_101 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_72;
      real_t tmp_102 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_103 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_104 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_105 = tmp_102*tmp_5 + tmp_103*tmp_21 + tmp_104*tmp_27;
      real_t tmp_106 = tmp_105 - 1.0/4.0;
      real_t tmp_107 = tmp_102*tmp_35 + tmp_103*tmp_36 + tmp_104*tmp_37;
      real_t tmp_108 = tmp_107 - 1.0/4.0;
      real_t tmp_109 = tmp_102*tmp_40 + tmp_103*tmp_41 + tmp_104*tmp_42;
      real_t tmp_110 = tmp_109 - 1.0/4.0;
      real_t tmp_111 = tmp_0*tmp_106 + tmp_1*tmp_108 + tmp_110*tmp_3;
      real_t tmp_112 = tmp_106*tmp_9 + tmp_108*tmp_4 + tmp_110*tmp_2;
      real_t tmp_113 = tmp_10*tmp_108 + tmp_106*tmp_8 + tmp_110*tmp_6;
      real_t tmp_114 = -tmp_105 - tmp_107 - tmp_109 + 1;
      real_t tmp_115 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_72;
      real_t tmp_116 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_117 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_118 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_119 = tmp_116*tmp_5 + tmp_117*tmp_21 + tmp_118*tmp_27;
      real_t tmp_120 = tmp_119 - 1.0/4.0;
      real_t tmp_121 = tmp_116*tmp_35 + tmp_117*tmp_36 + tmp_118*tmp_37;
      real_t tmp_122 = tmp_121 - 1.0/4.0;
      real_t tmp_123 = tmp_116*tmp_40 + tmp_117*tmp_41 + tmp_118*tmp_42;
      real_t tmp_124 = tmp_123 - 1.0/4.0;
      real_t tmp_125 = tmp_0*tmp_120 + tmp_1*tmp_122 + tmp_124*tmp_3;
      real_t tmp_126 = tmp_120*tmp_9 + tmp_122*tmp_4 + tmp_124*tmp_2;
      real_t tmp_127 = tmp_10*tmp_122 + tmp_120*tmp_8 + tmp_124*tmp_6;
      real_t tmp_128 = -tmp_119 - tmp_121 - tmp_123 + 1;
      real_t tmp_129 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_72;
      real_t tmp_130 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_131 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_132 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_133 = tmp_130*tmp_5 + tmp_131*tmp_21 + tmp_132*tmp_27;
      real_t tmp_134 = tmp_133 - 1.0/4.0;
      real_t tmp_135 = tmp_130*tmp_35 + tmp_131*tmp_36 + tmp_132*tmp_37;
      real_t tmp_136 = tmp_135 - 1.0/4.0;
      real_t tmp_137 = tmp_130*tmp_40 + tmp_131*tmp_41 + tmp_132*tmp_42;
      real_t tmp_138 = tmp_137 - 1.0/4.0;
      real_t tmp_139 = tmp_0*tmp_134 + tmp_1*tmp_136 + tmp_138*tmp_3;
      real_t tmp_140 = tmp_134*tmp_9 + tmp_136*tmp_4 + tmp_138*tmp_2;
      real_t tmp_141 = tmp_10*tmp_136 + tmp_134*tmp_8 + tmp_138*tmp_6;
      real_t tmp_142 = -tmp_133 - tmp_135 - tmp_137 + 1;
      real_t tmp_143 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_72;
      real_t tmp_144 = 0.25*p_affine_13_2*tmp_15;
      real_t tmp_145 = tmp_144*tmp_21;
      real_t tmp_146 = tmp_144*tmp_27;
      real_t tmp_147 = 0.5*p_affine_13_0*tmp_49 + 0.5*p_affine_13_1*tmp_56 + 0.5*p_affine_13_2*tmp_63;
      real_t tmp_148 = tmp_144*tmp_36;
      real_t tmp_149 = tmp_144*tmp_37;
      real_t tmp_150 = 0.5*p_affine_13_0*tmp_48 + 0.5*p_affine_13_1*tmp_55 + 0.5*p_affine_13_2*tmp_62;
      real_t tmp_151 = tmp_144*tmp_41;
      real_t tmp_152 = tmp_144*tmp_42;
      real_t tmp_153 = 0.5*p_affine_13_0*tmp_47 + 0.5*p_affine_13_1*tmp_54 + 0.5*p_affine_13_2*tmp_61;
      real_t a_0_0 = tmp_101*(1.0*tmp_100*tmp_66*tmp_99 - tmp_100*tmp_71 - tmp_52*tmp_97 - tmp_58*tmp_98 - tmp_64*tmp_99) + tmp_115*(-tmp_111*tmp_52 - tmp_112*tmp_58 + 1.0*tmp_113*tmp_114*tmp_66 - tmp_113*tmp_64 - tmp_114*tmp_71) + tmp_129*(-tmp_125*tmp_52 - tmp_126*tmp_58 + 1.0*tmp_127*tmp_128*tmp_66 - tmp_127*tmp_64 - tmp_128*tmp_71) + tmp_143*(-tmp_139*tmp_52 - tmp_140*tmp_58 + 1.0*tmp_141*tmp_142*tmp_66 - tmp_141*tmp_64 - tmp_142*tmp_71) + tmp_73*(-tmp_45*tmp_52 - tmp_53*tmp_58 - tmp_59*tmp_64 + 1.0*tmp_59*tmp_66*tmp_67 - tmp_67*tmp_71) + tmp_87*(-tmp_52*tmp_83 - tmp_58*tmp_84 - tmp_64*tmp_85 + 1.0*tmp_66*tmp_85*tmp_86 - tmp_71*tmp_86);
      real_t a_0_1 = tmp_101*(-tmp_145*tmp_98 - tmp_146*tmp_97 - tmp_147*tmp_99 + 1.0*tmp_66*tmp_91*tmp_99 - tmp_71*tmp_91) + tmp_115*(1.0*tmp_105*tmp_113*tmp_66 - tmp_105*tmp_71 - tmp_111*tmp_146 - tmp_112*tmp_145 - tmp_113*tmp_147) + tmp_129*(1.0*tmp_119*tmp_127*tmp_66 - tmp_119*tmp_71 - tmp_125*tmp_146 - tmp_126*tmp_145 - tmp_127*tmp_147) + tmp_143*(1.0*tmp_133*tmp_141*tmp_66 - tmp_133*tmp_71 - tmp_139*tmp_146 - tmp_140*tmp_145 - tmp_141*tmp_147) + tmp_73*(-tmp_145*tmp_53 - tmp_146*tmp_45 - tmp_147*tmp_59 + 1.0*tmp_33*tmp_59*tmp_66 - tmp_33*tmp_71) + tmp_87*(-tmp_145*tmp_84 - tmp_146*tmp_83 - tmp_147*tmp_85 + 1.0*tmp_66*tmp_77*tmp_85 - tmp_71*tmp_77);
      real_t a_0_2 = tmp_101*(-tmp_148*tmp_98 - tmp_149*tmp_97 - tmp_150*tmp_99 + 1.0*tmp_66*tmp_93*tmp_99 - tmp_71*tmp_93) + tmp_115*(1.0*tmp_107*tmp_113*tmp_66 - tmp_107*tmp_71 - tmp_111*tmp_149 - tmp_112*tmp_148 - tmp_113*tmp_150) + tmp_129*(1.0*tmp_121*tmp_127*tmp_66 - tmp_121*tmp_71 - tmp_125*tmp_149 - tmp_126*tmp_148 - tmp_127*tmp_150) + tmp_143*(1.0*tmp_135*tmp_141*tmp_66 - tmp_135*tmp_71 - tmp_139*tmp_149 - tmp_140*tmp_148 - tmp_141*tmp_150) + tmp_73*(-tmp_148*tmp_53 - tmp_149*tmp_45 - tmp_150*tmp_59 + 1.0*tmp_38*tmp_59*tmp_66 - tmp_38*tmp_71) + tmp_87*(-tmp_148*tmp_84 - tmp_149*tmp_83 - tmp_150*tmp_85 + 1.0*tmp_66*tmp_79*tmp_85 - tmp_71*tmp_79);
      real_t a_0_3 = tmp_101*(-tmp_151*tmp_98 - tmp_152*tmp_97 - tmp_153*tmp_99 + 1.0*tmp_66*tmp_95*tmp_99 - tmp_71*tmp_95) + tmp_115*(1.0*tmp_109*tmp_113*tmp_66 - tmp_109*tmp_71 - tmp_111*tmp_152 - tmp_112*tmp_151 - tmp_113*tmp_153) + tmp_129*(1.0*tmp_123*tmp_127*tmp_66 - tmp_123*tmp_71 - tmp_125*tmp_152 - tmp_126*tmp_151 - tmp_127*tmp_153) + tmp_143*(1.0*tmp_137*tmp_141*tmp_66 - tmp_137*tmp_71 - tmp_139*tmp_152 - tmp_140*tmp_151 - tmp_141*tmp_153) + tmp_73*(-tmp_151*tmp_53 - tmp_152*tmp_45 - tmp_153*tmp_59 + 1.0*tmp_43*tmp_59*tmp_66 - tmp_43*tmp_71) + tmp_87*(-tmp_151*tmp_84 - tmp_152*tmp_83 - tmp_153*tmp_85 + 1.0*tmp_66*tmp_81*tmp_85 - tmp_71*tmp_81);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_11 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_12 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_13 = tmp_11*tmp_2;
      real_t tmp_14 = tmp_12*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_13 + tmp_0*tmp_9 - tmp_1*tmp_14 + tmp_10*tmp_3 - tmp_10*tmp_6 + tmp_11*tmp_12*tmp_4);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = 0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18;
      real_t tmp_20 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_21 = tmp_15*(tmp_19 + tmp_20);
      real_t tmp_22 = -tmp_1*tmp_8 + tmp_11*tmp_4;
      real_t tmp_23 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_24 = -tmp_23;
      real_t tmp_25 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_26 = 0.091576213509770743*tmp_24 + 0.81684757298045851*tmp_25;
      real_t tmp_27 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_28 = tmp_15*(tmp_26 + tmp_27);
      real_t tmp_29 = -tmp_13 + tmp_9;
      real_t tmp_30 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_31 = -tmp_30;
      real_t tmp_32 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_33 = 0.091576213509770743*tmp_31 + 0.81684757298045851*tmp_32;
      real_t tmp_34 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_35 = tmp_15*(tmp_33 + tmp_34);
      real_t tmp_36 = tmp_21*tmp_7 + tmp_22*tmp_28 + tmp_29*tmp_35 - 1.0/4.0;
      real_t tmp_37 = -tmp_0*tmp_2 + tmp_12*tmp_4;
      real_t tmp_38 = tmp_0*tmp_8 - tmp_10*tmp_4;
      real_t tmp_39 = tmp_10*tmp_2 - tmp_14;
      real_t tmp_40 = tmp_21*tmp_37 + tmp_28*tmp_38 + tmp_35*tmp_39 - 1.0/4.0;
      real_t tmp_41 = tmp_0*tmp_5 - tmp_1*tmp_12;
      real_t tmp_42 = -tmp_0*tmp_11 + tmp_1*tmp_10;
      real_t tmp_43 = -tmp_10*tmp_5 + tmp_11*tmp_12;
      real_t tmp_44 = tmp_21*tmp_41 + tmp_28*tmp_42 + tmp_35*tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_36 + tmp_1*tmp_40 + tmp_4*tmp_44;
      real_t tmp_46 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_47 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_48 = tmp_46*tmp_47;
      real_t tmp_49 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_50 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_51 = tmp_48 - tmp_49*tmp_50;
      real_t tmp_52 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_53 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_54 = tmp_49*tmp_53;
      real_t tmp_55 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_56 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_57 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_58 = tmp_53*tmp_55;
      real_t tmp_59 = tmp_46*tmp_56;
      real_t tmp_60 = tmp_50*tmp_57;
      real_t tmp_61 = 1.0 / (-tmp_47*tmp_58 + tmp_48*tmp_57 - tmp_49*tmp_60 + tmp_50*tmp_55*tmp_56 + tmp_52*tmp_54 - tmp_52*tmp_59);
      real_t tmp_62 = 0.5*tmp_61;
      real_t tmp_63 = tmp_51*tmp_62;
      real_t tmp_64 = -tmp_46*tmp_52 + tmp_50*tmp_55;
      real_t tmp_65 = tmp_62*tmp_64;
      real_t tmp_66 = -tmp_47*tmp_55 + tmp_49*tmp_52;
      real_t tmp_67 = tmp_62*tmp_66;
      real_t tmp_68 = -tmp_63 - tmp_65 - tmp_67;
      real_t tmp_69 = 0.5*p_affine_13_2;
      real_t tmp_70 = tmp_68*tmp_69;
      real_t tmp_71 = tmp_12*tmp_36 + tmp_2*tmp_44 + tmp_40*tmp_5;
      real_t tmp_72 = -tmp_47*tmp_53 + tmp_50*tmp_56;
      real_t tmp_73 = tmp_62*tmp_72;
      real_t tmp_74 = tmp_52*tmp_53 - tmp_60;
      real_t tmp_75 = tmp_62*tmp_74;
      real_t tmp_76 = tmp_47*tmp_57 - tmp_52*tmp_56;
      real_t tmp_77 = tmp_62*tmp_76;
      real_t tmp_78 = -tmp_73 - tmp_75 - tmp_77;
      real_t tmp_79 = tmp_69*tmp_78;
      real_t tmp_80 = tmp_10*tmp_36 + tmp_11*tmp_40 + tmp_44*tmp_8;
      real_t tmp_81 = tmp_54 - tmp_59;
      real_t tmp_82 = 1.0*tmp_61;
      real_t tmp_83 = tmp_81*tmp_82;
      real_t tmp_84 = tmp_46*tmp_57 - tmp_58;
      real_t tmp_85 = tmp_82*tmp_84;
      real_t tmp_86 = -tmp_49*tmp_57 + tmp_55*tmp_56;
      real_t tmp_87 = tmp_82*tmp_86;
      real_t tmp_88 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_78 + 0.5*p_affine_13_2*(-tmp_83 - tmp_85 - tmp_87);
      real_t tmp_89 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_90 = tmp_61*(tmp_19 + tmp_89);
      real_t tmp_91 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_92 = tmp_61*(tmp_26 + tmp_91);
      real_t tmp_93 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_94 = tmp_61*(tmp_33 + tmp_93);
      real_t tmp_95 = tmp_66*tmp_94 + tmp_76*tmp_92 + tmp_86*tmp_90;
      real_t tmp_96 = tmp_64*tmp_94 + tmp_74*tmp_92 + tmp_84*tmp_90;
      real_t tmp_97 = tmp_51*tmp_94 + tmp_72*tmp_92 + tmp_81*tmp_90;
      real_t tmp_98 = -tmp_95 - tmp_96 - tmp_97 + 1;
      real_t tmp_99 = (std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)*std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)) + (std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)*std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)) + (std::abs(tmp_23*tmp_32 - tmp_25*tmp_30)*std::abs(tmp_23*tmp_32 - tmp_25*tmp_30));
      real_t tmp_100 = 1.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_80;
      real_t tmp_102 = 1.0*tmp_15;
      real_t tmp_103 = 0.5*tmp_15;
      real_t tmp_104 = tmp_103*tmp_7;
      real_t tmp_105 = tmp_103*tmp_37;
      real_t tmp_106 = tmp_103*tmp_41;
      real_t tmp_107 = tmp_10*tmp_103;
      real_t tmp_108 = tmp_103*tmp_11;
      real_t tmp_109 = tmp_103*tmp_8;
      real_t tmp_110 = p_affine_13_0*(tmp_0*tmp_104 + tmp_1*tmp_105 + tmp_106*tmp_4 + tmp_107*tmp_29 + tmp_108*tmp_39 + tmp_109*tmp_43) + p_affine_13_1*(tmp_104*tmp_12 + tmp_105*tmp_5 + tmp_106*tmp_2 + tmp_107*tmp_22 + tmp_108*tmp_38 + tmp_109*tmp_42) + p_affine_13_2*(tmp_10*tmp_102*tmp_7 + tmp_102*tmp_11*tmp_37 + tmp_102*tmp_41*tmp_8);
      real_t tmp_111 = 2.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_112 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_111;
      real_t tmp_113 = 0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18;
      real_t tmp_114 = tmp_15*(tmp_113 + tmp_20);
      real_t tmp_115 = 0.44594849091596489*tmp_24 + 0.10810301816807022*tmp_25;
      real_t tmp_116 = tmp_15*(tmp_115 + tmp_27);
      real_t tmp_117 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32;
      real_t tmp_118 = tmp_15*(tmp_117 + tmp_34);
      real_t tmp_119 = tmp_114*tmp_7 + tmp_116*tmp_22 + tmp_118*tmp_29 - 1.0/4.0;
      real_t tmp_120 = tmp_114*tmp_37 + tmp_116*tmp_38 + tmp_118*tmp_39 - 1.0/4.0;
      real_t tmp_121 = tmp_114*tmp_41 + tmp_116*tmp_42 + tmp_118*tmp_43 - 1.0/4.0;
      real_t tmp_122 = tmp_0*tmp_119 + tmp_1*tmp_120 + tmp_121*tmp_4;
      real_t tmp_123 = tmp_119*tmp_12 + tmp_120*tmp_5 + tmp_121*tmp_2;
      real_t tmp_124 = tmp_10*tmp_119 + tmp_11*tmp_120 + tmp_121*tmp_8;
      real_t tmp_125 = tmp_61*(tmp_113 + tmp_89);
      real_t tmp_126 = tmp_61*(tmp_115 + tmp_91);
      real_t tmp_127 = tmp_61*(tmp_117 + tmp_93);
      real_t tmp_128 = tmp_125*tmp_86 + tmp_126*tmp_76 + tmp_127*tmp_66;
      real_t tmp_129 = tmp_125*tmp_84 + tmp_126*tmp_74 + tmp_127*tmp_64;
      real_t tmp_130 = tmp_125*tmp_81 + tmp_126*tmp_72 + tmp_127*tmp_51;
      real_t tmp_131 = -tmp_128 - tmp_129 - tmp_130 + 1;
      real_t tmp_132 = tmp_100*tmp_124;
      real_t tmp_133 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_111;
      real_t tmp_134 = 0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_135 = tmp_15*(tmp_134 + tmp_20);
      real_t tmp_136 = 0.81684757298045851*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_137 = tmp_15*(tmp_136 + tmp_27);
      real_t tmp_138 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_139 = tmp_15*(tmp_138 + tmp_34);
      real_t tmp_140 = tmp_135*tmp_7 + tmp_137*tmp_22 + tmp_139*tmp_29 - 1.0/4.0;
      real_t tmp_141 = tmp_135*tmp_37 + tmp_137*tmp_38 + tmp_139*tmp_39 - 1.0/4.0;
      real_t tmp_142 = tmp_135*tmp_41 + tmp_137*tmp_42 + tmp_139*tmp_43 - 1.0/4.0;
      real_t tmp_143 = tmp_0*tmp_140 + tmp_1*tmp_141 + tmp_142*tmp_4;
      real_t tmp_144 = tmp_12*tmp_140 + tmp_141*tmp_5 + tmp_142*tmp_2;
      real_t tmp_145 = tmp_10*tmp_140 + tmp_11*tmp_141 + tmp_142*tmp_8;
      real_t tmp_146 = tmp_61*(tmp_134 + tmp_89);
      real_t tmp_147 = tmp_61*(tmp_136 + tmp_91);
      real_t tmp_148 = tmp_61*(tmp_138 + tmp_93);
      real_t tmp_149 = tmp_146*tmp_86 + tmp_147*tmp_76 + tmp_148*tmp_66;
      real_t tmp_150 = tmp_146*tmp_84 + tmp_147*tmp_74 + tmp_148*tmp_64;
      real_t tmp_151 = tmp_146*tmp_81 + tmp_147*tmp_72 + tmp_148*tmp_51;
      real_t tmp_152 = -tmp_149 - tmp_150 - tmp_151 + 1;
      real_t tmp_153 = tmp_100*tmp_145;
      real_t tmp_154 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_111;
      real_t tmp_155 = 0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_156 = tmp_15*(tmp_155 + tmp_20);
      real_t tmp_157 = 0.10810301816807022*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_158 = tmp_15*(tmp_157 + tmp_27);
      real_t tmp_159 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_160 = tmp_15*(tmp_159 + tmp_34);
      real_t tmp_161 = tmp_156*tmp_7 + tmp_158*tmp_22 + tmp_160*tmp_29 - 1.0/4.0;
      real_t tmp_162 = tmp_156*tmp_37 + tmp_158*tmp_38 + tmp_160*tmp_39 - 1.0/4.0;
      real_t tmp_163 = tmp_156*tmp_41 + tmp_158*tmp_42 + tmp_160*tmp_43 - 1.0/4.0;
      real_t tmp_164 = tmp_0*tmp_161 + tmp_1*tmp_162 + tmp_163*tmp_4;
      real_t tmp_165 = tmp_12*tmp_161 + tmp_162*tmp_5 + tmp_163*tmp_2;
      real_t tmp_166 = tmp_10*tmp_161 + tmp_11*tmp_162 + tmp_163*tmp_8;
      real_t tmp_167 = tmp_61*(tmp_155 + tmp_89);
      real_t tmp_168 = tmp_61*(tmp_157 + tmp_91);
      real_t tmp_169 = tmp_61*(tmp_159 + tmp_93);
      real_t tmp_170 = tmp_167*tmp_86 + tmp_168*tmp_76 + tmp_169*tmp_66;
      real_t tmp_171 = tmp_167*tmp_84 + tmp_168*tmp_74 + tmp_169*tmp_64;
      real_t tmp_172 = tmp_167*tmp_81 + tmp_168*tmp_72 + tmp_169*tmp_51;
      real_t tmp_173 = -tmp_170 - tmp_171 - tmp_172 + 1;
      real_t tmp_174 = tmp_100*tmp_166;
      real_t tmp_175 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_111;
      real_t tmp_176 = 0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_177 = tmp_15*(tmp_176 + tmp_20);
      real_t tmp_178 = 0.091576213509770743*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_179 = tmp_15*(tmp_178 + tmp_27);
      real_t tmp_180 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_181 = tmp_15*(tmp_180 + tmp_34);
      real_t tmp_182 = tmp_177*tmp_7 + tmp_179*tmp_22 + tmp_181*tmp_29 - 1.0/4.0;
      real_t tmp_183 = tmp_177*tmp_37 + tmp_179*tmp_38 + tmp_181*tmp_39 - 1.0/4.0;
      real_t tmp_184 = tmp_177*tmp_41 + tmp_179*tmp_42 + tmp_181*tmp_43 - 1.0/4.0;
      real_t tmp_185 = tmp_0*tmp_182 + tmp_1*tmp_183 + tmp_184*tmp_4;
      real_t tmp_186 = tmp_12*tmp_182 + tmp_183*tmp_5 + tmp_184*tmp_2;
      real_t tmp_187 = tmp_10*tmp_182 + tmp_11*tmp_183 + tmp_184*tmp_8;
      real_t tmp_188 = tmp_61*(tmp_176 + tmp_89);
      real_t tmp_189 = tmp_61*(tmp_178 + tmp_91);
      real_t tmp_190 = tmp_61*(tmp_180 + tmp_93);
      real_t tmp_191 = tmp_188*tmp_86 + tmp_189*tmp_76 + tmp_190*tmp_66;
      real_t tmp_192 = tmp_188*tmp_84 + tmp_189*tmp_74 + tmp_190*tmp_64;
      real_t tmp_193 = tmp_188*tmp_81 + tmp_189*tmp_72 + tmp_190*tmp_51;
      real_t tmp_194 = -tmp_191 - tmp_192 - tmp_193 + 1;
      real_t tmp_195 = tmp_100*tmp_187;
      real_t tmp_196 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_111;
      real_t tmp_197 = 0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_198 = tmp_15*(tmp_197 + tmp_20);
      real_t tmp_199 = 0.44594849091596489*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_200 = tmp_15*(tmp_199 + tmp_27);
      real_t tmp_201 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_202 = tmp_15*(tmp_201 + tmp_34);
      real_t tmp_203 = tmp_198*tmp_7 + tmp_200*tmp_22 + tmp_202*tmp_29 - 1.0/4.0;
      real_t tmp_204 = tmp_198*tmp_37 + tmp_200*tmp_38 + tmp_202*tmp_39 - 1.0/4.0;
      real_t tmp_205 = tmp_198*tmp_41 + tmp_200*tmp_42 + tmp_202*tmp_43 - 1.0/4.0;
      real_t tmp_206 = tmp_0*tmp_203 + tmp_1*tmp_204 + tmp_205*tmp_4;
      real_t tmp_207 = tmp_12*tmp_203 + tmp_2*tmp_205 + tmp_204*tmp_5;
      real_t tmp_208 = tmp_10*tmp_203 + tmp_11*tmp_204 + tmp_205*tmp_8;
      real_t tmp_209 = tmp_61*(tmp_197 + tmp_89);
      real_t tmp_210 = tmp_61*(tmp_199 + tmp_91);
      real_t tmp_211 = tmp_61*(tmp_201 + tmp_93);
      real_t tmp_212 = tmp_209*tmp_86 + tmp_210*tmp_76 + tmp_211*tmp_66;
      real_t tmp_213 = tmp_209*tmp_84 + tmp_210*tmp_74 + tmp_211*tmp_64;
      real_t tmp_214 = tmp_209*tmp_81 + tmp_210*tmp_72 + tmp_211*tmp_51;
      real_t tmp_215 = -tmp_212 - tmp_213 - tmp_214 + 1;
      real_t tmp_216 = tmp_100*tmp_208;
      real_t tmp_217 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_111;
      real_t tmp_218 = 0.25*p_affine_13_2*tmp_61;
      real_t tmp_219 = tmp_218*tmp_76;
      real_t tmp_220 = tmp_218*tmp_66;
      real_t tmp_221 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*tmp_77 + 0.5*p_affine_13_2*tmp_87;
      real_t tmp_222 = tmp_218*tmp_74;
      real_t tmp_223 = tmp_218*tmp_64;
      real_t tmp_224 = 0.5*p_affine_13_0*tmp_65 + 0.5*p_affine_13_1*tmp_75 + 0.5*p_affine_13_2*tmp_85;
      real_t tmp_225 = tmp_218*tmp_72;
      real_t tmp_226 = tmp_218*tmp_51;
      real_t tmp_227 = 0.5*p_affine_13_0*tmp_63 + 0.5*p_affine_13_1*tmp_73 + 0.5*p_affine_13_2*tmp_83;
      real_t a_0_0 = tmp_112*(-tmp_101*tmp_98 + 0.5*tmp_110*tmp_98 - tmp_45*tmp_70 - tmp_71*tmp_79 - tmp_80*tmp_88) + tmp_133*(0.5*tmp_110*tmp_131 - tmp_122*tmp_70 - tmp_123*tmp_79 - tmp_124*tmp_88 - tmp_131*tmp_132) + tmp_154*(0.5*tmp_110*tmp_152 - tmp_143*tmp_70 - tmp_144*tmp_79 - tmp_145*tmp_88 - tmp_152*tmp_153) + tmp_175*(0.5*tmp_110*tmp_173 - tmp_164*tmp_70 - tmp_165*tmp_79 - tmp_166*tmp_88 - tmp_173*tmp_174) + tmp_196*(0.5*tmp_110*tmp_194 - tmp_185*tmp_70 - tmp_186*tmp_79 - tmp_187*tmp_88 - tmp_194*tmp_195) + tmp_217*(0.5*tmp_110*tmp_215 - tmp_206*tmp_70 - tmp_207*tmp_79 - tmp_208*tmp_88 - tmp_215*tmp_216);
      real_t a_0_1 = tmp_112*(-tmp_101*tmp_95 + 0.5*tmp_110*tmp_95 - tmp_219*tmp_71 - tmp_220*tmp_45 - tmp_221*tmp_80) + tmp_133*(0.5*tmp_110*tmp_128 - tmp_122*tmp_220 - tmp_123*tmp_219 - tmp_124*tmp_221 - tmp_128*tmp_132) + tmp_154*(0.5*tmp_110*tmp_149 - tmp_143*tmp_220 - tmp_144*tmp_219 - tmp_145*tmp_221 - tmp_149*tmp_153) + tmp_175*(0.5*tmp_110*tmp_170 - tmp_164*tmp_220 - tmp_165*tmp_219 - tmp_166*tmp_221 - tmp_170*tmp_174) + tmp_196*(0.5*tmp_110*tmp_191 - tmp_185*tmp_220 - tmp_186*tmp_219 - tmp_187*tmp_221 - tmp_191*tmp_195) + tmp_217*(0.5*tmp_110*tmp_212 - tmp_206*tmp_220 - tmp_207*tmp_219 - tmp_208*tmp_221 - tmp_212*tmp_216);
      real_t a_0_2 = tmp_112*(-tmp_101*tmp_96 + 0.5*tmp_110*tmp_96 - tmp_222*tmp_71 - tmp_223*tmp_45 - tmp_224*tmp_80) + tmp_133*(0.5*tmp_110*tmp_129 - tmp_122*tmp_223 - tmp_123*tmp_222 - tmp_124*tmp_224 - tmp_129*tmp_132) + tmp_154*(0.5*tmp_110*tmp_150 - tmp_143*tmp_223 - tmp_144*tmp_222 - tmp_145*tmp_224 - tmp_150*tmp_153) + tmp_175*(0.5*tmp_110*tmp_171 - tmp_164*tmp_223 - tmp_165*tmp_222 - tmp_166*tmp_224 - tmp_171*tmp_174) + tmp_196*(0.5*tmp_110*tmp_192 - tmp_185*tmp_223 - tmp_186*tmp_222 - tmp_187*tmp_224 - tmp_192*tmp_195) + tmp_217*(0.5*tmp_110*tmp_213 - tmp_206*tmp_223 - tmp_207*tmp_222 - tmp_208*tmp_224 - tmp_213*tmp_216);
      real_t a_0_3 = tmp_112*(-tmp_101*tmp_97 + 0.5*tmp_110*tmp_97 - tmp_225*tmp_71 - tmp_226*tmp_45 - tmp_227*tmp_80) + tmp_133*(0.5*tmp_110*tmp_130 - tmp_122*tmp_226 - tmp_123*tmp_225 - tmp_124*tmp_227 - tmp_130*tmp_132) + tmp_154*(0.5*tmp_110*tmp_151 - tmp_143*tmp_226 - tmp_144*tmp_225 - tmp_145*tmp_227 - tmp_151*tmp_153) + tmp_175*(0.5*tmp_110*tmp_172 - tmp_164*tmp_226 - tmp_165*tmp_225 - tmp_166*tmp_227 - tmp_172*tmp_174) + tmp_196*(0.5*tmp_110*tmp_193 - tmp_185*tmp_226 - tmp_186*tmp_225 - tmp_187*tmp_227 - tmp_193*tmp_195) + tmp_217*(0.5*tmp_110*tmp_214 - tmp_206*tmp_226 - tmp_207*tmp_225 - tmp_208*tmp_227 - tmp_214*tmp_216);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_6*tmp_9;
      real_t tmp_14 = tmp_4*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_12 + tmp_0*tmp_7 - tmp_1*tmp_13 + tmp_1*tmp_2*tmp_8 + tmp_11*tmp_3 - tmp_14*tmp_3);
      real_t tmp_16 = -p_affine_8_2 + p_affine_9_2;
      real_t tmp_17 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_18 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_19 = tmp_15*(0.091576213509770743*tmp_16 + 0.81684757298045851*tmp_17 + tmp_18);
      real_t tmp_20 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_21 = -p_affine_8_1 + p_affine_9_1;
      real_t tmp_22 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_23 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_24 = tmp_15*(0.091576213509770743*tmp_21 + 0.81684757298045851*tmp_22 + tmp_23);
      real_t tmp_25 = -tmp_12 + tmp_7;
      real_t tmp_26 = -p_affine_8_0 + p_affine_9_0;
      real_t tmp_27 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_28 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_29 = tmp_15*(0.091576213509770743*tmp_26 + 0.81684757298045851*tmp_27 + tmp_28);
      real_t tmp_30 = tmp_19*tmp_5 + tmp_20*tmp_24 + tmp_25*tmp_29 - 1.0/4.0;
      real_t tmp_31 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_32 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_33 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_34 = tmp_19*tmp_31 + tmp_24*tmp_32 + tmp_29*tmp_33 - 1.0/4.0;
      real_t tmp_35 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_36 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_37 = tmp_11 - tmp_14;
      real_t tmp_38 = tmp_19*tmp_35 + tmp_24*tmp_36 + tmp_29*tmp_37 - 1.0/4.0;
      real_t tmp_39 = tmp_0*tmp_30 + tmp_1*tmp_34 + tmp_3*tmp_38;
      real_t tmp_40 = 0.5*tmp_15;
      real_t tmp_41 = tmp_37*tmp_40;
      real_t tmp_42 = tmp_33*tmp_40;
      real_t tmp_43 = tmp_25*tmp_40;
      real_t tmp_44 = -tmp_41 - tmp_42 - tmp_43;
      real_t tmp_45 = p_affine_13_2*tmp_44;
      real_t tmp_46 = tmp_2*tmp_38 + tmp_30*tmp_9 + tmp_34*tmp_4;
      real_t tmp_47 = tmp_36*tmp_40;
      real_t tmp_48 = tmp_32*tmp_40;
      real_t tmp_49 = tmp_20*tmp_40;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = p_affine_13_2*tmp_50;
      real_t tmp_52 = 1.0*tmp_15;
      real_t tmp_53 = tmp_35*tmp_52;
      real_t tmp_54 = tmp_31*tmp_52;
      real_t tmp_55 = tmp_5*tmp_52;
      real_t tmp_56 = p_affine_13_0*tmp_44 + p_affine_13_1*tmp_50 + p_affine_13_2*(-tmp_53 - tmp_54 - tmp_55);
      real_t tmp_57 = tmp_10*tmp_34 + tmp_30*tmp_8 + tmp_38*tmp_6;
      real_t tmp_58 = tmp_15*(0.44594849091596489*tmp_16 + 0.10810301816807022*tmp_17 + tmp_18);
      real_t tmp_59 = tmp_15*(0.44594849091596489*tmp_21 + 0.10810301816807022*tmp_22 + tmp_23);
      real_t tmp_60 = tmp_15*(0.44594849091596489*tmp_26 + 0.10810301816807022*tmp_27 + tmp_28);
      real_t tmp_61 = tmp_20*tmp_59 + tmp_25*tmp_60 + tmp_5*tmp_58 - 1.0/4.0;
      real_t tmp_62 = tmp_31*tmp_58 + tmp_32*tmp_59 + tmp_33*tmp_60 - 1.0/4.0;
      real_t tmp_63 = tmp_35*tmp_58 + tmp_36*tmp_59 + tmp_37*tmp_60 - 1.0/4.0;
      real_t tmp_64 = tmp_0*tmp_61 + tmp_1*tmp_62 + tmp_3*tmp_63;
      real_t tmp_65 = tmp_2*tmp_63 + tmp_4*tmp_62 + tmp_61*tmp_9;
      real_t tmp_66 = tmp_10*tmp_62 + tmp_6*tmp_63 + tmp_61*tmp_8;
      real_t tmp_67 = tmp_15*(0.81684757298045851*tmp_16 + 0.091576213509770743*tmp_17 + tmp_18);
      real_t tmp_68 = tmp_15*(0.81684757298045851*tmp_21 + 0.091576213509770743*tmp_22 + tmp_23);
      real_t tmp_69 = tmp_15*(0.81684757298045851*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28);
      real_t tmp_70 = tmp_20*tmp_68 + tmp_25*tmp_69 + tmp_5*tmp_67 - 1.0/4.0;
      real_t tmp_71 = tmp_31*tmp_67 + tmp_32*tmp_68 + tmp_33*tmp_69 - 1.0/4.0;
      real_t tmp_72 = tmp_35*tmp_67 + tmp_36*tmp_68 + tmp_37*tmp_69 - 1.0/4.0;
      real_t tmp_73 = tmp_0*tmp_70 + tmp_1*tmp_71 + tmp_3*tmp_72;
      real_t tmp_74 = tmp_2*tmp_72 + tmp_4*tmp_71 + tmp_70*tmp_9;
      real_t tmp_75 = tmp_10*tmp_71 + tmp_6*tmp_72 + tmp_70*tmp_8;
      real_t tmp_76 = tmp_15*(0.10810301816807022*tmp_16 + 0.44594849091596489*tmp_17 + tmp_18);
      real_t tmp_77 = tmp_15*(0.10810301816807022*tmp_21 + 0.44594849091596489*tmp_22 + tmp_23);
      real_t tmp_78 = tmp_15*(0.10810301816807022*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28);
      real_t tmp_79 = tmp_20*tmp_77 + tmp_25*tmp_78 + tmp_5*tmp_76 - 1.0/4.0;
      real_t tmp_80 = tmp_31*tmp_76 + tmp_32*tmp_77 + tmp_33*tmp_78 - 1.0/4.0;
      real_t tmp_81 = tmp_35*tmp_76 + tmp_36*tmp_77 + tmp_37*tmp_78 - 1.0/4.0;
      real_t tmp_82 = tmp_0*tmp_79 + tmp_1*tmp_80 + tmp_3*tmp_81;
      real_t tmp_83 = tmp_2*tmp_81 + tmp_4*tmp_80 + tmp_79*tmp_9;
      real_t tmp_84 = tmp_10*tmp_80 + tmp_6*tmp_81 + tmp_79*tmp_8;
      real_t tmp_85 = tmp_15*(0.091576213509770743*tmp_16 + 0.091576213509770743*tmp_17 + tmp_18);
      real_t tmp_86 = tmp_15*(0.091576213509770743*tmp_21 + 0.091576213509770743*tmp_22 + tmp_23);
      real_t tmp_87 = tmp_15*(0.091576213509770743*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28);
      real_t tmp_88 = tmp_20*tmp_86 + tmp_25*tmp_87 + tmp_5*tmp_85 - 1.0/4.0;
      real_t tmp_89 = tmp_31*tmp_85 + tmp_32*tmp_86 + tmp_33*tmp_87 - 1.0/4.0;
      real_t tmp_90 = tmp_35*tmp_85 + tmp_36*tmp_86 + tmp_37*tmp_87 - 1.0/4.0;
      real_t tmp_91 = tmp_0*tmp_88 + tmp_1*tmp_89 + tmp_3*tmp_90;
      real_t tmp_92 = tmp_2*tmp_90 + tmp_4*tmp_89 + tmp_88*tmp_9;
      real_t tmp_93 = tmp_10*tmp_89 + tmp_6*tmp_90 + tmp_8*tmp_88;
      real_t tmp_94 = tmp_15*(0.44594849091596489*tmp_16 + 0.44594849091596489*tmp_17 + tmp_18);
      real_t tmp_95 = tmp_15*(0.44594849091596489*tmp_21 + 0.44594849091596489*tmp_22 + tmp_23);
      real_t tmp_96 = tmp_15*(0.44594849091596489*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28);
      real_t tmp_97 = tmp_20*tmp_95 + tmp_25*tmp_96 + tmp_5*tmp_94 - 1.0/4.0;
      real_t tmp_98 = tmp_31*tmp_94 + tmp_32*tmp_95 + tmp_33*tmp_96 - 1.0/4.0;
      real_t tmp_99 = tmp_35*tmp_94 + tmp_36*tmp_95 + tmp_37*tmp_96 - 1.0/4.0;
      real_t tmp_100 = tmp_0*tmp_97 + tmp_1*tmp_98 + tmp_3*tmp_99;
      real_t tmp_101 = tmp_2*tmp_99 + tmp_4*tmp_98 + tmp_9*tmp_97;
      real_t tmp_102 = tmp_10*tmp_98 + tmp_6*tmp_99 + tmp_8*tmp_97;
      real_t tmp_103 = p_affine_13_2*tmp_49;
      real_t tmp_104 = p_affine_13_2*tmp_43;
      real_t tmp_105 = p_affine_13_0*tmp_43 + p_affine_13_1*tmp_49 + p_affine_13_2*tmp_55;
      real_t tmp_106 = p_affine_13_2*tmp_48;
      real_t tmp_107 = p_affine_13_2*tmp_42;
      real_t tmp_108 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_48 + p_affine_13_2*tmp_54;
      real_t tmp_109 = p_affine_13_2*tmp_47;
      real_t tmp_110 = p_affine_13_2*tmp_41;
      real_t tmp_111 = p_affine_13_0*tmp_41 + p_affine_13_1*tmp_47 + p_affine_13_2*tmp_53;
      real_t a_0_0 = -0.11169079483900572*tmp_100*tmp_45 - 0.11169079483900572*tmp_101*tmp_51 - 0.11169079483900572*tmp_102*tmp_56 - 0.054975871827660928*tmp_39*tmp_45 - 0.11169079483900572*tmp_45*tmp_64 - 0.054975871827660928*tmp_45*tmp_73 - 0.11169079483900572*tmp_45*tmp_82 - 0.054975871827660928*tmp_45*tmp_91 - 0.054975871827660928*tmp_46*tmp_51 - 0.11169079483900572*tmp_51*tmp_65 - 0.054975871827660928*tmp_51*tmp_74 - 0.11169079483900572*tmp_51*tmp_83 - 0.054975871827660928*tmp_51*tmp_92 - 0.054975871827660928*tmp_56*tmp_57 - 0.11169079483900572*tmp_56*tmp_66 - 0.054975871827660928*tmp_56*tmp_75 - 0.11169079483900572*tmp_56*tmp_84 - 0.054975871827660928*tmp_56*tmp_93;
      real_t a_0_1 = -0.11169079483900572*tmp_100*tmp_104 - 0.11169079483900572*tmp_101*tmp_103 - 0.11169079483900572*tmp_102*tmp_105 - 0.054975871827660928*tmp_103*tmp_46 - 0.11169079483900572*tmp_103*tmp_65 - 0.054975871827660928*tmp_103*tmp_74 - 0.11169079483900572*tmp_103*tmp_83 - 0.054975871827660928*tmp_103*tmp_92 - 0.054975871827660928*tmp_104*tmp_39 - 0.11169079483900572*tmp_104*tmp_64 - 0.054975871827660928*tmp_104*tmp_73 - 0.11169079483900572*tmp_104*tmp_82 - 0.054975871827660928*tmp_104*tmp_91 - 0.054975871827660928*tmp_105*tmp_57 - 0.11169079483900572*tmp_105*tmp_66 - 0.054975871827660928*tmp_105*tmp_75 - 0.11169079483900572*tmp_105*tmp_84 - 0.054975871827660928*tmp_105*tmp_93;
      real_t a_0_2 = -0.11169079483900572*tmp_100*tmp_107 - 0.11169079483900572*tmp_101*tmp_106 - 0.11169079483900572*tmp_102*tmp_108 - 0.054975871827660928*tmp_106*tmp_46 - 0.11169079483900572*tmp_106*tmp_65 - 0.054975871827660928*tmp_106*tmp_74 - 0.11169079483900572*tmp_106*tmp_83 - 0.054975871827660928*tmp_106*tmp_92 - 0.054975871827660928*tmp_107*tmp_39 - 0.11169079483900572*tmp_107*tmp_64 - 0.054975871827660928*tmp_107*tmp_73 - 0.11169079483900572*tmp_107*tmp_82 - 0.054975871827660928*tmp_107*tmp_91 - 0.054975871827660928*tmp_108*tmp_57 - 0.11169079483900572*tmp_108*tmp_66 - 0.054975871827660928*tmp_108*tmp_75 - 0.11169079483900572*tmp_108*tmp_84 - 0.054975871827660928*tmp_108*tmp_93;
      real_t a_0_3 = -0.11169079483900572*tmp_100*tmp_110 - 0.11169079483900572*tmp_101*tmp_109 - 0.11169079483900572*tmp_102*tmp_111 - 0.054975871827660928*tmp_109*tmp_46 - 0.11169079483900572*tmp_109*tmp_65 - 0.054975871827660928*tmp_109*tmp_74 - 0.11169079483900572*tmp_109*tmp_83 - 0.054975871827660928*tmp_109*tmp_92 - 0.054975871827660928*tmp_110*tmp_39 - 0.11169079483900572*tmp_110*tmp_64 - 0.054975871827660928*tmp_110*tmp_73 - 0.11169079483900572*tmp_110*tmp_82 - 0.054975871827660928*tmp_110*tmp_91 - 0.054975871827660928*tmp_111*tmp_57 - 0.11169079483900572*tmp_111*tmp_66 - 0.054975871827660928*tmp_111*tmp_75 - 0.11169079483900572*tmp_111*tmp_84 - 0.054975871827660928*tmp_111*tmp_93;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }

};




class EGEpsilonFormP1E_2 : public hyteg::dg::DGForm
{

 public:
    EGEpsilonFormP1E_2(std::function< real_t ( const Point3D & ) > mu)
: callback_Scalar_Variable_Coefficient_3D_mu (mu)
    {}

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;



void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
}

      
 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t a_0_0 = 0;
      real_t a_1_0 = 0;
      real_t a_2_0 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t a_0_0 = 0;
      real_t a_1_0 = 0;
      real_t a_2_0 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElementInner[0]( 0 );
      const auto p_affine_0_1 = coordsElementInner[0]( 1 );

      const auto p_affine_1_0 = coordsElementInner[1]( 0 );
      const auto p_affine_1_1 = coordsElementInner[1]( 1 );

      const auto p_affine_2_0 = coordsElementInner[2]( 0 );
      const auto p_affine_2_1 = coordsElementInner[2]( 1 );

      const auto p_affine_3_0 = coordsElementOuter[0]( 0 );
      const auto p_affine_3_1 = coordsElementOuter[0]( 1 );

      const auto p_affine_4_0 = coordsElementOuter[1]( 0 );
      const auto p_affine_4_1 = coordsElementOuter[1]( 1 );

      const auto p_affine_5_0 = coordsElementOuter[2]( 0 );
      const auto p_affine_5_1 = coordsElementOuter[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t a_0_0 = 0;
      real_t a_1_0 = 0;
      real_t a_2_0 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t a_0_0 = 0;
      real_t a_1_0 = 0;
      real_t a_2_0 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

  void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

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

      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id6 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id7 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id8 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id9 );
      Scalar_Variable_Coefficient_3D_mu( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id10 );
      Scalar_Variable_Coefficient_3D_mu( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id11 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id12 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id13 );
      real_t tmp_0 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_9 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_10 = tmp_5*tmp_9;
      real_t tmp_11 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_12 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_13 = tmp_2*tmp_9;
      real_t tmp_14 = tmp_1*tmp_12;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_3 - tmp_0*tmp_6 + tmp_10*tmp_8 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_13 - tmp_14*tmp_8);
      real_t tmp_16 = 1.0*tmp_15;
      real_t tmp_17 = tmp_16*tmp_7;
      real_t tmp_18 = tmp_12*tmp_4 - tmp_13;
      real_t tmp_19 = tmp_16*tmp_18;
      real_t tmp_20 = tmp_10 - tmp_14;
      real_t tmp_21 = tmp_16*tmp_20;
      real_t tmp_22 = tmp_0*tmp_17 + tmp_11*tmp_19 + tmp_21*tmp_8;
      real_t tmp_23 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_22;
      real_t tmp_24 = -2*tmp_17 - 2*tmp_19 - 2*tmp_21;
      real_t tmp_25 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_26 = tmp_0*tmp_1 - tmp_11*tmp_9;
      real_t tmp_27 = 0.5*tmp_15;
      real_t tmp_28 = tmp_26*tmp_27;
      real_t tmp_29 = -tmp_0*tmp_4 + tmp_8*tmp_9;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = -tmp_1*tmp_8 + tmp_11*tmp_4;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_27*tmp_7;
      real_t tmp_34 = tmp_18*tmp_27;
      real_t tmp_35 = tmp_20*tmp_27;
      real_t tmp_36 = tmp_0*tmp_32 + tmp_11*tmp_30 + tmp_12*tmp_33 + tmp_2*tmp_35 + tmp_28*tmp_8 + tmp_34*tmp_5;
      real_t tmp_37 = tmp_36*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_38 = -tmp_0*tmp_5 + tmp_11*tmp_12;
      real_t tmp_39 = tmp_27*tmp_38;
      real_t tmp_40 = tmp_0*tmp_2 - tmp_12*tmp_8;
      real_t tmp_41 = tmp_27*tmp_40;
      real_t tmp_42 = -tmp_11*tmp_2 + tmp_5*tmp_8;
      real_t tmp_43 = tmp_27*tmp_42;
      real_t tmp_44 = tmp_0*tmp_43 + tmp_1*tmp_34 + tmp_11*tmp_41 + tmp_33*tmp_9 + tmp_35*tmp_4 + tmp_39*tmp_8;
      real_t tmp_45 = tmp_44*(-tmp_39 - tmp_41 - tmp_43);
      real_t tmp_46 = p_affine_0_0*p_affine_1_1;
      real_t tmp_47 = p_affine_0_0*p_affine_1_2;
      real_t tmp_48 = p_affine_2_1*p_affine_3_2;
      real_t tmp_49 = p_affine_0_1*p_affine_1_0;
      real_t tmp_50 = p_affine_0_1*p_affine_1_2;
      real_t tmp_51 = p_affine_2_2*p_affine_3_0;
      real_t tmp_52 = p_affine_0_2*p_affine_1_0;
      real_t tmp_53 = p_affine_0_2*p_affine_1_1;
      real_t tmp_54 = p_affine_2_0*p_affine_3_1;
      real_t tmp_55 = p_affine_2_2*p_affine_3_1;
      real_t tmp_56 = p_affine_2_0*p_affine_3_2;
      real_t tmp_57 = p_affine_2_1*p_affine_3_0;
      real_t tmp_58 = std::abs(p_affine_0_0*tmp_48 - p_affine_0_0*tmp_55 + p_affine_0_1*tmp_51 - p_affine_0_1*tmp_56 + p_affine_0_2*tmp_54 - p_affine_0_2*tmp_57 - p_affine_1_0*tmp_48 + p_affine_1_0*tmp_55 - p_affine_1_1*tmp_51 + p_affine_1_1*tmp_56 - p_affine_1_2*tmp_54 + p_affine_1_2*tmp_57 + p_affine_2_0*tmp_50 - p_affine_2_0*tmp_53 - p_affine_2_1*tmp_47 + p_affine_2_1*tmp_52 + p_affine_2_2*tmp_46 - p_affine_2_2*tmp_49 - p_affine_3_0*tmp_50 + p_affine_3_0*tmp_53 + p_affine_3_1*tmp_47 - p_affine_3_1*tmp_52 - p_affine_3_2*tmp_46 + p_affine_3_2*tmp_49);
      real_t tmp_59 = 0.018781320953002646*tmp_58;
      real_t tmp_60 = tmp_22*tmp_24;
      real_t tmp_61 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id1;
      real_t tmp_62 = 0.012248840519393657*tmp_58;
      real_t tmp_63 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id2;
      real_t tmp_64 = 0.0070910034628469103*tmp_58;
      real_t tmp_65 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id3;
      real_t tmp_66 = 0.0070910034628469103*tmp_58;
      real_t tmp_67 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id4;
      real_t tmp_68 = 0.0070910034628469103*tmp_58;
      real_t tmp_69 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id5;
      real_t tmp_70 = 0.0070910034628469103*tmp_58;
      real_t tmp_71 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id6;
      real_t tmp_72 = 0.018781320953002646*tmp_58;
      real_t tmp_73 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id7;
      real_t tmp_74 = 0.012248840519393657*tmp_58;
      real_t tmp_75 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id8;
      real_t tmp_76 = 0.018781320953002646*tmp_58;
      real_t tmp_77 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id9;
      real_t tmp_78 = 0.012248840519393657*tmp_58;
      real_t tmp_79 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id10;
      real_t tmp_80 = 0.018781320953002646*tmp_58;
      real_t tmp_81 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id11;
      real_t tmp_82 = 0.012248840519393657*tmp_58;
      real_t tmp_83 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id12;
      real_t tmp_84 = 0.0070910034628469103*tmp_58;
      real_t tmp_85 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id13;
      real_t tmp_86 = 0.0070910034628469103*tmp_58;
      real_t tmp_87 = 2.0*tmp_15;
      real_t tmp_88 = tmp_7*tmp_87;
      real_t tmp_89 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_87;
      real_t tmp_90 = tmp_31*tmp_36;
      real_t tmp_91 = tmp_42*tmp_44;
      real_t tmp_92 = tmp_22*tmp_88;
      real_t tmp_93 = Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_87;
      real_t tmp_94 = Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_87;
      real_t tmp_95 = Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_87;
      real_t tmp_96 = Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_87;
      real_t tmp_97 = Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_87;
      real_t tmp_98 = Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_87;
      real_t tmp_99 = Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_87;
      real_t tmp_100 = Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_87;
      real_t tmp_101 = Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_87;
      real_t tmp_102 = Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_87;
      real_t tmp_103 = Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_87;
      real_t tmp_104 = Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_87;
      real_t tmp_105 = Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_87;
      real_t tmp_106 = tmp_23*tmp_87;
      real_t tmp_107 = tmp_29*tmp_36;
      real_t tmp_108 = tmp_40*tmp_44;
      real_t tmp_109 = tmp_18*tmp_22;
      real_t tmp_110 = tmp_26*tmp_36;
      real_t tmp_111 = tmp_38*tmp_44;
      real_t tmp_112 = tmp_20*tmp_22;
      real_t a_0_0 = tmp_59*(tmp_23*tmp_24 + tmp_25*tmp_37 + tmp_25*tmp_45) + tmp_62*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_60 + tmp_37*tmp_61 + tmp_45*tmp_61) + tmp_64*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_60 + tmp_37*tmp_63 + tmp_45*tmp_63) + tmp_66*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_60 + tmp_37*tmp_65 + tmp_45*tmp_65) + tmp_68*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_60 + tmp_37*tmp_67 + tmp_45*tmp_67) + tmp_70*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_60 + tmp_37*tmp_69 + tmp_45*tmp_69) + tmp_72*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_60 + tmp_37*tmp_71 + tmp_45*tmp_71) + tmp_74*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_60 + tmp_37*tmp_73 + tmp_45*tmp_73) + tmp_76*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_60 + tmp_37*tmp_75 + tmp_45*tmp_75) + tmp_78*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_60 + tmp_37*tmp_77 + tmp_45*tmp_77) + tmp_80*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_60 + tmp_37*tmp_79 + tmp_45*tmp_79) + tmp_82*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_60 + tmp_37*tmp_81 + tmp_45*tmp_81) + tmp_84*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_60 + tmp_37*tmp_83 + tmp_45*tmp_83) + tmp_86*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_60 + tmp_37*tmp_85 + tmp_45*tmp_85);
      real_t a_1_0 = tmp_59*(tmp_23*tmp_88 + tmp_89*tmp_90 + tmp_89*tmp_91) + tmp_62*(Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_92 + tmp_90*tmp_93 + tmp_91*tmp_93) + tmp_64*(Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_92 + tmp_90*tmp_94 + tmp_91*tmp_94) + tmp_66*(Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_92 + tmp_90*tmp_95 + tmp_91*tmp_95) + tmp_68*(Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_92 + tmp_90*tmp_96 + tmp_91*tmp_96) + tmp_70*(Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_92 + tmp_90*tmp_97 + tmp_91*tmp_97) + tmp_72*(Scalar_Variable_Coefficient_3D_mu_out0_id6*tmp_92 + tmp_90*tmp_98 + tmp_91*tmp_98) + tmp_74*(Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_92 + tmp_90*tmp_99 + tmp_91*tmp_99) + tmp_76*(Scalar_Variable_Coefficient_3D_mu_out0_id8*tmp_92 + tmp_100*tmp_90 + tmp_100*tmp_91) + tmp_78*(Scalar_Variable_Coefficient_3D_mu_out0_id9*tmp_92 + tmp_101*tmp_90 + tmp_101*tmp_91) + tmp_80*(Scalar_Variable_Coefficient_3D_mu_out0_id10*tmp_92 + tmp_102*tmp_90 + tmp_102*tmp_91) + tmp_82*(Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_92 + tmp_103*tmp_90 + tmp_103*tmp_91) + tmp_84*(Scalar_Variable_Coefficient_3D_mu_out0_id12*tmp_92 + tmp_104*tmp_90 + tmp_104*tmp_91) + tmp_86*(Scalar_Variable_Coefficient_3D_mu_out0_id13*tmp_92 + tmp_105*tmp_90 + tmp_105*tmp_91);
      real_t a_2_0 = tmp_59*(tmp_106*tmp_18 + tmp_107*tmp_89 + tmp_108*tmp_89) + tmp_62*(tmp_107*tmp_93 + tmp_108*tmp_93 + tmp_109*tmp_93) + tmp_64*(tmp_107*tmp_94 + tmp_108*tmp_94 + tmp_109*tmp_94) + tmp_66*(tmp_107*tmp_95 + tmp_108*tmp_95 + tmp_109*tmp_95) + tmp_68*(tmp_107*tmp_96 + tmp_108*tmp_96 + tmp_109*tmp_96) + tmp_70*(tmp_107*tmp_97 + tmp_108*tmp_97 + tmp_109*tmp_97) + tmp_72*(tmp_107*tmp_98 + tmp_108*tmp_98 + tmp_109*tmp_98) + tmp_74*(tmp_107*tmp_99 + tmp_108*tmp_99 + tmp_109*tmp_99) + tmp_76*(tmp_100*tmp_107 + tmp_100*tmp_108 + tmp_100*tmp_109) + tmp_78*(tmp_101*tmp_107 + tmp_101*tmp_108 + tmp_101*tmp_109) + tmp_80*(tmp_102*tmp_107 + tmp_102*tmp_108 + tmp_102*tmp_109) + tmp_82*(tmp_103*tmp_107 + tmp_103*tmp_108 + tmp_103*tmp_109) + tmp_84*(tmp_104*tmp_107 + tmp_104*tmp_108 + tmp_104*tmp_109) + tmp_86*(tmp_105*tmp_107 + tmp_105*tmp_108 + tmp_105*tmp_109);
      real_t a_3_0 = tmp_59*(tmp_106*tmp_20 + tmp_110*tmp_89 + tmp_111*tmp_89) + tmp_62*(tmp_110*tmp_93 + tmp_111*tmp_93 + tmp_112*tmp_93) + tmp_64*(tmp_110*tmp_94 + tmp_111*tmp_94 + tmp_112*tmp_94) + tmp_66*(tmp_110*tmp_95 + tmp_111*tmp_95 + tmp_112*tmp_95) + tmp_68*(tmp_110*tmp_96 + tmp_111*tmp_96 + tmp_112*tmp_96) + tmp_70*(tmp_110*tmp_97 + tmp_111*tmp_97 + tmp_112*tmp_97) + tmp_72*(tmp_110*tmp_98 + tmp_111*tmp_98 + tmp_112*tmp_98) + tmp_74*(tmp_110*tmp_99 + tmp_111*tmp_99 + tmp_112*tmp_99) + tmp_76*(tmp_100*tmp_110 + tmp_100*tmp_111 + tmp_100*tmp_112) + tmp_78*(tmp_101*tmp_110 + tmp_101*tmp_111 + tmp_101*tmp_112) + tmp_80*(tmp_102*tmp_110 + tmp_102*tmp_111 + tmp_102*tmp_112) + tmp_82*(tmp_103*tmp_110 + tmp_103*tmp_111 + tmp_103*tmp_112) + tmp_84*(tmp_104*tmp_110 + tmp_104*tmp_111 + tmp_104*tmp_112) + tmp_86*(tmp_105*tmp_110 + tmp_105*tmp_111 + tmp_105*tmp_112);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_6*tmp_9;
      real_t tmp_14 = tmp_4*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_12 + tmp_0*tmp_7 - tmp_1*tmp_13 + tmp_1*tmp_2*tmp_8 + tmp_11*tmp_3 - tmp_14*tmp_3);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_7;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = tmp_33 - 1.0/4.0;
      real_t tmp_35 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_36 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_37 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_38 = tmp_20*tmp_35 + tmp_26*tmp_36 + tmp_32*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/4.0;
      real_t tmp_40 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_41 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_42 = tmp_11 - tmp_14;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_26*tmp_41 + tmp_32*tmp_42;
      real_t tmp_44 = tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_34 + tmp_1*tmp_39 + tmp_3*tmp_44;
      real_t tmp_46 = 0.5*tmp_15;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = tmp_37*tmp_46;
      real_t tmp_49 = tmp_27*tmp_46;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = 0.5*p_affine_13_2;
      real_t tmp_52 = tmp_50*tmp_51;
      real_t tmp_53 = tmp_2*tmp_44 + tmp_34*tmp_9 + tmp_39*tmp_4;
      real_t tmp_54 = tmp_41*tmp_46;
      real_t tmp_55 = tmp_36*tmp_46;
      real_t tmp_56 = tmp_21*tmp_46;
      real_t tmp_57 = -tmp_54 - tmp_55 - tmp_56;
      real_t tmp_58 = tmp_51*tmp_57;
      real_t tmp_59 = tmp_10*tmp_39 + tmp_34*tmp_8 + tmp_44*tmp_6;
      real_t tmp_60 = 1.0*tmp_15;
      real_t tmp_61 = tmp_40*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_5*tmp_60;
      real_t tmp_64 = 0.5*p_affine_13_0*tmp_50 + 0.5*p_affine_13_1*tmp_57 + 0.5*p_affine_13_2*(-tmp_61 - tmp_62 - tmp_63);
      real_t tmp_65 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_66 = std::pow(tmp_65, -0.25);
      real_t tmp_67 = -tmp_33 - tmp_38 - tmp_43 + 1;
      real_t tmp_68 = tmp_46*tmp_5;
      real_t tmp_69 = tmp_35*tmp_46;
      real_t tmp_70 = tmp_40*tmp_46;
      real_t tmp_71 = 0.5*p_affine_13_0*(tmp_0*tmp_68 + tmp_1*tmp_69 + tmp_10*tmp_48 + tmp_3*tmp_70 + tmp_47*tmp_6 + tmp_49*tmp_8) + 0.5*p_affine_13_1*(tmp_10*tmp_55 + tmp_2*tmp_70 + tmp_4*tmp_69 + tmp_54*tmp_6 + tmp_56*tmp_8 + tmp_68*tmp_9) + 0.5*p_affine_13_2*(tmp_10*tmp_62 + tmp_6*tmp_61 + tmp_63*tmp_8);
      real_t tmp_72 = 2.0*std::pow(tmp_65, 1.0/2.0);
      real_t tmp_73 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_72;
      real_t tmp_74 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_75 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_76 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_77 = tmp_21*tmp_75 + tmp_27*tmp_76 + tmp_5*tmp_74;
      real_t tmp_78 = tmp_77 - 1.0/4.0;
      real_t tmp_79 = tmp_35*tmp_74 + tmp_36*tmp_75 + tmp_37*tmp_76;
      real_t tmp_80 = tmp_79 - 1.0/4.0;
      real_t tmp_81 = tmp_40*tmp_74 + tmp_41*tmp_75 + tmp_42*tmp_76;
      real_t tmp_82 = tmp_81 - 1.0/4.0;
      real_t tmp_83 = tmp_0*tmp_78 + tmp_1*tmp_80 + tmp_3*tmp_82;
      real_t tmp_84 = tmp_2*tmp_82 + tmp_4*tmp_80 + tmp_78*tmp_9;
      real_t tmp_85 = tmp_10*tmp_80 + tmp_6*tmp_82 + tmp_78*tmp_8;
      real_t tmp_86 = -tmp_77 - tmp_79 - tmp_81 + 1;
      real_t tmp_87 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_72;
      real_t tmp_88 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_89 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_90 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_91 = tmp_21*tmp_89 + tmp_27*tmp_90 + tmp_5*tmp_88;
      real_t tmp_92 = tmp_91 - 1.0/4.0;
      real_t tmp_93 = tmp_35*tmp_88 + tmp_36*tmp_89 + tmp_37*tmp_90;
      real_t tmp_94 = tmp_93 - 1.0/4.0;
      real_t tmp_95 = tmp_40*tmp_88 + tmp_41*tmp_89 + tmp_42*tmp_90;
      real_t tmp_96 = tmp_95 - 1.0/4.0;
      real_t tmp_97 = tmp_0*tmp_92 + tmp_1*tmp_94 + tmp_3*tmp_96;
      real_t tmp_98 = tmp_2*tmp_96 + tmp_4*tmp_94 + tmp_9*tmp_92;
      real_t tmp_99 = tmp_10*tmp_94 + tmp_6*tmp_96 + tmp_8*tmp_92;
      real_t tmp_100 = -tmp_91 - tmp_93 - tmp_95 + 1;
      real_t tmp_101 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_72;
      real_t tmp_102 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_103 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_104 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_105 = tmp_102*tmp_5 + tmp_103*tmp_21 + tmp_104*tmp_27;
      real_t tmp_106 = tmp_105 - 1.0/4.0;
      real_t tmp_107 = tmp_102*tmp_35 + tmp_103*tmp_36 + tmp_104*tmp_37;
      real_t tmp_108 = tmp_107 - 1.0/4.0;
      real_t tmp_109 = tmp_102*tmp_40 + tmp_103*tmp_41 + tmp_104*tmp_42;
      real_t tmp_110 = tmp_109 - 1.0/4.0;
      real_t tmp_111 = tmp_0*tmp_106 + tmp_1*tmp_108 + tmp_110*tmp_3;
      real_t tmp_112 = tmp_106*tmp_9 + tmp_108*tmp_4 + tmp_110*tmp_2;
      real_t tmp_113 = tmp_10*tmp_108 + tmp_106*tmp_8 + tmp_110*tmp_6;
      real_t tmp_114 = -tmp_105 - tmp_107 - tmp_109 + 1;
      real_t tmp_115 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_72;
      real_t tmp_116 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_117 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_118 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_119 = tmp_116*tmp_5 + tmp_117*tmp_21 + tmp_118*tmp_27;
      real_t tmp_120 = tmp_119 - 1.0/4.0;
      real_t tmp_121 = tmp_116*tmp_35 + tmp_117*tmp_36 + tmp_118*tmp_37;
      real_t tmp_122 = tmp_121 - 1.0/4.0;
      real_t tmp_123 = tmp_116*tmp_40 + tmp_117*tmp_41 + tmp_118*tmp_42;
      real_t tmp_124 = tmp_123 - 1.0/4.0;
      real_t tmp_125 = tmp_0*tmp_120 + tmp_1*tmp_122 + tmp_124*tmp_3;
      real_t tmp_126 = tmp_120*tmp_9 + tmp_122*tmp_4 + tmp_124*tmp_2;
      real_t tmp_127 = tmp_10*tmp_122 + tmp_120*tmp_8 + tmp_124*tmp_6;
      real_t tmp_128 = -tmp_119 - tmp_121 - tmp_123 + 1;
      real_t tmp_129 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_72;
      real_t tmp_130 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_131 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_132 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_133 = tmp_130*tmp_5 + tmp_131*tmp_21 + tmp_132*tmp_27;
      real_t tmp_134 = tmp_133 - 1.0/4.0;
      real_t tmp_135 = tmp_130*tmp_35 + tmp_131*tmp_36 + tmp_132*tmp_37;
      real_t tmp_136 = tmp_135 - 1.0/4.0;
      real_t tmp_137 = tmp_130*tmp_40 + tmp_131*tmp_41 + tmp_132*tmp_42;
      real_t tmp_138 = tmp_137 - 1.0/4.0;
      real_t tmp_139 = tmp_0*tmp_134 + tmp_1*tmp_136 + tmp_138*tmp_3;
      real_t tmp_140 = tmp_134*tmp_9 + tmp_136*tmp_4 + tmp_138*tmp_2;
      real_t tmp_141 = tmp_10*tmp_136 + tmp_134*tmp_8 + tmp_138*tmp_6;
      real_t tmp_142 = -tmp_133 - tmp_135 - tmp_137 + 1;
      real_t tmp_143 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_72;
      real_t tmp_144 = 0.25*p_affine_13_2*tmp_15;
      real_t tmp_145 = tmp_144*tmp_21;
      real_t tmp_146 = tmp_144*tmp_27;
      real_t tmp_147 = 0.5*p_affine_13_0*tmp_49 + 0.5*p_affine_13_1*tmp_56 + 0.5*p_affine_13_2*tmp_63;
      real_t tmp_148 = tmp_144*tmp_36;
      real_t tmp_149 = tmp_144*tmp_37;
      real_t tmp_150 = 0.5*p_affine_13_0*tmp_48 + 0.5*p_affine_13_1*tmp_55 + 0.5*p_affine_13_2*tmp_62;
      real_t tmp_151 = tmp_144*tmp_41;
      real_t tmp_152 = tmp_144*tmp_42;
      real_t tmp_153 = 0.5*p_affine_13_0*tmp_47 + 0.5*p_affine_13_1*tmp_54 + 0.5*p_affine_13_2*tmp_61;
      real_t a_0_0 = tmp_101*(1.0*tmp_100*tmp_66*tmp_99 - tmp_100*tmp_71 - tmp_52*tmp_97 - tmp_58*tmp_98 - tmp_64*tmp_99) + tmp_115*(-tmp_111*tmp_52 - tmp_112*tmp_58 + 1.0*tmp_113*tmp_114*tmp_66 - tmp_113*tmp_64 - tmp_114*tmp_71) + tmp_129*(-tmp_125*tmp_52 - tmp_126*tmp_58 + 1.0*tmp_127*tmp_128*tmp_66 - tmp_127*tmp_64 - tmp_128*tmp_71) + tmp_143*(-tmp_139*tmp_52 - tmp_140*tmp_58 + 1.0*tmp_141*tmp_142*tmp_66 - tmp_141*tmp_64 - tmp_142*tmp_71) + tmp_73*(-tmp_45*tmp_52 - tmp_53*tmp_58 - tmp_59*tmp_64 + 1.0*tmp_59*tmp_66*tmp_67 - tmp_67*tmp_71) + tmp_87*(-tmp_52*tmp_83 - tmp_58*tmp_84 - tmp_64*tmp_85 + 1.0*tmp_66*tmp_85*tmp_86 - tmp_71*tmp_86);
      real_t a_1_0 = tmp_101*(-tmp_145*tmp_98 - tmp_146*tmp_97 - tmp_147*tmp_99 + 1.0*tmp_66*tmp_91*tmp_99 - tmp_71*tmp_91) + tmp_115*(1.0*tmp_105*tmp_113*tmp_66 - tmp_105*tmp_71 - tmp_111*tmp_146 - tmp_112*tmp_145 - tmp_113*tmp_147) + tmp_129*(1.0*tmp_119*tmp_127*tmp_66 - tmp_119*tmp_71 - tmp_125*tmp_146 - tmp_126*tmp_145 - tmp_127*tmp_147) + tmp_143*(1.0*tmp_133*tmp_141*tmp_66 - tmp_133*tmp_71 - tmp_139*tmp_146 - tmp_140*tmp_145 - tmp_141*tmp_147) + tmp_73*(-tmp_145*tmp_53 - tmp_146*tmp_45 - tmp_147*tmp_59 + 1.0*tmp_33*tmp_59*tmp_66 - tmp_33*tmp_71) + tmp_87*(-tmp_145*tmp_84 - tmp_146*tmp_83 - tmp_147*tmp_85 + 1.0*tmp_66*tmp_77*tmp_85 - tmp_71*tmp_77);
      real_t a_2_0 = tmp_101*(-tmp_148*tmp_98 - tmp_149*tmp_97 - tmp_150*tmp_99 + 1.0*tmp_66*tmp_93*tmp_99 - tmp_71*tmp_93) + tmp_115*(1.0*tmp_107*tmp_113*tmp_66 - tmp_107*tmp_71 - tmp_111*tmp_149 - tmp_112*tmp_148 - tmp_113*tmp_150) + tmp_129*(1.0*tmp_121*tmp_127*tmp_66 - tmp_121*tmp_71 - tmp_125*tmp_149 - tmp_126*tmp_148 - tmp_127*tmp_150) + tmp_143*(1.0*tmp_135*tmp_141*tmp_66 - tmp_135*tmp_71 - tmp_139*tmp_149 - tmp_140*tmp_148 - tmp_141*tmp_150) + tmp_73*(-tmp_148*tmp_53 - tmp_149*tmp_45 - tmp_150*tmp_59 + 1.0*tmp_38*tmp_59*tmp_66 - tmp_38*tmp_71) + tmp_87*(-tmp_148*tmp_84 - tmp_149*tmp_83 - tmp_150*tmp_85 + 1.0*tmp_66*tmp_79*tmp_85 - tmp_71*tmp_79);
      real_t a_3_0 = tmp_101*(-tmp_151*tmp_98 - tmp_152*tmp_97 - tmp_153*tmp_99 + 1.0*tmp_66*tmp_95*tmp_99 - tmp_71*tmp_95) + tmp_115*(1.0*tmp_109*tmp_113*tmp_66 - tmp_109*tmp_71 - tmp_111*tmp_152 - tmp_112*tmp_151 - tmp_113*tmp_153) + tmp_129*(1.0*tmp_123*tmp_127*tmp_66 - tmp_123*tmp_71 - tmp_125*tmp_152 - tmp_126*tmp_151 - tmp_127*tmp_153) + tmp_143*(1.0*tmp_137*tmp_141*tmp_66 - tmp_137*tmp_71 - tmp_139*tmp_152 - tmp_140*tmp_151 - tmp_141*tmp_153) + tmp_73*(-tmp_151*tmp_53 - tmp_152*tmp_45 - tmp_153*tmp_59 + 1.0*tmp_43*tmp_59*tmp_66 - tmp_43*tmp_71) + tmp_87*(-tmp_151*tmp_84 - tmp_152*tmp_83 - tmp_153*tmp_85 + 1.0*tmp_66*tmp_81*tmp_85 - tmp_71*tmp_81);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_1 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_2 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_5 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_11 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_12 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_13 = tmp_12*tmp_2;
      real_t tmp_14 = tmp_1*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_13 + tmp_0*tmp_9 + tmp_10*tmp_3 - tmp_10*tmp_6 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_14);
      real_t tmp_16 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_17 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_18 = -tmp_17;
      real_t tmp_19 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_20 = 0.091576213509770743*tmp_18 + 0.81684757298045851*tmp_19;
      real_t tmp_21 = tmp_15*(tmp_16 + tmp_20);
      real_t tmp_22 = tmp_12*tmp_4 - tmp_14;
      real_t tmp_23 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_24 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_25 = -tmp_24;
      real_t tmp_26 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_27 = 0.091576213509770743*tmp_25 + 0.81684757298045851*tmp_26;
      real_t tmp_28 = tmp_15*(tmp_23 + tmp_27);
      real_t tmp_29 = -tmp_13 + tmp_9;
      real_t tmp_30 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_31 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_32 = -tmp_31;
      real_t tmp_33 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_34 = 0.091576213509770743*tmp_32 + 0.81684757298045851*tmp_33;
      real_t tmp_35 = tmp_15*(tmp_30 + tmp_34);
      real_t tmp_36 = tmp_21*tmp_7 + tmp_22*tmp_28 + tmp_29*tmp_35 - 1.0/4.0;
      real_t tmp_37 = -tmp_0*tmp_2 + tmp_11*tmp_4;
      real_t tmp_38 = tmp_0*tmp_8 - tmp_10*tmp_4;
      real_t tmp_39 = tmp_10*tmp_2 - tmp_11*tmp_8;
      real_t tmp_40 = tmp_21*tmp_37 + tmp_28*tmp_38 + tmp_35*tmp_39 - 1.0/4.0;
      real_t tmp_41 = tmp_0*tmp_5 - tmp_1*tmp_11;
      real_t tmp_42 = -tmp_0*tmp_12 + tmp_1*tmp_10;
      real_t tmp_43 = -tmp_10*tmp_5 + tmp_11*tmp_12;
      real_t tmp_44 = tmp_21*tmp_41 + tmp_28*tmp_42 + tmp_35*tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_36 + tmp_1*tmp_40 + tmp_4*tmp_44;
      real_t tmp_46 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_47 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_48 = tmp_46*tmp_47;
      real_t tmp_49 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_50 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_51 = tmp_49*tmp_50;
      real_t tmp_52 = tmp_48 - tmp_51;
      real_t tmp_53 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_54 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_55 = tmp_49*tmp_54;
      real_t tmp_56 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_57 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_58 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_59 = tmp_47*tmp_57;
      real_t tmp_60 = tmp_46*tmp_54;
      real_t tmp_61 = 1.0 / (tmp_48*tmp_58 + tmp_50*tmp_56*tmp_57 - tmp_51*tmp_58 + tmp_53*tmp_55 - tmp_53*tmp_59 - tmp_56*tmp_60);
      real_t tmp_62 = 0.5*tmp_61;
      real_t tmp_63 = tmp_52*tmp_62;
      real_t tmp_64 = tmp_50*tmp_57 - tmp_60;
      real_t tmp_65 = tmp_62*tmp_64;
      real_t tmp_66 = tmp_55 - tmp_59;
      real_t tmp_67 = tmp_62*tmp_66;
      real_t tmp_68 = -tmp_63 - tmp_65 - tmp_67;
      real_t tmp_69 = 0.5*p_affine_13_2;
      real_t tmp_70 = tmp_68*tmp_69;
      real_t tmp_71 = tmp_11*tmp_36 + tmp_2*tmp_44 + tmp_40*tmp_5;
      real_t tmp_72 = -tmp_47*tmp_53 + tmp_50*tmp_56;
      real_t tmp_73 = tmp_62*tmp_72;
      real_t tmp_74 = -tmp_50*tmp_58 + tmp_53*tmp_54;
      real_t tmp_75 = tmp_62*tmp_74;
      real_t tmp_76 = tmp_47*tmp_58 - tmp_54*tmp_56;
      real_t tmp_77 = tmp_62*tmp_76;
      real_t tmp_78 = -tmp_73 - tmp_75 - tmp_77;
      real_t tmp_79 = tmp_69*tmp_78;
      real_t tmp_80 = tmp_10*tmp_36 + tmp_12*tmp_40 + tmp_44*tmp_8;
      real_t tmp_81 = -tmp_46*tmp_56 + tmp_49*tmp_53;
      real_t tmp_82 = 1.0*tmp_61;
      real_t tmp_83 = tmp_81*tmp_82;
      real_t tmp_84 = tmp_46*tmp_58 - tmp_53*tmp_57;
      real_t tmp_85 = tmp_82*tmp_84;
      real_t tmp_86 = -tmp_49*tmp_58 + tmp_56*tmp_57;
      real_t tmp_87 = tmp_82*tmp_86;
      real_t tmp_88 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_78 + 0.5*p_affine_13_2*(-tmp_83 - tmp_85 - tmp_87);
      real_t tmp_89 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_90 = tmp_61*(tmp_20 + tmp_89);
      real_t tmp_91 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_92 = tmp_61*(tmp_27 + tmp_91);
      real_t tmp_93 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_94 = tmp_61*(tmp_34 + tmp_93);
      real_t tmp_95 = tmp_66*tmp_94 + tmp_76*tmp_92 + tmp_86*tmp_90;
      real_t tmp_96 = tmp_64*tmp_94 + tmp_74*tmp_92 + tmp_84*tmp_90;
      real_t tmp_97 = tmp_52*tmp_94 + tmp_72*tmp_92 + tmp_81*tmp_90;
      real_t tmp_98 = -tmp_95 - tmp_96 - tmp_97 + 1;
      real_t tmp_99 = (std::abs(tmp_17*tmp_26 - tmp_19*tmp_24)*std::abs(tmp_17*tmp_26 - tmp_19*tmp_24)) + (std::abs(tmp_17*tmp_33 - tmp_19*tmp_31)*std::abs(tmp_17*tmp_33 - tmp_19*tmp_31)) + (std::abs(tmp_24*tmp_33 - tmp_26*tmp_31)*std::abs(tmp_24*tmp_33 - tmp_26*tmp_31));
      real_t tmp_100 = 1.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_80;
      real_t tmp_102 = 1.0*tmp_15;
      real_t tmp_103 = 0.5*tmp_15;
      real_t tmp_104 = tmp_103*tmp_7;
      real_t tmp_105 = tmp_103*tmp_37;
      real_t tmp_106 = tmp_103*tmp_41;
      real_t tmp_107 = tmp_10*tmp_103;
      real_t tmp_108 = tmp_103*tmp_12;
      real_t tmp_109 = tmp_103*tmp_8;
      real_t tmp_110 = 0.5*p_affine_13_0*(tmp_0*tmp_104 + tmp_1*tmp_105 + tmp_106*tmp_4 + tmp_107*tmp_29 + tmp_108*tmp_39 + tmp_109*tmp_43) + 0.5*p_affine_13_1*(tmp_104*tmp_11 + tmp_105*tmp_5 + tmp_106*tmp_2 + tmp_107*tmp_22 + tmp_108*tmp_38 + tmp_109*tmp_42) + 0.5*p_affine_13_2*(tmp_10*tmp_102*tmp_7 + tmp_102*tmp_12*tmp_37 + tmp_102*tmp_41*tmp_8);
      real_t tmp_111 = 2.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_112 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_111;
      real_t tmp_113 = 0.44594849091596489*tmp_18 + 0.10810301816807022*tmp_19;
      real_t tmp_114 = tmp_15*(tmp_113 + tmp_16);
      real_t tmp_115 = 0.44594849091596489*tmp_25 + 0.10810301816807022*tmp_26;
      real_t tmp_116 = tmp_15*(tmp_115 + tmp_23);
      real_t tmp_117 = 0.44594849091596489*tmp_32 + 0.10810301816807022*tmp_33;
      real_t tmp_118 = tmp_15*(tmp_117 + tmp_30);
      real_t tmp_119 = tmp_114*tmp_7 + tmp_116*tmp_22 + tmp_118*tmp_29 - 1.0/4.0;
      real_t tmp_120 = tmp_114*tmp_37 + tmp_116*tmp_38 + tmp_118*tmp_39 - 1.0/4.0;
      real_t tmp_121 = tmp_114*tmp_41 + tmp_116*tmp_42 + tmp_118*tmp_43 - 1.0/4.0;
      real_t tmp_122 = tmp_0*tmp_119 + tmp_1*tmp_120 + tmp_121*tmp_4;
      real_t tmp_123 = tmp_11*tmp_119 + tmp_120*tmp_5 + tmp_121*tmp_2;
      real_t tmp_124 = tmp_10*tmp_119 + tmp_12*tmp_120 + tmp_121*tmp_8;
      real_t tmp_125 = tmp_61*(tmp_113 + tmp_89);
      real_t tmp_126 = tmp_61*(tmp_115 + tmp_91);
      real_t tmp_127 = tmp_61*(tmp_117 + tmp_93);
      real_t tmp_128 = tmp_125*tmp_86 + tmp_126*tmp_76 + tmp_127*tmp_66;
      real_t tmp_129 = tmp_125*tmp_84 + tmp_126*tmp_74 + tmp_127*tmp_64;
      real_t tmp_130 = tmp_125*tmp_81 + tmp_126*tmp_72 + tmp_127*tmp_52;
      real_t tmp_131 = -tmp_128 - tmp_129 - tmp_130 + 1;
      real_t tmp_132 = tmp_100*tmp_124;
      real_t tmp_133 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_111;
      real_t tmp_134 = 0.81684757298045851*tmp_18 + 0.091576213509770743*tmp_19;
      real_t tmp_135 = tmp_15*(tmp_134 + tmp_16);
      real_t tmp_136 = 0.81684757298045851*tmp_25 + 0.091576213509770743*tmp_26;
      real_t tmp_137 = tmp_15*(tmp_136 + tmp_23);
      real_t tmp_138 = 0.81684757298045851*tmp_32 + 0.091576213509770743*tmp_33;
      real_t tmp_139 = tmp_15*(tmp_138 + tmp_30);
      real_t tmp_140 = tmp_135*tmp_7 + tmp_137*tmp_22 + tmp_139*tmp_29 - 1.0/4.0;
      real_t tmp_141 = tmp_135*tmp_37 + tmp_137*tmp_38 + tmp_139*tmp_39 - 1.0/4.0;
      real_t tmp_142 = tmp_135*tmp_41 + tmp_137*tmp_42 + tmp_139*tmp_43 - 1.0/4.0;
      real_t tmp_143 = tmp_0*tmp_140 + tmp_1*tmp_141 + tmp_142*tmp_4;
      real_t tmp_144 = tmp_11*tmp_140 + tmp_141*tmp_5 + tmp_142*tmp_2;
      real_t tmp_145 = tmp_10*tmp_140 + tmp_12*tmp_141 + tmp_142*tmp_8;
      real_t tmp_146 = tmp_61*(tmp_134 + tmp_89);
      real_t tmp_147 = tmp_61*(tmp_136 + tmp_91);
      real_t tmp_148 = tmp_61*(tmp_138 + tmp_93);
      real_t tmp_149 = tmp_146*tmp_86 + tmp_147*tmp_76 + tmp_148*tmp_66;
      real_t tmp_150 = tmp_146*tmp_84 + tmp_147*tmp_74 + tmp_148*tmp_64;
      real_t tmp_151 = tmp_146*tmp_81 + tmp_147*tmp_72 + tmp_148*tmp_52;
      real_t tmp_152 = -tmp_149 - tmp_150 - tmp_151 + 1;
      real_t tmp_153 = tmp_100*tmp_145;
      real_t tmp_154 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_111;
      real_t tmp_155 = 0.10810301816807022*tmp_18 + 0.44594849091596489*tmp_19;
      real_t tmp_156 = tmp_15*(tmp_155 + tmp_16);
      real_t tmp_157 = 0.10810301816807022*tmp_25 + 0.44594849091596489*tmp_26;
      real_t tmp_158 = tmp_15*(tmp_157 + tmp_23);
      real_t tmp_159 = 0.10810301816807022*tmp_32 + 0.44594849091596489*tmp_33;
      real_t tmp_160 = tmp_15*(tmp_159 + tmp_30);
      real_t tmp_161 = tmp_156*tmp_7 + tmp_158*tmp_22 + tmp_160*tmp_29 - 1.0/4.0;
      real_t tmp_162 = tmp_156*tmp_37 + tmp_158*tmp_38 + tmp_160*tmp_39 - 1.0/4.0;
      real_t tmp_163 = tmp_156*tmp_41 + tmp_158*tmp_42 + tmp_160*tmp_43 - 1.0/4.0;
      real_t tmp_164 = tmp_0*tmp_161 + tmp_1*tmp_162 + tmp_163*tmp_4;
      real_t tmp_165 = tmp_11*tmp_161 + tmp_162*tmp_5 + tmp_163*tmp_2;
      real_t tmp_166 = tmp_10*tmp_161 + tmp_12*tmp_162 + tmp_163*tmp_8;
      real_t tmp_167 = tmp_61*(tmp_155 + tmp_89);
      real_t tmp_168 = tmp_61*(tmp_157 + tmp_91);
      real_t tmp_169 = tmp_61*(tmp_159 + tmp_93);
      real_t tmp_170 = tmp_167*tmp_86 + tmp_168*tmp_76 + tmp_169*tmp_66;
      real_t tmp_171 = tmp_167*tmp_84 + tmp_168*tmp_74 + tmp_169*tmp_64;
      real_t tmp_172 = tmp_167*tmp_81 + tmp_168*tmp_72 + tmp_169*tmp_52;
      real_t tmp_173 = -tmp_170 - tmp_171 - tmp_172 + 1;
      real_t tmp_174 = tmp_100*tmp_166;
      real_t tmp_175 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_111;
      real_t tmp_176 = 0.091576213509770743*tmp_18 + 0.091576213509770743*tmp_19;
      real_t tmp_177 = tmp_15*(tmp_16 + tmp_176);
      real_t tmp_178 = 0.091576213509770743*tmp_25 + 0.091576213509770743*tmp_26;
      real_t tmp_179 = tmp_15*(tmp_178 + tmp_23);
      real_t tmp_180 = 0.091576213509770743*tmp_32 + 0.091576213509770743*tmp_33;
      real_t tmp_181 = tmp_15*(tmp_180 + tmp_30);
      real_t tmp_182 = tmp_177*tmp_7 + tmp_179*tmp_22 + tmp_181*tmp_29 - 1.0/4.0;
      real_t tmp_183 = tmp_177*tmp_37 + tmp_179*tmp_38 + tmp_181*tmp_39 - 1.0/4.0;
      real_t tmp_184 = tmp_177*tmp_41 + tmp_179*tmp_42 + tmp_181*tmp_43 - 1.0/4.0;
      real_t tmp_185 = tmp_0*tmp_182 + tmp_1*tmp_183 + tmp_184*tmp_4;
      real_t tmp_186 = tmp_11*tmp_182 + tmp_183*tmp_5 + tmp_184*tmp_2;
      real_t tmp_187 = tmp_10*tmp_182 + tmp_12*tmp_183 + tmp_184*tmp_8;
      real_t tmp_188 = tmp_61*(tmp_176 + tmp_89);
      real_t tmp_189 = tmp_61*(tmp_178 + tmp_91);
      real_t tmp_190 = tmp_61*(tmp_180 + tmp_93);
      real_t tmp_191 = tmp_188*tmp_86 + tmp_189*tmp_76 + tmp_190*tmp_66;
      real_t tmp_192 = tmp_188*tmp_84 + tmp_189*tmp_74 + tmp_190*tmp_64;
      real_t tmp_193 = tmp_188*tmp_81 + tmp_189*tmp_72 + tmp_190*tmp_52;
      real_t tmp_194 = -tmp_191 - tmp_192 - tmp_193 + 1;
      real_t tmp_195 = tmp_100*tmp_187;
      real_t tmp_196 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_111;
      real_t tmp_197 = 0.44594849091596489*tmp_18 + 0.44594849091596489*tmp_19;
      real_t tmp_198 = tmp_15*(tmp_16 + tmp_197);
      real_t tmp_199 = 0.44594849091596489*tmp_25 + 0.44594849091596489*tmp_26;
      real_t tmp_200 = tmp_15*(tmp_199 + tmp_23);
      real_t tmp_201 = 0.44594849091596489*tmp_32 + 0.44594849091596489*tmp_33;
      real_t tmp_202 = tmp_15*(tmp_201 + tmp_30);
      real_t tmp_203 = tmp_198*tmp_7 + tmp_200*tmp_22 + tmp_202*tmp_29 - 1.0/4.0;
      real_t tmp_204 = tmp_198*tmp_37 + tmp_200*tmp_38 + tmp_202*tmp_39 - 1.0/4.0;
      real_t tmp_205 = tmp_198*tmp_41 + tmp_200*tmp_42 + tmp_202*tmp_43 - 1.0/4.0;
      real_t tmp_206 = tmp_0*tmp_203 + tmp_1*tmp_204 + tmp_205*tmp_4;
      real_t tmp_207 = tmp_11*tmp_203 + tmp_2*tmp_205 + tmp_204*tmp_5;
      real_t tmp_208 = tmp_10*tmp_203 + tmp_12*tmp_204 + tmp_205*tmp_8;
      real_t tmp_209 = tmp_61*(tmp_197 + tmp_89);
      real_t tmp_210 = tmp_61*(tmp_199 + tmp_91);
      real_t tmp_211 = tmp_61*(tmp_201 + tmp_93);
      real_t tmp_212 = tmp_209*tmp_86 + tmp_210*tmp_76 + tmp_211*tmp_66;
      real_t tmp_213 = tmp_209*tmp_84 + tmp_210*tmp_74 + tmp_211*tmp_64;
      real_t tmp_214 = tmp_209*tmp_81 + tmp_210*tmp_72 + tmp_211*tmp_52;
      real_t tmp_215 = -tmp_212 - tmp_213 - tmp_214 + 1;
      real_t tmp_216 = tmp_100*tmp_208;
      real_t tmp_217 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_111;
      real_t tmp_218 = 0.25*p_affine_13_2*tmp_61;
      real_t tmp_219 = tmp_218*tmp_76;
      real_t tmp_220 = tmp_218*tmp_66;
      real_t tmp_221 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*tmp_77 + 0.5*p_affine_13_2*tmp_87;
      real_t tmp_222 = tmp_218*tmp_74;
      real_t tmp_223 = tmp_218*tmp_64;
      real_t tmp_224 = 0.5*p_affine_13_0*tmp_65 + 0.5*p_affine_13_1*tmp_75 + 0.5*p_affine_13_2*tmp_85;
      real_t tmp_225 = tmp_218*tmp_72;
      real_t tmp_226 = tmp_218*tmp_52;
      real_t tmp_227 = 0.5*p_affine_13_0*tmp_63 + 0.5*p_affine_13_1*tmp_73 + 0.5*p_affine_13_2*tmp_83;
      real_t a_0_0 = tmp_112*(-tmp_101*tmp_98 - tmp_110*tmp_98 + tmp_45*tmp_70 + tmp_71*tmp_79 + tmp_80*tmp_88) + tmp_133*(-tmp_110*tmp_131 + tmp_122*tmp_70 + tmp_123*tmp_79 + tmp_124*tmp_88 - tmp_131*tmp_132) + tmp_154*(-tmp_110*tmp_152 + tmp_143*tmp_70 + tmp_144*tmp_79 + tmp_145*tmp_88 - tmp_152*tmp_153) + tmp_175*(-tmp_110*tmp_173 + tmp_164*tmp_70 + tmp_165*tmp_79 + tmp_166*tmp_88 - tmp_173*tmp_174) + tmp_196*(-tmp_110*tmp_194 + tmp_185*tmp_70 + tmp_186*tmp_79 + tmp_187*tmp_88 - tmp_194*tmp_195) + tmp_217*(-tmp_110*tmp_215 + tmp_206*tmp_70 + tmp_207*tmp_79 + tmp_208*tmp_88 - tmp_215*tmp_216);
      real_t a_1_0 = tmp_112*(-tmp_101*tmp_95 - tmp_110*tmp_95 + tmp_219*tmp_71 + tmp_220*tmp_45 + tmp_221*tmp_80) + tmp_133*(-tmp_110*tmp_128 + tmp_122*tmp_220 + tmp_123*tmp_219 + tmp_124*tmp_221 - tmp_128*tmp_132) + tmp_154*(-tmp_110*tmp_149 + tmp_143*tmp_220 + tmp_144*tmp_219 + tmp_145*tmp_221 - tmp_149*tmp_153) + tmp_175*(-tmp_110*tmp_170 + tmp_164*tmp_220 + tmp_165*tmp_219 + tmp_166*tmp_221 - tmp_170*tmp_174) + tmp_196*(-tmp_110*tmp_191 + tmp_185*tmp_220 + tmp_186*tmp_219 + tmp_187*tmp_221 - tmp_191*tmp_195) + tmp_217*(-tmp_110*tmp_212 + tmp_206*tmp_220 + tmp_207*tmp_219 + tmp_208*tmp_221 - tmp_212*tmp_216);
      real_t a_2_0 = tmp_112*(-tmp_101*tmp_96 - tmp_110*tmp_96 + tmp_222*tmp_71 + tmp_223*tmp_45 + tmp_224*tmp_80) + tmp_133*(-tmp_110*tmp_129 + tmp_122*tmp_223 + tmp_123*tmp_222 + tmp_124*tmp_224 - tmp_129*tmp_132) + tmp_154*(-tmp_110*tmp_150 + tmp_143*tmp_223 + tmp_144*tmp_222 + tmp_145*tmp_224 - tmp_150*tmp_153) + tmp_175*(-tmp_110*tmp_171 + tmp_164*tmp_223 + tmp_165*tmp_222 + tmp_166*tmp_224 - tmp_171*tmp_174) + tmp_196*(-tmp_110*tmp_192 + tmp_185*tmp_223 + tmp_186*tmp_222 + tmp_187*tmp_224 - tmp_192*tmp_195) + tmp_217*(-tmp_110*tmp_213 + tmp_206*tmp_223 + tmp_207*tmp_222 + tmp_208*tmp_224 - tmp_213*tmp_216);
      real_t a_3_0 = tmp_112*(-tmp_101*tmp_97 - tmp_110*tmp_97 + tmp_225*tmp_71 + tmp_226*tmp_45 + tmp_227*tmp_80) + tmp_133*(-tmp_110*tmp_130 + tmp_122*tmp_226 + tmp_123*tmp_225 + tmp_124*tmp_227 - tmp_130*tmp_132) + tmp_154*(-tmp_110*tmp_151 + tmp_143*tmp_226 + tmp_144*tmp_225 + tmp_145*tmp_227 - tmp_151*tmp_153) + tmp_175*(-tmp_110*tmp_172 + tmp_164*tmp_226 + tmp_165*tmp_225 + tmp_166*tmp_227 - tmp_172*tmp_174) + tmp_196*(-tmp_110*tmp_193 + tmp_185*tmp_226 + tmp_186*tmp_225 + tmp_187*tmp_227 - tmp_193*tmp_195) + tmp_217*(-tmp_110*tmp_214 + tmp_206*tmp_226 + tmp_207*tmp_225 + tmp_208*tmp_227 - tmp_214*tmp_216);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_6*tmp_9;
      real_t tmp_14 = tmp_4*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_12 + tmp_0*tmp_7 - tmp_1*tmp_13 + tmp_1*tmp_2*tmp_8 + tmp_11*tmp_3 - tmp_14*tmp_3);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_7;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = tmp_33 - 1.0/4.0;
      real_t tmp_35 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_36 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_37 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_38 = tmp_20*tmp_35 + tmp_26*tmp_36 + tmp_32*tmp_37;
      real_t tmp_39 = tmp_38 - 1.0/4.0;
      real_t tmp_40 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_41 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_42 = tmp_11 - tmp_14;
      real_t tmp_43 = tmp_20*tmp_40 + tmp_26*tmp_41 + tmp_32*tmp_42;
      real_t tmp_44 = tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_34 + tmp_1*tmp_39 + tmp_3*tmp_44;
      real_t tmp_46 = 0.5*tmp_15;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = tmp_37*tmp_46;
      real_t tmp_49 = tmp_27*tmp_46;
      real_t tmp_50 = -tmp_47 - tmp_48 - tmp_49;
      real_t tmp_51 = p_affine_13_2*tmp_50;
      real_t tmp_52 = tmp_2*tmp_44 + tmp_34*tmp_9 + tmp_39*tmp_4;
      real_t tmp_53 = tmp_41*tmp_46;
      real_t tmp_54 = tmp_36*tmp_46;
      real_t tmp_55 = tmp_21*tmp_46;
      real_t tmp_56 = -tmp_53 - tmp_54 - tmp_55;
      real_t tmp_57 = p_affine_13_2*tmp_56;
      real_t tmp_58 = 1.0*tmp_15;
      real_t tmp_59 = tmp_40*tmp_58;
      real_t tmp_60 = tmp_35*tmp_58;
      real_t tmp_61 = tmp_5*tmp_58;
      real_t tmp_62 = p_affine_13_0*tmp_50 + p_affine_13_1*tmp_56 + p_affine_13_2*(-tmp_59 - tmp_60 - tmp_61);
      real_t tmp_63 = tmp_10*tmp_39 + tmp_34*tmp_8 + tmp_44*tmp_6;
      real_t tmp_64 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_65 = std::pow(tmp_64, -0.25);
      real_t tmp_66 = -tmp_33 - tmp_38 - tmp_43 + 1;
      real_t tmp_67 = tmp_46*tmp_5;
      real_t tmp_68 = tmp_35*tmp_46;
      real_t tmp_69 = tmp_40*tmp_46;
      real_t tmp_70 = p_affine_13_0*(tmp_0*tmp_67 + tmp_1*tmp_68 + tmp_10*tmp_48 + tmp_3*tmp_69 + tmp_47*tmp_6 + tmp_49*tmp_8) + p_affine_13_1*(tmp_10*tmp_54 + tmp_2*tmp_69 + tmp_4*tmp_68 + tmp_53*tmp_6 + tmp_55*tmp_8 + tmp_67*tmp_9) + p_affine_13_2*(tmp_10*tmp_60 + tmp_59*tmp_6 + tmp_61*tmp_8);
      real_t tmp_71 = 2.0*std::pow(tmp_64, 1.0/2.0);
      real_t tmp_72 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_71;
      real_t tmp_73 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_74 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_75 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_76 = tmp_21*tmp_74 + tmp_27*tmp_75 + tmp_5*tmp_73;
      real_t tmp_77 = tmp_76 - 1.0/4.0;
      real_t tmp_78 = tmp_35*tmp_73 + tmp_36*tmp_74 + tmp_37*tmp_75;
      real_t tmp_79 = tmp_78 - 1.0/4.0;
      real_t tmp_80 = tmp_40*tmp_73 + tmp_41*tmp_74 + tmp_42*tmp_75;
      real_t tmp_81 = tmp_80 - 1.0/4.0;
      real_t tmp_82 = tmp_0*tmp_77 + tmp_1*tmp_79 + tmp_3*tmp_81;
      real_t tmp_83 = tmp_2*tmp_81 + tmp_4*tmp_79 + tmp_77*tmp_9;
      real_t tmp_84 = tmp_10*tmp_79 + tmp_6*tmp_81 + tmp_77*tmp_8;
      real_t tmp_85 = -tmp_76 - tmp_78 - tmp_80 + 1;
      real_t tmp_86 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_71;
      real_t tmp_87 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_88 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_89 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_90 = tmp_21*tmp_88 + tmp_27*tmp_89 + tmp_5*tmp_87;
      real_t tmp_91 = tmp_90 - 1.0/4.0;
      real_t tmp_92 = tmp_35*tmp_87 + tmp_36*tmp_88 + tmp_37*tmp_89;
      real_t tmp_93 = tmp_92 - 1.0/4.0;
      real_t tmp_94 = tmp_40*tmp_87 + tmp_41*tmp_88 + tmp_42*tmp_89;
      real_t tmp_95 = tmp_94 - 1.0/4.0;
      real_t tmp_96 = tmp_0*tmp_91 + tmp_1*tmp_93 + tmp_3*tmp_95;
      real_t tmp_97 = tmp_2*tmp_95 + tmp_4*tmp_93 + tmp_9*tmp_91;
      real_t tmp_98 = tmp_10*tmp_93 + tmp_6*tmp_95 + tmp_8*tmp_91;
      real_t tmp_99 = -tmp_90 - tmp_92 - tmp_94 + 1;
      real_t tmp_100 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_71;
      real_t tmp_101 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_102 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_103 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_104 = tmp_101*tmp_5 + tmp_102*tmp_21 + tmp_103*tmp_27;
      real_t tmp_105 = tmp_104 - 1.0/4.0;
      real_t tmp_106 = tmp_101*tmp_35 + tmp_102*tmp_36 + tmp_103*tmp_37;
      real_t tmp_107 = tmp_106 - 1.0/4.0;
      real_t tmp_108 = tmp_101*tmp_40 + tmp_102*tmp_41 + tmp_103*tmp_42;
      real_t tmp_109 = tmp_108 - 1.0/4.0;
      real_t tmp_110 = tmp_0*tmp_105 + tmp_1*tmp_107 + tmp_109*tmp_3;
      real_t tmp_111 = tmp_105*tmp_9 + tmp_107*tmp_4 + tmp_109*tmp_2;
      real_t tmp_112 = tmp_10*tmp_107 + tmp_105*tmp_8 + tmp_109*tmp_6;
      real_t tmp_113 = -tmp_104 - tmp_106 - tmp_108 + 1;
      real_t tmp_114 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_71;
      real_t tmp_115 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_116 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_117 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_118 = tmp_115*tmp_5 + tmp_116*tmp_21 + tmp_117*tmp_27;
      real_t tmp_119 = tmp_118 - 1.0/4.0;
      real_t tmp_120 = tmp_115*tmp_35 + tmp_116*tmp_36 + tmp_117*tmp_37;
      real_t tmp_121 = tmp_120 - 1.0/4.0;
      real_t tmp_122 = tmp_115*tmp_40 + tmp_116*tmp_41 + tmp_117*tmp_42;
      real_t tmp_123 = tmp_122 - 1.0/4.0;
      real_t tmp_124 = tmp_0*tmp_119 + tmp_1*tmp_121 + tmp_123*tmp_3;
      real_t tmp_125 = tmp_119*tmp_9 + tmp_121*tmp_4 + tmp_123*tmp_2;
      real_t tmp_126 = tmp_10*tmp_121 + tmp_119*tmp_8 + tmp_123*tmp_6;
      real_t tmp_127 = -tmp_118 - tmp_120 - tmp_122 + 1;
      real_t tmp_128 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_71;
      real_t tmp_129 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_130 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_131 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_132 = tmp_129*tmp_5 + tmp_130*tmp_21 + tmp_131*tmp_27;
      real_t tmp_133 = tmp_132 - 1.0/4.0;
      real_t tmp_134 = tmp_129*tmp_35 + tmp_130*tmp_36 + tmp_131*tmp_37;
      real_t tmp_135 = tmp_134 - 1.0/4.0;
      real_t tmp_136 = tmp_129*tmp_40 + tmp_130*tmp_41 + tmp_131*tmp_42;
      real_t tmp_137 = tmp_136 - 1.0/4.0;
      real_t tmp_138 = tmp_0*tmp_133 + tmp_1*tmp_135 + tmp_137*tmp_3;
      real_t tmp_139 = tmp_133*tmp_9 + tmp_135*tmp_4 + tmp_137*tmp_2;
      real_t tmp_140 = tmp_10*tmp_135 + tmp_133*tmp_8 + tmp_137*tmp_6;
      real_t tmp_141 = -tmp_132 - tmp_134 - tmp_136 + 1;
      real_t tmp_142 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_71;
      real_t tmp_143 = p_affine_13_2*tmp_55;
      real_t tmp_144 = p_affine_13_2*tmp_49;
      real_t tmp_145 = p_affine_13_0*tmp_49 + p_affine_13_1*tmp_55 + p_affine_13_2*tmp_61;
      real_t tmp_146 = p_affine_13_2*tmp_54;
      real_t tmp_147 = p_affine_13_2*tmp_48;
      real_t tmp_148 = p_affine_13_0*tmp_48 + p_affine_13_1*tmp_54 + p_affine_13_2*tmp_60;
      real_t tmp_149 = p_affine_13_2*tmp_53;
      real_t tmp_150 = p_affine_13_2*tmp_47;
      real_t tmp_151 = p_affine_13_0*tmp_47 + p_affine_13_1*tmp_53 + p_affine_13_2*tmp_59;
      real_t a_0_0 = tmp_100*(-tmp_51*tmp_96 - tmp_57*tmp_97 - tmp_62*tmp_98 + 1.0*tmp_65*tmp_98*tmp_99 - tmp_70*tmp_99) + tmp_114*(-tmp_110*tmp_51 - tmp_111*tmp_57 + 1.0*tmp_112*tmp_113*tmp_65 - tmp_112*tmp_62 - tmp_113*tmp_70) + tmp_128*(-tmp_124*tmp_51 - tmp_125*tmp_57 + 1.0*tmp_126*tmp_127*tmp_65 - tmp_126*tmp_62 - tmp_127*tmp_70) + tmp_142*(-tmp_138*tmp_51 - tmp_139*tmp_57 + 1.0*tmp_140*tmp_141*tmp_65 - tmp_140*tmp_62 - tmp_141*tmp_70) + tmp_72*(-tmp_45*tmp_51 - tmp_52*tmp_57 - tmp_62*tmp_63 + 1.0*tmp_63*tmp_65*tmp_66 - tmp_66*tmp_70) + tmp_86*(-tmp_51*tmp_82 - tmp_57*tmp_83 - tmp_62*tmp_84 + 1.0*tmp_65*tmp_84*tmp_85 - tmp_70*tmp_85);
      real_t a_1_0 = tmp_100*(-tmp_143*tmp_97 - tmp_144*tmp_96 - tmp_145*tmp_98 + 1.0*tmp_65*tmp_90*tmp_98 - tmp_70*tmp_90) + tmp_114*(1.0*tmp_104*tmp_112*tmp_65 - tmp_104*tmp_70 - tmp_110*tmp_144 - tmp_111*tmp_143 - tmp_112*tmp_145) + tmp_128*(1.0*tmp_118*tmp_126*tmp_65 - tmp_118*tmp_70 - tmp_124*tmp_144 - tmp_125*tmp_143 - tmp_126*tmp_145) + tmp_142*(1.0*tmp_132*tmp_140*tmp_65 - tmp_132*tmp_70 - tmp_138*tmp_144 - tmp_139*tmp_143 - tmp_140*tmp_145) + tmp_72*(-tmp_143*tmp_52 - tmp_144*tmp_45 - tmp_145*tmp_63 + 1.0*tmp_33*tmp_63*tmp_65 - tmp_33*tmp_70) + tmp_86*(-tmp_143*tmp_83 - tmp_144*tmp_82 - tmp_145*tmp_84 + 1.0*tmp_65*tmp_76*tmp_84 - tmp_70*tmp_76);
      real_t a_2_0 = tmp_100*(-tmp_146*tmp_97 - tmp_147*tmp_96 - tmp_148*tmp_98 + 1.0*tmp_65*tmp_92*tmp_98 - tmp_70*tmp_92) + tmp_114*(1.0*tmp_106*tmp_112*tmp_65 - tmp_106*tmp_70 - tmp_110*tmp_147 - tmp_111*tmp_146 - tmp_112*tmp_148) + tmp_128*(1.0*tmp_120*tmp_126*tmp_65 - tmp_120*tmp_70 - tmp_124*tmp_147 - tmp_125*tmp_146 - tmp_126*tmp_148) + tmp_142*(1.0*tmp_134*tmp_140*tmp_65 - tmp_134*tmp_70 - tmp_138*tmp_147 - tmp_139*tmp_146 - tmp_140*tmp_148) + tmp_72*(-tmp_146*tmp_52 - tmp_147*tmp_45 - tmp_148*tmp_63 + 1.0*tmp_38*tmp_63*tmp_65 - tmp_38*tmp_70) + tmp_86*(-tmp_146*tmp_83 - tmp_147*tmp_82 - tmp_148*tmp_84 + 1.0*tmp_65*tmp_78*tmp_84 - tmp_70*tmp_78);
      real_t a_3_0 = tmp_100*(-tmp_149*tmp_97 - tmp_150*tmp_96 - tmp_151*tmp_98 + 1.0*tmp_65*tmp_94*tmp_98 - tmp_70*tmp_94) + tmp_114*(1.0*tmp_108*tmp_112*tmp_65 - tmp_108*tmp_70 - tmp_110*tmp_150 - tmp_111*tmp_149 - tmp_112*tmp_151) + tmp_128*(1.0*tmp_122*tmp_126*tmp_65 - tmp_122*tmp_70 - tmp_124*tmp_150 - tmp_125*tmp_149 - tmp_126*tmp_151) + tmp_142*(1.0*tmp_136*tmp_140*tmp_65 - tmp_136*tmp_70 - tmp_138*tmp_150 - tmp_139*tmp_149 - tmp_140*tmp_151) + tmp_72*(-tmp_149*tmp_52 - tmp_150*tmp_45 - tmp_151*tmp_63 + 1.0*tmp_43*tmp_63*tmp_65 - tmp_43*tmp_70) + tmp_86*(-tmp_149*tmp_83 - tmp_150*tmp_82 - tmp_151*tmp_84 + 1.0*tmp_65*tmp_80*tmp_84 - tmp_70*tmp_80);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }

};




class EGEpsilonFormEE : public hyteg::dg::DGForm
{

 public:
    EGEpsilonFormEE(std::function< real_t ( const Point3D & ) > mu)
: callback_Scalar_Variable_Coefficient_3D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_mu (mu)
, callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g2 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_mu (mu)
    {}

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g2;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_mu;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_g0( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g0( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g2( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g2( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g1( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g1( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
}

      
 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      real_t tmp_0 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_5 = -tmp_4;
      real_t tmp_6 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_7 = -tmp_6;
      real_t tmp_8 = 1.0 / (tmp_3 - tmp_5*tmp_7);
      real_t tmp_9 = 1.0*tmp_8;
      real_t tmp_10 = tmp_3*tmp_9;
      real_t tmp_11 = ((tmp_10 + tmp_5*tmp_6*tmp_9)*(tmp_10 + tmp_5*tmp_6*tmp_9));
      real_t tmp_12 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_13 = ((tmp_10 + tmp_4*tmp_7*tmp_9)*(tmp_10 + tmp_4*tmp_7*tmp_9));
      real_t tmp_14 = tmp_1*tmp_8;
      real_t tmp_15 = tmp_2*tmp_8;
      real_t tmp_16 = 1.0*((tmp_14*tmp_4 + tmp_14*tmp_5 + tmp_15*tmp_6 + tmp_15*tmp_7)*(tmp_14*tmp_4 + tmp_14*tmp_5 + tmp_15*tmp_6 + tmp_15*tmp_7));
      real_t tmp_17 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id1;
      real_t tmp_18 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id2;
      real_t tmp_19 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id3;
      real_t tmp_20 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id4;
      real_t tmp_21 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id5;
      real_t a_0_0 = 0.054975871827660928*tmp_0*(Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_16 + tmp_11*tmp_12 + tmp_12*tmp_13) + 0.11169079483900572*tmp_0*(Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_16 + tmp_11*tmp_17 + tmp_13*tmp_17) + 0.054975871827660928*tmp_0*(Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_16 + tmp_11*tmp_18 + tmp_13*tmp_18) + 0.11169079483900572*tmp_0*(Scalar_Variable_Coefficient_2D_mu_out0_id3*tmp_16 + tmp_11*tmp_19 + tmp_13*tmp_19) + 0.054975871827660928*tmp_0*(Scalar_Variable_Coefficient_2D_mu_out0_id4*tmp_16 + tmp_11*tmp_20 + tmp_13*tmp_20) + 0.11169079483900572*tmp_0*(Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_16 + tmp_11*tmp_21 + tmp_13*tmp_21);
      elMat( 0, 0) = a_0_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -tmp_5;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = -tmp_9;
      real_t tmp_11 = 1.0 / (-tmp_10*tmp_8 + tmp_7);
      real_t tmp_12 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_13 = tmp_11*(0.21132486540518713*tmp_1 + tmp_12);
      real_t tmp_14 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_15 = tmp_11*(0.21132486540518713*tmp_0 + tmp_14);
      real_t tmp_16 = tmp_13*tmp_5 + tmp_15*tmp_6 - 1.0/3.0;
      real_t tmp_17 = tmp_13*tmp_4 + tmp_15*tmp_9 - 1.0/3.0;
      real_t tmp_18 = tmp_16*tmp_4 + tmp_17*tmp_8;
      real_t tmp_19 = tmp_10*tmp_16 + tmp_17*tmp_6;
      real_t tmp_20 = 1.0*tmp_11;
      real_t tmp_21 = tmp_20*tmp_7;
      real_t tmp_22 = 0.5*tmp_11;
      real_t tmp_23 = tmp_22*tmp_4;
      real_t tmp_24 = tmp_22*tmp_6;
      real_t tmp_25 = tmp_10*tmp_24 + tmp_23*tmp_5 + tmp_23*tmp_8 + tmp_24*tmp_9;
      real_t tmp_26 = 1.0*p_affine_10_0*(tmp_20*tmp_8*tmp_9 + tmp_21) + 1.0*p_affine_10_1*tmp_25;
      real_t tmp_27 = 1.0*p_affine_10_0*tmp_25 + 1.0*p_affine_10_1*(tmp_10*tmp_20*tmp_5 + tmp_21);
      real_t tmp_28 = 2*tmp_2;
      real_t tmp_29 = tmp_11*(0.78867513459481287*tmp_1 + tmp_12);
      real_t tmp_30 = tmp_11*(0.78867513459481287*tmp_0 + tmp_14);
      real_t tmp_31 = tmp_29*tmp_5 + tmp_30*tmp_6 - 1.0/3.0;
      real_t tmp_32 = tmp_29*tmp_4 + tmp_30*tmp_9 - 1.0/3.0;
      real_t tmp_33 = tmp_31*tmp_4 + tmp_32*tmp_8;
      real_t tmp_34 = tmp_10*tmp_31 + tmp_32*tmp_6;
      real_t a_0_0 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_28*(-tmp_18*tmp_26 - tmp_19*tmp_27 + tmp_3*((tmp_18*tmp_18) + (tmp_19*tmp_19))) + 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_28*(-tmp_26*tmp_33 - tmp_27*tmp_34 + tmp_3*((tmp_33*tmp_33) + (tmp_34*tmp_34)));
      elMat( 0, 0) = a_0_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElementInner[0]( 0 );
      const auto p_affine_0_1 = coordsElementInner[0]( 1 );

      const auto p_affine_1_0 = coordsElementInner[1]( 0 );
      const auto p_affine_1_1 = coordsElementInner[1]( 1 );

      const auto p_affine_2_0 = coordsElementInner[2]( 0 );
      const auto p_affine_2_1 = coordsElementInner[2]( 1 );

      const auto p_affine_3_0 = coordsElementOuter[0]( 0 );
      const auto p_affine_3_1 = coordsElementOuter[0]( 1 );

      const auto p_affine_4_0 = coordsElementOuter[1]( 0 );
      const auto p_affine_4_1 = coordsElementOuter[1]( 1 );

      const auto p_affine_5_0 = coordsElementOuter[2]( 0 );
      const auto p_affine_5_1 = coordsElementOuter[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_1 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_2 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_3 = tmp_0*tmp_2;
      real_t tmp_4 = -tmp_1;
      real_t tmp_5 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (tmp_3 - tmp_4*tmp_6);
      real_t tmp_8 = -p_affine_3_1;
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.21132486540518713*tmp_9;
      real_t tmp_11 = tmp_7*(tmp_10 + tmp_8);
      real_t tmp_12 = -p_affine_3_0;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + 0.21132486540518713*tmp_13;
      real_t tmp_15 = tmp_7*(tmp_12 + tmp_14);
      real_t tmp_16 = tmp_1*tmp_11 + tmp_15*tmp_2 - 1.0/3.0;
      real_t tmp_17 = tmp_0*tmp_11 + tmp_15*tmp_5 - 1.0/3.0;
      real_t tmp_18 = tmp_0*tmp_16 + tmp_17*tmp_4;
      real_t tmp_19 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_20 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_25 = -tmp_24;
      real_t tmp_26 = 1.0 / (tmp_21 - tmp_23*tmp_25);
      real_t tmp_27 = 1.0*tmp_26;
      real_t tmp_28 = tmp_21*tmp_27;
      real_t tmp_29 = 0.5*tmp_26;
      real_t tmp_30 = tmp_19*tmp_29;
      real_t tmp_31 = tmp_20*tmp_29;
      real_t tmp_32 = tmp_22*tmp_30 + tmp_23*tmp_30 + tmp_24*tmp_31 + tmp_25*tmp_31;
      real_t tmp_33 = p_affine_10_0*(tmp_23*tmp_24*tmp_27 + tmp_28) + p_affine_10_1*tmp_32;
      real_t tmp_34 = -p_affine_0_1;
      real_t tmp_35 = tmp_26*(tmp_10 + tmp_34);
      real_t tmp_36 = -p_affine_0_0;
      real_t tmp_37 = tmp_26*(tmp_14 + tmp_36);
      real_t tmp_38 = tmp_20*tmp_37 + tmp_22*tmp_35 - 1.0/3.0;
      real_t tmp_39 = tmp_19*tmp_35 + tmp_24*tmp_37 - 1.0/3.0;
      real_t tmp_40 = tmp_19*tmp_38 + tmp_23*tmp_39;
      real_t tmp_41 = 1.0*tmp_7;
      real_t tmp_42 = tmp_3*tmp_41;
      real_t tmp_43 = 0.5*tmp_7;
      real_t tmp_44 = tmp_0*tmp_43;
      real_t tmp_45 = tmp_2*tmp_43;
      real_t tmp_46 = tmp_1*tmp_44 + tmp_4*tmp_44 + tmp_45*tmp_5 + tmp_45*tmp_6;
      real_t tmp_47 = 0.5*p_affine_10_0*(tmp_4*tmp_41*tmp_5 + tmp_42) + 0.5*p_affine_10_1*tmp_46;
      real_t tmp_48 = tmp_16*tmp_6 + tmp_17*tmp_2;
      real_t tmp_49 = p_affine_10_0*tmp_32 + p_affine_10_1*(tmp_22*tmp_25*tmp_27 + tmp_28);
      real_t tmp_50 = tmp_20*tmp_39 + tmp_25*tmp_38;
      real_t tmp_51 = 0.5*p_affine_10_0*tmp_46 + 0.5*p_affine_10_1*(tmp_1*tmp_41*tmp_6 + tmp_42);
      real_t tmp_52 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_53 = 1.0 / (tmp_52);
      real_t tmp_54 = 2*tmp_52;
      real_t tmp_55 = p_affine_6_1 + 0.78867513459481287*tmp_9;
      real_t tmp_56 = tmp_7*(tmp_55 + tmp_8);
      real_t tmp_57 = p_affine_6_0 + 0.78867513459481287*tmp_13;
      real_t tmp_58 = tmp_7*(tmp_12 + tmp_57);
      real_t tmp_59 = tmp_1*tmp_56 + tmp_2*tmp_58 - 1.0/3.0;
      real_t tmp_60 = tmp_0*tmp_56 + tmp_5*tmp_58 - 1.0/3.0;
      real_t tmp_61 = tmp_0*tmp_59 + tmp_4*tmp_60;
      real_t tmp_62 = tmp_26*(tmp_34 + tmp_55);
      real_t tmp_63 = tmp_26*(tmp_36 + tmp_57);
      real_t tmp_64 = tmp_20*tmp_63 + tmp_22*tmp_62 - 1.0/3.0;
      real_t tmp_65 = tmp_19*tmp_62 + tmp_24*tmp_63 - 1.0/3.0;
      real_t tmp_66 = tmp_19*tmp_64 + tmp_23*tmp_65;
      real_t tmp_67 = tmp_2*tmp_60 + tmp_59*tmp_6;
      real_t tmp_68 = tmp_20*tmp_65 + tmp_25*tmp_64;
      real_t a_0_0 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_54*(0.5*tmp_18*tmp_33 - tmp_40*tmp_47 + 0.5*tmp_48*tmp_49 - tmp_50*tmp_51 - tmp_53*(tmp_18*tmp_40 + tmp_48*tmp_50)) + 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_54*(0.5*tmp_33*tmp_61 - tmp_47*tmp_66 + 0.5*tmp_49*tmp_67 - tmp_51*tmp_68 - tmp_53*(tmp_61*tmp_66 + tmp_67*tmp_68));
      elMat( 0, 0) = a_0_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -tmp_5;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = -tmp_9;
      real_t tmp_11 = 1.0 / (-tmp_10*tmp_8 + tmp_7);
      real_t tmp_12 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_13 = tmp_11*(0.21132486540518713*tmp_1 + tmp_12);
      real_t tmp_14 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_15 = tmp_11*(0.21132486540518713*tmp_0 + tmp_14);
      real_t tmp_16 = tmp_13*tmp_5 + tmp_15*tmp_6 - 1.0/3.0;
      real_t tmp_17 = tmp_13*tmp_4 + tmp_15*tmp_9 - 1.0/3.0;
      real_t tmp_18 = tmp_16*tmp_4 + tmp_17*tmp_8;
      real_t tmp_19 = tmp_10*tmp_16 + tmp_17*tmp_6;
      real_t tmp_20 = 1.0*tmp_11;
      real_t tmp_21 = tmp_20*tmp_7;
      real_t tmp_22 = 0.5*tmp_11;
      real_t tmp_23 = tmp_22*tmp_4;
      real_t tmp_24 = tmp_22*tmp_6;
      real_t tmp_25 = tmp_10*tmp_24 + tmp_23*tmp_5 + tmp_23*tmp_8 + tmp_24*tmp_9;
      real_t tmp_26 = 2*p_affine_10_0*(tmp_20*tmp_8*tmp_9 + tmp_21) + 2*p_affine_10_1*tmp_25;
      real_t tmp_27 = 2*p_affine_10_0*tmp_25 + 2*p_affine_10_1*(tmp_10*tmp_20*tmp_5 + tmp_21);
      real_t tmp_28 = 2*tmp_2;
      real_t tmp_29 = tmp_11*(0.78867513459481287*tmp_1 + tmp_12);
      real_t tmp_30 = tmp_11*(0.78867513459481287*tmp_0 + tmp_14);
      real_t tmp_31 = tmp_29*tmp_5 + tmp_30*tmp_6 - 1.0/3.0;
      real_t tmp_32 = tmp_29*tmp_4 + tmp_30*tmp_9 - 1.0/3.0;
      real_t tmp_33 = tmp_31*tmp_4 + tmp_32*tmp_8;
      real_t tmp_34 = tmp_10*tmp_31 + tmp_32*tmp_6;
      real_t a_0_0 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_28*(-tmp_18*tmp_26 - tmp_19*tmp_27 + tmp_3*((tmp_18*tmp_18) + (tmp_19*tmp_19))) + 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id1*tmp_28*(-tmp_26*tmp_33 - tmp_27*tmp_34 + tmp_3*((tmp_33*tmp_33) + (tmp_34*tmp_34)));
      elMat( 0, 0) = a_0_0;
   }

  void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_g0( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id0 );
      Scalar_Variable_Coefficient_2D_g1( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_g0( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id3 );
      Scalar_Variable_Coefficient_2D_g1( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -tmp_5;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = -tmp_9;
      real_t tmp_11 = 1.0 / (-tmp_10*tmp_8 + tmp_7);
      real_t tmp_12 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_13 = tmp_11*(0.21132486540518713*tmp_1 + tmp_12);
      real_t tmp_14 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_15 = tmp_11*(0.21132486540518713*tmp_0 + tmp_14);
      real_t tmp_16 = tmp_13*tmp_5 + tmp_15*tmp_6 - 1.0/3.0;
      real_t tmp_17 = tmp_13*tmp_4 + tmp_15*tmp_9 - 1.0/3.0;
      real_t tmp_18 = 1.0*tmp_11;
      real_t tmp_19 = tmp_18*tmp_7;
      real_t tmp_20 = 0.5*tmp_11;
      real_t tmp_21 = tmp_20*tmp_4;
      real_t tmp_22 = tmp_20*tmp_6;
      real_t tmp_23 = tmp_10*tmp_22 + tmp_21*tmp_5 + tmp_21*tmp_8 + tmp_22*tmp_9;
      real_t tmp_24 = p_affine_10_0*(tmp_18*tmp_8*tmp_9 + tmp_19) + p_affine_10_1*tmp_23;
      real_t tmp_25 = p_affine_10_0*tmp_23 + p_affine_10_1*(tmp_10*tmp_18*tmp_5 + tmp_19);
      real_t tmp_26 = 2*tmp_2;
      real_t tmp_27 = tmp_11*(0.78867513459481287*tmp_1 + tmp_12);
      real_t tmp_28 = tmp_11*(0.78867513459481287*tmp_0 + tmp_14);
      real_t tmp_29 = tmp_27*tmp_5 + tmp_28*tmp_6 - 1.0/3.0;
      real_t tmp_30 = tmp_27*tmp_4 + tmp_28*tmp_9 - 1.0/3.0;
      real_t a_0_0 = 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id2*tmp_26*(Scalar_Variable_Coefficient_2D_g0_out0_id0*(-tmp_24 + tmp_3*(tmp_16*tmp_4 + tmp_17*tmp_8)) + Scalar_Variable_Coefficient_2D_g1_out0_id1*(-tmp_25 + tmp_3*(tmp_10*tmp_16 + tmp_17*tmp_6))) + 0.5*Scalar_Variable_Coefficient_2D_mu_out0_id5*tmp_26*(Scalar_Variable_Coefficient_2D_g0_out0_id3*(-tmp_24 + tmp_3*(tmp_29*tmp_4 + tmp_30*tmp_8)) + Scalar_Variable_Coefficient_2D_g1_out0_id4*(-tmp_25 + tmp_3*(tmp_10*tmp_29 + tmp_30*tmp_6)));
      elMat( 0, 0) = a_0_0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      real_t Scalar_Variable_Coefficient_3D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id13 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id14 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id15 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id16 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id17 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id18 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id19 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id20 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id21 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id22 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id23 = 0;
      Scalar_Variable_Coefficient_3D_g0( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id0 );
      Scalar_Variable_Coefficient_3D_g2( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id1 );
      Scalar_Variable_Coefficient_3D_g1( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_g0( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id4 );
      Scalar_Variable_Coefficient_3D_g2( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id5 );
      Scalar_Variable_Coefficient_3D_g1( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id6 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id7 );
      Scalar_Variable_Coefficient_3D_g0( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id8 );
      Scalar_Variable_Coefficient_3D_g2( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id9 );
      Scalar_Variable_Coefficient_3D_g1( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id10 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id11 );
      Scalar_Variable_Coefficient_3D_g0( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id12 );
      Scalar_Variable_Coefficient_3D_g2( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id13 );
      Scalar_Variable_Coefficient_3D_g1( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id14 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id15 );
      Scalar_Variable_Coefficient_3D_g0( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id16 );
      Scalar_Variable_Coefficient_3D_g2( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id17 );
      Scalar_Variable_Coefficient_3D_g1( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id18 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id19 );
      Scalar_Variable_Coefficient_3D_g0( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id20 );
      Scalar_Variable_Coefficient_3D_g2( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id21 );
      Scalar_Variable_Coefficient_3D_g1( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id22 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id23 );
      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_13 = tmp_10*tmp_9 - tmp_11*tmp_12;
      real_t tmp_14 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_15 = tmp_12*tmp_14;
      real_t tmp_16 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_17 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_18 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_19 = tmp_17*tmp_18;
      real_t tmp_20 = tmp_10*tmp_18;
      real_t tmp_21 = tmp_14*tmp_17;
      real_t tmp_22 = tmp_12*tmp_16;
      real_t tmp_23 = 1.0 / (tmp_10*tmp_16*tmp_9 + tmp_11*tmp_19 - tmp_11*tmp_22 + tmp_15*tmp_8 - tmp_20*tmp_8 - tmp_21*tmp_9);
      real_t tmp_24 = -tmp_4;
      real_t tmp_25 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_26 = tmp_23*(0.81684757298045851*tmp_24 + tmp_25 + 0.091576213509770743*tmp_5);
      real_t tmp_27 = tmp_11*tmp_18 - tmp_14*tmp_9;
      real_t tmp_28 = -tmp_1;
      real_t tmp_29 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_30 = tmp_23*(0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_28 + tmp_29);
      real_t tmp_31 = tmp_15 - tmp_20;
      real_t tmp_32 = -tmp_3;
      real_t tmp_33 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_34 = tmp_23*(0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_32 + tmp_33);
      real_t tmp_35 = tmp_13*tmp_26 + tmp_27*tmp_30 + tmp_31*tmp_34 - 1.0/4.0;
      real_t tmp_36 = -tmp_10*tmp_8 + tmp_11*tmp_17;
      real_t tmp_37 = -tmp_11*tmp_16 + tmp_14*tmp_8;
      real_t tmp_38 = tmp_10*tmp_16 - tmp_21;
      real_t tmp_39 = tmp_26*tmp_36 + tmp_30*tmp_37 + tmp_34*tmp_38 - 1.0/4.0;
      real_t tmp_40 = tmp_12*tmp_8 - tmp_17*tmp_9;
      real_t tmp_41 = tmp_16*tmp_9 - tmp_18*tmp_8;
      real_t tmp_42 = tmp_19 - tmp_22;
      real_t tmp_43 = tmp_26*tmp_40 + tmp_30*tmp_41 + tmp_34*tmp_42 - 1.0/4.0;
      real_t tmp_44 = 1.0*tmp_23;
      real_t tmp_45 = 0.5*tmp_23;
      real_t tmp_46 = tmp_45*tmp_8;
      real_t tmp_47 = tmp_45*tmp_9;
      real_t tmp_48 = tmp_11*tmp_45;
      real_t tmp_49 = tmp_31*tmp_45;
      real_t tmp_50 = tmp_38*tmp_45;
      real_t tmp_51 = tmp_42*tmp_45;
      real_t tmp_52 = tmp_10*tmp_51 + tmp_12*tmp_50 + tmp_17*tmp_49 + tmp_27*tmp_46 + tmp_37*tmp_47 + tmp_41*tmp_48;
      real_t tmp_53 = tmp_13*tmp_46 + tmp_14*tmp_51 + tmp_16*tmp_49 + tmp_18*tmp_50 + tmp_36*tmp_47 + tmp_40*tmp_48;
      real_t tmp_54 = p_affine_13_0*(tmp_11*tmp_42*tmp_44 + tmp_31*tmp_44*tmp_8 + tmp_38*tmp_44*tmp_9) + p_affine_13_1*tmp_52 + p_affine_13_2*tmp_53;
      real_t tmp_55 = tmp_10*tmp_40*tmp_45 + tmp_12*tmp_36*tmp_45 + tmp_13*tmp_17*tmp_45 + tmp_14*tmp_41*tmp_45 + tmp_16*tmp_27*tmp_45 + tmp_18*tmp_37*tmp_45;
      real_t tmp_56 = p_affine_13_0*tmp_52 + p_affine_13_1*(tmp_10*tmp_41*tmp_44 + tmp_12*tmp_37*tmp_44 + tmp_17*tmp_27*tmp_44) + p_affine_13_2*tmp_55;
      real_t tmp_57 = p_affine_13_0*tmp_53 + p_affine_13_1*tmp_55 + p_affine_13_2*(tmp_13*tmp_16*tmp_44 + tmp_14*tmp_40*tmp_44 + tmp_18*tmp_36*tmp_44);
      real_t tmp_58 = 2.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_59 = tmp_23*(0.10810301816807022*tmp_24 + tmp_25 + 0.44594849091596489*tmp_5);
      real_t tmp_60 = tmp_23*(0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_28 + tmp_29);
      real_t tmp_61 = tmp_23*(0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_32 + tmp_33);
      real_t tmp_62 = tmp_13*tmp_59 + tmp_27*tmp_60 + tmp_31*tmp_61 - 1.0/4.0;
      real_t tmp_63 = tmp_36*tmp_59 + tmp_37*tmp_60 + tmp_38*tmp_61 - 1.0/4.0;
      real_t tmp_64 = tmp_40*tmp_59 + tmp_41*tmp_60 + tmp_42*tmp_61 - 1.0/4.0;
      real_t tmp_65 = tmp_23*(0.091576213509770743*tmp_24 + tmp_25 + 0.091576213509770743*tmp_5);
      real_t tmp_66 = tmp_23*(0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_28 + tmp_29);
      real_t tmp_67 = tmp_23*(0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_32 + tmp_33);
      real_t tmp_68 = tmp_13*tmp_65 + tmp_27*tmp_66 + tmp_31*tmp_67 - 1.0/4.0;
      real_t tmp_69 = tmp_36*tmp_65 + tmp_37*tmp_66 + tmp_38*tmp_67 - 1.0/4.0;
      real_t tmp_70 = tmp_40*tmp_65 + tmp_41*tmp_66 + tmp_42*tmp_67 - 1.0/4.0;
      real_t tmp_71 = tmp_23*(0.44594849091596489*tmp_24 + tmp_25 + 0.44594849091596489*tmp_5);
      real_t tmp_72 = tmp_23*(0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_28 + tmp_29);
      real_t tmp_73 = tmp_23*(0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_32 + tmp_33);
      real_t tmp_74 = tmp_13*tmp_71 + tmp_27*tmp_72 + tmp_31*tmp_73 - 1.0/4.0;
      real_t tmp_75 = tmp_36*tmp_71 + tmp_37*tmp_72 + tmp_38*tmp_73 - 1.0/4.0;
      real_t tmp_76 = tmp_40*tmp_71 + tmp_41*tmp_72 + tmp_42*tmp_73 - 1.0/4.0;
      real_t tmp_77 = tmp_23*(0.091576213509770743*tmp_24 + tmp_25 + 0.81684757298045851*tmp_5);
      real_t tmp_78 = tmp_23*(0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_28 + tmp_29);
      real_t tmp_79 = tmp_23*(0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_32 + tmp_33);
      real_t tmp_80 = tmp_13*tmp_77 + tmp_27*tmp_78 + tmp_31*tmp_79 - 1.0/4.0;
      real_t tmp_81 = tmp_36*tmp_77 + tmp_37*tmp_78 + tmp_38*tmp_79 - 1.0/4.0;
      real_t tmp_82 = tmp_40*tmp_77 + tmp_41*tmp_78 + tmp_42*tmp_79 - 1.0/4.0;
      real_t tmp_83 = tmp_23*(0.44594849091596489*tmp_24 + tmp_25 + 0.10810301816807022*tmp_5);
      real_t tmp_84 = tmp_23*(0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_28 + tmp_29);
      real_t tmp_85 = tmp_23*(0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_32 + tmp_33);
      real_t tmp_86 = tmp_13*tmp_83 + tmp_27*tmp_84 + tmp_31*tmp_85 - 1.0/4.0;
      real_t tmp_87 = tmp_36*tmp_83 + tmp_37*tmp_84 + tmp_38*tmp_85 - 1.0/4.0;
      real_t tmp_88 = tmp_40*tmp_83 + tmp_41*tmp_84 + tmp_42*tmp_85 - 1.0/4.0;
      real_t a_0_0 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id11*tmp_58*(Scalar_Variable_Coefficient_3D_g0_out0_id8*(-tmp_54 + 1.0*tmp_7*(tmp_11*tmp_43 + tmp_35*tmp_8 + tmp_39*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id10*(-tmp_56 + 1.0*tmp_7*(tmp_10*tmp_43 + tmp_12*tmp_39 + tmp_17*tmp_35)) + Scalar_Variable_Coefficient_3D_g2_out0_id9*(-tmp_57 + 1.0*tmp_7*(tmp_14*tmp_43 + tmp_16*tmp_35 + tmp_18*tmp_39))) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id15*tmp_58*(Scalar_Variable_Coefficient_3D_g0_out0_id12*(-tmp_54 + 1.0*tmp_7*(tmp_11*tmp_64 + tmp_62*tmp_8 + tmp_63*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id14*(-tmp_56 + 1.0*tmp_7*(tmp_10*tmp_64 + tmp_12*tmp_63 + tmp_17*tmp_62)) + Scalar_Variable_Coefficient_3D_g2_out0_id13*(-tmp_57 + 1.0*tmp_7*(tmp_14*tmp_64 + tmp_16*tmp_62 + tmp_18*tmp_63))) + 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id19*tmp_58*(Scalar_Variable_Coefficient_3D_g0_out0_id16*(-tmp_54 + 1.0*tmp_7*(tmp_11*tmp_70 + tmp_68*tmp_8 + tmp_69*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id18*(-tmp_56 + 1.0*tmp_7*(tmp_10*tmp_70 + tmp_12*tmp_69 + tmp_17*tmp_68)) + Scalar_Variable_Coefficient_3D_g2_out0_id17*(-tmp_57 + 1.0*tmp_7*(tmp_14*tmp_70 + tmp_16*tmp_68 + tmp_18*tmp_69))) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id23*tmp_58*(Scalar_Variable_Coefficient_3D_g0_out0_id20*(-tmp_54 + 1.0*tmp_7*(tmp_11*tmp_76 + tmp_74*tmp_8 + tmp_75*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id22*(-tmp_56 + 1.0*tmp_7*(tmp_10*tmp_76 + tmp_12*tmp_75 + tmp_17*tmp_74)) + Scalar_Variable_Coefficient_3D_g2_out0_id21*(-tmp_57 + 1.0*tmp_7*(tmp_14*tmp_76 + tmp_16*tmp_74 + tmp_18*tmp_75))) + 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_58*(Scalar_Variable_Coefficient_3D_g0_out0_id0*(-tmp_54 + 1.0*tmp_7*(tmp_11*tmp_82 + tmp_8*tmp_80 + tmp_81*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id2*(-tmp_56 + 1.0*tmp_7*(tmp_10*tmp_82 + tmp_12*tmp_81 + tmp_17*tmp_80)) + Scalar_Variable_Coefficient_3D_g2_out0_id1*(-tmp_57 + 1.0*tmp_7*(tmp_14*tmp_82 + tmp_16*tmp_80 + tmp_18*tmp_81))) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id7*tmp_58*(Scalar_Variable_Coefficient_3D_g0_out0_id4*(-tmp_54 + 1.0*tmp_7*(tmp_11*tmp_88 + tmp_8*tmp_86 + tmp_87*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id6*(-tmp_56 + 1.0*tmp_7*(tmp_10*tmp_88 + tmp_12*tmp_87 + tmp_17*tmp_86)) + Scalar_Variable_Coefficient_3D_g2_out0_id5*(-tmp_57 + 1.0*tmp_7*(tmp_14*tmp_88 + tmp_16*tmp_86 + tmp_18*tmp_87)));
      elMat( 0, 0) = a_0_0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

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

      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id6 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id7 );
      Scalar_Variable_Coefficient_3D_mu( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id8 );
      Scalar_Variable_Coefficient_3D_mu( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id9 );
      Scalar_Variable_Coefficient_3D_mu( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id10 );
      Scalar_Variable_Coefficient_3D_mu( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id11 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id12 );
      Scalar_Variable_Coefficient_3D_mu( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id13 );
      real_t tmp_0 = p_affine_0_0*p_affine_1_1;
      real_t tmp_1 = p_affine_0_0*p_affine_1_2;
      real_t tmp_2 = p_affine_2_1*p_affine_3_2;
      real_t tmp_3 = p_affine_0_1*p_affine_1_0;
      real_t tmp_4 = p_affine_0_1*p_affine_1_2;
      real_t tmp_5 = p_affine_2_2*p_affine_3_0;
      real_t tmp_6 = p_affine_0_2*p_affine_1_0;
      real_t tmp_7 = p_affine_0_2*p_affine_1_1;
      real_t tmp_8 = p_affine_2_0*p_affine_3_1;
      real_t tmp_9 = p_affine_2_2*p_affine_3_1;
      real_t tmp_10 = p_affine_2_0*p_affine_3_2;
      real_t tmp_11 = p_affine_2_1*p_affine_3_0;
      real_t tmp_12 = std::abs(p_affine_0_0*tmp_2 - p_affine_0_0*tmp_9 - p_affine_0_1*tmp_10 + p_affine_0_1*tmp_5 - p_affine_0_2*tmp_11 + p_affine_0_2*tmp_8 - p_affine_1_0*tmp_2 + p_affine_1_0*tmp_9 + p_affine_1_1*tmp_10 - p_affine_1_1*tmp_5 + p_affine_1_2*tmp_11 - p_affine_1_2*tmp_8 + p_affine_2_0*tmp_4 - p_affine_2_0*tmp_7 - p_affine_2_1*tmp_1 + p_affine_2_1*tmp_6 + p_affine_2_2*tmp_0 - p_affine_2_2*tmp_3 - p_affine_3_0*tmp_4 + p_affine_3_0*tmp_7 + p_affine_3_1*tmp_1 - p_affine_3_1*tmp_6 - p_affine_3_2*tmp_0 + p_affine_3_2*tmp_3);
      real_t tmp_13 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_14 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_17 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_18 = tmp_16*tmp_17;
      real_t tmp_19 = tmp_15 - tmp_18;
      real_t tmp_20 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_21 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_22 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_23 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_24 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_25 = tmp_17*tmp_24;
      real_t tmp_26 = tmp_14*tmp_24;
      real_t tmp_27 = tmp_13*tmp_22;
      real_t tmp_28 = 1.0 / (tmp_15*tmp_20 + tmp_16*tmp_21*tmp_22 - tmp_18*tmp_20 - tmp_21*tmp_26 + tmp_23*tmp_25 - tmp_23*tmp_27);
      real_t tmp_29 = tmp_20*tmp_28;
      real_t tmp_30 = tmp_16*tmp_22 - tmp_26;
      real_t tmp_31 = tmp_21*tmp_28;
      real_t tmp_32 = tmp_25 - tmp_27;
      real_t tmp_33 = tmp_23*tmp_28;
      real_t tmp_34 = ((1.0*tmp_19*tmp_29 + 1.0*tmp_30*tmp_31 + 1.0*tmp_32*tmp_33)*(1.0*tmp_19*tmp_29 + 1.0*tmp_30*tmp_31 + 1.0*tmp_32*tmp_33));
      real_t tmp_35 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_36 = -tmp_14*tmp_21 + tmp_17*tmp_23;
      real_t tmp_37 = tmp_24*tmp_28;
      real_t tmp_38 = tmp_14*tmp_20 - tmp_22*tmp_23;
      real_t tmp_39 = tmp_13*tmp_28;
      real_t tmp_40 = -tmp_17*tmp_20 + tmp_21*tmp_22;
      real_t tmp_41 = tmp_16*tmp_28;
      real_t tmp_42 = ((1.0*tmp_36*tmp_37 + 1.0*tmp_38*tmp_39 + 1.0*tmp_40*tmp_41)*(1.0*tmp_36*tmp_37 + 1.0*tmp_38*tmp_39 + 1.0*tmp_40*tmp_41));
      real_t tmp_43 = -tmp_13*tmp_23 + tmp_16*tmp_21;
      real_t tmp_44 = tmp_22*tmp_28;
      real_t tmp_45 = -tmp_16*tmp_20 + tmp_23*tmp_24;
      real_t tmp_46 = tmp_17*tmp_28;
      real_t tmp_47 = tmp_13*tmp_20 - tmp_21*tmp_24;
      real_t tmp_48 = tmp_14*tmp_28;
      real_t tmp_49 = ((1.0*tmp_43*tmp_44 + 1.0*tmp_45*tmp_46 + 1.0*tmp_47*tmp_48)*(1.0*tmp_43*tmp_44 + 1.0*tmp_45*tmp_46 + 1.0*tmp_47*tmp_48));
      real_t tmp_50 = tmp_19*tmp_28;
      real_t tmp_51 = tmp_28*tmp_30;
      real_t tmp_52 = tmp_28*tmp_32;
      real_t tmp_53 = ((tmp_14*tmp_52 + tmp_17*tmp_51 + tmp_22*tmp_50 + tmp_29*tmp_43 + tmp_31*tmp_45 + tmp_33*tmp_47)*(tmp_14*tmp_52 + tmp_17*tmp_51 + tmp_22*tmp_50 + tmp_29*tmp_43 + tmp_31*tmp_45 + tmp_33*tmp_47));
      real_t tmp_54 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_55 = ((tmp_13*tmp_51 + tmp_16*tmp_52 + tmp_24*tmp_50 + tmp_29*tmp_36 + tmp_31*tmp_38 + tmp_33*tmp_40)*(tmp_13*tmp_51 + tmp_16*tmp_52 + tmp_24*tmp_50 + tmp_29*tmp_36 + tmp_31*tmp_38 + tmp_33*tmp_40));
      real_t tmp_56 = ((tmp_36*tmp_44 + tmp_37*tmp_43 + tmp_38*tmp_46 + tmp_39*tmp_45 + tmp_40*tmp_48 + tmp_41*tmp_47)*(tmp_36*tmp_44 + tmp_37*tmp_43 + tmp_38*tmp_46 + tmp_39*tmp_45 + tmp_40*tmp_48 + tmp_41*tmp_47));
      real_t tmp_57 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id1;
      real_t tmp_58 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id1;
      real_t tmp_59 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id2;
      real_t tmp_60 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id2;
      real_t tmp_61 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id3;
      real_t tmp_62 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id3;
      real_t tmp_63 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id4;
      real_t tmp_64 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id4;
      real_t tmp_65 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id5;
      real_t tmp_66 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id5;
      real_t tmp_67 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id6;
      real_t tmp_68 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id6;
      real_t tmp_69 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id7;
      real_t tmp_70 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id7;
      real_t tmp_71 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id8;
      real_t tmp_72 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id8;
      real_t tmp_73 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id9;
      real_t tmp_74 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id9;
      real_t tmp_75 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id10;
      real_t tmp_76 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id10;
      real_t tmp_77 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id11;
      real_t tmp_78 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id11;
      real_t tmp_79 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id12;
      real_t tmp_80 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id12;
      real_t tmp_81 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id13;
      real_t tmp_82 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id13;
      real_t a_0_0 = 0.018781320953002646*tmp_12*(tmp_34*tmp_35 + tmp_35*tmp_42 + tmp_35*tmp_49 + tmp_53*tmp_54 + tmp_54*tmp_55 + tmp_54*tmp_56) + 0.012248840519393657*tmp_12*(tmp_34*tmp_57 + tmp_42*tmp_57 + tmp_49*tmp_57 + tmp_53*tmp_58 + tmp_55*tmp_58 + tmp_56*tmp_58) + 0.0070910034628469103*tmp_12*(tmp_34*tmp_59 + tmp_42*tmp_59 + tmp_49*tmp_59 + tmp_53*tmp_60 + tmp_55*tmp_60 + tmp_56*tmp_60) + 0.0070910034628469103*tmp_12*(tmp_34*tmp_61 + tmp_42*tmp_61 + tmp_49*tmp_61 + tmp_53*tmp_62 + tmp_55*tmp_62 + tmp_56*tmp_62) + 0.0070910034628469103*tmp_12*(tmp_34*tmp_63 + tmp_42*tmp_63 + tmp_49*tmp_63 + tmp_53*tmp_64 + tmp_55*tmp_64 + tmp_56*tmp_64) + 0.0070910034628469103*tmp_12*(tmp_34*tmp_65 + tmp_42*tmp_65 + tmp_49*tmp_65 + tmp_53*tmp_66 + tmp_55*tmp_66 + tmp_56*tmp_66) + 0.018781320953002646*tmp_12*(tmp_34*tmp_67 + tmp_42*tmp_67 + tmp_49*tmp_67 + tmp_53*tmp_68 + tmp_55*tmp_68 + tmp_56*tmp_68) + 0.012248840519393657*tmp_12*(tmp_34*tmp_69 + tmp_42*tmp_69 + tmp_49*tmp_69 + tmp_53*tmp_70 + tmp_55*tmp_70 + tmp_56*tmp_70) + 0.018781320953002646*tmp_12*(tmp_34*tmp_71 + tmp_42*tmp_71 + tmp_49*tmp_71 + tmp_53*tmp_72 + tmp_55*tmp_72 + tmp_56*tmp_72) + 0.012248840519393657*tmp_12*(tmp_34*tmp_73 + tmp_42*tmp_73 + tmp_49*tmp_73 + tmp_53*tmp_74 + tmp_55*tmp_74 + tmp_56*tmp_74) + 0.018781320953002646*tmp_12*(tmp_34*tmp_75 + tmp_42*tmp_75 + tmp_49*tmp_75 + tmp_53*tmp_76 + tmp_55*tmp_76 + tmp_56*tmp_76) + 0.012248840519393657*tmp_12*(tmp_34*tmp_77 + tmp_42*tmp_77 + tmp_49*tmp_77 + tmp_53*tmp_78 + tmp_55*tmp_78 + tmp_56*tmp_78) + 0.0070910034628469103*tmp_12*(tmp_34*tmp_79 + tmp_42*tmp_79 + tmp_49*tmp_79 + tmp_53*tmp_80 + tmp_55*tmp_80 + tmp_56*tmp_80) + 0.0070910034628469103*tmp_12*(tmp_34*tmp_81 + tmp_42*tmp_81 + tmp_49*tmp_81 + tmp_53*tmp_82 + tmp_55*tmp_82 + tmp_56*tmp_82);
      elMat( 0, 0) = a_0_0;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_11 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_12 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_13 = tmp_12*tmp_2;
      real_t tmp_14 = tmp_1*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_13 + tmp_0*tmp_9 + tmp_10*tmp_3 - tmp_10*tmp_6 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_14);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = tmp_12*tmp_4 - tmp_14;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_13 + tmp_9;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_7 + tmp_21*tmp_26 + tmp_27*tmp_32 - 1.0/4.0;
      real_t tmp_34 = -tmp_0*tmp_2 + tmp_11*tmp_4;
      real_t tmp_35 = tmp_0*tmp_8 - tmp_10*tmp_4;
      real_t tmp_36 = tmp_10*tmp_2 - tmp_11*tmp_8;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36 - 1.0/4.0;
      real_t tmp_38 = tmp_0*tmp_5 - tmp_1*tmp_11;
      real_t tmp_39 = -tmp_0*tmp_12 + tmp_1*tmp_10;
      real_t tmp_40 = -tmp_10*tmp_5 + tmp_11*tmp_12;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40 - 1.0/4.0;
      real_t tmp_42 = tmp_0*tmp_33 + tmp_1*tmp_37 + tmp_4*tmp_41;
      real_t tmp_43 = 1.0*tmp_15;
      real_t tmp_44 = 0.5*tmp_15;
      real_t tmp_45 = tmp_0*tmp_44;
      real_t tmp_46 = tmp_1*tmp_44;
      real_t tmp_47 = tmp_4*tmp_44;
      real_t tmp_48 = tmp_27*tmp_44;
      real_t tmp_49 = tmp_36*tmp_44;
      real_t tmp_50 = tmp_40*tmp_44;
      real_t tmp_51 = tmp_11*tmp_48 + tmp_2*tmp_50 + tmp_21*tmp_45 + tmp_35*tmp_46 + tmp_39*tmp_47 + tmp_49*tmp_5;
      real_t tmp_52 = tmp_10*tmp_48 + tmp_12*tmp_49 + tmp_34*tmp_46 + tmp_38*tmp_47 + tmp_45*tmp_7 + tmp_50*tmp_8;
      real_t tmp_53 = 1.0*p_affine_13_0*(tmp_0*tmp_27*tmp_43 + tmp_1*tmp_36*tmp_43 + tmp_4*tmp_40*tmp_43) + 1.0*p_affine_13_1*tmp_51 + 1.0*p_affine_13_2*tmp_52;
      real_t tmp_54 = tmp_10*tmp_33 + tmp_12*tmp_37 + tmp_41*tmp_8;
      real_t tmp_55 = tmp_10*tmp_21*tmp_44 + tmp_11*tmp_44*tmp_7 + tmp_12*tmp_35*tmp_44 + tmp_2*tmp_38*tmp_44 + tmp_34*tmp_44*tmp_5 + tmp_39*tmp_44*tmp_8;
      real_t tmp_56 = 1.0*p_affine_13_0*tmp_52 + 1.0*p_affine_13_1*tmp_55 + 1.0*p_affine_13_2*(tmp_10*tmp_43*tmp_7 + tmp_12*tmp_34*tmp_43 + tmp_38*tmp_43*tmp_8);
      real_t tmp_57 = tmp_11*tmp_33 + tmp_2*tmp_41 + tmp_37*tmp_5;
      real_t tmp_58 = 1.0*p_affine_13_0*tmp_51 + 1.0*p_affine_13_1*(tmp_11*tmp_21*tmp_43 + tmp_2*tmp_39*tmp_43 + tmp_35*tmp_43*tmp_5) + 1.0*p_affine_13_2*tmp_55;
      real_t tmp_59 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_60 = std::pow(tmp_59, -0.25);
      real_t tmp_61 = 2.0*std::pow(tmp_59, 1.0/2.0);
      real_t tmp_62 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_63 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_64 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_65 = tmp_21*tmp_63 + tmp_27*tmp_64 + tmp_62*tmp_7 - 1.0/4.0;
      real_t tmp_66 = tmp_34*tmp_62 + tmp_35*tmp_63 + tmp_36*tmp_64 - 1.0/4.0;
      real_t tmp_67 = tmp_38*tmp_62 + tmp_39*tmp_63 + tmp_40*tmp_64 - 1.0/4.0;
      real_t tmp_68 = tmp_0*tmp_65 + tmp_1*tmp_66 + tmp_4*tmp_67;
      real_t tmp_69 = tmp_10*tmp_65 + tmp_12*tmp_66 + tmp_67*tmp_8;
      real_t tmp_70 = tmp_11*tmp_65 + tmp_2*tmp_67 + tmp_5*tmp_66;
      real_t tmp_71 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_72 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_73 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_74 = tmp_21*tmp_72 + tmp_27*tmp_73 + tmp_7*tmp_71 - 1.0/4.0;
      real_t tmp_75 = tmp_34*tmp_71 + tmp_35*tmp_72 + tmp_36*tmp_73 - 1.0/4.0;
      real_t tmp_76 = tmp_38*tmp_71 + tmp_39*tmp_72 + tmp_40*tmp_73 - 1.0/4.0;
      real_t tmp_77 = tmp_0*tmp_74 + tmp_1*tmp_75 + tmp_4*tmp_76;
      real_t tmp_78 = tmp_10*tmp_74 + tmp_12*tmp_75 + tmp_76*tmp_8;
      real_t tmp_79 = tmp_11*tmp_74 + tmp_2*tmp_76 + tmp_5*tmp_75;
      real_t tmp_80 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_81 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_82 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_83 = tmp_21*tmp_81 + tmp_27*tmp_82 + tmp_7*tmp_80 - 1.0/4.0;
      real_t tmp_84 = tmp_34*tmp_80 + tmp_35*tmp_81 + tmp_36*tmp_82 - 1.0/4.0;
      real_t tmp_85 = tmp_38*tmp_80 + tmp_39*tmp_81 + tmp_40*tmp_82 - 1.0/4.0;
      real_t tmp_86 = tmp_0*tmp_83 + tmp_1*tmp_84 + tmp_4*tmp_85;
      real_t tmp_87 = tmp_10*tmp_83 + tmp_12*tmp_84 + tmp_8*tmp_85;
      real_t tmp_88 = tmp_11*tmp_83 + tmp_2*tmp_85 + tmp_5*tmp_84;
      real_t tmp_89 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_90 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_91 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_92 = tmp_21*tmp_90 + tmp_27*tmp_91 + tmp_7*tmp_89 - 1.0/4.0;
      real_t tmp_93 = tmp_34*tmp_89 + tmp_35*tmp_90 + tmp_36*tmp_91 - 1.0/4.0;
      real_t tmp_94 = tmp_38*tmp_89 + tmp_39*tmp_90 + tmp_40*tmp_91 - 1.0/4.0;
      real_t tmp_95 = tmp_0*tmp_92 + tmp_1*tmp_93 + tmp_4*tmp_94;
      real_t tmp_96 = tmp_10*tmp_92 + tmp_12*tmp_93 + tmp_8*tmp_94;
      real_t tmp_97 = tmp_11*tmp_92 + tmp_2*tmp_94 + tmp_5*tmp_93;
      real_t tmp_98 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_99 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_100 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_101 = tmp_100*tmp_27 + tmp_21*tmp_99 + tmp_7*tmp_98 - 1.0/4.0;
      real_t tmp_102 = tmp_100*tmp_36 + tmp_34*tmp_98 + tmp_35*tmp_99 - 1.0/4.0;
      real_t tmp_103 = tmp_100*tmp_40 + tmp_38*tmp_98 + tmp_39*tmp_99 - 1.0/4.0;
      real_t tmp_104 = tmp_0*tmp_101 + tmp_1*tmp_102 + tmp_103*tmp_4;
      real_t tmp_105 = tmp_10*tmp_101 + tmp_102*tmp_12 + tmp_103*tmp_8;
      real_t tmp_106 = tmp_101*tmp_11 + tmp_102*tmp_5 + tmp_103*tmp_2;
      real_t a_0_0 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_61*(-tmp_42*tmp_53 - tmp_54*tmp_56 - tmp_57*tmp_58 + 1.0*tmp_60*((tmp_42*tmp_42) + (tmp_54*tmp_54) + (tmp_57*tmp_57))) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_61*(-tmp_53*tmp_68 - tmp_56*tmp_69 - tmp_58*tmp_70 + 1.0*tmp_60*((tmp_68*tmp_68) + (tmp_69*tmp_69) + (tmp_70*tmp_70))) + 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_61*(-tmp_53*tmp_77 - tmp_56*tmp_78 - tmp_58*tmp_79 + 1.0*tmp_60*((tmp_77*tmp_77) + (tmp_78*tmp_78) + (tmp_79*tmp_79))) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_61*(-tmp_53*tmp_86 - tmp_56*tmp_87 - tmp_58*tmp_88 + 1.0*tmp_60*((tmp_86*tmp_86) + (tmp_87*tmp_87) + (tmp_88*tmp_88))) + 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_61*(-tmp_53*tmp_95 - tmp_56*tmp_96 - tmp_58*tmp_97 + 1.0*tmp_60*((tmp_95*tmp_95) + (tmp_96*tmp_96) + (tmp_97*tmp_97))) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_61*(-tmp_104*tmp_53 - tmp_105*tmp_56 - tmp_106*tmp_58 + 1.0*tmp_60*((tmp_104*tmp_104) + (tmp_105*tmp_105) + (tmp_106*tmp_106)));
      elMat( 0, 0) = a_0_0;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_1 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_2 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_5 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_11 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_12 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_13 = tmp_12*tmp_2;
      real_t tmp_14 = tmp_1*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_13 + tmp_0*tmp_9 + tmp_10*tmp_3 - tmp_10*tmp_6 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_14);
      real_t tmp_16 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_17 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_18 = -tmp_17;
      real_t tmp_19 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_20 = 0.091576213509770743*tmp_18 + 0.81684757298045851*tmp_19;
      real_t tmp_21 = tmp_15*(tmp_16 + tmp_20);
      real_t tmp_22 = tmp_12*tmp_4 - tmp_14;
      real_t tmp_23 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_24 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_25 = -tmp_24;
      real_t tmp_26 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_27 = 0.091576213509770743*tmp_25 + 0.81684757298045851*tmp_26;
      real_t tmp_28 = tmp_15*(tmp_23 + tmp_27);
      real_t tmp_29 = -tmp_13 + tmp_9;
      real_t tmp_30 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_31 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_32 = -tmp_31;
      real_t tmp_33 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_34 = 0.091576213509770743*tmp_32 + 0.81684757298045851*tmp_33;
      real_t tmp_35 = tmp_15*(tmp_30 + tmp_34);
      real_t tmp_36 = tmp_21*tmp_7 + tmp_22*tmp_28 + tmp_29*tmp_35 - 1.0/4.0;
      real_t tmp_37 = -tmp_0*tmp_2 + tmp_11*tmp_4;
      real_t tmp_38 = tmp_0*tmp_8 - tmp_10*tmp_4;
      real_t tmp_39 = tmp_10*tmp_2 - tmp_11*tmp_8;
      real_t tmp_40 = tmp_21*tmp_37 + tmp_28*tmp_38 + tmp_35*tmp_39 - 1.0/4.0;
      real_t tmp_41 = tmp_0*tmp_5 - tmp_1*tmp_11;
      real_t tmp_42 = -tmp_0*tmp_12 + tmp_1*tmp_10;
      real_t tmp_43 = -tmp_10*tmp_5 + tmp_11*tmp_12;
      real_t tmp_44 = tmp_21*tmp_41 + tmp_28*tmp_42 + tmp_35*tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_0*tmp_36 + tmp_1*tmp_40 + tmp_4*tmp_44;
      real_t tmp_46 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_47 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_48 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_49 = tmp_47*tmp_48;
      real_t tmp_50 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_51 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_52 = tmp_50*tmp_51;
      real_t tmp_53 = tmp_49 - tmp_52;
      real_t tmp_54 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_55 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_56 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_57 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_58 = tmp_51*tmp_57;
      real_t tmp_59 = tmp_48*tmp_57;
      real_t tmp_60 = tmp_47*tmp_55;
      real_t tmp_61 = 1.0 / (tmp_46*tmp_49 - tmp_46*tmp_52 + tmp_50*tmp_54*tmp_55 - tmp_54*tmp_59 + tmp_56*tmp_58 - tmp_56*tmp_60);
      real_t tmp_62 = 1.0*tmp_61;
      real_t tmp_63 = tmp_50*tmp_55 - tmp_59;
      real_t tmp_64 = tmp_58 - tmp_60;
      real_t tmp_65 = -tmp_48*tmp_54 + tmp_51*tmp_56;
      real_t tmp_66 = 0.5*tmp_61;
      real_t tmp_67 = tmp_46*tmp_66;
      real_t tmp_68 = tmp_46*tmp_48 - tmp_55*tmp_56;
      real_t tmp_69 = tmp_54*tmp_66;
      real_t tmp_70 = -tmp_46*tmp_51 + tmp_54*tmp_55;
      real_t tmp_71 = tmp_56*tmp_66;
      real_t tmp_72 = tmp_53*tmp_66;
      real_t tmp_73 = tmp_63*tmp_66;
      real_t tmp_74 = tmp_64*tmp_66;
      real_t tmp_75 = tmp_47*tmp_73 + tmp_50*tmp_74 + tmp_57*tmp_72 + tmp_65*tmp_67 + tmp_68*tmp_69 + tmp_70*tmp_71;
      real_t tmp_76 = -tmp_47*tmp_56 + tmp_50*tmp_54;
      real_t tmp_77 = -tmp_46*tmp_50 + tmp_56*tmp_57;
      real_t tmp_78 = tmp_46*tmp_47 - tmp_54*tmp_57;
      real_t tmp_79 = tmp_48*tmp_74 + tmp_51*tmp_73 + tmp_55*tmp_72 + tmp_67*tmp_76 + tmp_69*tmp_77 + tmp_71*tmp_78;
      real_t tmp_80 = p_affine_13_0*(tmp_46*tmp_53*tmp_62 + tmp_54*tmp_62*tmp_63 + tmp_56*tmp_62*tmp_64) + p_affine_13_1*tmp_75 + p_affine_13_2*tmp_79;
      real_t tmp_81 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_82 = tmp_61*(tmp_20 + tmp_81);
      real_t tmp_83 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_84 = tmp_61*(tmp_27 + tmp_83);
      real_t tmp_85 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_86 = tmp_61*(tmp_34 + tmp_85);
      real_t tmp_87 = tmp_53*tmp_86 + tmp_65*tmp_84 + tmp_76*tmp_82 - 1.0/4.0;
      real_t tmp_88 = tmp_63*tmp_86 + tmp_68*tmp_84 + tmp_77*tmp_82 - 1.0/4.0;
      real_t tmp_89 = tmp_64*tmp_86 + tmp_70*tmp_84 + tmp_78*tmp_82 - 1.0/4.0;
      real_t tmp_90 = tmp_46*tmp_87 + tmp_54*tmp_88 + tmp_56*tmp_89;
      real_t tmp_91 = 1.0*tmp_15;
      real_t tmp_92 = 0.5*tmp_15;
      real_t tmp_93 = tmp_0*tmp_92;
      real_t tmp_94 = tmp_1*tmp_92;
      real_t tmp_95 = tmp_4*tmp_92;
      real_t tmp_96 = tmp_29*tmp_92;
      real_t tmp_97 = tmp_39*tmp_92;
      real_t tmp_98 = tmp_43*tmp_92;
      real_t tmp_99 = tmp_11*tmp_96 + tmp_2*tmp_98 + tmp_22*tmp_93 + tmp_38*tmp_94 + tmp_42*tmp_95 + tmp_5*tmp_97;
      real_t tmp_100 = tmp_10*tmp_96 + tmp_12*tmp_97 + tmp_37*tmp_94 + tmp_41*tmp_95 + tmp_7*tmp_93 + tmp_8*tmp_98;
      real_t tmp_101 = 0.5*p_affine_13_0*(tmp_0*tmp_29*tmp_91 + tmp_1*tmp_39*tmp_91 + tmp_4*tmp_43*tmp_91) + 0.5*p_affine_13_1*tmp_99 + 0.5*p_affine_13_2*tmp_100;
      real_t tmp_102 = tmp_10*tmp_36 + tmp_12*tmp_40 + tmp_44*tmp_8;
      real_t tmp_103 = tmp_47*tmp_66*tmp_77 + tmp_48*tmp_66*tmp_70 + tmp_50*tmp_66*tmp_78 + tmp_51*tmp_66*tmp_68 + tmp_55*tmp_65*tmp_66 + tmp_57*tmp_66*tmp_76;
      real_t tmp_104 = p_affine_13_0*tmp_79 + p_affine_13_1*tmp_103 + p_affine_13_2*(tmp_48*tmp_62*tmp_78 + tmp_51*tmp_62*tmp_77 + tmp_55*tmp_62*tmp_76);
      real_t tmp_105 = tmp_11*tmp_36 + tmp_2*tmp_44 + tmp_40*tmp_5;
      real_t tmp_106 = p_affine_13_0*tmp_75 + p_affine_13_1*(tmp_47*tmp_62*tmp_68 + tmp_50*tmp_62*tmp_70 + tmp_57*tmp_62*tmp_65) + p_affine_13_2*tmp_103;
      real_t tmp_107 = tmp_48*tmp_89 + tmp_51*tmp_88 + tmp_55*tmp_87;
      real_t tmp_108 = tmp_10*tmp_22*tmp_92 + tmp_11*tmp_7*tmp_92 + tmp_12*tmp_38*tmp_92 + tmp_2*tmp_41*tmp_92 + tmp_37*tmp_5*tmp_92 + tmp_42*tmp_8*tmp_92;
      real_t tmp_109 = 0.5*p_affine_13_0*tmp_100 + 0.5*p_affine_13_1*tmp_108 + 0.5*p_affine_13_2*(tmp_10*tmp_7*tmp_91 + tmp_12*tmp_37*tmp_91 + tmp_41*tmp_8*tmp_91);
      real_t tmp_110 = tmp_47*tmp_88 + tmp_50*tmp_89 + tmp_57*tmp_87;
      real_t tmp_111 = 0.5*p_affine_13_0*tmp_99 + 0.5*p_affine_13_1*(tmp_11*tmp_22*tmp_91 + tmp_2*tmp_42*tmp_91 + tmp_38*tmp_5*tmp_91) + 0.5*p_affine_13_2*tmp_108;
      real_t tmp_112 = (std::abs(tmp_17*tmp_26 - tmp_19*tmp_24)*std::abs(tmp_17*tmp_26 - tmp_19*tmp_24)) + (std::abs(tmp_17*tmp_33 - tmp_19*tmp_31)*std::abs(tmp_17*tmp_33 - tmp_19*tmp_31)) + (std::abs(tmp_24*tmp_33 - tmp_26*tmp_31)*std::abs(tmp_24*tmp_33 - tmp_26*tmp_31));
      real_t tmp_113 = 1.0*std::pow(tmp_112, -0.25);
      real_t tmp_114 = 2.0*std::pow(tmp_112, 1.0/2.0);
      real_t tmp_115 = 0.44594849091596489*tmp_18 + 0.10810301816807022*tmp_19;
      real_t tmp_116 = tmp_15*(tmp_115 + tmp_16);
      real_t tmp_117 = 0.44594849091596489*tmp_25 + 0.10810301816807022*tmp_26;
      real_t tmp_118 = tmp_15*(tmp_117 + tmp_23);
      real_t tmp_119 = 0.44594849091596489*tmp_32 + 0.10810301816807022*tmp_33;
      real_t tmp_120 = tmp_15*(tmp_119 + tmp_30);
      real_t tmp_121 = tmp_116*tmp_7 + tmp_118*tmp_22 + tmp_120*tmp_29 - 1.0/4.0;
      real_t tmp_122 = tmp_116*tmp_37 + tmp_118*tmp_38 + tmp_120*tmp_39 - 1.0/4.0;
      real_t tmp_123 = tmp_116*tmp_41 + tmp_118*tmp_42 + tmp_120*tmp_43 - 1.0/4.0;
      real_t tmp_124 = tmp_0*tmp_121 + tmp_1*tmp_122 + tmp_123*tmp_4;
      real_t tmp_125 = tmp_61*(tmp_115 + tmp_81);
      real_t tmp_126 = tmp_61*(tmp_117 + tmp_83);
      real_t tmp_127 = tmp_61*(tmp_119 + tmp_85);
      real_t tmp_128 = tmp_125*tmp_76 + tmp_126*tmp_65 + tmp_127*tmp_53 - 1.0/4.0;
      real_t tmp_129 = tmp_125*tmp_77 + tmp_126*tmp_68 + tmp_127*tmp_63 - 1.0/4.0;
      real_t tmp_130 = tmp_125*tmp_78 + tmp_126*tmp_70 + tmp_127*tmp_64 - 1.0/4.0;
      real_t tmp_131 = tmp_128*tmp_46 + tmp_129*tmp_54 + tmp_130*tmp_56;
      real_t tmp_132 = tmp_10*tmp_121 + tmp_12*tmp_122 + tmp_123*tmp_8;
      real_t tmp_133 = tmp_11*tmp_121 + tmp_122*tmp_5 + tmp_123*tmp_2;
      real_t tmp_134 = tmp_128*tmp_55 + tmp_129*tmp_51 + tmp_130*tmp_48;
      real_t tmp_135 = tmp_128*tmp_57 + tmp_129*tmp_47 + tmp_130*tmp_50;
      real_t tmp_136 = 0.81684757298045851*tmp_18 + 0.091576213509770743*tmp_19;
      real_t tmp_137 = tmp_15*(tmp_136 + tmp_16);
      real_t tmp_138 = 0.81684757298045851*tmp_25 + 0.091576213509770743*tmp_26;
      real_t tmp_139 = tmp_15*(tmp_138 + tmp_23);
      real_t tmp_140 = 0.81684757298045851*tmp_32 + 0.091576213509770743*tmp_33;
      real_t tmp_141 = tmp_15*(tmp_140 + tmp_30);
      real_t tmp_142 = tmp_137*tmp_7 + tmp_139*tmp_22 + tmp_141*tmp_29 - 1.0/4.0;
      real_t tmp_143 = tmp_137*tmp_37 + tmp_139*tmp_38 + tmp_141*tmp_39 - 1.0/4.0;
      real_t tmp_144 = tmp_137*tmp_41 + tmp_139*tmp_42 + tmp_141*tmp_43 - 1.0/4.0;
      real_t tmp_145 = tmp_0*tmp_142 + tmp_1*tmp_143 + tmp_144*tmp_4;
      real_t tmp_146 = tmp_61*(tmp_136 + tmp_81);
      real_t tmp_147 = tmp_61*(tmp_138 + tmp_83);
      real_t tmp_148 = tmp_61*(tmp_140 + tmp_85);
      real_t tmp_149 = tmp_146*tmp_76 + tmp_147*tmp_65 + tmp_148*tmp_53 - 1.0/4.0;
      real_t tmp_150 = tmp_146*tmp_77 + tmp_147*tmp_68 + tmp_148*tmp_63 - 1.0/4.0;
      real_t tmp_151 = tmp_146*tmp_78 + tmp_147*tmp_70 + tmp_148*tmp_64 - 1.0/4.0;
      real_t tmp_152 = tmp_149*tmp_46 + tmp_150*tmp_54 + tmp_151*tmp_56;
      real_t tmp_153 = tmp_10*tmp_142 + tmp_12*tmp_143 + tmp_144*tmp_8;
      real_t tmp_154 = tmp_11*tmp_142 + tmp_143*tmp_5 + tmp_144*tmp_2;
      real_t tmp_155 = tmp_149*tmp_55 + tmp_150*tmp_51 + tmp_151*tmp_48;
      real_t tmp_156 = tmp_149*tmp_57 + tmp_150*tmp_47 + tmp_151*tmp_50;
      real_t tmp_157 = 0.10810301816807022*tmp_18 + 0.44594849091596489*tmp_19;
      real_t tmp_158 = tmp_15*(tmp_157 + tmp_16);
      real_t tmp_159 = 0.10810301816807022*tmp_25 + 0.44594849091596489*tmp_26;
      real_t tmp_160 = tmp_15*(tmp_159 + tmp_23);
      real_t tmp_161 = 0.10810301816807022*tmp_32 + 0.44594849091596489*tmp_33;
      real_t tmp_162 = tmp_15*(tmp_161 + tmp_30);
      real_t tmp_163 = tmp_158*tmp_7 + tmp_160*tmp_22 + tmp_162*tmp_29 - 1.0/4.0;
      real_t tmp_164 = tmp_158*tmp_37 + tmp_160*tmp_38 + tmp_162*tmp_39 - 1.0/4.0;
      real_t tmp_165 = tmp_158*tmp_41 + tmp_160*tmp_42 + tmp_162*tmp_43 - 1.0/4.0;
      real_t tmp_166 = tmp_0*tmp_163 + tmp_1*tmp_164 + tmp_165*tmp_4;
      real_t tmp_167 = tmp_61*(tmp_157 + tmp_81);
      real_t tmp_168 = tmp_61*(tmp_159 + tmp_83);
      real_t tmp_169 = tmp_61*(tmp_161 + tmp_85);
      real_t tmp_170 = tmp_167*tmp_76 + tmp_168*tmp_65 + tmp_169*tmp_53 - 1.0/4.0;
      real_t tmp_171 = tmp_167*tmp_77 + tmp_168*tmp_68 + tmp_169*tmp_63 - 1.0/4.0;
      real_t tmp_172 = tmp_167*tmp_78 + tmp_168*tmp_70 + tmp_169*tmp_64 - 1.0/4.0;
      real_t tmp_173 = tmp_170*tmp_46 + tmp_171*tmp_54 + tmp_172*tmp_56;
      real_t tmp_174 = tmp_10*tmp_163 + tmp_12*tmp_164 + tmp_165*tmp_8;
      real_t tmp_175 = tmp_11*tmp_163 + tmp_164*tmp_5 + tmp_165*tmp_2;
      real_t tmp_176 = tmp_170*tmp_55 + tmp_171*tmp_51 + tmp_172*tmp_48;
      real_t tmp_177 = tmp_170*tmp_57 + tmp_171*tmp_47 + tmp_172*tmp_50;
      real_t tmp_178 = 0.091576213509770743*tmp_18 + 0.091576213509770743*tmp_19;
      real_t tmp_179 = tmp_15*(tmp_16 + tmp_178);
      real_t tmp_180 = 0.091576213509770743*tmp_25 + 0.091576213509770743*tmp_26;
      real_t tmp_181 = tmp_15*(tmp_180 + tmp_23);
      real_t tmp_182 = 0.091576213509770743*tmp_32 + 0.091576213509770743*tmp_33;
      real_t tmp_183 = tmp_15*(tmp_182 + tmp_30);
      real_t tmp_184 = tmp_179*tmp_7 + tmp_181*tmp_22 + tmp_183*tmp_29 - 1.0/4.0;
      real_t tmp_185 = tmp_179*tmp_37 + tmp_181*tmp_38 + tmp_183*tmp_39 - 1.0/4.0;
      real_t tmp_186 = tmp_179*tmp_41 + tmp_181*tmp_42 + tmp_183*tmp_43 - 1.0/4.0;
      real_t tmp_187 = tmp_0*tmp_184 + tmp_1*tmp_185 + tmp_186*tmp_4;
      real_t tmp_188 = tmp_61*(tmp_178 + tmp_81);
      real_t tmp_189 = tmp_61*(tmp_180 + tmp_83);
      real_t tmp_190 = tmp_61*(tmp_182 + tmp_85);
      real_t tmp_191 = tmp_188*tmp_76 + tmp_189*tmp_65 + tmp_190*tmp_53 - 1.0/4.0;
      real_t tmp_192 = tmp_188*tmp_77 + tmp_189*tmp_68 + tmp_190*tmp_63 - 1.0/4.0;
      real_t tmp_193 = tmp_188*tmp_78 + tmp_189*tmp_70 + tmp_190*tmp_64 - 1.0/4.0;
      real_t tmp_194 = tmp_191*tmp_46 + tmp_192*tmp_54 + tmp_193*tmp_56;
      real_t tmp_195 = tmp_10*tmp_184 + tmp_12*tmp_185 + tmp_186*tmp_8;
      real_t tmp_196 = tmp_11*tmp_184 + tmp_185*tmp_5 + tmp_186*tmp_2;
      real_t tmp_197 = tmp_191*tmp_55 + tmp_192*tmp_51 + tmp_193*tmp_48;
      real_t tmp_198 = tmp_191*tmp_57 + tmp_192*tmp_47 + tmp_193*tmp_50;
      real_t tmp_199 = 0.44594849091596489*tmp_18 + 0.44594849091596489*tmp_19;
      real_t tmp_200 = tmp_15*(tmp_16 + tmp_199);
      real_t tmp_201 = 0.44594849091596489*tmp_25 + 0.44594849091596489*tmp_26;
      real_t tmp_202 = tmp_15*(tmp_201 + tmp_23);
      real_t tmp_203 = 0.44594849091596489*tmp_32 + 0.44594849091596489*tmp_33;
      real_t tmp_204 = tmp_15*(tmp_203 + tmp_30);
      real_t tmp_205 = tmp_200*tmp_7 + tmp_202*tmp_22 + tmp_204*tmp_29 - 1.0/4.0;
      real_t tmp_206 = tmp_200*tmp_37 + tmp_202*tmp_38 + tmp_204*tmp_39 - 1.0/4.0;
      real_t tmp_207 = tmp_200*tmp_41 + tmp_202*tmp_42 + tmp_204*tmp_43 - 1.0/4.0;
      real_t tmp_208 = tmp_0*tmp_205 + tmp_1*tmp_206 + tmp_207*tmp_4;
      real_t tmp_209 = tmp_61*(tmp_199 + tmp_81);
      real_t tmp_210 = tmp_61*(tmp_201 + tmp_83);
      real_t tmp_211 = tmp_61*(tmp_203 + tmp_85);
      real_t tmp_212 = tmp_209*tmp_76 + tmp_210*tmp_65 + tmp_211*tmp_53 - 1.0/4.0;
      real_t tmp_213 = tmp_209*tmp_77 + tmp_210*tmp_68 + tmp_211*tmp_63 - 1.0/4.0;
      real_t tmp_214 = tmp_209*tmp_78 + tmp_210*tmp_70 + tmp_211*tmp_64 - 1.0/4.0;
      real_t tmp_215 = tmp_212*tmp_46 + tmp_213*tmp_54 + tmp_214*tmp_56;
      real_t tmp_216 = tmp_10*tmp_205 + tmp_12*tmp_206 + tmp_207*tmp_8;
      real_t tmp_217 = tmp_11*tmp_205 + tmp_2*tmp_207 + tmp_206*tmp_5;
      real_t tmp_218 = tmp_212*tmp_55 + tmp_213*tmp_51 + tmp_214*tmp_48;
      real_t tmp_219 = tmp_212*tmp_57 + tmp_213*tmp_47 + tmp_214*tmp_50;
      real_t a_0_0 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_114*(-tmp_101*tmp_90 + 0.5*tmp_102*tmp_104 + 0.5*tmp_105*tmp_106 - tmp_107*tmp_109 - tmp_110*tmp_111 - tmp_113*(tmp_102*tmp_107 + tmp_105*tmp_110 + tmp_45*tmp_90) + 0.5*tmp_45*tmp_80) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_114*(-tmp_101*tmp_131 + 0.5*tmp_104*tmp_132 + 0.5*tmp_106*tmp_133 - tmp_109*tmp_134 - tmp_111*tmp_135 - tmp_113*(tmp_124*tmp_131 + tmp_132*tmp_134 + tmp_133*tmp_135) + 0.5*tmp_124*tmp_80) + 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_114*(-tmp_101*tmp_152 + 0.5*tmp_104*tmp_153 + 0.5*tmp_106*tmp_154 - tmp_109*tmp_155 - tmp_111*tmp_156 - tmp_113*(tmp_145*tmp_152 + tmp_153*tmp_155 + tmp_154*tmp_156) + 0.5*tmp_145*tmp_80) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_114*(-tmp_101*tmp_173 + 0.5*tmp_104*tmp_174 + 0.5*tmp_106*tmp_175 - tmp_109*tmp_176 - tmp_111*tmp_177 - tmp_113*(tmp_166*tmp_173 + tmp_174*tmp_176 + tmp_175*tmp_177) + 0.5*tmp_166*tmp_80) + 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_114*(-tmp_101*tmp_194 + 0.5*tmp_104*tmp_195 + 0.5*tmp_106*tmp_196 - tmp_109*tmp_197 - tmp_111*tmp_198 - tmp_113*(tmp_187*tmp_194 + tmp_195*tmp_197 + tmp_196*tmp_198) + 0.5*tmp_187*tmp_80) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_114*(-tmp_101*tmp_215 + 0.5*tmp_104*tmp_216 + 0.5*tmp_106*tmp_217 - tmp_109*tmp_218 - tmp_111*tmp_219 - tmp_113*(tmp_208*tmp_215 + tmp_216*tmp_218 + tmp_217*tmp_219) + 0.5*tmp_208*tmp_80);
      elMat( 0, 0) = a_0_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      Scalar_Variable_Coefficient_3D_mu( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id1 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id2 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id3 );
      Scalar_Variable_Coefficient_3D_mu( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id4 );
      Scalar_Variable_Coefficient_3D_mu( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_mu_out0_id5 );
      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_4*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_6*tmp_9;
      real_t tmp_14 = tmp_4*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_12 + tmp_0*tmp_7 - tmp_1*tmp_13 + tmp_1*tmp_2*tmp_8 + tmp_11*tmp_3 - tmp_14*tmp_3);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.091576213509770743*tmp_23 + 0.81684757298045851*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_12 + tmp_7;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.091576213509770743*tmp_29 + 0.81684757298045851*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_21*tmp_26 + tmp_27*tmp_32 - 1.0/4.0;
      real_t tmp_34 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_35 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_36 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36 - 1.0/4.0;
      real_t tmp_38 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_39 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_40 = tmp_11 - tmp_14;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40 - 1.0/4.0;
      real_t tmp_42 = tmp_0*tmp_33 + tmp_1*tmp_37 + tmp_3*tmp_41;
      real_t tmp_43 = 1.0*tmp_15;
      real_t tmp_44 = 0.5*tmp_15;
      real_t tmp_45 = tmp_0*tmp_44;
      real_t tmp_46 = tmp_1*tmp_44;
      real_t tmp_47 = tmp_3*tmp_44;
      real_t tmp_48 = tmp_27*tmp_44;
      real_t tmp_49 = tmp_36*tmp_44;
      real_t tmp_50 = tmp_40*tmp_44;
      real_t tmp_51 = tmp_2*tmp_50 + tmp_21*tmp_45 + tmp_35*tmp_46 + tmp_39*tmp_47 + tmp_4*tmp_49 + tmp_48*tmp_9;
      real_t tmp_52 = tmp_10*tmp_49 + tmp_34*tmp_46 + tmp_38*tmp_47 + tmp_45*tmp_5 + tmp_48*tmp_8 + tmp_50*tmp_6;
      real_t tmp_53 = 2*p_affine_13_0*(tmp_0*tmp_27*tmp_43 + tmp_1*tmp_36*tmp_43 + tmp_3*tmp_40*tmp_43) + 2*p_affine_13_1*tmp_51 + 2*p_affine_13_2*tmp_52;
      real_t tmp_54 = tmp_10*tmp_37 + tmp_33*tmp_8 + tmp_41*tmp_6;
      real_t tmp_55 = tmp_10*tmp_35*tmp_44 + tmp_2*tmp_38*tmp_44 + tmp_21*tmp_44*tmp_8 + tmp_34*tmp_4*tmp_44 + tmp_39*tmp_44*tmp_6 + tmp_44*tmp_5*tmp_9;
      real_t tmp_56 = 2*p_affine_13_0*tmp_52 + 2*p_affine_13_1*tmp_55 + 2*p_affine_13_2*(tmp_10*tmp_34*tmp_43 + tmp_38*tmp_43*tmp_6 + tmp_43*tmp_5*tmp_8);
      real_t tmp_57 = tmp_2*tmp_41 + tmp_33*tmp_9 + tmp_37*tmp_4;
      real_t tmp_58 = 2*p_affine_13_0*tmp_51 + 2*p_affine_13_1*(tmp_2*tmp_39*tmp_43 + tmp_21*tmp_43*tmp_9 + tmp_35*tmp_4*tmp_43) + 2*p_affine_13_2*tmp_55;
      real_t tmp_59 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_60 = std::pow(tmp_59, -0.25);
      real_t tmp_61 = 2.0*std::pow(tmp_59, 1.0/2.0);
      real_t tmp_62 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_63 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_64 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_65 = tmp_21*tmp_63 + tmp_27*tmp_64 + tmp_5*tmp_62 - 1.0/4.0;
      real_t tmp_66 = tmp_34*tmp_62 + tmp_35*tmp_63 + tmp_36*tmp_64 - 1.0/4.0;
      real_t tmp_67 = tmp_38*tmp_62 + tmp_39*tmp_63 + tmp_40*tmp_64 - 1.0/4.0;
      real_t tmp_68 = tmp_0*tmp_65 + tmp_1*tmp_66 + tmp_3*tmp_67;
      real_t tmp_69 = tmp_10*tmp_66 + tmp_6*tmp_67 + tmp_65*tmp_8;
      real_t tmp_70 = tmp_2*tmp_67 + tmp_4*tmp_66 + tmp_65*tmp_9;
      real_t tmp_71 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_72 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_73 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_74 = tmp_21*tmp_72 + tmp_27*tmp_73 + tmp_5*tmp_71 - 1.0/4.0;
      real_t tmp_75 = tmp_34*tmp_71 + tmp_35*tmp_72 + tmp_36*tmp_73 - 1.0/4.0;
      real_t tmp_76 = tmp_38*tmp_71 + tmp_39*tmp_72 + tmp_40*tmp_73 - 1.0/4.0;
      real_t tmp_77 = tmp_0*tmp_74 + tmp_1*tmp_75 + tmp_3*tmp_76;
      real_t tmp_78 = tmp_10*tmp_75 + tmp_6*tmp_76 + tmp_74*tmp_8;
      real_t tmp_79 = tmp_2*tmp_76 + tmp_4*tmp_75 + tmp_74*tmp_9;
      real_t tmp_80 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_81 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_82 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_83 = tmp_21*tmp_81 + tmp_27*tmp_82 + tmp_5*tmp_80 - 1.0/4.0;
      real_t tmp_84 = tmp_34*tmp_80 + tmp_35*tmp_81 + tmp_36*tmp_82 - 1.0/4.0;
      real_t tmp_85 = tmp_38*tmp_80 + tmp_39*tmp_81 + tmp_40*tmp_82 - 1.0/4.0;
      real_t tmp_86 = tmp_0*tmp_83 + tmp_1*tmp_84 + tmp_3*tmp_85;
      real_t tmp_87 = tmp_10*tmp_84 + tmp_6*tmp_85 + tmp_8*tmp_83;
      real_t tmp_88 = tmp_2*tmp_85 + tmp_4*tmp_84 + tmp_83*tmp_9;
      real_t tmp_89 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_90 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_91 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_92 = tmp_21*tmp_90 + tmp_27*tmp_91 + tmp_5*tmp_89 - 1.0/4.0;
      real_t tmp_93 = tmp_34*tmp_89 + tmp_35*tmp_90 + tmp_36*tmp_91 - 1.0/4.0;
      real_t tmp_94 = tmp_38*tmp_89 + tmp_39*tmp_90 + tmp_40*tmp_91 - 1.0/4.0;
      real_t tmp_95 = tmp_0*tmp_92 + tmp_1*tmp_93 + tmp_3*tmp_94;
      real_t tmp_96 = tmp_10*tmp_93 + tmp_6*tmp_94 + tmp_8*tmp_92;
      real_t tmp_97 = tmp_2*tmp_94 + tmp_4*tmp_93 + tmp_9*tmp_92;
      real_t tmp_98 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_99 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_100 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_101 = tmp_100*tmp_27 + tmp_21*tmp_99 + tmp_5*tmp_98 - 1.0/4.0;
      real_t tmp_102 = tmp_100*tmp_36 + tmp_34*tmp_98 + tmp_35*tmp_99 - 1.0/4.0;
      real_t tmp_103 = tmp_100*tmp_40 + tmp_38*tmp_98 + tmp_39*tmp_99 - 1.0/4.0;
      real_t tmp_104 = tmp_0*tmp_101 + tmp_1*tmp_102 + tmp_103*tmp_3;
      real_t tmp_105 = tmp_10*tmp_102 + tmp_101*tmp_8 + tmp_103*tmp_6;
      real_t tmp_106 = tmp_101*tmp_9 + tmp_102*tmp_4 + tmp_103*tmp_2;
      real_t a_0_0 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_61*(-tmp_42*tmp_53 - tmp_54*tmp_56 - tmp_57*tmp_58 + 1.0*tmp_60*((tmp_42*tmp_42) + (tmp_54*tmp_54) + (tmp_57*tmp_57))) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id1*tmp_61*(-tmp_53*tmp_68 - tmp_56*tmp_69 - tmp_58*tmp_70 + 1.0*tmp_60*((tmp_68*tmp_68) + (tmp_69*tmp_69) + (tmp_70*tmp_70))) + 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id2*tmp_61*(-tmp_53*tmp_77 - tmp_56*tmp_78 - tmp_58*tmp_79 + 1.0*tmp_60*((tmp_77*tmp_77) + (tmp_78*tmp_78) + (tmp_79*tmp_79))) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id3*tmp_61*(-tmp_53*tmp_86 - tmp_56*tmp_87 - tmp_58*tmp_88 + 1.0*tmp_60*((tmp_86*tmp_86) + (tmp_87*tmp_87) + (tmp_88*tmp_88))) + 0.054975871827660928*Scalar_Variable_Coefficient_3D_mu_out0_id4*tmp_61*(-tmp_53*tmp_95 - tmp_56*tmp_96 - tmp_58*tmp_97 + 1.0*tmp_60*((tmp_95*tmp_95) + (tmp_96*tmp_96) + (tmp_97*tmp_97))) + 0.11169079483900572*Scalar_Variable_Coefficient_3D_mu_out0_id5*tmp_61*(-tmp_104*tmp_53 - tmp_105*tmp_56 - tmp_106*tmp_58 + 1.0*tmp_60*((tmp_104*tmp_104) + (tmp_105*tmp_105) + (tmp_106*tmp_106)));
      elMat( 0, 0) = a_0_0;
   }

};


} //eg
} // dg
} // hyteg
