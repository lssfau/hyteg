
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

class EGVectorLaplaceFormNitscheBC_P1P1_00 : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_P1P1_00()
: callback_Scalar_Variable_Coefficient_3D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_g0( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g0( Point3D( {in_0, in_1, in_2} ) );
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = 1.0 / (tmp_4);
      real_t tmp_6 = tmp_0*tmp_5;
      real_t tmp_7 = tmp_2*tmp_5;
      real_t tmp_8 = -tmp_6 - tmp_7;
      real_t tmp_9 = tmp_1*tmp_5;
      real_t tmp_10 = tmp_3*tmp_5;
      real_t tmp_11 = -tmp_10 - tmp_9;
      real_t tmp_12 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_13 = tmp_12*((tmp_11*tmp_11) + (tmp_8*tmp_8));
      real_t tmp_14 = tmp_12*(tmp_11*tmp_9 + tmp_7*tmp_8);
      real_t tmp_15 = 0.5*tmp_14;
      real_t tmp_16 = tmp_12*(tmp_10*tmp_11 + tmp_6*tmp_8);
      real_t tmp_17 = 0.5*tmp_16;
      real_t tmp_18 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_19 = tmp_12*((tmp_1*tmp_1)*tmp_18 + tmp_18*(tmp_2*tmp_2));
      real_t tmp_20 = tmp_12*(tmp_0*tmp_18*tmp_2 + tmp_1*tmp_18*tmp_3);
      real_t tmp_21 = 0.5*tmp_20;
      real_t tmp_22 = tmp_12*((tmp_0*tmp_0)*tmp_18 + tmp_18*(tmp_3*tmp_3));
      real_t a_0_0 = 0.5*tmp_13;
      real_t a_0_1 = tmp_15;
      real_t a_0_2 = tmp_17;
      real_t a_1_0 = tmp_15;
      real_t a_1_1 = 0.5*tmp_19;
      real_t a_1_2 = tmp_21;
      real_t a_2_0 = tmp_17;
      real_t a_2_1 = tmp_21;
      real_t a_2_2 = 0.5*tmp_22;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_5 = 0.21132486540518713*tmp_1 + tmp_4;
      real_t tmp_6 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_7 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_8 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = 1.0 / (-tmp_6*tmp_9 + tmp_7*tmp_8);
      real_t tmp_11 = tmp_10*tmp_6;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_0 + tmp_12;
      real_t tmp_14 = tmp_10*tmp_8;
      real_t tmp_15 = tmp_11*tmp_5 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_10*tmp_7;
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_5;
      real_t tmp_19 = -tmp_15 - tmp_18 + 1;
      real_t tmp_20 = p_affine_10_0*(-tmp_14 - tmp_17) + p_affine_10_1*(-tmp_11 - tmp_16);
      real_t tmp_21 = 1.0*tmp_20;
      real_t tmp_22 = 0.5*tmp_2;
      real_t tmp_23 = 0.78867513459481287*tmp_1 + tmp_4;
      real_t tmp_24 = 0.78867513459481287*tmp_0 + tmp_12;
      real_t tmp_25 = tmp_11*tmp_23 + tmp_14*tmp_24;
      real_t tmp_26 = tmp_16*tmp_23 + tmp_17*tmp_24;
      real_t tmp_27 = -tmp_25 - tmp_26 + 1;
      real_t tmp_28 = 0.5*tmp_2;
      real_t tmp_29 = 0.5*tmp_20;
      real_t tmp_30 = p_affine_10_0*tmp_14 + p_affine_10_1*tmp_11;
      real_t tmp_31 = 0.5*tmp_30;
      real_t tmp_32 = tmp_22*(3*tmp_15*tmp_19*tmp_3 - tmp_15*tmp_29 - tmp_19*tmp_31) + tmp_28*(3*tmp_25*tmp_27*tmp_3 - tmp_25*tmp_29 - tmp_27*tmp_31);
      real_t tmp_33 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_16;
      real_t tmp_34 = 0.5*tmp_33;
      real_t tmp_35 = tmp_22*(3*tmp_18*tmp_19*tmp_3 - tmp_18*tmp_29 - tmp_19*tmp_34) + tmp_28*(3*tmp_26*tmp_27*tmp_3 - tmp_26*tmp_29 - tmp_27*tmp_34);
      real_t tmp_36 = 1.0*tmp_30;
      real_t tmp_37 = tmp_22*(3*tmp_15*tmp_18*tmp_3 - tmp_15*tmp_34 - tmp_18*tmp_31) + tmp_28*(3*tmp_25*tmp_26*tmp_3 - tmp_25*tmp_34 - tmp_26*tmp_31);
      real_t tmp_38 = 1.0*tmp_33;
      real_t a_0_0 = tmp_22*(3*(tmp_19*tmp_19)*tmp_3 - tmp_19*tmp_21) + tmp_28*(-tmp_21*tmp_27 + 3*(tmp_27*tmp_27)*tmp_3);
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_22*(3*(tmp_15*tmp_15)*tmp_3 - tmp_15*tmp_36) + tmp_28*(3*(tmp_25*tmp_25)*tmp_3 - tmp_25*tmp_36);
      real_t a_1_2 = tmp_37;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_37;
      real_t a_2_2 = tmp_22*(3*(tmp_18*tmp_18)*tmp_3 - tmp_18*tmp_38) + tmp_28*(3*(tmp_26*tmp_26)*tmp_3 - tmp_26*tmp_38);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
                                          MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = 1.0 / (tmp_0*tmp_1 - tmp_2*tmp_3);
      real_t tmp_5 = tmp_0*tmp_4;
      real_t tmp_6 = tmp_3*tmp_4;
      real_t tmp_7 = tmp_1*tmp_4;
      real_t tmp_8 = tmp_2*tmp_4;
      real_t tmp_9 = p_affine_10_0*(-tmp_5 - tmp_6) + p_affine_10_1*(-tmp_7 - tmp_8);
      real_t tmp_10 = -p_affine_3_1;
      real_t tmp_11 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_12 = p_affine_6_1 + 0.21132486540518713*tmp_11;
      real_t tmp_13 = tmp_10 + tmp_12;
      real_t tmp_14 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_15 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_16 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_17 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_18 = 1.0 / (-tmp_14*tmp_17 + tmp_15*tmp_16);
      real_t tmp_19 = tmp_14*tmp_18;
      real_t tmp_20 = -p_affine_3_0;
      real_t tmp_21 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_22 = p_affine_6_0 + 0.21132486540518713*tmp_21;
      real_t tmp_23 = tmp_20 + tmp_22;
      real_t tmp_24 = tmp_16*tmp_18;
      real_t tmp_25 = tmp_13*tmp_19 + tmp_23*tmp_24;
      real_t tmp_26 = tmp_15*tmp_18;
      real_t tmp_27 = tmp_17*tmp_18;
      real_t tmp_28 = tmp_13*tmp_26 + tmp_23*tmp_27;
      real_t tmp_29 = -tmp_25 - tmp_28 + 1;
      real_t tmp_30 = -p_affine_0_1;
      real_t tmp_31 = tmp_12 + tmp_30;
      real_t tmp_32 = -p_affine_0_0;
      real_t tmp_33 = tmp_22 + tmp_32;
      real_t tmp_34 = tmp_31*tmp_8 + tmp_33*tmp_5;
      real_t tmp_35 = tmp_31*tmp_7 + tmp_33*tmp_6;
      real_t tmp_36 = -tmp_34 - tmp_35 + 1;
      real_t tmp_37 = 0.5*p_affine_10_0*(-tmp_24 - tmp_27) + 0.5*p_affine_10_1*(-tmp_19 - tmp_26);
      real_t tmp_38 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_21*tmp_21), 1.0/2.0));
      real_t tmp_39 = 3/tmp_38;
      real_t tmp_40 = tmp_36*tmp_39;
      real_t tmp_41 = 0.5*tmp_38;
      real_t tmp_42 = p_affine_6_1 + 0.78867513459481287*tmp_11;
      real_t tmp_43 = tmp_10 + tmp_42;
      real_t tmp_44 = p_affine_6_0 + 0.78867513459481287*tmp_21;
      real_t tmp_45 = tmp_20 + tmp_44;
      real_t tmp_46 = tmp_19*tmp_43 + tmp_24*tmp_45;
      real_t tmp_47 = tmp_26*tmp_43 + tmp_27*tmp_45;
      real_t tmp_48 = -tmp_46 - tmp_47 + 1;
      real_t tmp_49 = tmp_30 + tmp_42;
      real_t tmp_50 = tmp_32 + tmp_44;
      real_t tmp_51 = tmp_49*tmp_8 + tmp_5*tmp_50;
      real_t tmp_52 = tmp_49*tmp_7 + tmp_50*tmp_6;
      real_t tmp_53 = -tmp_51 - tmp_52 + 1;
      real_t tmp_54 = tmp_39*tmp_53;
      real_t tmp_55 = 0.5*tmp_38;
      real_t tmp_56 = 0.5*p_affine_10_0*tmp_24 + 0.5*p_affine_10_1*tmp_19;
      real_t tmp_57 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*tmp_26;
      real_t tmp_58 = p_affine_10_0*tmp_5 + p_affine_10_1*tmp_8;
      real_t tmp_59 = tmp_34*tmp_39;
      real_t tmp_60 = tmp_39*tmp_51;
      real_t tmp_61 = p_affine_10_0*tmp_6 + p_affine_10_1*tmp_7;
      real_t tmp_62 = tmp_35*tmp_39;
      real_t tmp_63 = tmp_39*tmp_52;
      real_t a_0_0 = tmp_41*(-tmp_29*tmp_40 + 0.5*tmp_29*tmp_9 - tmp_36*tmp_37) + tmp_55*(-tmp_37*tmp_53 - tmp_48*tmp_54 + 0.5*tmp_48*tmp_9);
      real_t a_0_1 = tmp_41*(-tmp_25*tmp_40 + 0.5*tmp_25*tmp_9 - tmp_36*tmp_56) + tmp_55*(-tmp_46*tmp_54 + 0.5*tmp_46*tmp_9 - tmp_53*tmp_56);
      real_t a_0_2 = tmp_41*(-tmp_28*tmp_40 + 0.5*tmp_28*tmp_9 - tmp_36*tmp_57) + tmp_55*(-tmp_47*tmp_54 + 0.5*tmp_47*tmp_9 - tmp_53*tmp_57);
      real_t a_1_0 = tmp_41*(0.5*tmp_29*tmp_58 - tmp_29*tmp_59 - tmp_34*tmp_37) + tmp_55*(-tmp_37*tmp_51 + 0.5*tmp_48*tmp_58 - tmp_48*tmp_60);
      real_t a_1_1 = tmp_41*(0.5*tmp_25*tmp_58 - tmp_25*tmp_59 - tmp_34*tmp_56) + tmp_55*(0.5*tmp_46*tmp_58 - tmp_46*tmp_60 - tmp_51*tmp_56);
      real_t a_1_2 = tmp_41*(0.5*tmp_28*tmp_58 - tmp_28*tmp_59 - tmp_34*tmp_57) + tmp_55*(0.5*tmp_47*tmp_58 - tmp_47*tmp_60 - tmp_51*tmp_57);
      real_t a_2_0 = tmp_41*(0.5*tmp_29*tmp_61 - tmp_29*tmp_62 - tmp_35*tmp_37) + tmp_55*(-tmp_37*tmp_52 + 0.5*tmp_48*tmp_61 - tmp_48*tmp_63);
      real_t a_2_1 = tmp_41*(0.5*tmp_25*tmp_61 - tmp_25*tmp_62 - tmp_35*tmp_56) + tmp_55*(0.5*tmp_46*tmp_61 - tmp_46*tmp_63 - tmp_52*tmp_56);
      real_t a_2_2 = tmp_41*(0.5*tmp_28*tmp_61 - tmp_28*tmp_62 - tmp_35*tmp_57) + tmp_55*(0.5*tmp_47*tmp_61 - tmp_47*tmp_63 - tmp_52*tmp_57);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_5 = 0.21132486540518713*tmp_1 + tmp_4;
      real_t tmp_6 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_7 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_8 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = 1.0 / (-tmp_6*tmp_9 + tmp_7*tmp_8);
      real_t tmp_11 = tmp_10*tmp_6;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_0 + tmp_12;
      real_t tmp_14 = tmp_10*tmp_8;
      real_t tmp_15 = tmp_11*tmp_5 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_10*tmp_7;
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_5;
      real_t tmp_19 = -tmp_15 - tmp_18 + 1;
      real_t tmp_20 = p_affine_10_0*(-tmp_14 - tmp_17) + p_affine_10_1*(-tmp_11 - tmp_16);
      real_t tmp_21 = 2*tmp_20;
      real_t tmp_22 = 0.5*tmp_2;
      real_t tmp_23 = 0.78867513459481287*tmp_1 + tmp_4;
      real_t tmp_24 = 0.78867513459481287*tmp_0 + tmp_12;
      real_t tmp_25 = tmp_11*tmp_23 + tmp_14*tmp_24;
      real_t tmp_26 = tmp_16*tmp_23 + tmp_17*tmp_24;
      real_t tmp_27 = -tmp_25 - tmp_26 + 1;
      real_t tmp_28 = 0.5*tmp_2;
      real_t tmp_29 = p_affine_10_0*tmp_14 + p_affine_10_1*tmp_11;
      real_t tmp_30 = tmp_22*(3*tmp_15*tmp_19*tmp_3 - tmp_15*tmp_20 - tmp_19*tmp_29) + tmp_28*(-tmp_20*tmp_25 + 3*tmp_25*tmp_27*tmp_3 - tmp_27*tmp_29);
      real_t tmp_31 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_16;
      real_t tmp_32 = tmp_22*(3*tmp_18*tmp_19*tmp_3 - tmp_18*tmp_20 - tmp_19*tmp_31) + tmp_28*(-tmp_20*tmp_26 + 3*tmp_26*tmp_27*tmp_3 - tmp_27*tmp_31);
      real_t tmp_33 = 2*tmp_29;
      real_t tmp_34 = tmp_22*(3*tmp_15*tmp_18*tmp_3 - tmp_15*tmp_31 - tmp_18*tmp_29) + tmp_28*(3*tmp_25*tmp_26*tmp_3 - tmp_25*tmp_31 - tmp_26*tmp_29);
      real_t tmp_35 = 2*tmp_31;
      real_t a_0_0 = tmp_22*(3*(tmp_19*tmp_19)*tmp_3 - tmp_19*tmp_21) + tmp_28*(-tmp_21*tmp_27 + 3*(tmp_27*tmp_27)*tmp_3);
      real_t a_0_1 = tmp_30;
      real_t a_0_2 = tmp_32;
      real_t a_1_0 = tmp_30;
      real_t a_1_1 = tmp_22*(3*(tmp_15*tmp_15)*tmp_3 - tmp_15*tmp_33) + tmp_28*(3*(tmp_25*tmp_25)*tmp_3 - tmp_25*tmp_33);
      real_t a_1_2 = tmp_34;
      real_t a_2_0 = tmp_32;
      real_t a_2_1 = tmp_34;
      real_t a_2_2 = tmp_22*(3*(tmp_18*tmp_18)*tmp_3 - tmp_18*tmp_35) + tmp_28*(3*(tmp_26*tmp_26)*tmp_3 - tmp_26*tmp_35);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_g0( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id0 );
      Scalar_Variable_Coefficient_2D_g0( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id1 );
      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_5 = 0.21132486540518713*tmp_1 + tmp_4;
      real_t tmp_6 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_7 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_8 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = 1.0 / (-tmp_6*tmp_9 + tmp_7*tmp_8);
      real_t tmp_11 = tmp_10*tmp_6;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_0 + tmp_12;
      real_t tmp_14 = tmp_10*tmp_8;
      real_t tmp_15 = tmp_11*tmp_5 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_10*tmp_7;
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_5;
      real_t tmp_19 = p_affine_10_0*(-tmp_14 - tmp_17) + p_affine_10_1*(-tmp_11 - tmp_16);
      real_t tmp_20 = 0.5*Scalar_Variable_Coefficient_2D_g0_out0_id0*tmp_2;
      real_t tmp_21 = 0.78867513459481287*tmp_1 + tmp_4;
      real_t tmp_22 = 0.78867513459481287*tmp_0 + tmp_12;
      real_t tmp_23 = tmp_11*tmp_21 + tmp_14*tmp_22;
      real_t tmp_24 = tmp_16*tmp_21 + tmp_17*tmp_22;
      real_t tmp_25 = 0.5*Scalar_Variable_Coefficient_2D_g0_out0_id1*tmp_2;
      real_t tmp_26 = p_affine_10_0*tmp_14 + p_affine_10_1*tmp_11;
      real_t tmp_27 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_16;
      real_t a_0_0 = tmp_20*(-tmp_19 + 3*tmp_3*(-tmp_15 - tmp_18 + 1)) + tmp_25*(-tmp_19 + 3*tmp_3*(-tmp_23 - tmp_24 + 1));
      real_t a_1_0 = tmp_20*(3*tmp_15*tmp_3 - tmp_26) + tmp_25*(3*tmp_23*tmp_3 - tmp_26);
      real_t a_2_0 = tmp_20*(3*tmp_18*tmp_3 - tmp_27) + tmp_25*(3*tmp_24*tmp_3 - tmp_27);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_g0( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id0 );
      Scalar_Variable_Coefficient_3D_g0( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id1 );
      Scalar_Variable_Coefficient_3D_g0( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id2 );
      Scalar_Variable_Coefficient_3D_g0( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id3 );
      Scalar_Variable_Coefficient_3D_g0( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id4 );
      Scalar_Variable_Coefficient_3D_g0( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id5 );
      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -tmp_4;
      real_t tmp_9 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_10 = 0.81684757298045851*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_11 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_13 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_14 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_15 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_16 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_19 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_20 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_12*tmp_20;
      real_t tmp_23 = tmp_16*tmp_19;
      real_t tmp_24 = tmp_14*tmp_18;
      real_t tmp_25 = 1.0 / (tmp_11*tmp_12*tmp_18 - tmp_11*tmp_23 + tmp_13*tmp_21 - tmp_13*tmp_24 + tmp_15*tmp_17 - tmp_15*tmp_22);
      real_t tmp_26 = tmp_25*(tmp_11*tmp_12 - tmp_13*tmp_14);
      real_t tmp_27 = -tmp_1;
      real_t tmp_28 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_29 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_30 = tmp_25*(-tmp_11*tmp_16 + tmp_13*tmp_20);
      real_t tmp_31 = -tmp_3;
      real_t tmp_32 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_33 = 0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_34 = tmp_25*(tmp_17 - tmp_22);
      real_t tmp_35 = tmp_10*tmp_26 + tmp_29*tmp_30 + tmp_33*tmp_34;
      real_t tmp_36 = tmp_25*(-tmp_12*tmp_15 + tmp_13*tmp_19);
      real_t tmp_37 = tmp_25*(-tmp_13*tmp_18 + tmp_15*tmp_16);
      real_t tmp_38 = tmp_25*(tmp_12*tmp_18 - tmp_23);
      real_t tmp_39 = tmp_10*tmp_36 + tmp_29*tmp_37 + tmp_33*tmp_38;
      real_t tmp_40 = tmp_25*(-tmp_11*tmp_19 + tmp_14*tmp_15);
      real_t tmp_41 = tmp_25*(tmp_11*tmp_18 - tmp_15*tmp_20);
      real_t tmp_42 = tmp_25*(tmp_21 - tmp_24);
      real_t tmp_43 = tmp_10*tmp_40 + tmp_29*tmp_41 + tmp_33*tmp_42;
      real_t tmp_44 = p_affine_13_0*(-tmp_34 - tmp_38 - tmp_42) + p_affine_13_1*(-tmp_30 - tmp_37 - tmp_41) + p_affine_13_2*(-tmp_26 - tmp_36 - tmp_40);
      real_t tmp_45 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_46 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_g0_out0_id0*tmp_45;
      real_t tmp_47 = 0.10810301816807022*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_48 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_49 = 0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_50 = tmp_26*tmp_47 + tmp_30*tmp_48 + tmp_34*tmp_49;
      real_t tmp_51 = tmp_36*tmp_47 + tmp_37*tmp_48 + tmp_38*tmp_49;
      real_t tmp_52 = tmp_40*tmp_47 + tmp_41*tmp_48 + tmp_42*tmp_49;
      real_t tmp_53 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_g0_out0_id1*tmp_45;
      real_t tmp_54 = 0.091576213509770743*tmp_5 + 0.81684757298045851*tmp_8 + tmp_9;
      real_t tmp_55 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_56 = 0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_31 + tmp_32;
      real_t tmp_57 = tmp_26*tmp_54 + tmp_30*tmp_55 + tmp_34*tmp_56;
      real_t tmp_58 = tmp_36*tmp_54 + tmp_37*tmp_55 + tmp_38*tmp_56;
      real_t tmp_59 = tmp_40*tmp_54 + tmp_41*tmp_55 + tmp_42*tmp_56;
      real_t tmp_60 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_g0_out0_id2*tmp_45;
      real_t tmp_61 = 0.44594849091596489*tmp_5 + 0.10810301816807022*tmp_8 + tmp_9;
      real_t tmp_62 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_63 = 0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_31 + tmp_32;
      real_t tmp_64 = tmp_26*tmp_61 + tmp_30*tmp_62 + tmp_34*tmp_63;
      real_t tmp_65 = tmp_36*tmp_61 + tmp_37*tmp_62 + tmp_38*tmp_63;
      real_t tmp_66 = tmp_40*tmp_61 + tmp_41*tmp_62 + tmp_42*tmp_63;
      real_t tmp_67 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_g0_out0_id3*tmp_45;
      real_t tmp_68 = 0.091576213509770743*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_69 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_70 = 0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_71 = tmp_26*tmp_68 + tmp_30*tmp_69 + tmp_34*tmp_70;
      real_t tmp_72 = tmp_36*tmp_68 + tmp_37*tmp_69 + tmp_38*tmp_70;
      real_t tmp_73 = tmp_40*tmp_68 + tmp_41*tmp_69 + tmp_42*tmp_70;
      real_t tmp_74 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_g0_out0_id4*tmp_45;
      real_t tmp_75 = 0.44594849091596489*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_76 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_77 = 0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_78 = tmp_26*tmp_75 + tmp_30*tmp_76 + tmp_34*tmp_77;
      real_t tmp_79 = tmp_36*tmp_75 + tmp_37*tmp_76 + tmp_38*tmp_77;
      real_t tmp_80 = tmp_40*tmp_75 + tmp_41*tmp_76 + tmp_42*tmp_77;
      real_t tmp_81 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_g0_out0_id5*tmp_45;
      real_t tmp_82 = p_affine_13_0*tmp_34 + p_affine_13_1*tmp_30 + p_affine_13_2*tmp_26;
      real_t tmp_83 = p_affine_13_0*tmp_38 + p_affine_13_1*tmp_37 + p_affine_13_2*tmp_36;
      real_t tmp_84 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_40;
      real_t a_0_0 = tmp_46*(-tmp_44 + 3.0*tmp_7*(-tmp_35 - tmp_39 - tmp_43 + 1)) + tmp_53*(-tmp_44 + 3.0*tmp_7*(-tmp_50 - tmp_51 - tmp_52 + 1)) + tmp_60*(-tmp_44 + 3.0*tmp_7*(-tmp_57 - tmp_58 - tmp_59 + 1)) + tmp_67*(-tmp_44 + 3.0*tmp_7*(-tmp_64 - tmp_65 - tmp_66 + 1)) + tmp_74*(-tmp_44 + 3.0*tmp_7*(-tmp_71 - tmp_72 - tmp_73 + 1)) + tmp_81*(-tmp_44 + 3.0*tmp_7*(-tmp_78 - tmp_79 - tmp_80 + 1));
      real_t a_1_0 = tmp_46*(3.0*tmp_35*tmp_7 - tmp_82) + tmp_53*(3.0*tmp_50*tmp_7 - tmp_82) + tmp_60*(3.0*tmp_57*tmp_7 - tmp_82) + tmp_67*(3.0*tmp_64*tmp_7 - tmp_82) + tmp_74*(3.0*tmp_7*tmp_71 - tmp_82) + tmp_81*(3.0*tmp_7*tmp_78 - tmp_82);
      real_t a_2_0 = tmp_46*(3.0*tmp_39*tmp_7 - tmp_83) + tmp_53*(3.0*tmp_51*tmp_7 - tmp_83) + tmp_60*(3.0*tmp_58*tmp_7 - tmp_83) + tmp_67*(3.0*tmp_65*tmp_7 - tmp_83) + tmp_74*(3.0*tmp_7*tmp_72 - tmp_83) + tmp_81*(3.0*tmp_7*tmp_79 - tmp_83);
      real_t a_3_0 = tmp_46*(3.0*tmp_43*tmp_7 - tmp_84) + tmp_53*(3.0*tmp_52*tmp_7 - tmp_84) + tmp_60*(3.0*tmp_59*tmp_7 - tmp_84) + tmp_67*(3.0*tmp_66*tmp_7 - tmp_84) + tmp_74*(3.0*tmp_7*tmp_73 - tmp_84) + tmp_81*(3.0*tmp_7*tmp_80 - tmp_84);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = tmp_2 - tmp_5;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_10 = tmp_3*tmp_9;
      real_t tmp_11 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_12 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_1*tmp_12;
      real_t tmp_15 = tmp_10*tmp_8 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_13 - tmp_14*tmp_8 + tmp_2*tmp_7 - tmp_5*tmp_7;
      real_t tmp_16 = 1.0 / (tmp_15);
      real_t tmp_17 = tmp_16*tmp_6;
      real_t tmp_18 = tmp_12*tmp_4 - tmp_13;
      real_t tmp_19 = tmp_16*tmp_18;
      real_t tmp_20 = tmp_10 - tmp_14;
      real_t tmp_21 = tmp_16*tmp_20;
      real_t tmp_22 = -tmp_17 - tmp_19 - tmp_21;
      real_t tmp_23 = -tmp_0*tmp_11 + tmp_3*tmp_8;
      real_t tmp_24 = tmp_16*tmp_23;
      real_t tmp_25 = tmp_0*tmp_7 - tmp_12*tmp_8;
      real_t tmp_26 = tmp_16*tmp_25;
      real_t tmp_27 = tmp_11*tmp_12 - tmp_3*tmp_7;
      real_t tmp_28 = tmp_16*tmp_27;
      real_t tmp_29 = -tmp_24 - tmp_26 - tmp_28;
      real_t tmp_30 = -tmp_1*tmp_8 + tmp_11*tmp_4;
      real_t tmp_31 = tmp_16*tmp_30;
      real_t tmp_32 = -tmp_4*tmp_7 + tmp_8*tmp_9;
      real_t tmp_33 = tmp_16*tmp_32;
      real_t tmp_34 = tmp_1*tmp_7 - tmp_11*tmp_9;
      real_t tmp_35 = tmp_16*tmp_34;
      real_t tmp_36 = -tmp_31 - tmp_33 - tmp_35;
      real_t tmp_37 = p_affine_0_0*p_affine_1_1;
      real_t tmp_38 = p_affine_0_0*p_affine_1_2;
      real_t tmp_39 = p_affine_2_1*p_affine_3_2;
      real_t tmp_40 = p_affine_0_1*p_affine_1_0;
      real_t tmp_41 = p_affine_0_1*p_affine_1_2;
      real_t tmp_42 = p_affine_2_2*p_affine_3_0;
      real_t tmp_43 = p_affine_0_2*p_affine_1_0;
      real_t tmp_44 = p_affine_0_2*p_affine_1_1;
      real_t tmp_45 = p_affine_2_0*p_affine_3_1;
      real_t tmp_46 = p_affine_2_2*p_affine_3_1;
      real_t tmp_47 = p_affine_2_0*p_affine_3_2;
      real_t tmp_48 = p_affine_2_1*p_affine_3_0;
      real_t tmp_49 = std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = tmp_49*((tmp_22*tmp_22) + (tmp_29*tmp_29) + (tmp_36*tmp_36));
      real_t tmp_51 = tmp_49*(tmp_21*tmp_22 + tmp_28*tmp_29 + tmp_35*tmp_36);
      real_t tmp_52 = 0.16666666666666666*tmp_51;
      real_t tmp_53 = tmp_49*(tmp_19*tmp_22 + tmp_26*tmp_29 + tmp_33*tmp_36);
      real_t tmp_54 = 0.16666666666666666*tmp_53;
      real_t tmp_55 = tmp_49*(tmp_17*tmp_22 + tmp_24*tmp_29 + tmp_31*tmp_36);
      real_t tmp_56 = 0.16666666666666666*tmp_55;
      real_t tmp_57 = 1.0 / (tmp_15*tmp_15);
      real_t tmp_58 = tmp_49*((tmp_20*tmp_20)*tmp_57 + (tmp_27*tmp_27)*tmp_57 + (tmp_34*tmp_34)*tmp_57);
      real_t tmp_59 = tmp_20*tmp_57;
      real_t tmp_60 = tmp_27*tmp_57;
      real_t tmp_61 = tmp_34*tmp_57;
      real_t tmp_62 = tmp_49*(tmp_18*tmp_59 + tmp_25*tmp_60 + tmp_32*tmp_61);
      real_t tmp_63 = 0.16666666666666666*tmp_62;
      real_t tmp_64 = tmp_49*(tmp_23*tmp_60 + tmp_30*tmp_61 + tmp_59*tmp_6);
      real_t tmp_65 = 0.16666666666666666*tmp_64;
      real_t tmp_66 = tmp_49*((tmp_18*tmp_18)*tmp_57 + (tmp_25*tmp_25)*tmp_57 + (tmp_32*tmp_32)*tmp_57);
      real_t tmp_67 = tmp_49*(tmp_18*tmp_57*tmp_6 + tmp_23*tmp_25*tmp_57 + tmp_30*tmp_32*tmp_57);
      real_t tmp_68 = 0.16666666666666666*tmp_67;
      real_t tmp_69 = tmp_49*((tmp_23*tmp_23)*tmp_57 + (tmp_30*tmp_30)*tmp_57 + tmp_57*(tmp_6*tmp_6));
      real_t a_0_0 = 0.16666666666666666*tmp_50;
      real_t a_0_1 = tmp_52;
      real_t a_0_2 = tmp_54;
      real_t a_0_3 = tmp_56;
      real_t a_1_0 = tmp_52;
      real_t a_1_1 = 0.16666666666666666*tmp_58;
      real_t a_1_2 = tmp_63;
      real_t a_1_3 = tmp_65;
      real_t a_2_0 = tmp_54;
      real_t a_2_1 = tmp_63;
      real_t a_2_2 = 0.16666666666666666*tmp_66;
      real_t a_2_3 = tmp_68;
      real_t a_3_0 = tmp_56;
      real_t a_3_1 = tmp_65;
      real_t a_3_2 = tmp_68;
      real_t a_3_3 = 0.16666666666666666*tmp_69;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const override
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

         real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -tmp_4;
      real_t tmp_9 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_10 = 0.81684757298045851*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_11 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_13 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_14 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_15 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_16 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_19 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_20 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_12*tmp_16;
      real_t tmp_23 = tmp_11*tmp_19;
      real_t tmp_24 = tmp_13*tmp_18;
      real_t tmp_25 = 1.0 / (tmp_11*tmp_12*tmp_18 + tmp_13*tmp_21 - tmp_14*tmp_24 + tmp_15*tmp_17 - tmp_15*tmp_23 - tmp_20*tmp_22);
      real_t tmp_26 = tmp_25*(tmp_11*tmp_12 - tmp_13*tmp_14);
      real_t tmp_27 = -tmp_1;
      real_t tmp_28 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_29 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_30 = tmp_25*(-tmp_11*tmp_15 + tmp_13*tmp_20);
      real_t tmp_31 = -tmp_3;
      real_t tmp_32 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_33 = 0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_34 = tmp_25*(-tmp_12*tmp_20 + tmp_14*tmp_15);
      real_t tmp_35 = tmp_10*tmp_26 + tmp_29*tmp_30 + tmp_33*tmp_34;
      real_t tmp_36 = tmp_25*(tmp_13*tmp_19 - tmp_22);
      real_t tmp_37 = tmp_25*(tmp_15*tmp_16 - tmp_24);
      real_t tmp_38 = tmp_25*(tmp_12*tmp_18 - tmp_15*tmp_19);
      real_t tmp_39 = tmp_10*tmp_36 + tmp_29*tmp_37 + tmp_33*tmp_38;
      real_t tmp_40 = tmp_25*(tmp_17 - tmp_23);
      real_t tmp_41 = tmp_25*(tmp_11*tmp_18 - tmp_16*tmp_20);
      real_t tmp_42 = tmp_25*(-tmp_14*tmp_18 + tmp_21);
      real_t tmp_43 = tmp_10*tmp_40 + tmp_29*tmp_41 + tmp_33*tmp_42;
      real_t tmp_44 = -tmp_35 - tmp_39 - tmp_43 + 1;
      real_t tmp_45 = p_affine_13_0*(-tmp_34 - tmp_38 - tmp_42) + p_affine_13_1*(-tmp_30 - tmp_37 - tmp_41) + p_affine_13_2*(-tmp_26 - tmp_36 - tmp_40);
      real_t tmp_46 = 1.0*tmp_45;
      real_t tmp_47 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_48 = 0.054975871827660928*tmp_47;
      real_t tmp_49 = 0.10810301816807022*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_50 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_51 = 0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_52 = tmp_26*tmp_49 + tmp_30*tmp_50 + tmp_34*tmp_51;
      real_t tmp_53 = tmp_36*tmp_49 + tmp_37*tmp_50 + tmp_38*tmp_51;
      real_t tmp_54 = tmp_40*tmp_49 + tmp_41*tmp_50 + tmp_42*tmp_51;
      real_t tmp_55 = -tmp_52 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11169079483900572*tmp_47;
      real_t tmp_57 = 0.091576213509770743*tmp_5 + 0.81684757298045851*tmp_8 + tmp_9;
      real_t tmp_58 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_59 = 0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_31 + tmp_32;
      real_t tmp_60 = tmp_26*tmp_57 + tmp_30*tmp_58 + tmp_34*tmp_59;
      real_t tmp_61 = tmp_36*tmp_57 + tmp_37*tmp_58 + tmp_38*tmp_59;
      real_t tmp_62 = tmp_40*tmp_57 + tmp_41*tmp_58 + tmp_42*tmp_59;
      real_t tmp_63 = -tmp_60 - tmp_61 - tmp_62 + 1;
      real_t tmp_64 = 0.054975871827660928*tmp_47;
      real_t tmp_65 = 0.44594849091596489*tmp_5 + 0.10810301816807022*tmp_8 + tmp_9;
      real_t tmp_66 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_67 = 0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_31 + tmp_32;
      real_t tmp_68 = tmp_26*tmp_65 + tmp_30*tmp_66 + tmp_34*tmp_67;
      real_t tmp_69 = tmp_36*tmp_65 + tmp_37*tmp_66 + tmp_38*tmp_67;
      real_t tmp_70 = tmp_40*tmp_65 + tmp_41*tmp_66 + tmp_42*tmp_67;
      real_t tmp_71 = -tmp_68 - tmp_69 - tmp_70 + 1;
      real_t tmp_72 = 0.11169079483900572*tmp_47;
      real_t tmp_73 = 0.091576213509770743*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_74 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_75 = 0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_76 = tmp_26*tmp_73 + tmp_30*tmp_74 + tmp_34*tmp_75;
      real_t tmp_77 = tmp_36*tmp_73 + tmp_37*tmp_74 + tmp_38*tmp_75;
      real_t tmp_78 = tmp_40*tmp_73 + tmp_41*tmp_74 + tmp_42*tmp_75;
      real_t tmp_79 = -tmp_76 - tmp_77 - tmp_78 + 1;
      real_t tmp_80 = 0.054975871827660928*tmp_47;
      real_t tmp_81 = 0.44594849091596489*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_82 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_83 = 0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_84 = tmp_26*tmp_81 + tmp_30*tmp_82 + tmp_34*tmp_83;
      real_t tmp_85 = tmp_36*tmp_81 + tmp_37*tmp_82 + tmp_38*tmp_83;
      real_t tmp_86 = tmp_40*tmp_81 + tmp_41*tmp_82 + tmp_42*tmp_83;
      real_t tmp_87 = -tmp_84 - tmp_85 - tmp_86 + 1;
      real_t tmp_88 = 0.11169079483900572*tmp_47;
      real_t tmp_89 = 0.5*tmp_45;
      real_t tmp_90 = p_affine_13_0*tmp_34 + p_affine_13_1*tmp_30 + p_affine_13_2*tmp_26;
      real_t tmp_91 = 0.5*tmp_90;
      real_t tmp_92 = tmp_48*(3.0*tmp_35*tmp_44*tmp_7 - tmp_35*tmp_89 - tmp_44*tmp_91) + tmp_56*(3.0*tmp_52*tmp_55*tmp_7 - tmp_52*tmp_89 - tmp_55*tmp_91) + tmp_64*(3.0*tmp_60*tmp_63*tmp_7 - tmp_60*tmp_89 - tmp_63*tmp_91) + tmp_72*(3.0*tmp_68*tmp_7*tmp_71 - tmp_68*tmp_89 - tmp_71*tmp_91) + tmp_80*(3.0*tmp_7*tmp_76*tmp_79 - tmp_76*tmp_89 - tmp_79*tmp_91) + tmp_88*(3.0*tmp_7*tmp_84*tmp_87 - tmp_84*tmp_89 - tmp_87*tmp_91);
      real_t tmp_93 = p_affine_13_0*tmp_38 + p_affine_13_1*tmp_37 + p_affine_13_2*tmp_36;
      real_t tmp_94 = 0.5*tmp_93;
      real_t tmp_95 = tmp_48*(3.0*tmp_39*tmp_44*tmp_7 - tmp_39*tmp_89 - tmp_44*tmp_94) + tmp_56*(3.0*tmp_53*tmp_55*tmp_7 - tmp_53*tmp_89 - tmp_55*tmp_94) + tmp_64*(3.0*tmp_61*tmp_63*tmp_7 - tmp_61*tmp_89 - tmp_63*tmp_94) + tmp_72*(3.0*tmp_69*tmp_7*tmp_71 - tmp_69*tmp_89 - tmp_71*tmp_94) + tmp_80*(3.0*tmp_7*tmp_77*tmp_79 - tmp_77*tmp_89 - tmp_79*tmp_94) + tmp_88*(3.0*tmp_7*tmp_85*tmp_87 - tmp_85*tmp_89 - tmp_87*tmp_94);
      real_t tmp_96 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_40;
      real_t tmp_97 = 0.5*tmp_96;
      real_t tmp_98 = tmp_48*(3.0*tmp_43*tmp_44*tmp_7 - tmp_43*tmp_89 - tmp_44*tmp_97) + tmp_56*(3.0*tmp_54*tmp_55*tmp_7 - tmp_54*tmp_89 - tmp_55*tmp_97) + tmp_64*(3.0*tmp_62*tmp_63*tmp_7 - tmp_62*tmp_89 - tmp_63*tmp_97) + tmp_72*(3.0*tmp_7*tmp_70*tmp_71 - tmp_70*tmp_89 - tmp_71*tmp_97) + tmp_80*(3.0*tmp_7*tmp_78*tmp_79 - tmp_78*tmp_89 - tmp_79*tmp_97) + tmp_88*(3.0*tmp_7*tmp_86*tmp_87 - tmp_86*tmp_89 - tmp_87*tmp_97);
      real_t tmp_99 = 1.0*tmp_90;
      real_t tmp_100 = tmp_48*(3.0*tmp_35*tmp_39*tmp_7 - tmp_35*tmp_94 - tmp_39*tmp_91) + tmp_56*(3.0*tmp_52*tmp_53*tmp_7 - tmp_52*tmp_94 - tmp_53*tmp_91) + tmp_64*(3.0*tmp_60*tmp_61*tmp_7 - tmp_60*tmp_94 - tmp_61*tmp_91) + tmp_72*(3.0*tmp_68*tmp_69*tmp_7 - tmp_68*tmp_94 - tmp_69*tmp_91) + tmp_80*(3.0*tmp_7*tmp_76*tmp_77 - tmp_76*tmp_94 - tmp_77*tmp_91) + tmp_88*(3.0*tmp_7*tmp_84*tmp_85 - tmp_84*tmp_94 - tmp_85*tmp_91);
      real_t tmp_101 = tmp_48*(3.0*tmp_35*tmp_43*tmp_7 - tmp_35*tmp_97 - tmp_43*tmp_91) + tmp_56*(3.0*tmp_52*tmp_54*tmp_7 - tmp_52*tmp_97 - tmp_54*tmp_91) + tmp_64*(3.0*tmp_60*tmp_62*tmp_7 - tmp_60*tmp_97 - tmp_62*tmp_91) + tmp_72*(3.0*tmp_68*tmp_7*tmp_70 - tmp_68*tmp_97 - tmp_70*tmp_91) + tmp_80*(3.0*tmp_7*tmp_76*tmp_78 - tmp_76*tmp_97 - tmp_78*tmp_91) + tmp_88*(3.0*tmp_7*tmp_84*tmp_86 - tmp_84*tmp_97 - tmp_86*tmp_91);
      real_t tmp_102 = 1.0*tmp_93;
      real_t tmp_103 = tmp_48*(3.0*tmp_39*tmp_43*tmp_7 - tmp_39*tmp_97 - tmp_43*tmp_94) + tmp_56*(3.0*tmp_53*tmp_54*tmp_7 - tmp_53*tmp_97 - tmp_54*tmp_94) + tmp_64*(3.0*tmp_61*tmp_62*tmp_7 - tmp_61*tmp_97 - tmp_62*tmp_94) + tmp_72*(3.0*tmp_69*tmp_7*tmp_70 - tmp_69*tmp_97 - tmp_70*tmp_94) + tmp_80*(3.0*tmp_7*tmp_77*tmp_78 - tmp_77*tmp_97 - tmp_78*tmp_94) + tmp_88*(3.0*tmp_7*tmp_85*tmp_86 - tmp_85*tmp_97 - tmp_86*tmp_94);
      real_t tmp_104 = 1.0*tmp_96;
      real_t a_0_0 = tmp_48*(3.0*(tmp_44*tmp_44)*tmp_7 - tmp_44*tmp_46) + tmp_56*(-tmp_46*tmp_55 + 3.0*(tmp_55*tmp_55)*tmp_7) + tmp_64*(-tmp_46*tmp_63 + 3.0*(tmp_63*tmp_63)*tmp_7) + tmp_72*(-tmp_46*tmp_71 + 3.0*tmp_7*(tmp_71*tmp_71)) + tmp_80*(-tmp_46*tmp_79 + 3.0*tmp_7*(tmp_79*tmp_79)) + tmp_88*(-tmp_46*tmp_87 + 3.0*tmp_7*(tmp_87*tmp_87));
      real_t a_0_1 = tmp_92;
      real_t a_0_2 = tmp_95;
      real_t a_0_3 = tmp_98;
      real_t a_1_0 = tmp_92;
      real_t a_1_1 = tmp_48*(3.0*(tmp_35*tmp_35)*tmp_7 - tmp_35*tmp_99) + tmp_56*(3.0*(tmp_52*tmp_52)*tmp_7 - tmp_52*tmp_99) + tmp_64*(3.0*(tmp_60*tmp_60)*tmp_7 - tmp_60*tmp_99) + tmp_72*(3.0*(tmp_68*tmp_68)*tmp_7 - tmp_68*tmp_99) + tmp_80*(3.0*tmp_7*(tmp_76*tmp_76) - tmp_76*tmp_99) + tmp_88*(3.0*tmp_7*(tmp_84*tmp_84) - tmp_84*tmp_99);
      real_t a_1_2 = tmp_100;
      real_t a_1_3 = tmp_101;
      real_t a_2_0 = tmp_95;
      real_t a_2_1 = tmp_100;
      real_t a_2_2 = tmp_48*(-tmp_102*tmp_39 + 3.0*(tmp_39*tmp_39)*tmp_7) + tmp_56*(-tmp_102*tmp_53 + 3.0*(tmp_53*tmp_53)*tmp_7) + tmp_64*(-tmp_102*tmp_61 + 3.0*(tmp_61*tmp_61)*tmp_7) + tmp_72*(-tmp_102*tmp_69 + 3.0*(tmp_69*tmp_69)*tmp_7) + tmp_80*(-tmp_102*tmp_77 + 3.0*tmp_7*(tmp_77*tmp_77)) + tmp_88*(-tmp_102*tmp_85 + 3.0*tmp_7*(tmp_85*tmp_85));
      real_t a_2_3 = tmp_103;
      real_t a_3_0 = tmp_98;
      real_t a_3_1 = tmp_101;
      real_t a_3_2 = tmp_103;
      real_t a_3_3 = tmp_48*(-tmp_104*tmp_43 + 3.0*(tmp_43*tmp_43)*tmp_7) + tmp_56*(-tmp_104*tmp_54 + 3.0*(tmp_54*tmp_54)*tmp_7) + tmp_64*(-tmp_104*tmp_62 + 3.0*(tmp_62*tmp_62)*tmp_7) + tmp_72*(-tmp_104*tmp_70 + 3.0*tmp_7*(tmp_70*tmp_70)) + tmp_80*(-tmp_104*tmp_78 + 3.0*tmp_7*(tmp_78*tmp_78)) + tmp_88*(-tmp_104*tmp_86 + 3.0*tmp_7*(tmp_86*tmp_86));
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
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
                                  MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_26 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_27 = -tmp_26;
      real_t tmp_28 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_29 = 0.091576213509770743*tmp_27 + 0.81684757298045851*tmp_28;
      real_t tmp_30 = tmp_25 + tmp_29;
      real_t tmp_31 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_32 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_33 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_34 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_35 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_36 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_39 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_40 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_41 = tmp_39*tmp_40;
      real_t tmp_42 = tmp_32*tmp_36;
      real_t tmp_43 = tmp_31*tmp_39;
      real_t tmp_44 = tmp_33*tmp_38;
      real_t tmp_45 = 1.0 / (tmp_31*tmp_32*tmp_38 + tmp_33*tmp_41 - tmp_34*tmp_44 + tmp_35*tmp_37 - tmp_35*tmp_43 - tmp_40*tmp_42);
      real_t tmp_46 = tmp_45*(tmp_31*tmp_32 - tmp_33*tmp_34);
      real_t tmp_47 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_48 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_49 = -tmp_48;
      real_t tmp_50 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_51 = 0.091576213509770743*tmp_49 + 0.81684757298045851*tmp_50;
      real_t tmp_52 = tmp_47 + tmp_51;
      real_t tmp_53 = tmp_45*(-tmp_31*tmp_35 + tmp_33*tmp_40);
      real_t tmp_54 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_55 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_56 = -tmp_55;
      real_t tmp_57 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_58 = 0.091576213509770743*tmp_56 + 0.81684757298045851*tmp_57;
      real_t tmp_59 = tmp_54 + tmp_58;
      real_t tmp_60 = tmp_45*(-tmp_32*tmp_40 + tmp_34*tmp_35);
      real_t tmp_61 = tmp_30*tmp_46 + tmp_52*tmp_53 + tmp_59*tmp_60;
      real_t tmp_62 = tmp_45*(tmp_33*tmp_39 - tmp_42);
      real_t tmp_63 = tmp_45*(tmp_35*tmp_36 - tmp_44);
      real_t tmp_64 = tmp_45*(tmp_32*tmp_38 - tmp_35*tmp_39);
      real_t tmp_65 = tmp_30*tmp_62 + tmp_52*tmp_63 + tmp_59*tmp_64;
      real_t tmp_66 = tmp_45*(tmp_37 - tmp_43);
      real_t tmp_67 = tmp_45*(tmp_31*tmp_38 - tmp_36*tmp_40);
      real_t tmp_68 = tmp_45*(-tmp_34*tmp_38 + tmp_41);
      real_t tmp_69 = tmp_30*tmp_66 + tmp_52*tmp_67 + tmp_59*tmp_68;
      real_t tmp_70 = -tmp_61 - tmp_65 - tmp_69 + 1;
      real_t tmp_71 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_72 = tmp_29 + tmp_71;
      real_t tmp_73 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_74 = tmp_51 + tmp_73;
      real_t tmp_75 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_76 = tmp_58 + tmp_75;
      real_t tmp_77 = tmp_17*tmp_76 + tmp_20*tmp_74 + tmp_23*tmp_72;
      real_t tmp_78 = tmp_16*tmp_76 + tmp_19*tmp_74 + tmp_22*tmp_72;
      real_t tmp_79 = tmp_15*tmp_76 + tmp_18*tmp_74 + tmp_21*tmp_72;
      real_t tmp_80 = -tmp_77 - tmp_78 - tmp_79 + 1;
      real_t tmp_81 = 0.5*p_affine_13_0*(-tmp_60 - tmp_64 - tmp_68) + 0.5*p_affine_13_1*(-tmp_53 - tmp_63 - tmp_67) + 0.5*p_affine_13_2*(-tmp_46 - tmp_62 - tmp_66);
      real_t tmp_82 = (std::abs(tmp_26*tmp_50 - tmp_28*tmp_48)*std::abs(tmp_26*tmp_50 - tmp_28*tmp_48)) + (std::abs(tmp_26*tmp_57 - tmp_28*tmp_55)*std::abs(tmp_26*tmp_57 - tmp_28*tmp_55)) + (std::abs(tmp_48*tmp_57 - tmp_50*tmp_55)*std::abs(tmp_48*tmp_57 - tmp_50*tmp_55));
      real_t tmp_83 = 3.0*std::pow(tmp_82, -0.25);
      real_t tmp_84 = tmp_80*tmp_83;
      real_t tmp_85 = 1.0*std::pow(tmp_82, 1.0/2.0);
      real_t tmp_86 = 0.054975871827660928*tmp_85;
      real_t tmp_87 = 0.44594849091596489*tmp_27 + 0.10810301816807022*tmp_28;
      real_t tmp_88 = tmp_25 + tmp_87;
      real_t tmp_89 = 0.44594849091596489*tmp_49 + 0.10810301816807022*tmp_50;
      real_t tmp_90 = tmp_47 + tmp_89;
      real_t tmp_91 = 0.44594849091596489*tmp_56 + 0.10810301816807022*tmp_57;
      real_t tmp_92 = tmp_54 + tmp_91;
      real_t tmp_93 = tmp_46*tmp_88 + tmp_53*tmp_90 + tmp_60*tmp_92;
      real_t tmp_94 = tmp_62*tmp_88 + tmp_63*tmp_90 + tmp_64*tmp_92;
      real_t tmp_95 = tmp_66*tmp_88 + tmp_67*tmp_90 + tmp_68*tmp_92;
      real_t tmp_96 = -tmp_93 - tmp_94 - tmp_95 + 1;
      real_t tmp_97 = tmp_71 + tmp_87;
      real_t tmp_98 = tmp_73 + tmp_89;
      real_t tmp_99 = tmp_75 + tmp_91;
      real_t tmp_100 = tmp_17*tmp_99 + tmp_20*tmp_98 + tmp_23*tmp_97;
      real_t tmp_101 = tmp_16*tmp_99 + tmp_19*tmp_98 + tmp_22*tmp_97;
      real_t tmp_102 = tmp_15*tmp_99 + tmp_18*tmp_98 + tmp_21*tmp_97;
      real_t tmp_103 = -tmp_100 - tmp_101 - tmp_102 + 1;
      real_t tmp_104 = tmp_103*tmp_83;
      real_t tmp_105 = 0.11169079483900572*tmp_85;
      real_t tmp_106 = 0.81684757298045851*tmp_27 + 0.091576213509770743*tmp_28;
      real_t tmp_107 = tmp_106 + tmp_25;
      real_t tmp_108 = 0.81684757298045851*tmp_49 + 0.091576213509770743*tmp_50;
      real_t tmp_109 = tmp_108 + tmp_47;
      real_t tmp_110 = 0.81684757298045851*tmp_56 + 0.091576213509770743*tmp_57;
      real_t tmp_111 = tmp_110 + tmp_54;
      real_t tmp_112 = tmp_107*tmp_46 + tmp_109*tmp_53 + tmp_111*tmp_60;
      real_t tmp_113 = tmp_107*tmp_62 + tmp_109*tmp_63 + tmp_111*tmp_64;
      real_t tmp_114 = tmp_107*tmp_66 + tmp_109*tmp_67 + tmp_111*tmp_68;
      real_t tmp_115 = -tmp_112 - tmp_113 - tmp_114 + 1;
      real_t tmp_116 = tmp_106 + tmp_71;
      real_t tmp_117 = tmp_108 + tmp_73;
      real_t tmp_118 = tmp_110 + tmp_75;
      real_t tmp_119 = tmp_116*tmp_23 + tmp_117*tmp_20 + tmp_118*tmp_17;
      real_t tmp_120 = tmp_116*tmp_22 + tmp_117*tmp_19 + tmp_118*tmp_16;
      real_t tmp_121 = tmp_116*tmp_21 + tmp_117*tmp_18 + tmp_118*tmp_15;
      real_t tmp_122 = -tmp_119 - tmp_120 - tmp_121 + 1;
      real_t tmp_123 = tmp_122*tmp_83;
      real_t tmp_124 = 0.054975871827660928*tmp_85;
      real_t tmp_125 = 0.10810301816807022*tmp_27 + 0.44594849091596489*tmp_28;
      real_t tmp_126 = tmp_125 + tmp_25;
      real_t tmp_127 = 0.10810301816807022*tmp_49 + 0.44594849091596489*tmp_50;
      real_t tmp_128 = tmp_127 + tmp_47;
      real_t tmp_129 = 0.10810301816807022*tmp_56 + 0.44594849091596489*tmp_57;
      real_t tmp_130 = tmp_129 + tmp_54;
      real_t tmp_131 = tmp_126*tmp_46 + tmp_128*tmp_53 + tmp_130*tmp_60;
      real_t tmp_132 = tmp_126*tmp_62 + tmp_128*tmp_63 + tmp_130*tmp_64;
      real_t tmp_133 = tmp_126*tmp_66 + tmp_128*tmp_67 + tmp_130*tmp_68;
      real_t tmp_134 = -tmp_131 - tmp_132 - tmp_133 + 1;
      real_t tmp_135 = tmp_125 + tmp_71;
      real_t tmp_136 = tmp_127 + tmp_73;
      real_t tmp_137 = tmp_129 + tmp_75;
      real_t tmp_138 = tmp_135*tmp_23 + tmp_136*tmp_20 + tmp_137*tmp_17;
      real_t tmp_139 = tmp_135*tmp_22 + tmp_136*tmp_19 + tmp_137*tmp_16;
      real_t tmp_140 = tmp_135*tmp_21 + tmp_136*tmp_18 + tmp_137*tmp_15;
      real_t tmp_141 = -tmp_138 - tmp_139 - tmp_140 + 1;
      real_t tmp_142 = tmp_141*tmp_83;
      real_t tmp_143 = 0.11169079483900572*tmp_85;
      real_t tmp_144 = 0.091576213509770743*tmp_27 + 0.091576213509770743*tmp_28;
      real_t tmp_145 = tmp_144 + tmp_25;
      real_t tmp_146 = 0.091576213509770743*tmp_49 + 0.091576213509770743*tmp_50;
      real_t tmp_147 = tmp_146 + tmp_47;
      real_t tmp_148 = 0.091576213509770743*tmp_56 + 0.091576213509770743*tmp_57;
      real_t tmp_149 = tmp_148 + tmp_54;
      real_t tmp_150 = tmp_145*tmp_46 + tmp_147*tmp_53 + tmp_149*tmp_60;
      real_t tmp_151 = tmp_145*tmp_62 + tmp_147*tmp_63 + tmp_149*tmp_64;
      real_t tmp_152 = tmp_145*tmp_66 + tmp_147*tmp_67 + tmp_149*tmp_68;
      real_t tmp_153 = -tmp_150 - tmp_151 - tmp_152 + 1;
      real_t tmp_154 = tmp_144 + tmp_71;
      real_t tmp_155 = tmp_146 + tmp_73;
      real_t tmp_156 = tmp_148 + tmp_75;
      real_t tmp_157 = tmp_154*tmp_23 + tmp_155*tmp_20 + tmp_156*tmp_17;
      real_t tmp_158 = tmp_154*tmp_22 + tmp_155*tmp_19 + tmp_156*tmp_16;
      real_t tmp_159 = tmp_15*tmp_156 + tmp_154*tmp_21 + tmp_155*tmp_18;
      real_t tmp_160 = -tmp_157 - tmp_158 - tmp_159 + 1;
      real_t tmp_161 = tmp_160*tmp_83;
      real_t tmp_162 = 0.054975871827660928*tmp_85;
      real_t tmp_163 = 0.44594849091596489*tmp_27 + 0.44594849091596489*tmp_28;
      real_t tmp_164 = tmp_163 + tmp_25;
      real_t tmp_165 = 0.44594849091596489*tmp_49 + 0.44594849091596489*tmp_50;
      real_t tmp_166 = tmp_165 + tmp_47;
      real_t tmp_167 = 0.44594849091596489*tmp_56 + 0.44594849091596489*tmp_57;
      real_t tmp_168 = tmp_167 + tmp_54;
      real_t tmp_169 = tmp_164*tmp_46 + tmp_166*tmp_53 + tmp_168*tmp_60;
      real_t tmp_170 = tmp_164*tmp_62 + tmp_166*tmp_63 + tmp_168*tmp_64;
      real_t tmp_171 = tmp_164*tmp_66 + tmp_166*tmp_67 + tmp_168*tmp_68;
      real_t tmp_172 = -tmp_169 - tmp_170 - tmp_171 + 1;
      real_t tmp_173 = tmp_163 + tmp_71;
      real_t tmp_174 = tmp_165 + tmp_73;
      real_t tmp_175 = tmp_167 + tmp_75;
      real_t tmp_176 = tmp_17*tmp_175 + tmp_173*tmp_23 + tmp_174*tmp_20;
      real_t tmp_177 = tmp_16*tmp_175 + tmp_173*tmp_22 + tmp_174*tmp_19;
      real_t tmp_178 = tmp_15*tmp_175 + tmp_173*tmp_21 + tmp_174*tmp_18;
      real_t tmp_179 = -tmp_176 - tmp_177 - tmp_178 + 1;
      real_t tmp_180 = tmp_179*tmp_83;
      real_t tmp_181 = 0.11169079483900572*tmp_85;
      real_t tmp_182 = 0.5*p_affine_13_0*tmp_60 + 0.5*p_affine_13_1*tmp_53 + 0.5*p_affine_13_2*tmp_46;
      real_t tmp_183 = 0.5*p_affine_13_0*tmp_64 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_62;
      real_t tmp_184 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_67 + 0.5*p_affine_13_2*tmp_66;
      real_t tmp_185 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_186 = tmp_77*tmp_83;
      real_t tmp_187 = tmp_100*tmp_83;
      real_t tmp_188 = tmp_119*tmp_83;
      real_t tmp_189 = tmp_138*tmp_83;
      real_t tmp_190 = tmp_157*tmp_83;
      real_t tmp_191 = tmp_176*tmp_83;
      real_t tmp_192 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_193 = tmp_78*tmp_83;
      real_t tmp_194 = tmp_101*tmp_83;
      real_t tmp_195 = tmp_120*tmp_83;
      real_t tmp_196 = tmp_139*tmp_83;
      real_t tmp_197 = tmp_158*tmp_83;
      real_t tmp_198 = tmp_177*tmp_83;
      real_t tmp_199 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t tmp_200 = tmp_79*tmp_83;
      real_t tmp_201 = tmp_102*tmp_83;
      real_t tmp_202 = tmp_121*tmp_83;
      real_t tmp_203 = tmp_140*tmp_83;
      real_t tmp_204 = tmp_159*tmp_83;
      real_t tmp_205 = tmp_178*tmp_83;
      real_t a_0_0 = tmp_105*(-tmp_103*tmp_81 - tmp_104*tmp_96 + 0.5*tmp_24*tmp_96) + tmp_124*(-tmp_115*tmp_123 + 0.5*tmp_115*tmp_24 - tmp_122*tmp_81) + tmp_143*(-tmp_134*tmp_142 + 0.5*tmp_134*tmp_24 - tmp_141*tmp_81) + tmp_162*(-tmp_153*tmp_161 + 0.5*tmp_153*tmp_24 - tmp_160*tmp_81) + tmp_181*(-tmp_172*tmp_180 + 0.5*tmp_172*tmp_24 - tmp_179*tmp_81) + tmp_86*(0.5*tmp_24*tmp_70 - tmp_70*tmp_84 - tmp_80*tmp_81);
      real_t a_0_1 = tmp_105*(-tmp_103*tmp_182 - tmp_104*tmp_93 + 0.5*tmp_24*tmp_93) + tmp_124*(-tmp_112*tmp_123 + 0.5*tmp_112*tmp_24 - tmp_122*tmp_182) + tmp_143*(-tmp_131*tmp_142 + 0.5*tmp_131*tmp_24 - tmp_141*tmp_182) + tmp_162*(-tmp_150*tmp_161 + 0.5*tmp_150*tmp_24 - tmp_160*tmp_182) + tmp_181*(-tmp_169*tmp_180 + 0.5*tmp_169*tmp_24 - tmp_179*tmp_182) + tmp_86*(-tmp_182*tmp_80 + 0.5*tmp_24*tmp_61 - tmp_61*tmp_84);
      real_t a_0_2 = tmp_105*(-tmp_103*tmp_183 - tmp_104*tmp_94 + 0.5*tmp_24*tmp_94) + tmp_124*(-tmp_113*tmp_123 + 0.5*tmp_113*tmp_24 - tmp_122*tmp_183) + tmp_143*(-tmp_132*tmp_142 + 0.5*tmp_132*tmp_24 - tmp_141*tmp_183) + tmp_162*(-tmp_151*tmp_161 + 0.5*tmp_151*tmp_24 - tmp_160*tmp_183) + tmp_181*(-tmp_170*tmp_180 + 0.5*tmp_170*tmp_24 - tmp_179*tmp_183) + tmp_86*(-tmp_183*tmp_80 + 0.5*tmp_24*tmp_65 - tmp_65*tmp_84);
      real_t a_0_3 = tmp_105*(-tmp_103*tmp_184 - tmp_104*tmp_95 + 0.5*tmp_24*tmp_95) + tmp_124*(-tmp_114*tmp_123 + 0.5*tmp_114*tmp_24 - tmp_122*tmp_184) + tmp_143*(-tmp_133*tmp_142 + 0.5*tmp_133*tmp_24 - tmp_141*tmp_184) + tmp_162*(-tmp_152*tmp_161 + 0.5*tmp_152*tmp_24 - tmp_160*tmp_184) + tmp_181*(-tmp_171*tmp_180 + 0.5*tmp_171*tmp_24 - tmp_179*tmp_184) + tmp_86*(-tmp_184*tmp_80 + 0.5*tmp_24*tmp_69 - tmp_69*tmp_84);
      real_t a_1_0 = tmp_105*(-tmp_100*tmp_81 + 0.5*tmp_185*tmp_96 - tmp_187*tmp_96) + tmp_124*(0.5*tmp_115*tmp_185 - tmp_115*tmp_188 - tmp_119*tmp_81) + tmp_143*(0.5*tmp_134*tmp_185 - tmp_134*tmp_189 - tmp_138*tmp_81) + tmp_162*(0.5*tmp_153*tmp_185 - tmp_153*tmp_190 - tmp_157*tmp_81) + tmp_181*(0.5*tmp_172*tmp_185 - tmp_172*tmp_191 - tmp_176*tmp_81) + tmp_86*(0.5*tmp_185*tmp_70 - tmp_186*tmp_70 - tmp_77*tmp_81);
      real_t a_1_1 = tmp_105*(-tmp_100*tmp_182 + 0.5*tmp_185*tmp_93 - tmp_187*tmp_93) + tmp_124*(0.5*tmp_112*tmp_185 - tmp_112*tmp_188 - tmp_119*tmp_182) + tmp_143*(0.5*tmp_131*tmp_185 - tmp_131*tmp_189 - tmp_138*tmp_182) + tmp_162*(0.5*tmp_150*tmp_185 - tmp_150*tmp_190 - tmp_157*tmp_182) + tmp_181*(0.5*tmp_169*tmp_185 - tmp_169*tmp_191 - tmp_176*tmp_182) + tmp_86*(-tmp_182*tmp_77 + 0.5*tmp_185*tmp_61 - tmp_186*tmp_61);
      real_t a_1_2 = tmp_105*(-tmp_100*tmp_183 + 0.5*tmp_185*tmp_94 - tmp_187*tmp_94) + tmp_124*(0.5*tmp_113*tmp_185 - tmp_113*tmp_188 - tmp_119*tmp_183) + tmp_143*(0.5*tmp_132*tmp_185 - tmp_132*tmp_189 - tmp_138*tmp_183) + tmp_162*(0.5*tmp_151*tmp_185 - tmp_151*tmp_190 - tmp_157*tmp_183) + tmp_181*(0.5*tmp_170*tmp_185 - tmp_170*tmp_191 - tmp_176*tmp_183) + tmp_86*(-tmp_183*tmp_77 + 0.5*tmp_185*tmp_65 - tmp_186*tmp_65);
      real_t a_1_3 = tmp_105*(-tmp_100*tmp_184 + 0.5*tmp_185*tmp_95 - tmp_187*tmp_95) + tmp_124*(0.5*tmp_114*tmp_185 - tmp_114*tmp_188 - tmp_119*tmp_184) + tmp_143*(0.5*tmp_133*tmp_185 - tmp_133*tmp_189 - tmp_138*tmp_184) + tmp_162*(0.5*tmp_152*tmp_185 - tmp_152*tmp_190 - tmp_157*tmp_184) + tmp_181*(0.5*tmp_171*tmp_185 - tmp_171*tmp_191 - tmp_176*tmp_184) + tmp_86*(-tmp_184*tmp_77 + 0.5*tmp_185*tmp_69 - tmp_186*tmp_69);
      real_t a_2_0 = tmp_105*(-tmp_101*tmp_81 + 0.5*tmp_192*tmp_96 - tmp_194*tmp_96) + tmp_124*(0.5*tmp_115*tmp_192 - tmp_115*tmp_195 - tmp_120*tmp_81) + tmp_143*(0.5*tmp_134*tmp_192 - tmp_134*tmp_196 - tmp_139*tmp_81) + tmp_162*(0.5*tmp_153*tmp_192 - tmp_153*tmp_197 - tmp_158*tmp_81) + tmp_181*(0.5*tmp_172*tmp_192 - tmp_172*tmp_198 - tmp_177*tmp_81) + tmp_86*(0.5*tmp_192*tmp_70 - tmp_193*tmp_70 - tmp_78*tmp_81);
      real_t a_2_1 = tmp_105*(-tmp_101*tmp_182 + 0.5*tmp_192*tmp_93 - tmp_194*tmp_93) + tmp_124*(0.5*tmp_112*tmp_192 - tmp_112*tmp_195 - tmp_120*tmp_182) + tmp_143*(0.5*tmp_131*tmp_192 - tmp_131*tmp_196 - tmp_139*tmp_182) + tmp_162*(0.5*tmp_150*tmp_192 - tmp_150*tmp_197 - tmp_158*tmp_182) + tmp_181*(0.5*tmp_169*tmp_192 - tmp_169*tmp_198 - tmp_177*tmp_182) + tmp_86*(-tmp_182*tmp_78 + 0.5*tmp_192*tmp_61 - tmp_193*tmp_61);
      real_t a_2_2 = tmp_105*(-tmp_101*tmp_183 + 0.5*tmp_192*tmp_94 - tmp_194*tmp_94) + tmp_124*(0.5*tmp_113*tmp_192 - tmp_113*tmp_195 - tmp_120*tmp_183) + tmp_143*(0.5*tmp_132*tmp_192 - tmp_132*tmp_196 - tmp_139*tmp_183) + tmp_162*(0.5*tmp_151*tmp_192 - tmp_151*tmp_197 - tmp_158*tmp_183) + tmp_181*(0.5*tmp_170*tmp_192 - tmp_170*tmp_198 - tmp_177*tmp_183) + tmp_86*(-tmp_183*tmp_78 + 0.5*tmp_192*tmp_65 - tmp_193*tmp_65);
      real_t a_2_3 = tmp_105*(-tmp_101*tmp_184 + 0.5*tmp_192*tmp_95 - tmp_194*tmp_95) + tmp_124*(0.5*tmp_114*tmp_192 - tmp_114*tmp_195 - tmp_120*tmp_184) + tmp_143*(0.5*tmp_133*tmp_192 - tmp_133*tmp_196 - tmp_139*tmp_184) + tmp_162*(0.5*tmp_152*tmp_192 - tmp_152*tmp_197 - tmp_158*tmp_184) + tmp_181*(0.5*tmp_171*tmp_192 - tmp_171*tmp_198 - tmp_177*tmp_184) + tmp_86*(-tmp_184*tmp_78 + 0.5*tmp_192*tmp_69 - tmp_193*tmp_69);
      real_t a_3_0 = tmp_105*(-tmp_102*tmp_81 + 0.5*tmp_199*tmp_96 - tmp_201*tmp_96) + tmp_124*(0.5*tmp_115*tmp_199 - tmp_115*tmp_202 - tmp_121*tmp_81) + tmp_143*(0.5*tmp_134*tmp_199 - tmp_134*tmp_203 - tmp_140*tmp_81) + tmp_162*(0.5*tmp_153*tmp_199 - tmp_153*tmp_204 - tmp_159*tmp_81) + tmp_181*(0.5*tmp_172*tmp_199 - tmp_172*tmp_205 - tmp_178*tmp_81) + tmp_86*(0.5*tmp_199*tmp_70 - tmp_200*tmp_70 - tmp_79*tmp_81);
      real_t a_3_1 = tmp_105*(-tmp_102*tmp_182 + 0.5*tmp_199*tmp_93 - tmp_201*tmp_93) + tmp_124*(0.5*tmp_112*tmp_199 - tmp_112*tmp_202 - tmp_121*tmp_182) + tmp_143*(0.5*tmp_131*tmp_199 - tmp_131*tmp_203 - tmp_140*tmp_182) + tmp_162*(0.5*tmp_150*tmp_199 - tmp_150*tmp_204 - tmp_159*tmp_182) + tmp_181*(0.5*tmp_169*tmp_199 - tmp_169*tmp_205 - tmp_178*tmp_182) + tmp_86*(-tmp_182*tmp_79 + 0.5*tmp_199*tmp_61 - tmp_200*tmp_61);
      real_t a_3_2 = tmp_105*(-tmp_102*tmp_183 + 0.5*tmp_199*tmp_94 - tmp_201*tmp_94) + tmp_124*(0.5*tmp_113*tmp_199 - tmp_113*tmp_202 - tmp_121*tmp_183) + tmp_143*(0.5*tmp_132*tmp_199 - tmp_132*tmp_203 - tmp_140*tmp_183) + tmp_162*(0.5*tmp_151*tmp_199 - tmp_151*tmp_204 - tmp_159*tmp_183) + tmp_181*(0.5*tmp_170*tmp_199 - tmp_170*tmp_205 - tmp_178*tmp_183) + tmp_86*(-tmp_183*tmp_79 + 0.5*tmp_199*tmp_65 - tmp_200*tmp_65);
      real_t a_3_3 = tmp_105*(-tmp_102*tmp_184 + 0.5*tmp_199*tmp_95 - tmp_201*tmp_95) + tmp_124*(0.5*tmp_114*tmp_199 - tmp_114*tmp_202 - tmp_121*tmp_184) + tmp_143*(0.5*tmp_133*tmp_199 - tmp_133*tmp_203 - tmp_140*tmp_184) + tmp_162*(0.5*tmp_152*tmp_199 - tmp_152*tmp_204 - tmp_159*tmp_184) + tmp_181*(0.5*tmp_171*tmp_199 - tmp_171*tmp_205 - tmp_178*tmp_184) + tmp_86*(-tmp_184*tmp_79 + 0.5*tmp_199*tmp_69 - tmp_200*tmp_69);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
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
                                        MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -tmp_4;
      real_t tmp_9 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_10 = 0.81684757298045851*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_11 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_13 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_14 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_15 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_16 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_19 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_20 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_12*tmp_20;
      real_t tmp_23 = tmp_16*tmp_19;
      real_t tmp_24 = tmp_14*tmp_18;
      real_t tmp_25 = 1.0 / (tmp_11*tmp_12*tmp_18 - tmp_11*tmp_23 + tmp_13*tmp_21 - tmp_13*tmp_24 + tmp_15*tmp_17 - tmp_15*tmp_22);
      real_t tmp_26 = tmp_25*(tmp_11*tmp_12 - tmp_13*tmp_14);
      real_t tmp_27 = -tmp_1;
      real_t tmp_28 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_29 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_30 = tmp_25*(-tmp_11*tmp_16 + tmp_13*tmp_20);
      real_t tmp_31 = -tmp_3;
      real_t tmp_32 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_33 = 0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_34 = tmp_25*(tmp_17 - tmp_22);
      real_t tmp_35 = tmp_10*tmp_26 + tmp_29*tmp_30 + tmp_33*tmp_34;
      real_t tmp_36 = tmp_25*(-tmp_12*tmp_15 + tmp_13*tmp_19);
      real_t tmp_37 = tmp_25*(-tmp_13*tmp_18 + tmp_15*tmp_16);
      real_t tmp_38 = tmp_25*(tmp_12*tmp_18 - tmp_23);
      real_t tmp_39 = tmp_10*tmp_36 + tmp_29*tmp_37 + tmp_33*tmp_38;
      real_t tmp_40 = tmp_25*(-tmp_11*tmp_19 + tmp_14*tmp_15);
      real_t tmp_41 = tmp_25*(tmp_11*tmp_18 - tmp_15*tmp_20);
      real_t tmp_42 = tmp_25*(tmp_21 - tmp_24);
      real_t tmp_43 = tmp_10*tmp_40 + tmp_29*tmp_41 + tmp_33*tmp_42;
      real_t tmp_44 = -tmp_35 - tmp_39 - tmp_43 + 1;
      real_t tmp_45 = p_affine_13_0*(-tmp_34 - tmp_38 - tmp_42) + p_affine_13_1*(-tmp_30 - tmp_37 - tmp_41) + p_affine_13_2*(-tmp_26 - tmp_36 - tmp_40);
      real_t tmp_46 = 2*tmp_45;
      real_t tmp_47 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_48 = 0.054975871827660928*tmp_47;
      real_t tmp_49 = 0.10810301816807022*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_50 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_51 = 0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_52 = tmp_26*tmp_49 + tmp_30*tmp_50 + tmp_34*tmp_51;
      real_t tmp_53 = tmp_36*tmp_49 + tmp_37*tmp_50 + tmp_38*tmp_51;
      real_t tmp_54 = tmp_40*tmp_49 + tmp_41*tmp_50 + tmp_42*tmp_51;
      real_t tmp_55 = -tmp_52 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11169079483900572*tmp_47;
      real_t tmp_57 = 0.091576213509770743*tmp_5 + 0.81684757298045851*tmp_8 + tmp_9;
      real_t tmp_58 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_59 = 0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_31 + tmp_32;
      real_t tmp_60 = tmp_26*tmp_57 + tmp_30*tmp_58 + tmp_34*tmp_59;
      real_t tmp_61 = tmp_36*tmp_57 + tmp_37*tmp_58 + tmp_38*tmp_59;
      real_t tmp_62 = tmp_40*tmp_57 + tmp_41*tmp_58 + tmp_42*tmp_59;
      real_t tmp_63 = -tmp_60 - tmp_61 - tmp_62 + 1;
      real_t tmp_64 = 0.054975871827660928*tmp_47;
      real_t tmp_65 = 0.44594849091596489*tmp_5 + 0.10810301816807022*tmp_8 + tmp_9;
      real_t tmp_66 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_67 = 0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_31 + tmp_32;
      real_t tmp_68 = tmp_26*tmp_65 + tmp_30*tmp_66 + tmp_34*tmp_67;
      real_t tmp_69 = tmp_36*tmp_65 + tmp_37*tmp_66 + tmp_38*tmp_67;
      real_t tmp_70 = tmp_40*tmp_65 + tmp_41*tmp_66 + tmp_42*tmp_67;
      real_t tmp_71 = -tmp_68 - tmp_69 - tmp_70 + 1;
      real_t tmp_72 = 0.11169079483900572*tmp_47;
      real_t tmp_73 = 0.091576213509770743*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_74 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_75 = 0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_76 = tmp_26*tmp_73 + tmp_30*tmp_74 + tmp_34*tmp_75;
      real_t tmp_77 = tmp_36*tmp_73 + tmp_37*tmp_74 + tmp_38*tmp_75;
      real_t tmp_78 = tmp_40*tmp_73 + tmp_41*tmp_74 + tmp_42*tmp_75;
      real_t tmp_79 = -tmp_76 - tmp_77 - tmp_78 + 1;
      real_t tmp_80 = 0.054975871827660928*tmp_47;
      real_t tmp_81 = 0.44594849091596489*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_82 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_83 = 0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_84 = tmp_26*tmp_81 + tmp_30*tmp_82 + tmp_34*tmp_83;
      real_t tmp_85 = tmp_36*tmp_81 + tmp_37*tmp_82 + tmp_38*tmp_83;
      real_t tmp_86 = tmp_40*tmp_81 + tmp_41*tmp_82 + tmp_42*tmp_83;
      real_t tmp_87 = -tmp_84 - tmp_85 - tmp_86 + 1;
      real_t tmp_88 = 0.11169079483900572*tmp_47;
      real_t tmp_89 = p_affine_13_0*tmp_34 + p_affine_13_1*tmp_30 + p_affine_13_2*tmp_26;
      real_t tmp_90 = tmp_48*(3.0*tmp_35*tmp_44*tmp_7 - tmp_35*tmp_45 - tmp_44*tmp_89) + tmp_56*(-tmp_45*tmp_52 + 3.0*tmp_52*tmp_55*tmp_7 - tmp_55*tmp_89) + tmp_64*(-tmp_45*tmp_60 + 3.0*tmp_60*tmp_63*tmp_7 - tmp_63*tmp_89) + tmp_72*(-tmp_45*tmp_68 + 3.0*tmp_68*tmp_7*tmp_71 - tmp_71*tmp_89) + tmp_80*(-tmp_45*tmp_76 + 3.0*tmp_7*tmp_76*tmp_79 - tmp_79*tmp_89) + tmp_88*(-tmp_45*tmp_84 + 3.0*tmp_7*tmp_84*tmp_87 - tmp_87*tmp_89);
      real_t tmp_91 = p_affine_13_0*tmp_38 + p_affine_13_1*tmp_37 + p_affine_13_2*tmp_36;
      real_t tmp_92 = tmp_48*(3.0*tmp_39*tmp_44*tmp_7 - tmp_39*tmp_45 - tmp_44*tmp_91) + tmp_56*(-tmp_45*tmp_53 + 3.0*tmp_53*tmp_55*tmp_7 - tmp_55*tmp_91) + tmp_64*(-tmp_45*tmp_61 + 3.0*tmp_61*tmp_63*tmp_7 - tmp_63*tmp_91) + tmp_72*(-tmp_45*tmp_69 + 3.0*tmp_69*tmp_7*tmp_71 - tmp_71*tmp_91) + tmp_80*(-tmp_45*tmp_77 + 3.0*tmp_7*tmp_77*tmp_79 - tmp_79*tmp_91) + tmp_88*(-tmp_45*tmp_85 + 3.0*tmp_7*tmp_85*tmp_87 - tmp_87*tmp_91);
      real_t tmp_93 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_40;
      real_t tmp_94 = tmp_48*(3.0*tmp_43*tmp_44*tmp_7 - tmp_43*tmp_45 - tmp_44*tmp_93) + tmp_56*(-tmp_45*tmp_54 + 3.0*tmp_54*tmp_55*tmp_7 - tmp_55*tmp_93) + tmp_64*(-tmp_45*tmp_62 + 3.0*tmp_62*tmp_63*tmp_7 - tmp_63*tmp_93) + tmp_72*(-tmp_45*tmp_70 + 3.0*tmp_7*tmp_70*tmp_71 - tmp_71*tmp_93) + tmp_80*(-tmp_45*tmp_78 + 3.0*tmp_7*tmp_78*tmp_79 - tmp_79*tmp_93) + tmp_88*(-tmp_45*tmp_86 + 3.0*tmp_7*tmp_86*tmp_87 - tmp_87*tmp_93);
      real_t tmp_95 = 2*tmp_89;
      real_t tmp_96 = tmp_48*(3.0*tmp_35*tmp_39*tmp_7 - tmp_35*tmp_91 - tmp_39*tmp_89) + tmp_56*(3.0*tmp_52*tmp_53*tmp_7 - tmp_52*tmp_91 - tmp_53*tmp_89) + tmp_64*(3.0*tmp_60*tmp_61*tmp_7 - tmp_60*tmp_91 - tmp_61*tmp_89) + tmp_72*(3.0*tmp_68*tmp_69*tmp_7 - tmp_68*tmp_91 - tmp_69*tmp_89) + tmp_80*(3.0*tmp_7*tmp_76*tmp_77 - tmp_76*tmp_91 - tmp_77*tmp_89) + tmp_88*(3.0*tmp_7*tmp_84*tmp_85 - tmp_84*tmp_91 - tmp_85*tmp_89);
      real_t tmp_97 = tmp_48*(3.0*tmp_35*tmp_43*tmp_7 - tmp_35*tmp_93 - tmp_43*tmp_89) + tmp_56*(3.0*tmp_52*tmp_54*tmp_7 - tmp_52*tmp_93 - tmp_54*tmp_89) + tmp_64*(3.0*tmp_60*tmp_62*tmp_7 - tmp_60*tmp_93 - tmp_62*tmp_89) + tmp_72*(3.0*tmp_68*tmp_7*tmp_70 - tmp_68*tmp_93 - tmp_70*tmp_89) + tmp_80*(3.0*tmp_7*tmp_76*tmp_78 - tmp_76*tmp_93 - tmp_78*tmp_89) + tmp_88*(3.0*tmp_7*tmp_84*tmp_86 - tmp_84*tmp_93 - tmp_86*tmp_89);
      real_t tmp_98 = 2*tmp_91;
      real_t tmp_99 = tmp_48*(3.0*tmp_39*tmp_43*tmp_7 - tmp_39*tmp_93 - tmp_43*tmp_91) + tmp_56*(3.0*tmp_53*tmp_54*tmp_7 - tmp_53*tmp_93 - tmp_54*tmp_91) + tmp_64*(3.0*tmp_61*tmp_62*tmp_7 - tmp_61*tmp_93 - tmp_62*tmp_91) + tmp_72*(3.0*tmp_69*tmp_7*tmp_70 - tmp_69*tmp_93 - tmp_70*tmp_91) + tmp_80*(3.0*tmp_7*tmp_77*tmp_78 - tmp_77*tmp_93 - tmp_78*tmp_91) + tmp_88*(3.0*tmp_7*tmp_85*tmp_86 - tmp_85*tmp_93 - tmp_86*tmp_91);
      real_t tmp_100 = 2*tmp_93;
      real_t a_0_0 = tmp_48*(3.0*(tmp_44*tmp_44)*tmp_7 - tmp_44*tmp_46) + tmp_56*(-tmp_46*tmp_55 + 3.0*(tmp_55*tmp_55)*tmp_7) + tmp_64*(-tmp_46*tmp_63 + 3.0*(tmp_63*tmp_63)*tmp_7) + tmp_72*(-tmp_46*tmp_71 + 3.0*tmp_7*(tmp_71*tmp_71)) + tmp_80*(-tmp_46*tmp_79 + 3.0*tmp_7*(tmp_79*tmp_79)) + tmp_88*(-tmp_46*tmp_87 + 3.0*tmp_7*(tmp_87*tmp_87));
      real_t a_0_1 = tmp_90;
      real_t a_0_2 = tmp_92;
      real_t a_0_3 = tmp_94;
      real_t a_1_0 = tmp_90;
      real_t a_1_1 = tmp_48*(3.0*(tmp_35*tmp_35)*tmp_7 - tmp_35*tmp_95) + tmp_56*(3.0*(tmp_52*tmp_52)*tmp_7 - tmp_52*tmp_95) + tmp_64*(3.0*(tmp_60*tmp_60)*tmp_7 - tmp_60*tmp_95) + tmp_72*(3.0*(tmp_68*tmp_68)*tmp_7 - tmp_68*tmp_95) + tmp_80*(3.0*tmp_7*(tmp_76*tmp_76) - tmp_76*tmp_95) + tmp_88*(3.0*tmp_7*(tmp_84*tmp_84) - tmp_84*tmp_95);
      real_t a_1_2 = tmp_96;
      real_t a_1_3 = tmp_97;
      real_t a_2_0 = tmp_92;
      real_t a_2_1 = tmp_96;
      real_t a_2_2 = tmp_48*(3.0*(tmp_39*tmp_39)*tmp_7 - tmp_39*tmp_98) + tmp_56*(3.0*(tmp_53*tmp_53)*tmp_7 - tmp_53*tmp_98) + tmp_64*(3.0*(tmp_61*tmp_61)*tmp_7 - tmp_61*tmp_98) + tmp_72*(3.0*(tmp_69*tmp_69)*tmp_7 - tmp_69*tmp_98) + tmp_80*(3.0*tmp_7*(tmp_77*tmp_77) - tmp_77*tmp_98) + tmp_88*(3.0*tmp_7*(tmp_85*tmp_85) - tmp_85*tmp_98);
      real_t a_2_3 = tmp_99;
      real_t a_3_0 = tmp_94;
      real_t a_3_1 = tmp_97;
      real_t a_3_2 = tmp_99;
      real_t a_3_3 = tmp_48*(-tmp_100*tmp_43 + 3.0*(tmp_43*tmp_43)*tmp_7) + tmp_56*(-tmp_100*tmp_54 + 3.0*(tmp_54*tmp_54)*tmp_7) + tmp_64*(-tmp_100*tmp_62 + 3.0*(tmp_62*tmp_62)*tmp_7) + tmp_72*(-tmp_100*tmp_70 + 3.0*tmp_7*(tmp_70*tmp_70)) + tmp_80*(-tmp_100*tmp_78 + 3.0*tmp_7*(tmp_78*tmp_78)) + tmp_88*(-tmp_100*tmp_86 + 3.0*tmp_7*(tmp_86*tmp_86));
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g0;

};




class EGVectorLaplaceFormNitscheBC_P1E_0 : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_P1E_0()

    {}





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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = -tmp_3;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = 1.0 / (tmp_2 + tmp_4*tmp_5);
      real_t tmp_7 = tmp_0*tmp_6;
      real_t tmp_8 = tmp_3*tmp_6;
      real_t tmp_9 = tmp_0*tmp_8 + tmp_4*tmp_7;
      real_t tmp_10 = tmp_1*tmp_6;
      real_t tmp_11 = tmp_5*tmp_6;
      real_t tmp_12 = tmp_11*tmp_4 + tmp_2*tmp_6;
      real_t tmp_13 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = tmp_13*(tmp_12*(-tmp_10 - tmp_11) + tmp_9*(-tmp_7 - tmp_8));
      real_t tmp_15 = tmp_13*(tmp_10*tmp_12 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_13*(tmp_11*tmp_12 + tmp_7*tmp_9);
      real_t a_0_0 = 0.5*tmp_14;
      real_t a_1_0 = 0.5*tmp_15;
      real_t a_2_0 = 0.5*tmp_16;
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
                                       MatrixXr&                                           elMat ) const override
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
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_3 = 0.21132486540518713*tmp_1 + tmp_2;
      real_t tmp_4 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_0*tmp_5;
      real_t tmp_7 = -tmp_4;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = 1.0 / (tmp_6 + tmp_7*tmp_8);
      real_t tmp_10 = tmp_4*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_11 + tmp_12;
      real_t tmp_14 = tmp_5*tmp_9;
      real_t tmp_15 = tmp_10*tmp_3 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_0*tmp_9;
      real_t tmp_17 = tmp_8*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_3;
      real_t tmp_19 = tmp_0*(tmp_15 - 1.0/3.0) + tmp_7*(tmp_18 - 1.0/3.0);
      real_t tmp_20 = 0.5*p_affine_10_0*(-tmp_14 - tmp_17) + 0.5*p_affine_10_1*(-tmp_10 - tmp_16);
      real_t tmp_21 = -tmp_15 - tmp_18 + 1;
      real_t tmp_22 = 0.5*p_affine_10_0*(tmp_17*tmp_7 + tmp_6*tmp_9) + 0.5*p_affine_10_1*(tmp_0*tmp_10 + tmp_16*tmp_7);
      real_t tmp_23 = std::abs(std::pow((tmp_1*tmp_1) + (tmp_11*tmp_11), 1.0/2.0));
      real_t tmp_24 = 1.0 / (tmp_23);
      real_t tmp_25 = 0.5*tmp_23;
      real_t tmp_26 = 0.78867513459481287*tmp_1 + tmp_2;
      real_t tmp_27 = 0.78867513459481287*tmp_11 + tmp_12;
      real_t tmp_28 = tmp_10*tmp_26 + tmp_14*tmp_27;
      real_t tmp_29 = tmp_16*tmp_26 + tmp_17*tmp_27;
      real_t tmp_30 = tmp_0*(tmp_28 - 1.0/3.0) + tmp_7*(tmp_29 - 1.0/3.0);
      real_t tmp_31 = -tmp_28 - tmp_29 + 1;
      real_t tmp_32 = 0.5*tmp_23;
      real_t tmp_33 = 0.5*p_affine_10_0*tmp_14 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_34 = 0.5*p_affine_10_0*tmp_17 + 0.5*p_affine_10_1*tmp_16;
      real_t a_0_0 = tmp_25*(-tmp_19*tmp_20 + 3*tmp_19*tmp_21*tmp_24 - tmp_21*tmp_22) + tmp_32*(-tmp_20*tmp_30 - tmp_22*tmp_31 + 3*tmp_24*tmp_30*tmp_31);
      real_t a_1_0 = tmp_25*(3*tmp_15*tmp_19*tmp_24 - tmp_15*tmp_22 - tmp_19*tmp_33) + tmp_32*(-tmp_22*tmp_28 + 3*tmp_24*tmp_28*tmp_30 - tmp_30*tmp_33);
      real_t a_2_0 = tmp_25*(3*tmp_18*tmp_19*tmp_24 - tmp_18*tmp_22 - tmp_19*tmp_34) + tmp_32*(-tmp_22*tmp_29 + 3*tmp_24*tmp_29*tmp_30 - tmp_30*tmp_34);
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
                                          MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = 1.0 / (tmp_0*tmp_1 - tmp_2*tmp_3);
      real_t tmp_5 = tmp_0*tmp_4;
      real_t tmp_6 = tmp_3*tmp_4;
      real_t tmp_7 = tmp_1*tmp_4;
      real_t tmp_8 = tmp_2*tmp_4;
      real_t tmp_9 = p_affine_10_0*(-tmp_5 - tmp_6) + p_affine_10_1*(-tmp_7 - tmp_8);
      real_t tmp_10 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_11 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_12 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = -tmp_11;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = 1.0 / (tmp_13 + tmp_14*tmp_15);
      real_t tmp_17 = -p_affine_3_1;
      real_t tmp_18 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_19 = p_affine_6_1 + 0.21132486540518713*tmp_18;
      real_t tmp_20 = tmp_16*(tmp_17 + tmp_19);
      real_t tmp_21 = -p_affine_3_0;
      real_t tmp_22 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_23 = p_affine_6_0 + 0.21132486540518713*tmp_22;
      real_t tmp_24 = tmp_16*(tmp_21 + tmp_23);
      real_t tmp_25 = tmp_10*(tmp_11*tmp_20 + tmp_12*tmp_24 - 1.0/3.0) + tmp_14*(tmp_10*tmp_20 + tmp_15*tmp_24 - 1.0/3.0);
      real_t tmp_26 = -p_affine_0_1;
      real_t tmp_27 = tmp_19 + tmp_26;
      real_t tmp_28 = -p_affine_0_0;
      real_t tmp_29 = tmp_23 + tmp_28;
      real_t tmp_30 = tmp_27*tmp_8 + tmp_29*tmp_5;
      real_t tmp_31 = tmp_27*tmp_7 + tmp_29*tmp_6;
      real_t tmp_32 = -tmp_30 - tmp_31 + 1;
      real_t tmp_33 = tmp_14*tmp_16;
      real_t tmp_34 = tmp_11*tmp_16;
      real_t tmp_35 = 0.5*p_affine_10_0*(tmp_13*tmp_16 + tmp_15*tmp_33) + 0.5*p_affine_10_1*(tmp_10*tmp_33 + tmp_10*tmp_34);
      real_t tmp_36 = std::abs(std::pow((tmp_18*tmp_18) + (tmp_22*tmp_22), 1.0/2.0));
      real_t tmp_37 = 3/tmp_36;
      real_t tmp_38 = tmp_25*tmp_37;
      real_t tmp_39 = 0.5*tmp_36;
      real_t tmp_40 = p_affine_6_1 + 0.78867513459481287*tmp_18;
      real_t tmp_41 = tmp_17 + tmp_40;
      real_t tmp_42 = p_affine_6_0 + 0.78867513459481287*tmp_22;
      real_t tmp_43 = tmp_16*(tmp_21 + tmp_42);
      real_t tmp_44 = tmp_10*(tmp_12*tmp_43 + tmp_34*tmp_41 - 1.0/3.0) + tmp_14*(tmp_10*tmp_16*tmp_41 + tmp_15*tmp_43 - 1.0/3.0);
      real_t tmp_45 = tmp_26 + tmp_40;
      real_t tmp_46 = tmp_28 + tmp_42;
      real_t tmp_47 = tmp_45*tmp_8 + tmp_46*tmp_5;
      real_t tmp_48 = tmp_45*tmp_7 + tmp_46*tmp_6;
      real_t tmp_49 = -tmp_47 - tmp_48 + 1;
      real_t tmp_50 = tmp_37*tmp_44;
      real_t tmp_51 = 0.5*tmp_36;
      real_t tmp_52 = p_affine_10_0*tmp_5 + p_affine_10_1*tmp_8;
      real_t tmp_53 = p_affine_10_0*tmp_6 + p_affine_10_1*tmp_7;
      real_t a_0_0 = tmp_39*(0.5*tmp_25*tmp_9 - tmp_32*tmp_35 - tmp_32*tmp_38) + tmp_51*(-tmp_35*tmp_49 + 0.5*tmp_44*tmp_9 - tmp_49*tmp_50);
      real_t a_1_0 = tmp_39*(0.5*tmp_25*tmp_52 - tmp_30*tmp_35 - tmp_30*tmp_38) + tmp_51*(-tmp_35*tmp_47 + 0.5*tmp_44*tmp_52 - tmp_47*tmp_50);
      real_t a_2_0 = tmp_39*(0.5*tmp_25*tmp_53 - tmp_31*tmp_35 - tmp_31*tmp_38) + tmp_51*(-tmp_35*tmp_48 + 0.5*tmp_44*tmp_53 - tmp_48*tmp_50);
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
                                                   MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = -tmp_3;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = 1.0 / (tmp_2 + tmp_4*tmp_5);
      real_t tmp_7 = tmp_0*tmp_6;
      real_t tmp_8 = tmp_5*tmp_6;
      real_t tmp_9 = tmp_1*tmp_6;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = p_affine_10_0*(-tmp_7 - tmp_8) + p_affine_10_1*(-tmp_10 - tmp_9);
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_14 = 0.21132486540518713*tmp_12 + tmp_13;
      real_t tmp_15 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_16 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_17 = 0.21132486540518713*tmp_15 + tmp_16;
      real_t tmp_18 = tmp_10*tmp_14 + tmp_17*tmp_7;
      real_t tmp_19 = tmp_14*tmp_9 + tmp_17*tmp_8;
      real_t tmp_20 = tmp_1*(tmp_18 - 1.0/3.0) + tmp_4*(tmp_19 - 1.0/3.0);
      real_t tmp_21 = p_affine_10_0*(tmp_2*tmp_6 + tmp_4*tmp_8) + p_affine_10_1*(tmp_1*tmp_10 + tmp_4*tmp_9);
      real_t tmp_22 = -tmp_18 - tmp_19 + 1;
      real_t tmp_23 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_15*tmp_15), 1.0/2.0));
      real_t tmp_24 = 1.0 / (tmp_23);
      real_t tmp_25 = 0.5*tmp_23;
      real_t tmp_26 = 0.78867513459481287*tmp_12 + tmp_13;
      real_t tmp_27 = 0.78867513459481287*tmp_15 + tmp_16;
      real_t tmp_28 = tmp_10*tmp_26 + tmp_27*tmp_7;
      real_t tmp_29 = tmp_26*tmp_9 + tmp_27*tmp_8;
      real_t tmp_30 = tmp_1*(tmp_28 - 1.0/3.0) + tmp_4*(tmp_29 - 1.0/3.0);
      real_t tmp_31 = -tmp_28 - tmp_29 + 1;
      real_t tmp_32 = 0.5*tmp_23;
      real_t tmp_33 = p_affine_10_0*tmp_7 + p_affine_10_1*tmp_10;
      real_t tmp_34 = p_affine_10_0*tmp_8 + p_affine_10_1*tmp_9;
      real_t a_0_0 = tmp_25*(-tmp_11*tmp_20 + 3*tmp_20*tmp_22*tmp_24 - tmp_21*tmp_22) + tmp_32*(-tmp_11*tmp_30 - tmp_21*tmp_31 + 3*tmp_24*tmp_30*tmp_31);
      real_t a_1_0 = tmp_25*(3*tmp_18*tmp_20*tmp_24 - tmp_18*tmp_21 - tmp_20*tmp_33) + tmp_32*(-tmp_21*tmp_28 + 3*tmp_24*tmp_28*tmp_30 - tmp_30*tmp_33);
      real_t a_2_0 = tmp_25*(3*tmp_19*tmp_20*tmp_24 - tmp_19*tmp_21 - tmp_20*tmp_34) + tmp_32*(-tmp_21*tmp_29 + 3*tmp_24*tmp_29*tmp_30 - tmp_30*tmp_34);
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
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_8 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_0*tmp_8;
      real_t tmp_13 = tmp_1*tmp_11;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_11*tmp_4 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_2*tmp_6 - tmp_5*tmp_6 + tmp_7*tmp_9);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_11*tmp_4 - tmp_12);
      real_t tmp_17 = tmp_14*(-tmp_13 + tmp_9);
      real_t tmp_18 = tmp_0*tmp_17 + tmp_11*tmp_15 + tmp_16*tmp_3;
      real_t tmp_19 = tmp_14*(-tmp_0*tmp_10 + tmp_3*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_0*tmp_6 - tmp_11*tmp_7);
      real_t tmp_21 = tmp_14*(tmp_10*tmp_11 - tmp_3*tmp_6);
      real_t tmp_22 = tmp_0*tmp_21 + tmp_11*tmp_19 + tmp_20*tmp_3;
      real_t tmp_23 = tmp_14*(-tmp_1*tmp_7 + tmp_10*tmp_4);
      real_t tmp_24 = tmp_14*(-tmp_4*tmp_6 + tmp_7*tmp_8);
      real_t tmp_25 = tmp_14*(tmp_1*tmp_6 - tmp_10*tmp_8);
      real_t tmp_26 = tmp_0*tmp_25 + tmp_11*tmp_23 + tmp_24*tmp_3;
      real_t tmp_27 = p_affine_0_0*p_affine_1_1;
      real_t tmp_28 = p_affine_0_0*p_affine_1_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_2;
      real_t tmp_30 = p_affine_0_1*p_affine_1_0;
      real_t tmp_31 = p_affine_0_1*p_affine_1_2;
      real_t tmp_32 = p_affine_2_2*p_affine_3_0;
      real_t tmp_33 = p_affine_0_2*p_affine_1_0;
      real_t tmp_34 = p_affine_0_2*p_affine_1_1;
      real_t tmp_35 = p_affine_2_0*p_affine_3_1;
      real_t tmp_36 = p_affine_2_2*p_affine_3_1;
      real_t tmp_37 = p_affine_2_0*p_affine_3_2;
      real_t tmp_38 = p_affine_2_1*p_affine_3_0;
      real_t tmp_39 = std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = tmp_39*(tmp_18*(-tmp_15 - tmp_16 - tmp_17) + tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25));
      real_t tmp_41 = tmp_39*(tmp_17*tmp_18 + tmp_21*tmp_22 + tmp_25*tmp_26);
      real_t tmp_42 = tmp_39*(tmp_16*tmp_18 + tmp_20*tmp_22 + tmp_24*tmp_26);
      real_t tmp_43 = tmp_39*(tmp_15*tmp_18 + tmp_19*tmp_22 + tmp_23*tmp_26);
      real_t a_0_0 = 0.16666666666666666*tmp_40;
      real_t a_1_0 = 0.16666666666666666*tmp_41;
      real_t a_2_0 = 0.16666666666666666*tmp_42;
      real_t a_3_0 = 0.16666666666666666*tmp_43;
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
                               MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_2 = -tmp_1;
      real_t tmp_3 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_4 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_5 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_3 + tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_7 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_8 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_9 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_13 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_14 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_14*tmp_7;
      real_t tmp_17 = tmp_10*tmp_13;
      real_t tmp_18 = tmp_12*tmp_9;
      real_t tmp_19 = 1.0 / (tmp_0*tmp_11 - tmp_0*tmp_16 + tmp_12*tmp_6*tmp_7 + tmp_15*tmp_8 - tmp_17*tmp_6 - tmp_18*tmp_8);
      real_t tmp_20 = tmp_19*(tmp_6*tmp_7 - tmp_8*tmp_9);
      real_t tmp_21 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_22 = -tmp_21;
      real_t tmp_23 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_24 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_25 = 0.091576213509770743*tmp_22 + 0.81684757298045851*tmp_23 + tmp_24;
      real_t tmp_26 = tmp_19*(-tmp_10*tmp_6 + tmp_14*tmp_8);
      real_t tmp_27 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_28 = -tmp_27;
      real_t tmp_29 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_30 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_31 = 0.091576213509770743*tmp_28 + 0.81684757298045851*tmp_29 + tmp_30;
      real_t tmp_32 = tmp_19*(tmp_11 - tmp_16);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_25*tmp_26 + tmp_31*tmp_32;
      real_t tmp_34 = tmp_19*(-tmp_0*tmp_7 + tmp_13*tmp_8);
      real_t tmp_35 = tmp_19*(tmp_0*tmp_10 - tmp_12*tmp_8);
      real_t tmp_36 = tmp_19*(tmp_12*tmp_7 - tmp_17);
      real_t tmp_37 = tmp_25*tmp_35 + tmp_31*tmp_36 + tmp_34*tmp_5;
      real_t tmp_38 = tmp_19*(tmp_0*tmp_9 - tmp_13*tmp_6);
      real_t tmp_39 = tmp_19*(-tmp_0*tmp_14 + tmp_12*tmp_6);
      real_t tmp_40 = tmp_19*(tmp_15 - tmp_18);
      real_t tmp_41 = tmp_25*tmp_39 + tmp_31*tmp_40 + tmp_38*tmp_5;
      real_t tmp_42 = tmp_0*(tmp_33 - 1.0/4.0) + tmp_6*(tmp_37 - 1.0/4.0) + tmp_8*(tmp_41 - 1.0/4.0);
      real_t tmp_43 = 0.5*p_affine_13_0*(-tmp_32 - tmp_36 - tmp_40) + 0.5*p_affine_13_1*(-tmp_26 - tmp_35 - tmp_39) + 0.5*p_affine_13_2*(-tmp_20 - tmp_34 - tmp_38);
      real_t tmp_44 = -tmp_33 - tmp_37 - tmp_41 + 1;
      real_t tmp_45 = 0.5*p_affine_13_0*(tmp_0*tmp_32 + tmp_36*tmp_6 + tmp_40*tmp_8) + 0.5*p_affine_13_1*(tmp_0*tmp_26 + tmp_35*tmp_6 + tmp_39*tmp_8) + 0.5*p_affine_13_2*(tmp_0*tmp_20 + tmp_34*tmp_6 + tmp_38*tmp_8);
      real_t tmp_46 = (std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)*std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)) + (std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)*std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)) + (std::abs(tmp_21*tmp_29 - tmp_23*tmp_27)*std::abs(tmp_21*tmp_29 - tmp_23*tmp_27));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_3 + tmp_4;
      real_t tmp_51 = 0.44594849091596489*tmp_22 + 0.10810301816807022*tmp_23 + tmp_24;
      real_t tmp_52 = 0.44594849091596489*tmp_28 + 0.10810301816807022*tmp_29 + tmp_30;
      real_t tmp_53 = tmp_20*tmp_50 + tmp_26*tmp_51 + tmp_32*tmp_52;
      real_t tmp_54 = tmp_34*tmp_50 + tmp_35*tmp_51 + tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_50 + tmp_39*tmp_51 + tmp_40*tmp_52;
      real_t tmp_56 = tmp_0*(tmp_53 - 1.0/4.0) + tmp_6*(tmp_54 - 1.0/4.0) + tmp_8*(tmp_55 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_60 = 0.81684757298045851*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_61 = 0.81684757298045851*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_62 = tmp_20*tmp_59 + tmp_26*tmp_60 + tmp_32*tmp_61;
      real_t tmp_63 = tmp_34*tmp_59 + tmp_35*tmp_60 + tmp_36*tmp_61;
      real_t tmp_64 = tmp_38*tmp_59 + tmp_39*tmp_60 + tmp_40*tmp_61;
      real_t tmp_65 = tmp_0*(tmp_62 - 1.0/4.0) + tmp_6*(tmp_63 - 1.0/4.0) + tmp_8*(tmp_64 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_69 = 0.10810301816807022*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_70 = 0.10810301816807022*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_71 = tmp_20*tmp_68 + tmp_26*tmp_69 + tmp_32*tmp_70;
      real_t tmp_72 = tmp_34*tmp_68 + tmp_35*tmp_69 + tmp_36*tmp_70;
      real_t tmp_73 = tmp_38*tmp_68 + tmp_39*tmp_69 + tmp_40*tmp_70;
      real_t tmp_74 = tmp_0*(tmp_71 - 1.0/4.0) + tmp_6*(tmp_72 - 1.0/4.0) + tmp_8*(tmp_73 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_78 = 0.091576213509770743*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_79 = 0.091576213509770743*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_80 = tmp_20*tmp_77 + tmp_26*tmp_78 + tmp_32*tmp_79;
      real_t tmp_81 = tmp_34*tmp_77 + tmp_35*tmp_78 + tmp_36*tmp_79;
      real_t tmp_82 = tmp_38*tmp_77 + tmp_39*tmp_78 + tmp_40*tmp_79;
      real_t tmp_83 = tmp_0*(tmp_80 - 1.0/4.0) + tmp_6*(tmp_81 - 1.0/4.0) + tmp_8*(tmp_82 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_87 = 0.44594849091596489*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_88 = 0.44594849091596489*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_89 = tmp_20*tmp_86 + tmp_26*tmp_87 + tmp_32*tmp_88;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = tmp_0*(tmp_89 - 1.0/4.0) + tmp_6*(tmp_90 - 1.0/4.0) + tmp_8*(tmp_91 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = 0.5*p_affine_13_0*tmp_32 + 0.5*p_affine_13_1*tmp_26 + 0.5*p_affine_13_2*tmp_20;
      real_t tmp_96 = 0.5*p_affine_13_0*tmp_36 + 0.5*p_affine_13_1*tmp_35 + 0.5*p_affine_13_2*tmp_34;
      real_t tmp_97 = 0.5*p_affine_13_0*tmp_40 + 0.5*p_affine_13_1*tmp_39 + 0.5*p_affine_13_2*tmp_38;
      real_t a_0_0 = tmp_49*(-tmp_42*tmp_43 + 3.0*tmp_42*tmp_44*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_43*tmp_56 - tmp_45*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_43*tmp_65 - tmp_45*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_43*tmp_74 - tmp_45*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_43*tmp_83 - tmp_45*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_43*tmp_92 - tmp_45*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_1_0 = tmp_49*(3.0*tmp_33*tmp_42*tmp_47 - tmp_33*tmp_45 - tmp_42*tmp_95) + tmp_58*(-tmp_45*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_45*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_45*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_45*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_45*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_2_0 = tmp_49*(3.0*tmp_37*tmp_42*tmp_47 - tmp_37*tmp_45 - tmp_42*tmp_96) + tmp_58*(-tmp_45*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_45*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_45*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_45*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_45*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_3_0 = tmp_49*(3.0*tmp_41*tmp_42*tmp_47 - tmp_41*tmp_45 - tmp_42*tmp_97) + tmp_58*(-tmp_45*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_45*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_45*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_45*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_45*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
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
                                  MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_26 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_27 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_28 = tmp_26*tmp_27;
      real_t tmp_29 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_30 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_31 = tmp_29*tmp_30;
      real_t tmp_32 = tmp_28 - tmp_31;
      real_t tmp_33 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_34 = tmp_30*tmp_33;
      real_t tmp_35 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_36 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_37 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_38 = tmp_27*tmp_37;
      real_t tmp_39 = tmp_26*tmp_33;
      real_t tmp_40 = 1.0 / (tmp_25*tmp_34 - tmp_25*tmp_38 + tmp_28*tmp_35 + tmp_29*tmp_36*tmp_37 - tmp_31*tmp_35 - tmp_36*tmp_39);
      real_t tmp_41 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_42 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_43 = -tmp_42;
      real_t tmp_44 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_45 = 0.091576213509770743*tmp_43 + 0.81684757298045851*tmp_44;
      real_t tmp_46 = tmp_40*(tmp_41 + tmp_45);
      real_t tmp_47 = tmp_29*tmp_37 - tmp_39;
      real_t tmp_48 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_49 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_50 = -tmp_49;
      real_t tmp_51 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_52 = 0.091576213509770743*tmp_50 + 0.81684757298045851*tmp_51;
      real_t tmp_53 = tmp_40*(tmp_48 + tmp_52);
      real_t tmp_54 = tmp_34 - tmp_38;
      real_t tmp_55 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = -tmp_56;
      real_t tmp_58 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_59 = 0.091576213509770743*tmp_57 + 0.81684757298045851*tmp_58;
      real_t tmp_60 = tmp_40*(tmp_55 + tmp_59);
      real_t tmp_61 = -tmp_25*tmp_27 + tmp_29*tmp_36;
      real_t tmp_62 = tmp_25*tmp_33 - tmp_29*tmp_35;
      real_t tmp_63 = tmp_27*tmp_35 - tmp_33*tmp_36;
      real_t tmp_64 = tmp_25*tmp_30 - tmp_26*tmp_36;
      real_t tmp_65 = -tmp_25*tmp_37 + tmp_26*tmp_35;
      real_t tmp_66 = -tmp_30*tmp_35 + tmp_36*tmp_37;
      real_t tmp_67 = tmp_25*(tmp_32*tmp_46 + tmp_47*tmp_53 + tmp_54*tmp_60 - 1.0/4.0) + tmp_26*(tmp_46*tmp_61 + tmp_53*tmp_62 + tmp_60*tmp_63 - 1.0/4.0) + tmp_29*(tmp_46*tmp_64 + tmp_53*tmp_65 + tmp_60*tmp_66 - 1.0/4.0);
      real_t tmp_68 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_69 = tmp_45 + tmp_68;
      real_t tmp_70 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_71 = tmp_52 + tmp_70;
      real_t tmp_72 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_73 = tmp_59 + tmp_72;
      real_t tmp_74 = tmp_17*tmp_73 + tmp_20*tmp_71 + tmp_23*tmp_69;
      real_t tmp_75 = tmp_16*tmp_73 + tmp_19*tmp_71 + tmp_22*tmp_69;
      real_t tmp_76 = tmp_15*tmp_73 + tmp_18*tmp_71 + tmp_21*tmp_69;
      real_t tmp_77 = -tmp_74 - tmp_75 - tmp_76 + 1;
      real_t tmp_78 = tmp_25*tmp_40;
      real_t tmp_79 = tmp_26*tmp_40;
      real_t tmp_80 = tmp_29*tmp_40;
      real_t tmp_81 = 0.5*p_affine_13_0*(tmp_54*tmp_78 + tmp_63*tmp_79 + tmp_66*tmp_80) + 0.5*p_affine_13_1*(tmp_47*tmp_78 + tmp_62*tmp_79 + tmp_65*tmp_80) + 0.5*p_affine_13_2*(tmp_32*tmp_78 + tmp_61*tmp_79 + tmp_64*tmp_80);
      real_t tmp_82 = (std::abs(tmp_42*tmp_51 - tmp_44*tmp_49)*std::abs(tmp_42*tmp_51 - tmp_44*tmp_49)) + (std::abs(tmp_42*tmp_58 - tmp_44*tmp_56)*std::abs(tmp_42*tmp_58 - tmp_44*tmp_56)) + (std::abs(tmp_49*tmp_58 - tmp_51*tmp_56)*std::abs(tmp_49*tmp_58 - tmp_51*tmp_56));
      real_t tmp_83 = 3.0*std::pow(tmp_82, -0.25);
      real_t tmp_84 = tmp_67*tmp_83;
      real_t tmp_85 = 1.0*std::pow(tmp_82, 1.0/2.0);
      real_t tmp_86 = 0.054975871827660928*tmp_85;
      real_t tmp_87 = 0.44594849091596489*tmp_43 + 0.10810301816807022*tmp_44;
      real_t tmp_88 = tmp_40*(tmp_41 + tmp_87);
      real_t tmp_89 = 0.44594849091596489*tmp_50 + 0.10810301816807022*tmp_51;
      real_t tmp_90 = tmp_40*(tmp_48 + tmp_89);
      real_t tmp_91 = 0.44594849091596489*tmp_57 + 0.10810301816807022*tmp_58;
      real_t tmp_92 = tmp_40*(tmp_55 + tmp_91);
      real_t tmp_93 = tmp_25*(tmp_32*tmp_88 + tmp_47*tmp_90 + tmp_54*tmp_92 - 1.0/4.0) + tmp_26*(tmp_61*tmp_88 + tmp_62*tmp_90 + tmp_63*tmp_92 - 1.0/4.0) + tmp_29*(tmp_64*tmp_88 + tmp_65*tmp_90 + tmp_66*tmp_92 - 1.0/4.0);
      real_t tmp_94 = tmp_68 + tmp_87;
      real_t tmp_95 = tmp_70 + tmp_89;
      real_t tmp_96 = tmp_72 + tmp_91;
      real_t tmp_97 = tmp_17*tmp_96 + tmp_20*tmp_95 + tmp_23*tmp_94;
      real_t tmp_98 = tmp_16*tmp_96 + tmp_19*tmp_95 + tmp_22*tmp_94;
      real_t tmp_99 = tmp_15*tmp_96 + tmp_18*tmp_95 + tmp_21*tmp_94;
      real_t tmp_100 = -tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_101 = tmp_83*tmp_93;
      real_t tmp_102 = 0.11169079483900572*tmp_85;
      real_t tmp_103 = 0.81684757298045851*tmp_43 + 0.091576213509770743*tmp_44;
      real_t tmp_104 = tmp_40*(tmp_103 + tmp_41);
      real_t tmp_105 = 0.81684757298045851*tmp_50 + 0.091576213509770743*tmp_51;
      real_t tmp_106 = tmp_40*(tmp_105 + tmp_48);
      real_t tmp_107 = 0.81684757298045851*tmp_57 + 0.091576213509770743*tmp_58;
      real_t tmp_108 = tmp_40*(tmp_107 + tmp_55);
      real_t tmp_109 = tmp_25*(tmp_104*tmp_32 + tmp_106*tmp_47 + tmp_108*tmp_54 - 1.0/4.0) + tmp_26*(tmp_104*tmp_61 + tmp_106*tmp_62 + tmp_108*tmp_63 - 1.0/4.0) + tmp_29*(tmp_104*tmp_64 + tmp_106*tmp_65 + tmp_108*tmp_66 - 1.0/4.0);
      real_t tmp_110 = tmp_103 + tmp_68;
      real_t tmp_111 = tmp_105 + tmp_70;
      real_t tmp_112 = tmp_107 + tmp_72;
      real_t tmp_113 = tmp_110*tmp_23 + tmp_111*tmp_20 + tmp_112*tmp_17;
      real_t tmp_114 = tmp_110*tmp_22 + tmp_111*tmp_19 + tmp_112*tmp_16;
      real_t tmp_115 = tmp_110*tmp_21 + tmp_111*tmp_18 + tmp_112*tmp_15;
      real_t tmp_116 = -tmp_113 - tmp_114 - tmp_115 + 1;
      real_t tmp_117 = tmp_109*tmp_83;
      real_t tmp_118 = 0.054975871827660928*tmp_85;
      real_t tmp_119 = 0.10810301816807022*tmp_43 + 0.44594849091596489*tmp_44;
      real_t tmp_120 = tmp_40*(tmp_119 + tmp_41);
      real_t tmp_121 = 0.10810301816807022*tmp_50 + 0.44594849091596489*tmp_51;
      real_t tmp_122 = tmp_40*(tmp_121 + tmp_48);
      real_t tmp_123 = 0.10810301816807022*tmp_57 + 0.44594849091596489*tmp_58;
      real_t tmp_124 = tmp_40*(tmp_123 + tmp_55);
      real_t tmp_125 = tmp_25*(tmp_120*tmp_32 + tmp_122*tmp_47 + tmp_124*tmp_54 - 1.0/4.0) + tmp_26*(tmp_120*tmp_61 + tmp_122*tmp_62 + tmp_124*tmp_63 - 1.0/4.0) + tmp_29*(tmp_120*tmp_64 + tmp_122*tmp_65 + tmp_124*tmp_66 - 1.0/4.0);
      real_t tmp_126 = tmp_119 + tmp_68;
      real_t tmp_127 = tmp_121 + tmp_70;
      real_t tmp_128 = tmp_123 + tmp_72;
      real_t tmp_129 = tmp_126*tmp_23 + tmp_127*tmp_20 + tmp_128*tmp_17;
      real_t tmp_130 = tmp_126*tmp_22 + tmp_127*tmp_19 + tmp_128*tmp_16;
      real_t tmp_131 = tmp_126*tmp_21 + tmp_127*tmp_18 + tmp_128*tmp_15;
      real_t tmp_132 = -tmp_129 - tmp_130 - tmp_131 + 1;
      real_t tmp_133 = tmp_125*tmp_83;
      real_t tmp_134 = 0.11169079483900572*tmp_85;
      real_t tmp_135 = 0.091576213509770743*tmp_43 + 0.091576213509770743*tmp_44;
      real_t tmp_136 = tmp_40*(tmp_135 + tmp_41);
      real_t tmp_137 = 0.091576213509770743*tmp_50 + 0.091576213509770743*tmp_51;
      real_t tmp_138 = tmp_40*(tmp_137 + tmp_48);
      real_t tmp_139 = 0.091576213509770743*tmp_57 + 0.091576213509770743*tmp_58;
      real_t tmp_140 = tmp_40*(tmp_139 + tmp_55);
      real_t tmp_141 = tmp_25*(tmp_136*tmp_32 + tmp_138*tmp_47 + tmp_140*tmp_54 - 1.0/4.0) + tmp_26*(tmp_136*tmp_61 + tmp_138*tmp_62 + tmp_140*tmp_63 - 1.0/4.0) + tmp_29*(tmp_136*tmp_64 + tmp_138*tmp_65 + tmp_140*tmp_66 - 1.0/4.0);
      real_t tmp_142 = tmp_135 + tmp_68;
      real_t tmp_143 = tmp_137 + tmp_70;
      real_t tmp_144 = tmp_139 + tmp_72;
      real_t tmp_145 = tmp_142*tmp_23 + tmp_143*tmp_20 + tmp_144*tmp_17;
      real_t tmp_146 = tmp_142*tmp_22 + tmp_143*tmp_19 + tmp_144*tmp_16;
      real_t tmp_147 = tmp_142*tmp_21 + tmp_143*tmp_18 + tmp_144*tmp_15;
      real_t tmp_148 = -tmp_145 - tmp_146 - tmp_147 + 1;
      real_t tmp_149 = tmp_141*tmp_83;
      real_t tmp_150 = 0.054975871827660928*tmp_85;
      real_t tmp_151 = 0.44594849091596489*tmp_43 + 0.44594849091596489*tmp_44;
      real_t tmp_152 = tmp_40*(tmp_151 + tmp_41);
      real_t tmp_153 = 0.44594849091596489*tmp_50 + 0.44594849091596489*tmp_51;
      real_t tmp_154 = tmp_40*(tmp_153 + tmp_48);
      real_t tmp_155 = 0.44594849091596489*tmp_57 + 0.44594849091596489*tmp_58;
      real_t tmp_156 = tmp_40*(tmp_155 + tmp_55);
      real_t tmp_157 = tmp_25*(tmp_152*tmp_32 + tmp_154*tmp_47 + tmp_156*tmp_54 - 1.0/4.0) + tmp_26*(tmp_152*tmp_61 + tmp_154*tmp_62 + tmp_156*tmp_63 - 1.0/4.0) + tmp_29*(tmp_152*tmp_64 + tmp_154*tmp_65 + tmp_156*tmp_66 - 1.0/4.0);
      real_t tmp_158 = tmp_151 + tmp_68;
      real_t tmp_159 = tmp_153 + tmp_70;
      real_t tmp_160 = tmp_155 + tmp_72;
      real_t tmp_161 = tmp_158*tmp_23 + tmp_159*tmp_20 + tmp_160*tmp_17;
      real_t tmp_162 = tmp_158*tmp_22 + tmp_159*tmp_19 + tmp_16*tmp_160;
      real_t tmp_163 = tmp_15*tmp_160 + tmp_158*tmp_21 + tmp_159*tmp_18;
      real_t tmp_164 = -tmp_161 - tmp_162 - tmp_163 + 1;
      real_t tmp_165 = tmp_157*tmp_83;
      real_t tmp_166 = 0.11169079483900572*tmp_85;
      real_t tmp_167 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_168 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_169 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t a_0_0 = tmp_102*(-tmp_100*tmp_101 - tmp_100*tmp_81 + 0.5*tmp_24*tmp_93) + tmp_118*(0.5*tmp_109*tmp_24 - tmp_116*tmp_117 - tmp_116*tmp_81) + tmp_134*(0.5*tmp_125*tmp_24 - tmp_132*tmp_133 - tmp_132*tmp_81) + tmp_150*(0.5*tmp_141*tmp_24 - tmp_148*tmp_149 - tmp_148*tmp_81) + tmp_166*(0.5*tmp_157*tmp_24 - tmp_164*tmp_165 - tmp_164*tmp_81) + tmp_86*(0.5*tmp_24*tmp_67 - tmp_77*tmp_81 - tmp_77*tmp_84);
      real_t a_1_0 = tmp_102*(-tmp_101*tmp_97 + 0.5*tmp_167*tmp_93 - tmp_81*tmp_97) + tmp_118*(0.5*tmp_109*tmp_167 - tmp_113*tmp_117 - tmp_113*tmp_81) + tmp_134*(0.5*tmp_125*tmp_167 - tmp_129*tmp_133 - tmp_129*tmp_81) + tmp_150*(0.5*tmp_141*tmp_167 - tmp_145*tmp_149 - tmp_145*tmp_81) + tmp_166*(0.5*tmp_157*tmp_167 - tmp_161*tmp_165 - tmp_161*tmp_81) + tmp_86*(0.5*tmp_167*tmp_67 - tmp_74*tmp_81 - tmp_74*tmp_84);
      real_t a_2_0 = tmp_102*(-tmp_101*tmp_98 + 0.5*tmp_168*tmp_93 - tmp_81*tmp_98) + tmp_118*(0.5*tmp_109*tmp_168 - tmp_114*tmp_117 - tmp_114*tmp_81) + tmp_134*(0.5*tmp_125*tmp_168 - tmp_130*tmp_133 - tmp_130*tmp_81) + tmp_150*(0.5*tmp_141*tmp_168 - tmp_146*tmp_149 - tmp_146*tmp_81) + tmp_166*(0.5*tmp_157*tmp_168 - tmp_162*tmp_165 - tmp_162*tmp_81) + tmp_86*(0.5*tmp_168*tmp_67 - tmp_75*tmp_81 - tmp_75*tmp_84);
      real_t a_3_0 = tmp_102*(-tmp_101*tmp_99 + 0.5*tmp_169*tmp_93 - tmp_81*tmp_99) + tmp_118*(0.5*tmp_109*tmp_169 - tmp_115*tmp_117 - tmp_115*tmp_81) + tmp_134*(0.5*tmp_125*tmp_169 - tmp_131*tmp_133 - tmp_131*tmp_81) + tmp_150*(0.5*tmp_141*tmp_169 - tmp_147*tmp_149 - tmp_147*tmp_81) + tmp_166*(0.5*tmp_157*tmp_169 - tmp_163*tmp_165 - tmp_163*tmp_81) + tmp_86*(0.5*tmp_169*tmp_67 - tmp_76*tmp_81 - tmp_76*tmp_84);
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
                                        MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_26 = -tmp_25;
      real_t tmp_27 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_28 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_29 = 0.091576213509770743*tmp_26 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -tmp_30;
      real_t tmp_32 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_33 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_34 = 0.091576213509770743*tmp_31 + 0.81684757298045851*tmp_32 + tmp_33;
      real_t tmp_35 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_36 = -tmp_35;
      real_t tmp_37 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_38 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_39 = 0.091576213509770743*tmp_36 + 0.81684757298045851*tmp_37 + tmp_38;
      real_t tmp_40 = tmp_17*tmp_39 + tmp_20*tmp_34 + tmp_23*tmp_29;
      real_t tmp_41 = tmp_16*tmp_39 + tmp_19*tmp_34 + tmp_22*tmp_29;
      real_t tmp_42 = tmp_15*tmp_39 + tmp_18*tmp_34 + tmp_21*tmp_29;
      real_t tmp_43 = tmp_11*(tmp_42 - 1.0/4.0) + tmp_6*(tmp_40 - 1.0/4.0) + tmp_9*(tmp_41 - 1.0/4.0);
      real_t tmp_44 = p_affine_13_0*(tmp_11*tmp_15 + tmp_16*tmp_9 + tmp_17*tmp_6) + p_affine_13_1*(tmp_11*tmp_18 + tmp_19*tmp_9 + tmp_20*tmp_6) + p_affine_13_2*(tmp_11*tmp_21 + tmp_22*tmp_9 + tmp_23*tmp_6);
      real_t tmp_45 = -tmp_40 - tmp_41 - tmp_42 + 1;
      real_t tmp_46 = (std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)*std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)) + (std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)*std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)) + (std::abs(tmp_30*tmp_37 - tmp_32*tmp_35)*std::abs(tmp_30*tmp_37 - tmp_32*tmp_35));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_26 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_51 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32 + tmp_33;
      real_t tmp_52 = 0.44594849091596489*tmp_36 + 0.10810301816807022*tmp_37 + tmp_38;
      real_t tmp_53 = tmp_17*tmp_52 + tmp_20*tmp_51 + tmp_23*tmp_50;
      real_t tmp_54 = tmp_16*tmp_52 + tmp_19*tmp_51 + tmp_22*tmp_50;
      real_t tmp_55 = tmp_15*tmp_52 + tmp_18*tmp_51 + tmp_21*tmp_50;
      real_t tmp_56 = tmp_11*(tmp_55 - 1.0/4.0) + tmp_6*(tmp_53 - 1.0/4.0) + tmp_9*(tmp_54 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_60 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_61 = 0.81684757298045851*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_62 = tmp_17*tmp_61 + tmp_20*tmp_60 + tmp_23*tmp_59;
      real_t tmp_63 = tmp_16*tmp_61 + tmp_19*tmp_60 + tmp_22*tmp_59;
      real_t tmp_64 = tmp_15*tmp_61 + tmp_18*tmp_60 + tmp_21*tmp_59;
      real_t tmp_65 = tmp_11*(tmp_64 - 1.0/4.0) + tmp_6*(tmp_62 - 1.0/4.0) + tmp_9*(tmp_63 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_69 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_70 = 0.10810301816807022*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_71 = tmp_17*tmp_70 + tmp_20*tmp_69 + tmp_23*tmp_68;
      real_t tmp_72 = tmp_16*tmp_70 + tmp_19*tmp_69 + tmp_22*tmp_68;
      real_t tmp_73 = tmp_15*tmp_70 + tmp_18*tmp_69 + tmp_21*tmp_68;
      real_t tmp_74 = tmp_11*(tmp_73 - 1.0/4.0) + tmp_6*(tmp_71 - 1.0/4.0) + tmp_9*(tmp_72 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_78 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_79 = 0.091576213509770743*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_80 = tmp_17*tmp_79 + tmp_20*tmp_78 + tmp_23*tmp_77;
      real_t tmp_81 = tmp_16*tmp_79 + tmp_19*tmp_78 + tmp_22*tmp_77;
      real_t tmp_82 = tmp_15*tmp_79 + tmp_18*tmp_78 + tmp_21*tmp_77;
      real_t tmp_83 = tmp_11*(tmp_82 - 1.0/4.0) + tmp_6*(tmp_80 - 1.0/4.0) + tmp_9*(tmp_81 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_87 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_88 = 0.44594849091596489*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_89 = tmp_17*tmp_88 + tmp_20*tmp_87 + tmp_23*tmp_86;
      real_t tmp_90 = tmp_16*tmp_88 + tmp_19*tmp_87 + tmp_22*tmp_86;
      real_t tmp_91 = tmp_15*tmp_88 + tmp_18*tmp_87 + tmp_21*tmp_86;
      real_t tmp_92 = tmp_11*(tmp_91 - 1.0/4.0) + tmp_6*(tmp_89 - 1.0/4.0) + tmp_9*(tmp_90 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_96 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_97 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t a_0_0 = tmp_49*(-tmp_24*tmp_43 + 3.0*tmp_43*tmp_45*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_24*tmp_56 - tmp_44*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_24*tmp_65 - tmp_44*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_24*tmp_74 - tmp_44*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_24*tmp_83 - tmp_44*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_24*tmp_92 - tmp_44*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_1_0 = tmp_49*(3.0*tmp_40*tmp_43*tmp_47 - tmp_40*tmp_44 - tmp_43*tmp_95) + tmp_58*(-tmp_44*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_44*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_44*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_44*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_44*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_2_0 = tmp_49*(3.0*tmp_41*tmp_43*tmp_47 - tmp_41*tmp_44 - tmp_43*tmp_96) + tmp_58*(-tmp_44*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_44*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_44*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_44*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_44*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_3_0 = tmp_49*(3.0*tmp_42*tmp_43*tmp_47 - tmp_42*tmp_44 - tmp_43*tmp_97) + tmp_58*(-tmp_44*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_44*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_44*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_44*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_44*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }

public:



};




class EGVectorLaplaceFormNitscheBC_P1P1_11 : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_P1P1_11()
: callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g1 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_g1( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g1( Point3D( {in_0, in_1, in_2} ) );
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = 1.0 / (tmp_4);
      real_t tmp_6 = tmp_0*tmp_5;
      real_t tmp_7 = tmp_2*tmp_5;
      real_t tmp_8 = -tmp_6 - tmp_7;
      real_t tmp_9 = tmp_1*tmp_5;
      real_t tmp_10 = tmp_3*tmp_5;
      real_t tmp_11 = -tmp_10 - tmp_9;
      real_t tmp_12 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_13 = tmp_12*((tmp_11*tmp_11) + (tmp_8*tmp_8));
      real_t tmp_14 = tmp_12*(tmp_11*tmp_9 + tmp_7*tmp_8);
      real_t tmp_15 = 0.5*tmp_14;
      real_t tmp_16 = tmp_12*(tmp_10*tmp_11 + tmp_6*tmp_8);
      real_t tmp_17 = 0.5*tmp_16;
      real_t tmp_18 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_19 = tmp_12*((tmp_1*tmp_1)*tmp_18 + tmp_18*(tmp_2*tmp_2));
      real_t tmp_20 = tmp_12*(tmp_0*tmp_18*tmp_2 + tmp_1*tmp_18*tmp_3);
      real_t tmp_21 = 0.5*tmp_20;
      real_t tmp_22 = tmp_12*((tmp_0*tmp_0)*tmp_18 + tmp_18*(tmp_3*tmp_3));
      real_t a_0_0 = 0.5*tmp_13;
      real_t a_0_1 = tmp_15;
      real_t a_0_2 = tmp_17;
      real_t a_1_0 = tmp_15;
      real_t a_1_1 = 0.5*tmp_19;
      real_t a_1_2 = tmp_21;
      real_t a_2_0 = tmp_17;
      real_t a_2_1 = tmp_21;
      real_t a_2_2 = 0.5*tmp_22;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_5 = 0.21132486540518713*tmp_1 + tmp_4;
      real_t tmp_6 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_7 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_8 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = 1.0 / (-tmp_6*tmp_9 + tmp_7*tmp_8);
      real_t tmp_11 = tmp_10*tmp_6;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_0 + tmp_12;
      real_t tmp_14 = tmp_10*tmp_8;
      real_t tmp_15 = tmp_11*tmp_5 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_10*tmp_7;
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_5;
      real_t tmp_19 = -tmp_15 - tmp_18 + 1;
      real_t tmp_20 = p_affine_10_0*(-tmp_14 - tmp_17) + p_affine_10_1*(-tmp_11 - tmp_16);
      real_t tmp_21 = 1.0*tmp_20;
      real_t tmp_22 = 0.5*tmp_2;
      real_t tmp_23 = 0.78867513459481287*tmp_1 + tmp_4;
      real_t tmp_24 = 0.78867513459481287*tmp_0 + tmp_12;
      real_t tmp_25 = tmp_11*tmp_23 + tmp_14*tmp_24;
      real_t tmp_26 = tmp_16*tmp_23 + tmp_17*tmp_24;
      real_t tmp_27 = -tmp_25 - tmp_26 + 1;
      real_t tmp_28 = 0.5*tmp_2;
      real_t tmp_29 = 0.5*tmp_20;
      real_t tmp_30 = p_affine_10_0*tmp_14 + p_affine_10_1*tmp_11;
      real_t tmp_31 = 0.5*tmp_30;
      real_t tmp_32 = tmp_22*(3*tmp_15*tmp_19*tmp_3 - tmp_15*tmp_29 - tmp_19*tmp_31) + tmp_28*(3*tmp_25*tmp_27*tmp_3 - tmp_25*tmp_29 - tmp_27*tmp_31);
      real_t tmp_33 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_16;
      real_t tmp_34 = 0.5*tmp_33;
      real_t tmp_35 = tmp_22*(3*tmp_18*tmp_19*tmp_3 - tmp_18*tmp_29 - tmp_19*tmp_34) + tmp_28*(3*tmp_26*tmp_27*tmp_3 - tmp_26*tmp_29 - tmp_27*tmp_34);
      real_t tmp_36 = 1.0*tmp_30;
      real_t tmp_37 = tmp_22*(3*tmp_15*tmp_18*tmp_3 - tmp_15*tmp_34 - tmp_18*tmp_31) + tmp_28*(3*tmp_25*tmp_26*tmp_3 - tmp_25*tmp_34 - tmp_26*tmp_31);
      real_t tmp_38 = 1.0*tmp_33;
      real_t a_0_0 = tmp_22*(3*(tmp_19*tmp_19)*tmp_3 - tmp_19*tmp_21) + tmp_28*(-tmp_21*tmp_27 + 3*(tmp_27*tmp_27)*tmp_3);
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_22*(3*(tmp_15*tmp_15)*tmp_3 - tmp_15*tmp_36) + tmp_28*(3*(tmp_25*tmp_25)*tmp_3 - tmp_25*tmp_36);
      real_t a_1_2 = tmp_37;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_37;
      real_t a_2_2 = tmp_22*(3*(tmp_18*tmp_18)*tmp_3 - tmp_18*tmp_38) + tmp_28*(3*(tmp_26*tmp_26)*tmp_3 - tmp_26*tmp_38);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
                                          MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = 1.0 / (tmp_0*tmp_1 - tmp_2*tmp_3);
      real_t tmp_5 = tmp_0*tmp_4;
      real_t tmp_6 = tmp_3*tmp_4;
      real_t tmp_7 = tmp_1*tmp_4;
      real_t tmp_8 = tmp_2*tmp_4;
      real_t tmp_9 = p_affine_10_0*(-tmp_5 - tmp_6) + p_affine_10_1*(-tmp_7 - tmp_8);
      real_t tmp_10 = -p_affine_3_1;
      real_t tmp_11 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_12 = p_affine_6_1 + 0.21132486540518713*tmp_11;
      real_t tmp_13 = tmp_10 + tmp_12;
      real_t tmp_14 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_15 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_16 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_17 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_18 = 1.0 / (-tmp_14*tmp_17 + tmp_15*tmp_16);
      real_t tmp_19 = tmp_14*tmp_18;
      real_t tmp_20 = -p_affine_3_0;
      real_t tmp_21 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_22 = p_affine_6_0 + 0.21132486540518713*tmp_21;
      real_t tmp_23 = tmp_20 + tmp_22;
      real_t tmp_24 = tmp_16*tmp_18;
      real_t tmp_25 = tmp_13*tmp_19 + tmp_23*tmp_24;
      real_t tmp_26 = tmp_15*tmp_18;
      real_t tmp_27 = tmp_17*tmp_18;
      real_t tmp_28 = tmp_13*tmp_26 + tmp_23*tmp_27;
      real_t tmp_29 = -tmp_25 - tmp_28 + 1;
      real_t tmp_30 = -p_affine_0_1;
      real_t tmp_31 = tmp_12 + tmp_30;
      real_t tmp_32 = -p_affine_0_0;
      real_t tmp_33 = tmp_22 + tmp_32;
      real_t tmp_34 = tmp_31*tmp_8 + tmp_33*tmp_5;
      real_t tmp_35 = tmp_31*tmp_7 + tmp_33*tmp_6;
      real_t tmp_36 = -tmp_34 - tmp_35 + 1;
      real_t tmp_37 = 0.5*p_affine_10_0*(-tmp_24 - tmp_27) + 0.5*p_affine_10_1*(-tmp_19 - tmp_26);
      real_t tmp_38 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_21*tmp_21), 1.0/2.0));
      real_t tmp_39 = 3/tmp_38;
      real_t tmp_40 = tmp_36*tmp_39;
      real_t tmp_41 = 0.5*tmp_38;
      real_t tmp_42 = p_affine_6_1 + 0.78867513459481287*tmp_11;
      real_t tmp_43 = tmp_10 + tmp_42;
      real_t tmp_44 = p_affine_6_0 + 0.78867513459481287*tmp_21;
      real_t tmp_45 = tmp_20 + tmp_44;
      real_t tmp_46 = tmp_19*tmp_43 + tmp_24*tmp_45;
      real_t tmp_47 = tmp_26*tmp_43 + tmp_27*tmp_45;
      real_t tmp_48 = -tmp_46 - tmp_47 + 1;
      real_t tmp_49 = tmp_30 + tmp_42;
      real_t tmp_50 = tmp_32 + tmp_44;
      real_t tmp_51 = tmp_49*tmp_8 + tmp_5*tmp_50;
      real_t tmp_52 = tmp_49*tmp_7 + tmp_50*tmp_6;
      real_t tmp_53 = -tmp_51 - tmp_52 + 1;
      real_t tmp_54 = tmp_39*tmp_53;
      real_t tmp_55 = 0.5*tmp_38;
      real_t tmp_56 = 0.5*p_affine_10_0*tmp_24 + 0.5*p_affine_10_1*tmp_19;
      real_t tmp_57 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*tmp_26;
      real_t tmp_58 = p_affine_10_0*tmp_5 + p_affine_10_1*tmp_8;
      real_t tmp_59 = tmp_34*tmp_39;
      real_t tmp_60 = tmp_39*tmp_51;
      real_t tmp_61 = p_affine_10_0*tmp_6 + p_affine_10_1*tmp_7;
      real_t tmp_62 = tmp_35*tmp_39;
      real_t tmp_63 = tmp_39*tmp_52;
      real_t a_0_0 = tmp_41*(-tmp_29*tmp_40 + 0.5*tmp_29*tmp_9 - tmp_36*tmp_37) + tmp_55*(-tmp_37*tmp_53 - tmp_48*tmp_54 + 0.5*tmp_48*tmp_9);
      real_t a_0_1 = tmp_41*(-tmp_25*tmp_40 + 0.5*tmp_25*tmp_9 - tmp_36*tmp_56) + tmp_55*(-tmp_46*tmp_54 + 0.5*tmp_46*tmp_9 - tmp_53*tmp_56);
      real_t a_0_2 = tmp_41*(-tmp_28*tmp_40 + 0.5*tmp_28*tmp_9 - tmp_36*tmp_57) + tmp_55*(-tmp_47*tmp_54 + 0.5*tmp_47*tmp_9 - tmp_53*tmp_57);
      real_t a_1_0 = tmp_41*(0.5*tmp_29*tmp_58 - tmp_29*tmp_59 - tmp_34*tmp_37) + tmp_55*(-tmp_37*tmp_51 + 0.5*tmp_48*tmp_58 - tmp_48*tmp_60);
      real_t a_1_1 = tmp_41*(0.5*tmp_25*tmp_58 - tmp_25*tmp_59 - tmp_34*tmp_56) + tmp_55*(0.5*tmp_46*tmp_58 - tmp_46*tmp_60 - tmp_51*tmp_56);
      real_t a_1_2 = tmp_41*(0.5*tmp_28*tmp_58 - tmp_28*tmp_59 - tmp_34*tmp_57) + tmp_55*(0.5*tmp_47*tmp_58 - tmp_47*tmp_60 - tmp_51*tmp_57);
      real_t a_2_0 = tmp_41*(0.5*tmp_29*tmp_61 - tmp_29*tmp_62 - tmp_35*tmp_37) + tmp_55*(-tmp_37*tmp_52 + 0.5*tmp_48*tmp_61 - tmp_48*tmp_63);
      real_t a_2_1 = tmp_41*(0.5*tmp_25*tmp_61 - tmp_25*tmp_62 - tmp_35*tmp_56) + tmp_55*(0.5*tmp_46*tmp_61 - tmp_46*tmp_63 - tmp_52*tmp_56);
      real_t a_2_2 = tmp_41*(0.5*tmp_28*tmp_61 - tmp_28*tmp_62 - tmp_35*tmp_57) + tmp_55*(0.5*tmp_47*tmp_61 - tmp_47*tmp_63 - tmp_52*tmp_57);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_5 = 0.21132486540518713*tmp_1 + tmp_4;
      real_t tmp_6 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_7 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_8 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = 1.0 / (-tmp_6*tmp_9 + tmp_7*tmp_8);
      real_t tmp_11 = tmp_10*tmp_6;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_0 + tmp_12;
      real_t tmp_14 = tmp_10*tmp_8;
      real_t tmp_15 = tmp_11*tmp_5 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_10*tmp_7;
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_5;
      real_t tmp_19 = -tmp_15 - tmp_18 + 1;
      real_t tmp_20 = p_affine_10_0*(-tmp_14 - tmp_17) + p_affine_10_1*(-tmp_11 - tmp_16);
      real_t tmp_21 = 2*tmp_20;
      real_t tmp_22 = 0.5*tmp_2;
      real_t tmp_23 = 0.78867513459481287*tmp_1 + tmp_4;
      real_t tmp_24 = 0.78867513459481287*tmp_0 + tmp_12;
      real_t tmp_25 = tmp_11*tmp_23 + tmp_14*tmp_24;
      real_t tmp_26 = tmp_16*tmp_23 + tmp_17*tmp_24;
      real_t tmp_27 = -tmp_25 - tmp_26 + 1;
      real_t tmp_28 = 0.5*tmp_2;
      real_t tmp_29 = p_affine_10_0*tmp_14 + p_affine_10_1*tmp_11;
      real_t tmp_30 = tmp_22*(3*tmp_15*tmp_19*tmp_3 - tmp_15*tmp_20 - tmp_19*tmp_29) + tmp_28*(-tmp_20*tmp_25 + 3*tmp_25*tmp_27*tmp_3 - tmp_27*tmp_29);
      real_t tmp_31 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_16;
      real_t tmp_32 = tmp_22*(3*tmp_18*tmp_19*tmp_3 - tmp_18*tmp_20 - tmp_19*tmp_31) + tmp_28*(-tmp_20*tmp_26 + 3*tmp_26*tmp_27*tmp_3 - tmp_27*tmp_31);
      real_t tmp_33 = 2*tmp_29;
      real_t tmp_34 = tmp_22*(3*tmp_15*tmp_18*tmp_3 - tmp_15*tmp_31 - tmp_18*tmp_29) + tmp_28*(3*tmp_25*tmp_26*tmp_3 - tmp_25*tmp_31 - tmp_26*tmp_29);
      real_t tmp_35 = 2*tmp_31;
      real_t a_0_0 = tmp_22*(3*(tmp_19*tmp_19)*tmp_3 - tmp_19*tmp_21) + tmp_28*(-tmp_21*tmp_27 + 3*(tmp_27*tmp_27)*tmp_3);
      real_t a_0_1 = tmp_30;
      real_t a_0_2 = tmp_32;
      real_t a_1_0 = tmp_30;
      real_t a_1_1 = tmp_22*(3*(tmp_15*tmp_15)*tmp_3 - tmp_15*tmp_33) + tmp_28*(3*(tmp_25*tmp_25)*tmp_3 - tmp_25*tmp_33);
      real_t a_1_2 = tmp_34;
      real_t a_2_0 = tmp_32;
      real_t a_2_1 = tmp_34;
      real_t a_2_2 = tmp_22*(3*(tmp_18*tmp_18)*tmp_3 - tmp_18*tmp_35) + tmp_28*(3*(tmp_26*tmp_26)*tmp_3 - tmp_26*tmp_35);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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

      real_t Scalar_Variable_Coefficient_2D_g1_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id1 = 0;
      Scalar_Variable_Coefficient_2D_g1( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id0 );
      Scalar_Variable_Coefficient_2D_g1( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id1 );
      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_5 = 0.21132486540518713*tmp_1 + tmp_4;
      real_t tmp_6 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_7 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_8 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = 1.0 / (-tmp_6*tmp_9 + tmp_7*tmp_8);
      real_t tmp_11 = tmp_10*tmp_6;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_0 + tmp_12;
      real_t tmp_14 = tmp_10*tmp_8;
      real_t tmp_15 = tmp_11*tmp_5 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_10*tmp_7;
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_5;
      real_t tmp_19 = p_affine_10_0*(-tmp_14 - tmp_17) + p_affine_10_1*(-tmp_11 - tmp_16);
      real_t tmp_20 = 0.5*Scalar_Variable_Coefficient_2D_g1_out0_id0*tmp_2;
      real_t tmp_21 = 0.78867513459481287*tmp_1 + tmp_4;
      real_t tmp_22 = 0.78867513459481287*tmp_0 + tmp_12;
      real_t tmp_23 = tmp_11*tmp_21 + tmp_14*tmp_22;
      real_t tmp_24 = tmp_16*tmp_21 + tmp_17*tmp_22;
      real_t tmp_25 = 0.5*Scalar_Variable_Coefficient_2D_g1_out0_id1*tmp_2;
      real_t tmp_26 = p_affine_10_0*tmp_14 + p_affine_10_1*tmp_11;
      real_t tmp_27 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_16;
      real_t a_0_0 = tmp_20*(-tmp_19 + 3*tmp_3*(-tmp_15 - tmp_18 + 1)) + tmp_25*(-tmp_19 + 3*tmp_3*(-tmp_23 - tmp_24 + 1));
      real_t a_1_0 = tmp_20*(3*tmp_15*tmp_3 - tmp_26) + tmp_25*(3*tmp_23*tmp_3 - tmp_26);
      real_t a_2_0 = tmp_20*(3*tmp_18*tmp_3 - tmp_27) + tmp_25*(3*tmp_24*tmp_3 - tmp_27);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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

      real_t Scalar_Variable_Coefficient_3D_g1_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_g1( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id0 );
      Scalar_Variable_Coefficient_3D_g1( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id1 );
      Scalar_Variable_Coefficient_3D_g1( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id2 );
      Scalar_Variable_Coefficient_3D_g1( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id3 );
      Scalar_Variable_Coefficient_3D_g1( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id4 );
      Scalar_Variable_Coefficient_3D_g1( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id5 );
      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -tmp_4;
      real_t tmp_9 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_10 = 0.81684757298045851*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_11 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_13 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_14 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_15 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_16 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_19 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_20 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_12*tmp_20;
      real_t tmp_23 = tmp_16*tmp_19;
      real_t tmp_24 = tmp_14*tmp_18;
      real_t tmp_25 = 1.0 / (tmp_11*tmp_12*tmp_18 - tmp_11*tmp_23 + tmp_13*tmp_21 - tmp_13*tmp_24 + tmp_15*tmp_17 - tmp_15*tmp_22);
      real_t tmp_26 = tmp_25*(tmp_11*tmp_12 - tmp_13*tmp_14);
      real_t tmp_27 = -tmp_1;
      real_t tmp_28 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_29 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_30 = tmp_25*(-tmp_11*tmp_16 + tmp_13*tmp_20);
      real_t tmp_31 = -tmp_3;
      real_t tmp_32 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_33 = 0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_34 = tmp_25*(tmp_17 - tmp_22);
      real_t tmp_35 = tmp_10*tmp_26 + tmp_29*tmp_30 + tmp_33*tmp_34;
      real_t tmp_36 = tmp_25*(-tmp_12*tmp_15 + tmp_13*tmp_19);
      real_t tmp_37 = tmp_25*(-tmp_13*tmp_18 + tmp_15*tmp_16);
      real_t tmp_38 = tmp_25*(tmp_12*tmp_18 - tmp_23);
      real_t tmp_39 = tmp_10*tmp_36 + tmp_29*tmp_37 + tmp_33*tmp_38;
      real_t tmp_40 = tmp_25*(-tmp_11*tmp_19 + tmp_14*tmp_15);
      real_t tmp_41 = tmp_25*(tmp_11*tmp_18 - tmp_15*tmp_20);
      real_t tmp_42 = tmp_25*(tmp_21 - tmp_24);
      real_t tmp_43 = tmp_10*tmp_40 + tmp_29*tmp_41 + tmp_33*tmp_42;
      real_t tmp_44 = p_affine_13_0*(-tmp_34 - tmp_38 - tmp_42) + p_affine_13_1*(-tmp_30 - tmp_37 - tmp_41) + p_affine_13_2*(-tmp_26 - tmp_36 - tmp_40);
      real_t tmp_45 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_46 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_g1_out0_id0*tmp_45;
      real_t tmp_47 = 0.10810301816807022*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_48 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_49 = 0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_50 = tmp_26*tmp_47 + tmp_30*tmp_48 + tmp_34*tmp_49;
      real_t tmp_51 = tmp_36*tmp_47 + tmp_37*tmp_48 + tmp_38*tmp_49;
      real_t tmp_52 = tmp_40*tmp_47 + tmp_41*tmp_48 + tmp_42*tmp_49;
      real_t tmp_53 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_g1_out0_id1*tmp_45;
      real_t tmp_54 = 0.091576213509770743*tmp_5 + 0.81684757298045851*tmp_8 + tmp_9;
      real_t tmp_55 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_56 = 0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_31 + tmp_32;
      real_t tmp_57 = tmp_26*tmp_54 + tmp_30*tmp_55 + tmp_34*tmp_56;
      real_t tmp_58 = tmp_36*tmp_54 + tmp_37*tmp_55 + tmp_38*tmp_56;
      real_t tmp_59 = tmp_40*tmp_54 + tmp_41*tmp_55 + tmp_42*tmp_56;
      real_t tmp_60 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_g1_out0_id2*tmp_45;
      real_t tmp_61 = 0.44594849091596489*tmp_5 + 0.10810301816807022*tmp_8 + tmp_9;
      real_t tmp_62 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_63 = 0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_31 + tmp_32;
      real_t tmp_64 = tmp_26*tmp_61 + tmp_30*tmp_62 + tmp_34*tmp_63;
      real_t tmp_65 = tmp_36*tmp_61 + tmp_37*tmp_62 + tmp_38*tmp_63;
      real_t tmp_66 = tmp_40*tmp_61 + tmp_41*tmp_62 + tmp_42*tmp_63;
      real_t tmp_67 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_g1_out0_id3*tmp_45;
      real_t tmp_68 = 0.091576213509770743*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_69 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_70 = 0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_71 = tmp_26*tmp_68 + tmp_30*tmp_69 + tmp_34*tmp_70;
      real_t tmp_72 = tmp_36*tmp_68 + tmp_37*tmp_69 + tmp_38*tmp_70;
      real_t tmp_73 = tmp_40*tmp_68 + tmp_41*tmp_69 + tmp_42*tmp_70;
      real_t tmp_74 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_g1_out0_id4*tmp_45;
      real_t tmp_75 = 0.44594849091596489*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_76 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_77 = 0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_78 = tmp_26*tmp_75 + tmp_30*tmp_76 + tmp_34*tmp_77;
      real_t tmp_79 = tmp_36*tmp_75 + tmp_37*tmp_76 + tmp_38*tmp_77;
      real_t tmp_80 = tmp_40*tmp_75 + tmp_41*tmp_76 + tmp_42*tmp_77;
      real_t tmp_81 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_g1_out0_id5*tmp_45;
      real_t tmp_82 = p_affine_13_0*tmp_34 + p_affine_13_1*tmp_30 + p_affine_13_2*tmp_26;
      real_t tmp_83 = p_affine_13_0*tmp_38 + p_affine_13_1*tmp_37 + p_affine_13_2*tmp_36;
      real_t tmp_84 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_40;
      real_t a_0_0 = tmp_46*(-tmp_44 + 3.0*tmp_7*(-tmp_35 - tmp_39 - tmp_43 + 1)) + tmp_53*(-tmp_44 + 3.0*tmp_7*(-tmp_50 - tmp_51 - tmp_52 + 1)) + tmp_60*(-tmp_44 + 3.0*tmp_7*(-tmp_57 - tmp_58 - tmp_59 + 1)) + tmp_67*(-tmp_44 + 3.0*tmp_7*(-tmp_64 - tmp_65 - tmp_66 + 1)) + tmp_74*(-tmp_44 + 3.0*tmp_7*(-tmp_71 - tmp_72 - tmp_73 + 1)) + tmp_81*(-tmp_44 + 3.0*tmp_7*(-tmp_78 - tmp_79 - tmp_80 + 1));
      real_t a_1_0 = tmp_46*(3.0*tmp_35*tmp_7 - tmp_82) + tmp_53*(3.0*tmp_50*tmp_7 - tmp_82) + tmp_60*(3.0*tmp_57*tmp_7 - tmp_82) + tmp_67*(3.0*tmp_64*tmp_7 - tmp_82) + tmp_74*(3.0*tmp_7*tmp_71 - tmp_82) + tmp_81*(3.0*tmp_7*tmp_78 - tmp_82);
      real_t a_2_0 = tmp_46*(3.0*tmp_39*tmp_7 - tmp_83) + tmp_53*(3.0*tmp_51*tmp_7 - tmp_83) + tmp_60*(3.0*tmp_58*tmp_7 - tmp_83) + tmp_67*(3.0*tmp_65*tmp_7 - tmp_83) + tmp_74*(3.0*tmp_7*tmp_72 - tmp_83) + tmp_81*(3.0*tmp_7*tmp_79 - tmp_83);
      real_t a_3_0 = tmp_46*(3.0*tmp_43*tmp_7 - tmp_84) + tmp_53*(3.0*tmp_52*tmp_7 - tmp_84) + tmp_60*(3.0*tmp_59*tmp_7 - tmp_84) + tmp_67*(3.0*tmp_66*tmp_7 - tmp_84) + tmp_74*(3.0*tmp_7*tmp_73 - tmp_84) + tmp_81*(3.0*tmp_7*tmp_80 - tmp_84);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = tmp_2 - tmp_5;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_10 = tmp_3*tmp_9;
      real_t tmp_11 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_12 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_1*tmp_12;
      real_t tmp_15 = tmp_10*tmp_8 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_13 - tmp_14*tmp_8 + tmp_2*tmp_7 - tmp_5*tmp_7;
      real_t tmp_16 = 1.0 / (tmp_15);
      real_t tmp_17 = tmp_16*tmp_6;
      real_t tmp_18 = tmp_12*tmp_4 - tmp_13;
      real_t tmp_19 = tmp_16*tmp_18;
      real_t tmp_20 = tmp_10 - tmp_14;
      real_t tmp_21 = tmp_16*tmp_20;
      real_t tmp_22 = -tmp_17 - tmp_19 - tmp_21;
      real_t tmp_23 = -tmp_0*tmp_11 + tmp_3*tmp_8;
      real_t tmp_24 = tmp_16*tmp_23;
      real_t tmp_25 = tmp_0*tmp_7 - tmp_12*tmp_8;
      real_t tmp_26 = tmp_16*tmp_25;
      real_t tmp_27 = tmp_11*tmp_12 - tmp_3*tmp_7;
      real_t tmp_28 = tmp_16*tmp_27;
      real_t tmp_29 = -tmp_24 - tmp_26 - tmp_28;
      real_t tmp_30 = -tmp_1*tmp_8 + tmp_11*tmp_4;
      real_t tmp_31 = tmp_16*tmp_30;
      real_t tmp_32 = -tmp_4*tmp_7 + tmp_8*tmp_9;
      real_t tmp_33 = tmp_16*tmp_32;
      real_t tmp_34 = tmp_1*tmp_7 - tmp_11*tmp_9;
      real_t tmp_35 = tmp_16*tmp_34;
      real_t tmp_36 = -tmp_31 - tmp_33 - tmp_35;
      real_t tmp_37 = p_affine_0_0*p_affine_1_1;
      real_t tmp_38 = p_affine_0_0*p_affine_1_2;
      real_t tmp_39 = p_affine_2_1*p_affine_3_2;
      real_t tmp_40 = p_affine_0_1*p_affine_1_0;
      real_t tmp_41 = p_affine_0_1*p_affine_1_2;
      real_t tmp_42 = p_affine_2_2*p_affine_3_0;
      real_t tmp_43 = p_affine_0_2*p_affine_1_0;
      real_t tmp_44 = p_affine_0_2*p_affine_1_1;
      real_t tmp_45 = p_affine_2_0*p_affine_3_1;
      real_t tmp_46 = p_affine_2_2*p_affine_3_1;
      real_t tmp_47 = p_affine_2_0*p_affine_3_2;
      real_t tmp_48 = p_affine_2_1*p_affine_3_0;
      real_t tmp_49 = std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = tmp_49*((tmp_22*tmp_22) + (tmp_29*tmp_29) + (tmp_36*tmp_36));
      real_t tmp_51 = tmp_49*(tmp_21*tmp_22 + tmp_28*tmp_29 + tmp_35*tmp_36);
      real_t tmp_52 = 0.16666666666666666*tmp_51;
      real_t tmp_53 = tmp_49*(tmp_19*tmp_22 + tmp_26*tmp_29 + tmp_33*tmp_36);
      real_t tmp_54 = 0.16666666666666666*tmp_53;
      real_t tmp_55 = tmp_49*(tmp_17*tmp_22 + tmp_24*tmp_29 + tmp_31*tmp_36);
      real_t tmp_56 = 0.16666666666666666*tmp_55;
      real_t tmp_57 = 1.0 / (tmp_15*tmp_15);
      real_t tmp_58 = tmp_49*((tmp_20*tmp_20)*tmp_57 + (tmp_27*tmp_27)*tmp_57 + (tmp_34*tmp_34)*tmp_57);
      real_t tmp_59 = tmp_20*tmp_57;
      real_t tmp_60 = tmp_27*tmp_57;
      real_t tmp_61 = tmp_34*tmp_57;
      real_t tmp_62 = tmp_49*(tmp_18*tmp_59 + tmp_25*tmp_60 + tmp_32*tmp_61);
      real_t tmp_63 = 0.16666666666666666*tmp_62;
      real_t tmp_64 = tmp_49*(tmp_23*tmp_60 + tmp_30*tmp_61 + tmp_59*tmp_6);
      real_t tmp_65 = 0.16666666666666666*tmp_64;
      real_t tmp_66 = tmp_49*((tmp_18*tmp_18)*tmp_57 + (tmp_25*tmp_25)*tmp_57 + (tmp_32*tmp_32)*tmp_57);
      real_t tmp_67 = tmp_49*(tmp_18*tmp_57*tmp_6 + tmp_23*tmp_25*tmp_57 + tmp_30*tmp_32*tmp_57);
      real_t tmp_68 = 0.16666666666666666*tmp_67;
      real_t tmp_69 = tmp_49*((tmp_23*tmp_23)*tmp_57 + (tmp_30*tmp_30)*tmp_57 + tmp_57*(tmp_6*tmp_6));
      real_t a_0_0 = 0.16666666666666666*tmp_50;
      real_t a_0_1 = tmp_52;
      real_t a_0_2 = tmp_54;
      real_t a_0_3 = tmp_56;
      real_t a_1_0 = tmp_52;
      real_t a_1_1 = 0.16666666666666666*tmp_58;
      real_t a_1_2 = tmp_63;
      real_t a_1_3 = tmp_65;
      real_t a_2_0 = tmp_54;
      real_t a_2_1 = tmp_63;
      real_t a_2_2 = 0.16666666666666666*tmp_66;
      real_t a_2_3 = tmp_68;
      real_t a_3_0 = tmp_56;
      real_t a_3_1 = tmp_65;
      real_t a_3_2 = tmp_68;
      real_t a_3_3 = 0.16666666666666666*tmp_69;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const override
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

         real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -tmp_4;
      real_t tmp_9 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_10 = 0.81684757298045851*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_11 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_13 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_14 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_15 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_16 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_19 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_20 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_12*tmp_16;
      real_t tmp_23 = tmp_11*tmp_19;
      real_t tmp_24 = tmp_13*tmp_18;
      real_t tmp_25 = 1.0 / (tmp_11*tmp_12*tmp_18 + tmp_13*tmp_21 - tmp_14*tmp_24 + tmp_15*tmp_17 - tmp_15*tmp_23 - tmp_20*tmp_22);
      real_t tmp_26 = tmp_25*(tmp_11*tmp_12 - tmp_13*tmp_14);
      real_t tmp_27 = -tmp_1;
      real_t tmp_28 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_29 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_30 = tmp_25*(-tmp_11*tmp_15 + tmp_13*tmp_20);
      real_t tmp_31 = -tmp_3;
      real_t tmp_32 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_33 = 0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_34 = tmp_25*(-tmp_12*tmp_20 + tmp_14*tmp_15);
      real_t tmp_35 = tmp_10*tmp_26 + tmp_29*tmp_30 + tmp_33*tmp_34;
      real_t tmp_36 = tmp_25*(tmp_13*tmp_19 - tmp_22);
      real_t tmp_37 = tmp_25*(tmp_15*tmp_16 - tmp_24);
      real_t tmp_38 = tmp_25*(tmp_12*tmp_18 - tmp_15*tmp_19);
      real_t tmp_39 = tmp_10*tmp_36 + tmp_29*tmp_37 + tmp_33*tmp_38;
      real_t tmp_40 = tmp_25*(tmp_17 - tmp_23);
      real_t tmp_41 = tmp_25*(tmp_11*tmp_18 - tmp_16*tmp_20);
      real_t tmp_42 = tmp_25*(-tmp_14*tmp_18 + tmp_21);
      real_t tmp_43 = tmp_10*tmp_40 + tmp_29*tmp_41 + tmp_33*tmp_42;
      real_t tmp_44 = -tmp_35 - tmp_39 - tmp_43 + 1;
      real_t tmp_45 = p_affine_13_0*(-tmp_34 - tmp_38 - tmp_42) + p_affine_13_1*(-tmp_30 - tmp_37 - tmp_41) + p_affine_13_2*(-tmp_26 - tmp_36 - tmp_40);
      real_t tmp_46 = 1.0*tmp_45;
      real_t tmp_47 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_48 = 0.054975871827660928*tmp_47;
      real_t tmp_49 = 0.10810301816807022*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_50 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_51 = 0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_52 = tmp_26*tmp_49 + tmp_30*tmp_50 + tmp_34*tmp_51;
      real_t tmp_53 = tmp_36*tmp_49 + tmp_37*tmp_50 + tmp_38*tmp_51;
      real_t tmp_54 = tmp_40*tmp_49 + tmp_41*tmp_50 + tmp_42*tmp_51;
      real_t tmp_55 = -tmp_52 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11169079483900572*tmp_47;
      real_t tmp_57 = 0.091576213509770743*tmp_5 + 0.81684757298045851*tmp_8 + tmp_9;
      real_t tmp_58 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_59 = 0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_31 + tmp_32;
      real_t tmp_60 = tmp_26*tmp_57 + tmp_30*tmp_58 + tmp_34*tmp_59;
      real_t tmp_61 = tmp_36*tmp_57 + tmp_37*tmp_58 + tmp_38*tmp_59;
      real_t tmp_62 = tmp_40*tmp_57 + tmp_41*tmp_58 + tmp_42*tmp_59;
      real_t tmp_63 = -tmp_60 - tmp_61 - tmp_62 + 1;
      real_t tmp_64 = 0.054975871827660928*tmp_47;
      real_t tmp_65 = 0.44594849091596489*tmp_5 + 0.10810301816807022*tmp_8 + tmp_9;
      real_t tmp_66 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_67 = 0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_31 + tmp_32;
      real_t tmp_68 = tmp_26*tmp_65 + tmp_30*tmp_66 + tmp_34*tmp_67;
      real_t tmp_69 = tmp_36*tmp_65 + tmp_37*tmp_66 + tmp_38*tmp_67;
      real_t tmp_70 = tmp_40*tmp_65 + tmp_41*tmp_66 + tmp_42*tmp_67;
      real_t tmp_71 = -tmp_68 - tmp_69 - tmp_70 + 1;
      real_t tmp_72 = 0.11169079483900572*tmp_47;
      real_t tmp_73 = 0.091576213509770743*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_74 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_75 = 0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_76 = tmp_26*tmp_73 + tmp_30*tmp_74 + tmp_34*tmp_75;
      real_t tmp_77 = tmp_36*tmp_73 + tmp_37*tmp_74 + tmp_38*tmp_75;
      real_t tmp_78 = tmp_40*tmp_73 + tmp_41*tmp_74 + tmp_42*tmp_75;
      real_t tmp_79 = -tmp_76 - tmp_77 - tmp_78 + 1;
      real_t tmp_80 = 0.054975871827660928*tmp_47;
      real_t tmp_81 = 0.44594849091596489*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_82 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_83 = 0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_84 = tmp_26*tmp_81 + tmp_30*tmp_82 + tmp_34*tmp_83;
      real_t tmp_85 = tmp_36*tmp_81 + tmp_37*tmp_82 + tmp_38*tmp_83;
      real_t tmp_86 = tmp_40*tmp_81 + tmp_41*tmp_82 + tmp_42*tmp_83;
      real_t tmp_87 = -tmp_84 - tmp_85 - tmp_86 + 1;
      real_t tmp_88 = 0.11169079483900572*tmp_47;
      real_t tmp_89 = 0.5*tmp_45;
      real_t tmp_90 = p_affine_13_0*tmp_34 + p_affine_13_1*tmp_30 + p_affine_13_2*tmp_26;
      real_t tmp_91 = 0.5*tmp_90;
      real_t tmp_92 = tmp_48*(3.0*tmp_35*tmp_44*tmp_7 - tmp_35*tmp_89 - tmp_44*tmp_91) + tmp_56*(3.0*tmp_52*tmp_55*tmp_7 - tmp_52*tmp_89 - tmp_55*tmp_91) + tmp_64*(3.0*tmp_60*tmp_63*tmp_7 - tmp_60*tmp_89 - tmp_63*tmp_91) + tmp_72*(3.0*tmp_68*tmp_7*tmp_71 - tmp_68*tmp_89 - tmp_71*tmp_91) + tmp_80*(3.0*tmp_7*tmp_76*tmp_79 - tmp_76*tmp_89 - tmp_79*tmp_91) + tmp_88*(3.0*tmp_7*tmp_84*tmp_87 - tmp_84*tmp_89 - tmp_87*tmp_91);
      real_t tmp_93 = p_affine_13_0*tmp_38 + p_affine_13_1*tmp_37 + p_affine_13_2*tmp_36;
      real_t tmp_94 = 0.5*tmp_93;
      real_t tmp_95 = tmp_48*(3.0*tmp_39*tmp_44*tmp_7 - tmp_39*tmp_89 - tmp_44*tmp_94) + tmp_56*(3.0*tmp_53*tmp_55*tmp_7 - tmp_53*tmp_89 - tmp_55*tmp_94) + tmp_64*(3.0*tmp_61*tmp_63*tmp_7 - tmp_61*tmp_89 - tmp_63*tmp_94) + tmp_72*(3.0*tmp_69*tmp_7*tmp_71 - tmp_69*tmp_89 - tmp_71*tmp_94) + tmp_80*(3.0*tmp_7*tmp_77*tmp_79 - tmp_77*tmp_89 - tmp_79*tmp_94) + tmp_88*(3.0*tmp_7*tmp_85*tmp_87 - tmp_85*tmp_89 - tmp_87*tmp_94);
      real_t tmp_96 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_40;
      real_t tmp_97 = 0.5*tmp_96;
      real_t tmp_98 = tmp_48*(3.0*tmp_43*tmp_44*tmp_7 - tmp_43*tmp_89 - tmp_44*tmp_97) + tmp_56*(3.0*tmp_54*tmp_55*tmp_7 - tmp_54*tmp_89 - tmp_55*tmp_97) + tmp_64*(3.0*tmp_62*tmp_63*tmp_7 - tmp_62*tmp_89 - tmp_63*tmp_97) + tmp_72*(3.0*tmp_7*tmp_70*tmp_71 - tmp_70*tmp_89 - tmp_71*tmp_97) + tmp_80*(3.0*tmp_7*tmp_78*tmp_79 - tmp_78*tmp_89 - tmp_79*tmp_97) + tmp_88*(3.0*tmp_7*tmp_86*tmp_87 - tmp_86*tmp_89 - tmp_87*tmp_97);
      real_t tmp_99 = 1.0*tmp_90;
      real_t tmp_100 = tmp_48*(3.0*tmp_35*tmp_39*tmp_7 - tmp_35*tmp_94 - tmp_39*tmp_91) + tmp_56*(3.0*tmp_52*tmp_53*tmp_7 - tmp_52*tmp_94 - tmp_53*tmp_91) + tmp_64*(3.0*tmp_60*tmp_61*tmp_7 - tmp_60*tmp_94 - tmp_61*tmp_91) + tmp_72*(3.0*tmp_68*tmp_69*tmp_7 - tmp_68*tmp_94 - tmp_69*tmp_91) + tmp_80*(3.0*tmp_7*tmp_76*tmp_77 - tmp_76*tmp_94 - tmp_77*tmp_91) + tmp_88*(3.0*tmp_7*tmp_84*tmp_85 - tmp_84*tmp_94 - tmp_85*tmp_91);
      real_t tmp_101 = tmp_48*(3.0*tmp_35*tmp_43*tmp_7 - tmp_35*tmp_97 - tmp_43*tmp_91) + tmp_56*(3.0*tmp_52*tmp_54*tmp_7 - tmp_52*tmp_97 - tmp_54*tmp_91) + tmp_64*(3.0*tmp_60*tmp_62*tmp_7 - tmp_60*tmp_97 - tmp_62*tmp_91) + tmp_72*(3.0*tmp_68*tmp_7*tmp_70 - tmp_68*tmp_97 - tmp_70*tmp_91) + tmp_80*(3.0*tmp_7*tmp_76*tmp_78 - tmp_76*tmp_97 - tmp_78*tmp_91) + tmp_88*(3.0*tmp_7*tmp_84*tmp_86 - tmp_84*tmp_97 - tmp_86*tmp_91);
      real_t tmp_102 = 1.0*tmp_93;
      real_t tmp_103 = tmp_48*(3.0*tmp_39*tmp_43*tmp_7 - tmp_39*tmp_97 - tmp_43*tmp_94) + tmp_56*(3.0*tmp_53*tmp_54*tmp_7 - tmp_53*tmp_97 - tmp_54*tmp_94) + tmp_64*(3.0*tmp_61*tmp_62*tmp_7 - tmp_61*tmp_97 - tmp_62*tmp_94) + tmp_72*(3.0*tmp_69*tmp_7*tmp_70 - tmp_69*tmp_97 - tmp_70*tmp_94) + tmp_80*(3.0*tmp_7*tmp_77*tmp_78 - tmp_77*tmp_97 - tmp_78*tmp_94) + tmp_88*(3.0*tmp_7*tmp_85*tmp_86 - tmp_85*tmp_97 - tmp_86*tmp_94);
      real_t tmp_104 = 1.0*tmp_96;
      real_t a_0_0 = tmp_48*(3.0*(tmp_44*tmp_44)*tmp_7 - tmp_44*tmp_46) + tmp_56*(-tmp_46*tmp_55 + 3.0*(tmp_55*tmp_55)*tmp_7) + tmp_64*(-tmp_46*tmp_63 + 3.0*(tmp_63*tmp_63)*tmp_7) + tmp_72*(-tmp_46*tmp_71 + 3.0*tmp_7*(tmp_71*tmp_71)) + tmp_80*(-tmp_46*tmp_79 + 3.0*tmp_7*(tmp_79*tmp_79)) + tmp_88*(-tmp_46*tmp_87 + 3.0*tmp_7*(tmp_87*tmp_87));
      real_t a_0_1 = tmp_92;
      real_t a_0_2 = tmp_95;
      real_t a_0_3 = tmp_98;
      real_t a_1_0 = tmp_92;
      real_t a_1_1 = tmp_48*(3.0*(tmp_35*tmp_35)*tmp_7 - tmp_35*tmp_99) + tmp_56*(3.0*(tmp_52*tmp_52)*tmp_7 - tmp_52*tmp_99) + tmp_64*(3.0*(tmp_60*tmp_60)*tmp_7 - tmp_60*tmp_99) + tmp_72*(3.0*(tmp_68*tmp_68)*tmp_7 - tmp_68*tmp_99) + tmp_80*(3.0*tmp_7*(tmp_76*tmp_76) - tmp_76*tmp_99) + tmp_88*(3.0*tmp_7*(tmp_84*tmp_84) - tmp_84*tmp_99);
      real_t a_1_2 = tmp_100;
      real_t a_1_3 = tmp_101;
      real_t a_2_0 = tmp_95;
      real_t a_2_1 = tmp_100;
      real_t a_2_2 = tmp_48*(-tmp_102*tmp_39 + 3.0*(tmp_39*tmp_39)*tmp_7) + tmp_56*(-tmp_102*tmp_53 + 3.0*(tmp_53*tmp_53)*tmp_7) + tmp_64*(-tmp_102*tmp_61 + 3.0*(tmp_61*tmp_61)*tmp_7) + tmp_72*(-tmp_102*tmp_69 + 3.0*(tmp_69*tmp_69)*tmp_7) + tmp_80*(-tmp_102*tmp_77 + 3.0*tmp_7*(tmp_77*tmp_77)) + tmp_88*(-tmp_102*tmp_85 + 3.0*tmp_7*(tmp_85*tmp_85));
      real_t a_2_3 = tmp_103;
      real_t a_3_0 = tmp_98;
      real_t a_3_1 = tmp_101;
      real_t a_3_2 = tmp_103;
      real_t a_3_3 = tmp_48*(-tmp_104*tmp_43 + 3.0*(tmp_43*tmp_43)*tmp_7) + tmp_56*(-tmp_104*tmp_54 + 3.0*(tmp_54*tmp_54)*tmp_7) + tmp_64*(-tmp_104*tmp_62 + 3.0*(tmp_62*tmp_62)*tmp_7) + tmp_72*(-tmp_104*tmp_70 + 3.0*tmp_7*(tmp_70*tmp_70)) + tmp_80*(-tmp_104*tmp_78 + 3.0*tmp_7*(tmp_78*tmp_78)) + tmp_88*(-tmp_104*tmp_86 + 3.0*tmp_7*(tmp_86*tmp_86));
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
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
                                  MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_26 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_27 = -tmp_26;
      real_t tmp_28 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_29 = 0.091576213509770743*tmp_27 + 0.81684757298045851*tmp_28;
      real_t tmp_30 = tmp_25 + tmp_29;
      real_t tmp_31 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_32 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_33 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_34 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_35 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_36 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_39 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_40 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_41 = tmp_39*tmp_40;
      real_t tmp_42 = tmp_32*tmp_36;
      real_t tmp_43 = tmp_31*tmp_39;
      real_t tmp_44 = tmp_33*tmp_38;
      real_t tmp_45 = 1.0 / (tmp_31*tmp_32*tmp_38 + tmp_33*tmp_41 - tmp_34*tmp_44 + tmp_35*tmp_37 - tmp_35*tmp_43 - tmp_40*tmp_42);
      real_t tmp_46 = tmp_45*(tmp_31*tmp_32 - tmp_33*tmp_34);
      real_t tmp_47 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_48 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_49 = -tmp_48;
      real_t tmp_50 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_51 = 0.091576213509770743*tmp_49 + 0.81684757298045851*tmp_50;
      real_t tmp_52 = tmp_47 + tmp_51;
      real_t tmp_53 = tmp_45*(-tmp_31*tmp_35 + tmp_33*tmp_40);
      real_t tmp_54 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_55 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_56 = -tmp_55;
      real_t tmp_57 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_58 = 0.091576213509770743*tmp_56 + 0.81684757298045851*tmp_57;
      real_t tmp_59 = tmp_54 + tmp_58;
      real_t tmp_60 = tmp_45*(-tmp_32*tmp_40 + tmp_34*tmp_35);
      real_t tmp_61 = tmp_30*tmp_46 + tmp_52*tmp_53 + tmp_59*tmp_60;
      real_t tmp_62 = tmp_45*(tmp_33*tmp_39 - tmp_42);
      real_t tmp_63 = tmp_45*(tmp_35*tmp_36 - tmp_44);
      real_t tmp_64 = tmp_45*(tmp_32*tmp_38 - tmp_35*tmp_39);
      real_t tmp_65 = tmp_30*tmp_62 + tmp_52*tmp_63 + tmp_59*tmp_64;
      real_t tmp_66 = tmp_45*(tmp_37 - tmp_43);
      real_t tmp_67 = tmp_45*(tmp_31*tmp_38 - tmp_36*tmp_40);
      real_t tmp_68 = tmp_45*(-tmp_34*tmp_38 + tmp_41);
      real_t tmp_69 = tmp_30*tmp_66 + tmp_52*tmp_67 + tmp_59*tmp_68;
      real_t tmp_70 = -tmp_61 - tmp_65 - tmp_69 + 1;
      real_t tmp_71 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_72 = tmp_29 + tmp_71;
      real_t tmp_73 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_74 = tmp_51 + tmp_73;
      real_t tmp_75 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_76 = tmp_58 + tmp_75;
      real_t tmp_77 = tmp_17*tmp_76 + tmp_20*tmp_74 + tmp_23*tmp_72;
      real_t tmp_78 = tmp_16*tmp_76 + tmp_19*tmp_74 + tmp_22*tmp_72;
      real_t tmp_79 = tmp_15*tmp_76 + tmp_18*tmp_74 + tmp_21*tmp_72;
      real_t tmp_80 = -tmp_77 - tmp_78 - tmp_79 + 1;
      real_t tmp_81 = 0.5*p_affine_13_0*(-tmp_60 - tmp_64 - tmp_68) + 0.5*p_affine_13_1*(-tmp_53 - tmp_63 - tmp_67) + 0.5*p_affine_13_2*(-tmp_46 - tmp_62 - tmp_66);
      real_t tmp_82 = (std::abs(tmp_26*tmp_50 - tmp_28*tmp_48)*std::abs(tmp_26*tmp_50 - tmp_28*tmp_48)) + (std::abs(tmp_26*tmp_57 - tmp_28*tmp_55)*std::abs(tmp_26*tmp_57 - tmp_28*tmp_55)) + (std::abs(tmp_48*tmp_57 - tmp_50*tmp_55)*std::abs(tmp_48*tmp_57 - tmp_50*tmp_55));
      real_t tmp_83 = 3.0*std::pow(tmp_82, -0.25);
      real_t tmp_84 = tmp_80*tmp_83;
      real_t tmp_85 = 1.0*std::pow(tmp_82, 1.0/2.0);
      real_t tmp_86 = 0.054975871827660928*tmp_85;
      real_t tmp_87 = 0.44594849091596489*tmp_27 + 0.10810301816807022*tmp_28;
      real_t tmp_88 = tmp_25 + tmp_87;
      real_t tmp_89 = 0.44594849091596489*tmp_49 + 0.10810301816807022*tmp_50;
      real_t tmp_90 = tmp_47 + tmp_89;
      real_t tmp_91 = 0.44594849091596489*tmp_56 + 0.10810301816807022*tmp_57;
      real_t tmp_92 = tmp_54 + tmp_91;
      real_t tmp_93 = tmp_46*tmp_88 + tmp_53*tmp_90 + tmp_60*tmp_92;
      real_t tmp_94 = tmp_62*tmp_88 + tmp_63*tmp_90 + tmp_64*tmp_92;
      real_t tmp_95 = tmp_66*tmp_88 + tmp_67*tmp_90 + tmp_68*tmp_92;
      real_t tmp_96 = -tmp_93 - tmp_94 - tmp_95 + 1;
      real_t tmp_97 = tmp_71 + tmp_87;
      real_t tmp_98 = tmp_73 + tmp_89;
      real_t tmp_99 = tmp_75 + tmp_91;
      real_t tmp_100 = tmp_17*tmp_99 + tmp_20*tmp_98 + tmp_23*tmp_97;
      real_t tmp_101 = tmp_16*tmp_99 + tmp_19*tmp_98 + tmp_22*tmp_97;
      real_t tmp_102 = tmp_15*tmp_99 + tmp_18*tmp_98 + tmp_21*tmp_97;
      real_t tmp_103 = -tmp_100 - tmp_101 - tmp_102 + 1;
      real_t tmp_104 = tmp_103*tmp_83;
      real_t tmp_105 = 0.11169079483900572*tmp_85;
      real_t tmp_106 = 0.81684757298045851*tmp_27 + 0.091576213509770743*tmp_28;
      real_t tmp_107 = tmp_106 + tmp_25;
      real_t tmp_108 = 0.81684757298045851*tmp_49 + 0.091576213509770743*tmp_50;
      real_t tmp_109 = tmp_108 + tmp_47;
      real_t tmp_110 = 0.81684757298045851*tmp_56 + 0.091576213509770743*tmp_57;
      real_t tmp_111 = tmp_110 + tmp_54;
      real_t tmp_112 = tmp_107*tmp_46 + tmp_109*tmp_53 + tmp_111*tmp_60;
      real_t tmp_113 = tmp_107*tmp_62 + tmp_109*tmp_63 + tmp_111*tmp_64;
      real_t tmp_114 = tmp_107*tmp_66 + tmp_109*tmp_67 + tmp_111*tmp_68;
      real_t tmp_115 = -tmp_112 - tmp_113 - tmp_114 + 1;
      real_t tmp_116 = tmp_106 + tmp_71;
      real_t tmp_117 = tmp_108 + tmp_73;
      real_t tmp_118 = tmp_110 + tmp_75;
      real_t tmp_119 = tmp_116*tmp_23 + tmp_117*tmp_20 + tmp_118*tmp_17;
      real_t tmp_120 = tmp_116*tmp_22 + tmp_117*tmp_19 + tmp_118*tmp_16;
      real_t tmp_121 = tmp_116*tmp_21 + tmp_117*tmp_18 + tmp_118*tmp_15;
      real_t tmp_122 = -tmp_119 - tmp_120 - tmp_121 + 1;
      real_t tmp_123 = tmp_122*tmp_83;
      real_t tmp_124 = 0.054975871827660928*tmp_85;
      real_t tmp_125 = 0.10810301816807022*tmp_27 + 0.44594849091596489*tmp_28;
      real_t tmp_126 = tmp_125 + tmp_25;
      real_t tmp_127 = 0.10810301816807022*tmp_49 + 0.44594849091596489*tmp_50;
      real_t tmp_128 = tmp_127 + tmp_47;
      real_t tmp_129 = 0.10810301816807022*tmp_56 + 0.44594849091596489*tmp_57;
      real_t tmp_130 = tmp_129 + tmp_54;
      real_t tmp_131 = tmp_126*tmp_46 + tmp_128*tmp_53 + tmp_130*tmp_60;
      real_t tmp_132 = tmp_126*tmp_62 + tmp_128*tmp_63 + tmp_130*tmp_64;
      real_t tmp_133 = tmp_126*tmp_66 + tmp_128*tmp_67 + tmp_130*tmp_68;
      real_t tmp_134 = -tmp_131 - tmp_132 - tmp_133 + 1;
      real_t tmp_135 = tmp_125 + tmp_71;
      real_t tmp_136 = tmp_127 + tmp_73;
      real_t tmp_137 = tmp_129 + tmp_75;
      real_t tmp_138 = tmp_135*tmp_23 + tmp_136*tmp_20 + tmp_137*tmp_17;
      real_t tmp_139 = tmp_135*tmp_22 + tmp_136*tmp_19 + tmp_137*tmp_16;
      real_t tmp_140 = tmp_135*tmp_21 + tmp_136*tmp_18 + tmp_137*tmp_15;
      real_t tmp_141 = -tmp_138 - tmp_139 - tmp_140 + 1;
      real_t tmp_142 = tmp_141*tmp_83;
      real_t tmp_143 = 0.11169079483900572*tmp_85;
      real_t tmp_144 = 0.091576213509770743*tmp_27 + 0.091576213509770743*tmp_28;
      real_t tmp_145 = tmp_144 + tmp_25;
      real_t tmp_146 = 0.091576213509770743*tmp_49 + 0.091576213509770743*tmp_50;
      real_t tmp_147 = tmp_146 + tmp_47;
      real_t tmp_148 = 0.091576213509770743*tmp_56 + 0.091576213509770743*tmp_57;
      real_t tmp_149 = tmp_148 + tmp_54;
      real_t tmp_150 = tmp_145*tmp_46 + tmp_147*tmp_53 + tmp_149*tmp_60;
      real_t tmp_151 = tmp_145*tmp_62 + tmp_147*tmp_63 + tmp_149*tmp_64;
      real_t tmp_152 = tmp_145*tmp_66 + tmp_147*tmp_67 + tmp_149*tmp_68;
      real_t tmp_153 = -tmp_150 - tmp_151 - tmp_152 + 1;
      real_t tmp_154 = tmp_144 + tmp_71;
      real_t tmp_155 = tmp_146 + tmp_73;
      real_t tmp_156 = tmp_148 + tmp_75;
      real_t tmp_157 = tmp_154*tmp_23 + tmp_155*tmp_20 + tmp_156*tmp_17;
      real_t tmp_158 = tmp_154*tmp_22 + tmp_155*tmp_19 + tmp_156*tmp_16;
      real_t tmp_159 = tmp_15*tmp_156 + tmp_154*tmp_21 + tmp_155*tmp_18;
      real_t tmp_160 = -tmp_157 - tmp_158 - tmp_159 + 1;
      real_t tmp_161 = tmp_160*tmp_83;
      real_t tmp_162 = 0.054975871827660928*tmp_85;
      real_t tmp_163 = 0.44594849091596489*tmp_27 + 0.44594849091596489*tmp_28;
      real_t tmp_164 = tmp_163 + tmp_25;
      real_t tmp_165 = 0.44594849091596489*tmp_49 + 0.44594849091596489*tmp_50;
      real_t tmp_166 = tmp_165 + tmp_47;
      real_t tmp_167 = 0.44594849091596489*tmp_56 + 0.44594849091596489*tmp_57;
      real_t tmp_168 = tmp_167 + tmp_54;
      real_t tmp_169 = tmp_164*tmp_46 + tmp_166*tmp_53 + tmp_168*tmp_60;
      real_t tmp_170 = tmp_164*tmp_62 + tmp_166*tmp_63 + tmp_168*tmp_64;
      real_t tmp_171 = tmp_164*tmp_66 + tmp_166*tmp_67 + tmp_168*tmp_68;
      real_t tmp_172 = -tmp_169 - tmp_170 - tmp_171 + 1;
      real_t tmp_173 = tmp_163 + tmp_71;
      real_t tmp_174 = tmp_165 + tmp_73;
      real_t tmp_175 = tmp_167 + tmp_75;
      real_t tmp_176 = tmp_17*tmp_175 + tmp_173*tmp_23 + tmp_174*tmp_20;
      real_t tmp_177 = tmp_16*tmp_175 + tmp_173*tmp_22 + tmp_174*tmp_19;
      real_t tmp_178 = tmp_15*tmp_175 + tmp_173*tmp_21 + tmp_174*tmp_18;
      real_t tmp_179 = -tmp_176 - tmp_177 - tmp_178 + 1;
      real_t tmp_180 = tmp_179*tmp_83;
      real_t tmp_181 = 0.11169079483900572*tmp_85;
      real_t tmp_182 = 0.5*p_affine_13_0*tmp_60 + 0.5*p_affine_13_1*tmp_53 + 0.5*p_affine_13_2*tmp_46;
      real_t tmp_183 = 0.5*p_affine_13_0*tmp_64 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_62;
      real_t tmp_184 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_67 + 0.5*p_affine_13_2*tmp_66;
      real_t tmp_185 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_186 = tmp_77*tmp_83;
      real_t tmp_187 = tmp_100*tmp_83;
      real_t tmp_188 = tmp_119*tmp_83;
      real_t tmp_189 = tmp_138*tmp_83;
      real_t tmp_190 = tmp_157*tmp_83;
      real_t tmp_191 = tmp_176*tmp_83;
      real_t tmp_192 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_193 = tmp_78*tmp_83;
      real_t tmp_194 = tmp_101*tmp_83;
      real_t tmp_195 = tmp_120*tmp_83;
      real_t tmp_196 = tmp_139*tmp_83;
      real_t tmp_197 = tmp_158*tmp_83;
      real_t tmp_198 = tmp_177*tmp_83;
      real_t tmp_199 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t tmp_200 = tmp_79*tmp_83;
      real_t tmp_201 = tmp_102*tmp_83;
      real_t tmp_202 = tmp_121*tmp_83;
      real_t tmp_203 = tmp_140*tmp_83;
      real_t tmp_204 = tmp_159*tmp_83;
      real_t tmp_205 = tmp_178*tmp_83;
      real_t a_0_0 = tmp_105*(-tmp_103*tmp_81 - tmp_104*tmp_96 + 0.5*tmp_24*tmp_96) + tmp_124*(-tmp_115*tmp_123 + 0.5*tmp_115*tmp_24 - tmp_122*tmp_81) + tmp_143*(-tmp_134*tmp_142 + 0.5*tmp_134*tmp_24 - tmp_141*tmp_81) + tmp_162*(-tmp_153*tmp_161 + 0.5*tmp_153*tmp_24 - tmp_160*tmp_81) + tmp_181*(-tmp_172*tmp_180 + 0.5*tmp_172*tmp_24 - tmp_179*tmp_81) + tmp_86*(0.5*tmp_24*tmp_70 - tmp_70*tmp_84 - tmp_80*tmp_81);
      real_t a_0_1 = tmp_105*(-tmp_103*tmp_182 - tmp_104*tmp_93 + 0.5*tmp_24*tmp_93) + tmp_124*(-tmp_112*tmp_123 + 0.5*tmp_112*tmp_24 - tmp_122*tmp_182) + tmp_143*(-tmp_131*tmp_142 + 0.5*tmp_131*tmp_24 - tmp_141*tmp_182) + tmp_162*(-tmp_150*tmp_161 + 0.5*tmp_150*tmp_24 - tmp_160*tmp_182) + tmp_181*(-tmp_169*tmp_180 + 0.5*tmp_169*tmp_24 - tmp_179*tmp_182) + tmp_86*(-tmp_182*tmp_80 + 0.5*tmp_24*tmp_61 - tmp_61*tmp_84);
      real_t a_0_2 = tmp_105*(-tmp_103*tmp_183 - tmp_104*tmp_94 + 0.5*tmp_24*tmp_94) + tmp_124*(-tmp_113*tmp_123 + 0.5*tmp_113*tmp_24 - tmp_122*tmp_183) + tmp_143*(-tmp_132*tmp_142 + 0.5*tmp_132*tmp_24 - tmp_141*tmp_183) + tmp_162*(-tmp_151*tmp_161 + 0.5*tmp_151*tmp_24 - tmp_160*tmp_183) + tmp_181*(-tmp_170*tmp_180 + 0.5*tmp_170*tmp_24 - tmp_179*tmp_183) + tmp_86*(-tmp_183*tmp_80 + 0.5*tmp_24*tmp_65 - tmp_65*tmp_84);
      real_t a_0_3 = tmp_105*(-tmp_103*tmp_184 - tmp_104*tmp_95 + 0.5*tmp_24*tmp_95) + tmp_124*(-tmp_114*tmp_123 + 0.5*tmp_114*tmp_24 - tmp_122*tmp_184) + tmp_143*(-tmp_133*tmp_142 + 0.5*tmp_133*tmp_24 - tmp_141*tmp_184) + tmp_162*(-tmp_152*tmp_161 + 0.5*tmp_152*tmp_24 - tmp_160*tmp_184) + tmp_181*(-tmp_171*tmp_180 + 0.5*tmp_171*tmp_24 - tmp_179*tmp_184) + tmp_86*(-tmp_184*tmp_80 + 0.5*tmp_24*tmp_69 - tmp_69*tmp_84);
      real_t a_1_0 = tmp_105*(-tmp_100*tmp_81 + 0.5*tmp_185*tmp_96 - tmp_187*tmp_96) + tmp_124*(0.5*tmp_115*tmp_185 - tmp_115*tmp_188 - tmp_119*tmp_81) + tmp_143*(0.5*tmp_134*tmp_185 - tmp_134*tmp_189 - tmp_138*tmp_81) + tmp_162*(0.5*tmp_153*tmp_185 - tmp_153*tmp_190 - tmp_157*tmp_81) + tmp_181*(0.5*tmp_172*tmp_185 - tmp_172*tmp_191 - tmp_176*tmp_81) + tmp_86*(0.5*tmp_185*tmp_70 - tmp_186*tmp_70 - tmp_77*tmp_81);
      real_t a_1_1 = tmp_105*(-tmp_100*tmp_182 + 0.5*tmp_185*tmp_93 - tmp_187*tmp_93) + tmp_124*(0.5*tmp_112*tmp_185 - tmp_112*tmp_188 - tmp_119*tmp_182) + tmp_143*(0.5*tmp_131*tmp_185 - tmp_131*tmp_189 - tmp_138*tmp_182) + tmp_162*(0.5*tmp_150*tmp_185 - tmp_150*tmp_190 - tmp_157*tmp_182) + tmp_181*(0.5*tmp_169*tmp_185 - tmp_169*tmp_191 - tmp_176*tmp_182) + tmp_86*(-tmp_182*tmp_77 + 0.5*tmp_185*tmp_61 - tmp_186*tmp_61);
      real_t a_1_2 = tmp_105*(-tmp_100*tmp_183 + 0.5*tmp_185*tmp_94 - tmp_187*tmp_94) + tmp_124*(0.5*tmp_113*tmp_185 - tmp_113*tmp_188 - tmp_119*tmp_183) + tmp_143*(0.5*tmp_132*tmp_185 - tmp_132*tmp_189 - tmp_138*tmp_183) + tmp_162*(0.5*tmp_151*tmp_185 - tmp_151*tmp_190 - tmp_157*tmp_183) + tmp_181*(0.5*tmp_170*tmp_185 - tmp_170*tmp_191 - tmp_176*tmp_183) + tmp_86*(-tmp_183*tmp_77 + 0.5*tmp_185*tmp_65 - tmp_186*tmp_65);
      real_t a_1_3 = tmp_105*(-tmp_100*tmp_184 + 0.5*tmp_185*tmp_95 - tmp_187*tmp_95) + tmp_124*(0.5*tmp_114*tmp_185 - tmp_114*tmp_188 - tmp_119*tmp_184) + tmp_143*(0.5*tmp_133*tmp_185 - tmp_133*tmp_189 - tmp_138*tmp_184) + tmp_162*(0.5*tmp_152*tmp_185 - tmp_152*tmp_190 - tmp_157*tmp_184) + tmp_181*(0.5*tmp_171*tmp_185 - tmp_171*tmp_191 - tmp_176*tmp_184) + tmp_86*(-tmp_184*tmp_77 + 0.5*tmp_185*tmp_69 - tmp_186*tmp_69);
      real_t a_2_0 = tmp_105*(-tmp_101*tmp_81 + 0.5*tmp_192*tmp_96 - tmp_194*tmp_96) + tmp_124*(0.5*tmp_115*tmp_192 - tmp_115*tmp_195 - tmp_120*tmp_81) + tmp_143*(0.5*tmp_134*tmp_192 - tmp_134*tmp_196 - tmp_139*tmp_81) + tmp_162*(0.5*tmp_153*tmp_192 - tmp_153*tmp_197 - tmp_158*tmp_81) + tmp_181*(0.5*tmp_172*tmp_192 - tmp_172*tmp_198 - tmp_177*tmp_81) + tmp_86*(0.5*tmp_192*tmp_70 - tmp_193*tmp_70 - tmp_78*tmp_81);
      real_t a_2_1 = tmp_105*(-tmp_101*tmp_182 + 0.5*tmp_192*tmp_93 - tmp_194*tmp_93) + tmp_124*(0.5*tmp_112*tmp_192 - tmp_112*tmp_195 - tmp_120*tmp_182) + tmp_143*(0.5*tmp_131*tmp_192 - tmp_131*tmp_196 - tmp_139*tmp_182) + tmp_162*(0.5*tmp_150*tmp_192 - tmp_150*tmp_197 - tmp_158*tmp_182) + tmp_181*(0.5*tmp_169*tmp_192 - tmp_169*tmp_198 - tmp_177*tmp_182) + tmp_86*(-tmp_182*tmp_78 + 0.5*tmp_192*tmp_61 - tmp_193*tmp_61);
      real_t a_2_2 = tmp_105*(-tmp_101*tmp_183 + 0.5*tmp_192*tmp_94 - tmp_194*tmp_94) + tmp_124*(0.5*tmp_113*tmp_192 - tmp_113*tmp_195 - tmp_120*tmp_183) + tmp_143*(0.5*tmp_132*tmp_192 - tmp_132*tmp_196 - tmp_139*tmp_183) + tmp_162*(0.5*tmp_151*tmp_192 - tmp_151*tmp_197 - tmp_158*tmp_183) + tmp_181*(0.5*tmp_170*tmp_192 - tmp_170*tmp_198 - tmp_177*tmp_183) + tmp_86*(-tmp_183*tmp_78 + 0.5*tmp_192*tmp_65 - tmp_193*tmp_65);
      real_t a_2_3 = tmp_105*(-tmp_101*tmp_184 + 0.5*tmp_192*tmp_95 - tmp_194*tmp_95) + tmp_124*(0.5*tmp_114*tmp_192 - tmp_114*tmp_195 - tmp_120*tmp_184) + tmp_143*(0.5*tmp_133*tmp_192 - tmp_133*tmp_196 - tmp_139*tmp_184) + tmp_162*(0.5*tmp_152*tmp_192 - tmp_152*tmp_197 - tmp_158*tmp_184) + tmp_181*(0.5*tmp_171*tmp_192 - tmp_171*tmp_198 - tmp_177*tmp_184) + tmp_86*(-tmp_184*tmp_78 + 0.5*tmp_192*tmp_69 - tmp_193*tmp_69);
      real_t a_3_0 = tmp_105*(-tmp_102*tmp_81 + 0.5*tmp_199*tmp_96 - tmp_201*tmp_96) + tmp_124*(0.5*tmp_115*tmp_199 - tmp_115*tmp_202 - tmp_121*tmp_81) + tmp_143*(0.5*tmp_134*tmp_199 - tmp_134*tmp_203 - tmp_140*tmp_81) + tmp_162*(0.5*tmp_153*tmp_199 - tmp_153*tmp_204 - tmp_159*tmp_81) + tmp_181*(0.5*tmp_172*tmp_199 - tmp_172*tmp_205 - tmp_178*tmp_81) + tmp_86*(0.5*tmp_199*tmp_70 - tmp_200*tmp_70 - tmp_79*tmp_81);
      real_t a_3_1 = tmp_105*(-tmp_102*tmp_182 + 0.5*tmp_199*tmp_93 - tmp_201*tmp_93) + tmp_124*(0.5*tmp_112*tmp_199 - tmp_112*tmp_202 - tmp_121*tmp_182) + tmp_143*(0.5*tmp_131*tmp_199 - tmp_131*tmp_203 - tmp_140*tmp_182) + tmp_162*(0.5*tmp_150*tmp_199 - tmp_150*tmp_204 - tmp_159*tmp_182) + tmp_181*(0.5*tmp_169*tmp_199 - tmp_169*tmp_205 - tmp_178*tmp_182) + tmp_86*(-tmp_182*tmp_79 + 0.5*tmp_199*tmp_61 - tmp_200*tmp_61);
      real_t a_3_2 = tmp_105*(-tmp_102*tmp_183 + 0.5*tmp_199*tmp_94 - tmp_201*tmp_94) + tmp_124*(0.5*tmp_113*tmp_199 - tmp_113*tmp_202 - tmp_121*tmp_183) + tmp_143*(0.5*tmp_132*tmp_199 - tmp_132*tmp_203 - tmp_140*tmp_183) + tmp_162*(0.5*tmp_151*tmp_199 - tmp_151*tmp_204 - tmp_159*tmp_183) + tmp_181*(0.5*tmp_170*tmp_199 - tmp_170*tmp_205 - tmp_178*tmp_183) + tmp_86*(-tmp_183*tmp_79 + 0.5*tmp_199*tmp_65 - tmp_200*tmp_65);
      real_t a_3_3 = tmp_105*(-tmp_102*tmp_184 + 0.5*tmp_199*tmp_95 - tmp_201*tmp_95) + tmp_124*(0.5*tmp_114*tmp_199 - tmp_114*tmp_202 - tmp_121*tmp_184) + tmp_143*(0.5*tmp_133*tmp_199 - tmp_133*tmp_203 - tmp_140*tmp_184) + tmp_162*(0.5*tmp_152*tmp_199 - tmp_152*tmp_204 - tmp_159*tmp_184) + tmp_181*(0.5*tmp_171*tmp_199 - tmp_171*tmp_205 - tmp_178*tmp_184) + tmp_86*(-tmp_184*tmp_79 + 0.5*tmp_199*tmp_69 - tmp_200*tmp_69);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
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
                                        MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -tmp_4;
      real_t tmp_9 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_10 = 0.81684757298045851*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_11 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_13 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_14 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_15 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_16 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_19 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_20 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_12*tmp_20;
      real_t tmp_23 = tmp_16*tmp_19;
      real_t tmp_24 = tmp_14*tmp_18;
      real_t tmp_25 = 1.0 / (tmp_11*tmp_12*tmp_18 - tmp_11*tmp_23 + tmp_13*tmp_21 - tmp_13*tmp_24 + tmp_15*tmp_17 - tmp_15*tmp_22);
      real_t tmp_26 = tmp_25*(tmp_11*tmp_12 - tmp_13*tmp_14);
      real_t tmp_27 = -tmp_1;
      real_t tmp_28 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_29 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_30 = tmp_25*(-tmp_11*tmp_16 + tmp_13*tmp_20);
      real_t tmp_31 = -tmp_3;
      real_t tmp_32 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_33 = 0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_34 = tmp_25*(tmp_17 - tmp_22);
      real_t tmp_35 = tmp_10*tmp_26 + tmp_29*tmp_30 + tmp_33*tmp_34;
      real_t tmp_36 = tmp_25*(-tmp_12*tmp_15 + tmp_13*tmp_19);
      real_t tmp_37 = tmp_25*(-tmp_13*tmp_18 + tmp_15*tmp_16);
      real_t tmp_38 = tmp_25*(tmp_12*tmp_18 - tmp_23);
      real_t tmp_39 = tmp_10*tmp_36 + tmp_29*tmp_37 + tmp_33*tmp_38;
      real_t tmp_40 = tmp_25*(-tmp_11*tmp_19 + tmp_14*tmp_15);
      real_t tmp_41 = tmp_25*(tmp_11*tmp_18 - tmp_15*tmp_20);
      real_t tmp_42 = tmp_25*(tmp_21 - tmp_24);
      real_t tmp_43 = tmp_10*tmp_40 + tmp_29*tmp_41 + tmp_33*tmp_42;
      real_t tmp_44 = -tmp_35 - tmp_39 - tmp_43 + 1;
      real_t tmp_45 = p_affine_13_0*(-tmp_34 - tmp_38 - tmp_42) + p_affine_13_1*(-tmp_30 - tmp_37 - tmp_41) + p_affine_13_2*(-tmp_26 - tmp_36 - tmp_40);
      real_t tmp_46 = 2*tmp_45;
      real_t tmp_47 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_48 = 0.054975871827660928*tmp_47;
      real_t tmp_49 = 0.10810301816807022*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_50 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_51 = 0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_52 = tmp_26*tmp_49 + tmp_30*tmp_50 + tmp_34*tmp_51;
      real_t tmp_53 = tmp_36*tmp_49 + tmp_37*tmp_50 + tmp_38*tmp_51;
      real_t tmp_54 = tmp_40*tmp_49 + tmp_41*tmp_50 + tmp_42*tmp_51;
      real_t tmp_55 = -tmp_52 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11169079483900572*tmp_47;
      real_t tmp_57 = 0.091576213509770743*tmp_5 + 0.81684757298045851*tmp_8 + tmp_9;
      real_t tmp_58 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_59 = 0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_31 + tmp_32;
      real_t tmp_60 = tmp_26*tmp_57 + tmp_30*tmp_58 + tmp_34*tmp_59;
      real_t tmp_61 = tmp_36*tmp_57 + tmp_37*tmp_58 + tmp_38*tmp_59;
      real_t tmp_62 = tmp_40*tmp_57 + tmp_41*tmp_58 + tmp_42*tmp_59;
      real_t tmp_63 = -tmp_60 - tmp_61 - tmp_62 + 1;
      real_t tmp_64 = 0.054975871827660928*tmp_47;
      real_t tmp_65 = 0.44594849091596489*tmp_5 + 0.10810301816807022*tmp_8 + tmp_9;
      real_t tmp_66 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_67 = 0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_31 + tmp_32;
      real_t tmp_68 = tmp_26*tmp_65 + tmp_30*tmp_66 + tmp_34*tmp_67;
      real_t tmp_69 = tmp_36*tmp_65 + tmp_37*tmp_66 + tmp_38*tmp_67;
      real_t tmp_70 = tmp_40*tmp_65 + tmp_41*tmp_66 + tmp_42*tmp_67;
      real_t tmp_71 = -tmp_68 - tmp_69 - tmp_70 + 1;
      real_t tmp_72 = 0.11169079483900572*tmp_47;
      real_t tmp_73 = 0.091576213509770743*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_74 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_75 = 0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_76 = tmp_26*tmp_73 + tmp_30*tmp_74 + tmp_34*tmp_75;
      real_t tmp_77 = tmp_36*tmp_73 + tmp_37*tmp_74 + tmp_38*tmp_75;
      real_t tmp_78 = tmp_40*tmp_73 + tmp_41*tmp_74 + tmp_42*tmp_75;
      real_t tmp_79 = -tmp_76 - tmp_77 - tmp_78 + 1;
      real_t tmp_80 = 0.054975871827660928*tmp_47;
      real_t tmp_81 = 0.44594849091596489*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_82 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_83 = 0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_84 = tmp_26*tmp_81 + tmp_30*tmp_82 + tmp_34*tmp_83;
      real_t tmp_85 = tmp_36*tmp_81 + tmp_37*tmp_82 + tmp_38*tmp_83;
      real_t tmp_86 = tmp_40*tmp_81 + tmp_41*tmp_82 + tmp_42*tmp_83;
      real_t tmp_87 = -tmp_84 - tmp_85 - tmp_86 + 1;
      real_t tmp_88 = 0.11169079483900572*tmp_47;
      real_t tmp_89 = p_affine_13_0*tmp_34 + p_affine_13_1*tmp_30 + p_affine_13_2*tmp_26;
      real_t tmp_90 = tmp_48*(3.0*tmp_35*tmp_44*tmp_7 - tmp_35*tmp_45 - tmp_44*tmp_89) + tmp_56*(-tmp_45*tmp_52 + 3.0*tmp_52*tmp_55*tmp_7 - tmp_55*tmp_89) + tmp_64*(-tmp_45*tmp_60 + 3.0*tmp_60*tmp_63*tmp_7 - tmp_63*tmp_89) + tmp_72*(-tmp_45*tmp_68 + 3.0*tmp_68*tmp_7*tmp_71 - tmp_71*tmp_89) + tmp_80*(-tmp_45*tmp_76 + 3.0*tmp_7*tmp_76*tmp_79 - tmp_79*tmp_89) + tmp_88*(-tmp_45*tmp_84 + 3.0*tmp_7*tmp_84*tmp_87 - tmp_87*tmp_89);
      real_t tmp_91 = p_affine_13_0*tmp_38 + p_affine_13_1*tmp_37 + p_affine_13_2*tmp_36;
      real_t tmp_92 = tmp_48*(3.0*tmp_39*tmp_44*tmp_7 - tmp_39*tmp_45 - tmp_44*tmp_91) + tmp_56*(-tmp_45*tmp_53 + 3.0*tmp_53*tmp_55*tmp_7 - tmp_55*tmp_91) + tmp_64*(-tmp_45*tmp_61 + 3.0*tmp_61*tmp_63*tmp_7 - tmp_63*tmp_91) + tmp_72*(-tmp_45*tmp_69 + 3.0*tmp_69*tmp_7*tmp_71 - tmp_71*tmp_91) + tmp_80*(-tmp_45*tmp_77 + 3.0*tmp_7*tmp_77*tmp_79 - tmp_79*tmp_91) + tmp_88*(-tmp_45*tmp_85 + 3.0*tmp_7*tmp_85*tmp_87 - tmp_87*tmp_91);
      real_t tmp_93 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_40;
      real_t tmp_94 = tmp_48*(3.0*tmp_43*tmp_44*tmp_7 - tmp_43*tmp_45 - tmp_44*tmp_93) + tmp_56*(-tmp_45*tmp_54 + 3.0*tmp_54*tmp_55*tmp_7 - tmp_55*tmp_93) + tmp_64*(-tmp_45*tmp_62 + 3.0*tmp_62*tmp_63*tmp_7 - tmp_63*tmp_93) + tmp_72*(-tmp_45*tmp_70 + 3.0*tmp_7*tmp_70*tmp_71 - tmp_71*tmp_93) + tmp_80*(-tmp_45*tmp_78 + 3.0*tmp_7*tmp_78*tmp_79 - tmp_79*tmp_93) + tmp_88*(-tmp_45*tmp_86 + 3.0*tmp_7*tmp_86*tmp_87 - tmp_87*tmp_93);
      real_t tmp_95 = 2*tmp_89;
      real_t tmp_96 = tmp_48*(3.0*tmp_35*tmp_39*tmp_7 - tmp_35*tmp_91 - tmp_39*tmp_89) + tmp_56*(3.0*tmp_52*tmp_53*tmp_7 - tmp_52*tmp_91 - tmp_53*tmp_89) + tmp_64*(3.0*tmp_60*tmp_61*tmp_7 - tmp_60*tmp_91 - tmp_61*tmp_89) + tmp_72*(3.0*tmp_68*tmp_69*tmp_7 - tmp_68*tmp_91 - tmp_69*tmp_89) + tmp_80*(3.0*tmp_7*tmp_76*tmp_77 - tmp_76*tmp_91 - tmp_77*tmp_89) + tmp_88*(3.0*tmp_7*tmp_84*tmp_85 - tmp_84*tmp_91 - tmp_85*tmp_89);
      real_t tmp_97 = tmp_48*(3.0*tmp_35*tmp_43*tmp_7 - tmp_35*tmp_93 - tmp_43*tmp_89) + tmp_56*(3.0*tmp_52*tmp_54*tmp_7 - tmp_52*tmp_93 - tmp_54*tmp_89) + tmp_64*(3.0*tmp_60*tmp_62*tmp_7 - tmp_60*tmp_93 - tmp_62*tmp_89) + tmp_72*(3.0*tmp_68*tmp_7*tmp_70 - tmp_68*tmp_93 - tmp_70*tmp_89) + tmp_80*(3.0*tmp_7*tmp_76*tmp_78 - tmp_76*tmp_93 - tmp_78*tmp_89) + tmp_88*(3.0*tmp_7*tmp_84*tmp_86 - tmp_84*tmp_93 - tmp_86*tmp_89);
      real_t tmp_98 = 2*tmp_91;
      real_t tmp_99 = tmp_48*(3.0*tmp_39*tmp_43*tmp_7 - tmp_39*tmp_93 - tmp_43*tmp_91) + tmp_56*(3.0*tmp_53*tmp_54*tmp_7 - tmp_53*tmp_93 - tmp_54*tmp_91) + tmp_64*(3.0*tmp_61*tmp_62*tmp_7 - tmp_61*tmp_93 - tmp_62*tmp_91) + tmp_72*(3.0*tmp_69*tmp_7*tmp_70 - tmp_69*tmp_93 - tmp_70*tmp_91) + tmp_80*(3.0*tmp_7*tmp_77*tmp_78 - tmp_77*tmp_93 - tmp_78*tmp_91) + tmp_88*(3.0*tmp_7*tmp_85*tmp_86 - tmp_85*tmp_93 - tmp_86*tmp_91);
      real_t tmp_100 = 2*tmp_93;
      real_t a_0_0 = tmp_48*(3.0*(tmp_44*tmp_44)*tmp_7 - tmp_44*tmp_46) + tmp_56*(-tmp_46*tmp_55 + 3.0*(tmp_55*tmp_55)*tmp_7) + tmp_64*(-tmp_46*tmp_63 + 3.0*(tmp_63*tmp_63)*tmp_7) + tmp_72*(-tmp_46*tmp_71 + 3.0*tmp_7*(tmp_71*tmp_71)) + tmp_80*(-tmp_46*tmp_79 + 3.0*tmp_7*(tmp_79*tmp_79)) + tmp_88*(-tmp_46*tmp_87 + 3.0*tmp_7*(tmp_87*tmp_87));
      real_t a_0_1 = tmp_90;
      real_t a_0_2 = tmp_92;
      real_t a_0_3 = tmp_94;
      real_t a_1_0 = tmp_90;
      real_t a_1_1 = tmp_48*(3.0*(tmp_35*tmp_35)*tmp_7 - tmp_35*tmp_95) + tmp_56*(3.0*(tmp_52*tmp_52)*tmp_7 - tmp_52*tmp_95) + tmp_64*(3.0*(tmp_60*tmp_60)*tmp_7 - tmp_60*tmp_95) + tmp_72*(3.0*(tmp_68*tmp_68)*tmp_7 - tmp_68*tmp_95) + tmp_80*(3.0*tmp_7*(tmp_76*tmp_76) - tmp_76*tmp_95) + tmp_88*(3.0*tmp_7*(tmp_84*tmp_84) - tmp_84*tmp_95);
      real_t a_1_2 = tmp_96;
      real_t a_1_3 = tmp_97;
      real_t a_2_0 = tmp_92;
      real_t a_2_1 = tmp_96;
      real_t a_2_2 = tmp_48*(3.0*(tmp_39*tmp_39)*tmp_7 - tmp_39*tmp_98) + tmp_56*(3.0*(tmp_53*tmp_53)*tmp_7 - tmp_53*tmp_98) + tmp_64*(3.0*(tmp_61*tmp_61)*tmp_7 - tmp_61*tmp_98) + tmp_72*(3.0*(tmp_69*tmp_69)*tmp_7 - tmp_69*tmp_98) + tmp_80*(3.0*tmp_7*(tmp_77*tmp_77) - tmp_77*tmp_98) + tmp_88*(3.0*tmp_7*(tmp_85*tmp_85) - tmp_85*tmp_98);
      real_t a_2_3 = tmp_99;
      real_t a_3_0 = tmp_94;
      real_t a_3_1 = tmp_97;
      real_t a_3_2 = tmp_99;
      real_t a_3_3 = tmp_48*(-tmp_100*tmp_43 + 3.0*(tmp_43*tmp_43)*tmp_7) + tmp_56*(-tmp_100*tmp_54 + 3.0*(tmp_54*tmp_54)*tmp_7) + tmp_64*(-tmp_100*tmp_62 + 3.0*(tmp_62*tmp_62)*tmp_7) + tmp_72*(-tmp_100*tmp_70 + 3.0*tmp_7*(tmp_70*tmp_70)) + tmp_80*(-tmp_100*tmp_78 + 3.0*tmp_7*(tmp_78*tmp_78)) + tmp_88*(-tmp_100*tmp_86 + 3.0*tmp_7*(tmp_86*tmp_86));
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;

};




class EGVectorLaplaceFormNitscheBC_P1E_1 : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_P1E_1()

    {}





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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_5 = -tmp_4;
      real_t tmp_6 = 1.0 / (tmp_2 + tmp_3*tmp_5);
      real_t tmp_7 = tmp_0*tmp_6;
      real_t tmp_8 = tmp_3*tmp_6;
      real_t tmp_9 = tmp_2*tmp_6 + tmp_5*tmp_8;
      real_t tmp_10 = tmp_1*tmp_6;
      real_t tmp_11 = tmp_4*tmp_6;
      real_t tmp_12 = tmp_1*tmp_11 + tmp_10*tmp_5;
      real_t tmp_13 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = tmp_13*(tmp_12*(-tmp_10 - tmp_11) + tmp_9*(-tmp_7 - tmp_8));
      real_t tmp_15 = tmp_13*(tmp_10*tmp_12 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_13*(tmp_11*tmp_12 + tmp_7*tmp_9);
      real_t a_0_0 = 0.5*tmp_14;
      real_t a_1_0 = 0.5*tmp_15;
      real_t a_2_0 = 0.5*tmp_16;
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
                                       MatrixXr&                                           elMat ) const override
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
      real_t tmp_2 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_4 = 0.21132486540518713*tmp_2 + tmp_3;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = 1.0 / (tmp_1*tmp_5 + tmp_8);
      real_t tmp_10 = tmp_5*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_11 + tmp_12;
      real_t tmp_14 = tmp_7*tmp_9;
      real_t tmp_15 = tmp_10*tmp_4 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_6*tmp_9;
      real_t tmp_17 = tmp_0*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_4;
      real_t tmp_19 = tmp_1*(tmp_15 - 1.0/3.0) + tmp_7*(tmp_18 - 1.0/3.0);
      real_t tmp_20 = 0.5*p_affine_10_0*(-tmp_14 - tmp_17) + 0.5*p_affine_10_1*(-tmp_10 - tmp_16);
      real_t tmp_21 = -tmp_15 - tmp_18 + 1;
      real_t tmp_22 = 0.5*p_affine_10_0*(tmp_1*tmp_14 + tmp_17*tmp_7) + 0.5*p_affine_10_1*(tmp_1*tmp_10 + tmp_8*tmp_9);
      real_t tmp_23 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_2*tmp_2), 1.0/2.0));
      real_t tmp_24 = 1.0 / (tmp_23);
      real_t tmp_25 = 0.5*tmp_23;
      real_t tmp_26 = 0.78867513459481287*tmp_2 + tmp_3;
      real_t tmp_27 = 0.78867513459481287*tmp_11 + tmp_12;
      real_t tmp_28 = tmp_10*tmp_26 + tmp_14*tmp_27;
      real_t tmp_29 = tmp_16*tmp_26 + tmp_17*tmp_27;
      real_t tmp_30 = tmp_1*(tmp_28 - 1.0/3.0) + tmp_7*(tmp_29 - 1.0/3.0);
      real_t tmp_31 = -tmp_28 - tmp_29 + 1;
      real_t tmp_32 = 0.5*tmp_23;
      real_t tmp_33 = 0.5*p_affine_10_0*tmp_14 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_34 = 0.5*p_affine_10_0*tmp_17 + 0.5*p_affine_10_1*tmp_16;
      real_t a_0_0 = tmp_25*(-tmp_19*tmp_20 + 3*tmp_19*tmp_21*tmp_24 - tmp_21*tmp_22) + tmp_32*(-tmp_20*tmp_30 - tmp_22*tmp_31 + 3*tmp_24*tmp_30*tmp_31);
      real_t a_1_0 = tmp_25*(3*tmp_15*tmp_19*tmp_24 - tmp_15*tmp_22 - tmp_19*tmp_33) + tmp_32*(-tmp_22*tmp_28 + 3*tmp_24*tmp_28*tmp_30 - tmp_30*tmp_33);
      real_t a_2_0 = tmp_25*(3*tmp_18*tmp_19*tmp_24 - tmp_18*tmp_22 - tmp_19*tmp_34) + tmp_32*(-tmp_22*tmp_29 + 3*tmp_24*tmp_29*tmp_30 - tmp_30*tmp_34);
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
                                          MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = 1.0 / (tmp_0*tmp_1 - tmp_2*tmp_3);
      real_t tmp_5 = tmp_0*tmp_4;
      real_t tmp_6 = tmp_3*tmp_4;
      real_t tmp_7 = tmp_1*tmp_4;
      real_t tmp_8 = tmp_2*tmp_4;
      real_t tmp_9 = p_affine_10_0*(-tmp_5 - tmp_6) + p_affine_10_1*(-tmp_7 - tmp_8);
      real_t tmp_10 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_11 = -tmp_10;
      real_t tmp_12 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_13 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_14 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = 1.0 / (tmp_11*tmp_12 + tmp_15);
      real_t tmp_17 = -p_affine_3_1;
      real_t tmp_18 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_19 = p_affine_6_1 + 0.21132486540518713*tmp_18;
      real_t tmp_20 = tmp_16*(tmp_17 + tmp_19);
      real_t tmp_21 = -p_affine_3_0;
      real_t tmp_22 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_23 = p_affine_6_0 + 0.21132486540518713*tmp_22;
      real_t tmp_24 = tmp_16*(tmp_21 + tmp_23);
      real_t tmp_25 = tmp_11*(tmp_12*tmp_20 + tmp_14*tmp_24 - 1.0/3.0) + tmp_14*(tmp_10*tmp_24 + tmp_13*tmp_20 - 1.0/3.0);
      real_t tmp_26 = -p_affine_0_1;
      real_t tmp_27 = tmp_19 + tmp_26;
      real_t tmp_28 = -p_affine_0_0;
      real_t tmp_29 = tmp_23 + tmp_28;
      real_t tmp_30 = tmp_27*tmp_8 + tmp_29*tmp_5;
      real_t tmp_31 = tmp_27*tmp_7 + tmp_29*tmp_6;
      real_t tmp_32 = -tmp_30 - tmp_31 + 1;
      real_t tmp_33 = tmp_14*tmp_16;
      real_t tmp_34 = tmp_12*tmp_16;
      real_t tmp_35 = 0.5*p_affine_10_0*(tmp_10*tmp_33 + tmp_11*tmp_33) + 0.5*p_affine_10_1*(tmp_11*tmp_34 + tmp_15*tmp_16);
      real_t tmp_36 = std::abs(std::pow((tmp_18*tmp_18) + (tmp_22*tmp_22), 1.0/2.0));
      real_t tmp_37 = 3/tmp_36;
      real_t tmp_38 = tmp_25*tmp_37;
      real_t tmp_39 = 0.5*tmp_36;
      real_t tmp_40 = p_affine_6_1 + 0.78867513459481287*tmp_18;
      real_t tmp_41 = tmp_17 + tmp_40;
      real_t tmp_42 = p_affine_6_0 + 0.78867513459481287*tmp_22;
      real_t tmp_43 = tmp_21 + tmp_42;
      real_t tmp_44 = tmp_11*(tmp_33*tmp_43 + tmp_34*tmp_41 - 1.0/3.0) + tmp_14*(tmp_10*tmp_16*tmp_43 + tmp_13*tmp_16*tmp_41 - 1.0/3.0);
      real_t tmp_45 = tmp_26 + tmp_40;
      real_t tmp_46 = tmp_28 + tmp_42;
      real_t tmp_47 = tmp_45*tmp_8 + tmp_46*tmp_5;
      real_t tmp_48 = tmp_45*tmp_7 + tmp_46*tmp_6;
      real_t tmp_49 = -tmp_47 - tmp_48 + 1;
      real_t tmp_50 = tmp_37*tmp_44;
      real_t tmp_51 = 0.5*tmp_36;
      real_t tmp_52 = p_affine_10_0*tmp_5 + p_affine_10_1*tmp_8;
      real_t tmp_53 = p_affine_10_0*tmp_6 + p_affine_10_1*tmp_7;
      real_t a_0_0 = tmp_39*(0.5*tmp_25*tmp_9 - tmp_32*tmp_35 - tmp_32*tmp_38) + tmp_51*(-tmp_35*tmp_49 + 0.5*tmp_44*tmp_9 - tmp_49*tmp_50);
      real_t a_1_0 = tmp_39*(0.5*tmp_25*tmp_52 - tmp_30*tmp_35 - tmp_30*tmp_38) + tmp_51*(-tmp_35*tmp_47 + 0.5*tmp_44*tmp_52 - tmp_47*tmp_50);
      real_t a_2_0 = tmp_39*(0.5*tmp_25*tmp_53 - tmp_31*tmp_35 - tmp_31*tmp_38) + tmp_51*(-tmp_35*tmp_48 + 0.5*tmp_44*tmp_53 - tmp_48*tmp_50);
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
                                                   MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_5 = -tmp_4;
      real_t tmp_6 = 1.0 / (tmp_2 + tmp_3*tmp_5);
      real_t tmp_7 = tmp_0*tmp_6;
      real_t tmp_8 = tmp_4*tmp_6;
      real_t tmp_9 = tmp_1*tmp_6;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = p_affine_10_0*(-tmp_7 - tmp_8) + p_affine_10_1*(-tmp_10 - tmp_9);
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_14 = 0.21132486540518713*tmp_12 + tmp_13;
      real_t tmp_15 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_16 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_17 = 0.21132486540518713*tmp_15 + tmp_16;
      real_t tmp_18 = tmp_10*tmp_14 + tmp_17*tmp_7;
      real_t tmp_19 = tmp_14*tmp_9 + tmp_17*tmp_8;
      real_t tmp_20 = tmp_0*(tmp_19 - 1.0/3.0) + tmp_5*(tmp_18 - 1.0/3.0);
      real_t tmp_21 = p_affine_10_0*(tmp_0*tmp_8 + tmp_5*tmp_7) + p_affine_10_1*(tmp_10*tmp_5 + tmp_2*tmp_6);
      real_t tmp_22 = -tmp_18 - tmp_19 + 1;
      real_t tmp_23 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_15*tmp_15), 1.0/2.0));
      real_t tmp_24 = 1.0 / (tmp_23);
      real_t tmp_25 = 0.5*tmp_23;
      real_t tmp_26 = 0.78867513459481287*tmp_12 + tmp_13;
      real_t tmp_27 = 0.78867513459481287*tmp_15 + tmp_16;
      real_t tmp_28 = tmp_10*tmp_26 + tmp_27*tmp_7;
      real_t tmp_29 = tmp_26*tmp_9 + tmp_27*tmp_8;
      real_t tmp_30 = tmp_0*(tmp_29 - 1.0/3.0) + tmp_5*(tmp_28 - 1.0/3.0);
      real_t tmp_31 = -tmp_28 - tmp_29 + 1;
      real_t tmp_32 = 0.5*tmp_23;
      real_t tmp_33 = p_affine_10_0*tmp_7 + p_affine_10_1*tmp_10;
      real_t tmp_34 = p_affine_10_0*tmp_8 + p_affine_10_1*tmp_9;
      real_t a_0_0 = tmp_25*(-tmp_11*tmp_20 + 3*tmp_20*tmp_22*tmp_24 - tmp_21*tmp_22) + tmp_32*(-tmp_11*tmp_30 - tmp_21*tmp_31 + 3*tmp_24*tmp_30*tmp_31);
      real_t a_1_0 = tmp_25*(3*tmp_18*tmp_20*tmp_24 - tmp_18*tmp_21 - tmp_20*tmp_33) + tmp_32*(-tmp_21*tmp_28 + 3*tmp_24*tmp_28*tmp_30 - tmp_30*tmp_33);
      real_t a_2_0 = tmp_25*(3*tmp_19*tmp_20*tmp_24 - tmp_19*tmp_21 - tmp_20*tmp_34) + tmp_32*(-tmp_21*tmp_29 + 3*tmp_24*tmp_29*tmp_30 - tmp_30*tmp_34);
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
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_8 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_0*tmp_8;
      real_t tmp_13 = tmp_1*tmp_11;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_11*tmp_4 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_2*tmp_6 - tmp_5*tmp_6 + tmp_7*tmp_9);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_11*tmp_4 - tmp_12);
      real_t tmp_17 = tmp_14*(-tmp_13 + tmp_9);
      real_t tmp_18 = tmp_1*tmp_16 + tmp_15*tmp_8 + tmp_17*tmp_4;
      real_t tmp_19 = tmp_14*(-tmp_0*tmp_10 + tmp_3*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_0*tmp_6 - tmp_11*tmp_7);
      real_t tmp_21 = tmp_14*(tmp_10*tmp_11 - tmp_3*tmp_6);
      real_t tmp_22 = tmp_1*tmp_20 + tmp_19*tmp_8 + tmp_21*tmp_4;
      real_t tmp_23 = tmp_14*(-tmp_1*tmp_7 + tmp_10*tmp_4);
      real_t tmp_24 = tmp_14*(-tmp_4*tmp_6 + tmp_7*tmp_8);
      real_t tmp_25 = tmp_14*(tmp_1*tmp_6 - tmp_10*tmp_8);
      real_t tmp_26 = tmp_1*tmp_24 + tmp_23*tmp_8 + tmp_25*tmp_4;
      real_t tmp_27 = p_affine_0_0*p_affine_1_1;
      real_t tmp_28 = p_affine_0_0*p_affine_1_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_2;
      real_t tmp_30 = p_affine_0_1*p_affine_1_0;
      real_t tmp_31 = p_affine_0_1*p_affine_1_2;
      real_t tmp_32 = p_affine_2_2*p_affine_3_0;
      real_t tmp_33 = p_affine_0_2*p_affine_1_0;
      real_t tmp_34 = p_affine_0_2*p_affine_1_1;
      real_t tmp_35 = p_affine_2_0*p_affine_3_1;
      real_t tmp_36 = p_affine_2_2*p_affine_3_1;
      real_t tmp_37 = p_affine_2_0*p_affine_3_2;
      real_t tmp_38 = p_affine_2_1*p_affine_3_0;
      real_t tmp_39 = std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = tmp_39*(tmp_18*(-tmp_15 - tmp_16 - tmp_17) + tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25));
      real_t tmp_41 = tmp_39*(tmp_17*tmp_18 + tmp_21*tmp_22 + tmp_25*tmp_26);
      real_t tmp_42 = tmp_39*(tmp_16*tmp_18 + tmp_20*tmp_22 + tmp_24*tmp_26);
      real_t tmp_43 = tmp_39*(tmp_15*tmp_18 + tmp_19*tmp_22 + tmp_23*tmp_26);
      real_t a_0_0 = 0.16666666666666666*tmp_40;
      real_t a_1_0 = 0.16666666666666666*tmp_41;
      real_t a_2_0 = 0.16666666666666666*tmp_42;
      real_t a_3_0 = 0.16666666666666666*tmp_43;
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
                               MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_2 = -tmp_1;
      real_t tmp_3 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_4 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_5 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_3 + tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_7 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_8 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_9 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_10 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_11 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_12 = tmp_11*tmp_9;
      real_t tmp_13 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_14 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_15 = tmp_0*tmp_14;
      real_t tmp_16 = tmp_14*tmp_7;
      real_t tmp_17 = tmp_0*tmp_11;
      real_t tmp_18 = tmp_13*tmp_9;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_16 + tmp_13*tmp_6*tmp_7 + tmp_15*tmp_8 - tmp_17*tmp_6 - tmp_18*tmp_8);
      real_t tmp_20 = tmp_19*(tmp_6*tmp_7 - tmp_8*tmp_9);
      real_t tmp_21 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_22 = -tmp_21;
      real_t tmp_23 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_24 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_25 = 0.091576213509770743*tmp_22 + 0.81684757298045851*tmp_23 + tmp_24;
      real_t tmp_26 = tmp_19*(-tmp_11*tmp_6 + tmp_14*tmp_8);
      real_t tmp_27 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_28 = -tmp_27;
      real_t tmp_29 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_30 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_31 = 0.091576213509770743*tmp_28 + 0.81684757298045851*tmp_29 + tmp_30;
      real_t tmp_32 = tmp_19*(tmp_12 - tmp_16);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_25*tmp_26 + tmp_31*tmp_32;
      real_t tmp_34 = tmp_19*(tmp_0*tmp_8 - tmp_10*tmp_7);
      real_t tmp_35 = tmp_19*(tmp_10*tmp_11 - tmp_13*tmp_8);
      real_t tmp_36 = tmp_19*(tmp_13*tmp_7 - tmp_17);
      real_t tmp_37 = tmp_25*tmp_35 + tmp_31*tmp_36 + tmp_34*tmp_5;
      real_t tmp_38 = tmp_19*(-tmp_0*tmp_6 + tmp_10*tmp_9);
      real_t tmp_39 = tmp_19*(-tmp_10*tmp_14 + tmp_13*tmp_6);
      real_t tmp_40 = tmp_19*(tmp_15 - tmp_18);
      real_t tmp_41 = tmp_25*tmp_39 + tmp_31*tmp_40 + tmp_38*tmp_5;
      real_t tmp_42 = tmp_0*(tmp_33 - 1.0/4.0) + tmp_7*(tmp_41 - 1.0/4.0) + tmp_9*(tmp_37 - 1.0/4.0);
      real_t tmp_43 = 0.5*p_affine_13_0*(-tmp_32 - tmp_36 - tmp_40) + 0.5*p_affine_13_1*(-tmp_26 - tmp_35 - tmp_39) + 0.5*p_affine_13_2*(-tmp_20 - tmp_34 - tmp_38);
      real_t tmp_44 = -tmp_33 - tmp_37 - tmp_41 + 1;
      real_t tmp_45 = 0.5*p_affine_13_0*(tmp_0*tmp_32 + tmp_36*tmp_9 + tmp_40*tmp_7) + 0.5*p_affine_13_1*(tmp_0*tmp_26 + tmp_35*tmp_9 + tmp_39*tmp_7) + 0.5*p_affine_13_2*(tmp_0*tmp_20 + tmp_34*tmp_9 + tmp_38*tmp_7);
      real_t tmp_46 = (std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)*std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)) + (std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)*std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)) + (std::abs(tmp_21*tmp_29 - tmp_23*tmp_27)*std::abs(tmp_21*tmp_29 - tmp_23*tmp_27));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_3 + tmp_4;
      real_t tmp_51 = 0.44594849091596489*tmp_22 + 0.10810301816807022*tmp_23 + tmp_24;
      real_t tmp_52 = 0.44594849091596489*tmp_28 + 0.10810301816807022*tmp_29 + tmp_30;
      real_t tmp_53 = tmp_20*tmp_50 + tmp_26*tmp_51 + tmp_32*tmp_52;
      real_t tmp_54 = tmp_34*tmp_50 + tmp_35*tmp_51 + tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_50 + tmp_39*tmp_51 + tmp_40*tmp_52;
      real_t tmp_56 = tmp_0*(tmp_53 - 1.0/4.0) + tmp_7*(tmp_55 - 1.0/4.0) + tmp_9*(tmp_54 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_60 = 0.81684757298045851*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_61 = 0.81684757298045851*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_62 = tmp_20*tmp_59 + tmp_26*tmp_60 + tmp_32*tmp_61;
      real_t tmp_63 = tmp_34*tmp_59 + tmp_35*tmp_60 + tmp_36*tmp_61;
      real_t tmp_64 = tmp_38*tmp_59 + tmp_39*tmp_60 + tmp_40*tmp_61;
      real_t tmp_65 = tmp_0*(tmp_62 - 1.0/4.0) + tmp_7*(tmp_64 - 1.0/4.0) + tmp_9*(tmp_63 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_69 = 0.10810301816807022*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_70 = 0.10810301816807022*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_71 = tmp_20*tmp_68 + tmp_26*tmp_69 + tmp_32*tmp_70;
      real_t tmp_72 = tmp_34*tmp_68 + tmp_35*tmp_69 + tmp_36*tmp_70;
      real_t tmp_73 = tmp_38*tmp_68 + tmp_39*tmp_69 + tmp_40*tmp_70;
      real_t tmp_74 = tmp_0*(tmp_71 - 1.0/4.0) + tmp_7*(tmp_73 - 1.0/4.0) + tmp_9*(tmp_72 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_78 = 0.091576213509770743*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_79 = 0.091576213509770743*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_80 = tmp_20*tmp_77 + tmp_26*tmp_78 + tmp_32*tmp_79;
      real_t tmp_81 = tmp_34*tmp_77 + tmp_35*tmp_78 + tmp_36*tmp_79;
      real_t tmp_82 = tmp_38*tmp_77 + tmp_39*tmp_78 + tmp_40*tmp_79;
      real_t tmp_83 = tmp_0*(tmp_80 - 1.0/4.0) + tmp_7*(tmp_82 - 1.0/4.0) + tmp_9*(tmp_81 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_87 = 0.44594849091596489*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_88 = 0.44594849091596489*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_89 = tmp_20*tmp_86 + tmp_26*tmp_87 + tmp_32*tmp_88;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = tmp_0*(tmp_89 - 1.0/4.0) + tmp_7*(tmp_91 - 1.0/4.0) + tmp_9*(tmp_90 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = 0.5*p_affine_13_0*tmp_32 + 0.5*p_affine_13_1*tmp_26 + 0.5*p_affine_13_2*tmp_20;
      real_t tmp_96 = 0.5*p_affine_13_0*tmp_36 + 0.5*p_affine_13_1*tmp_35 + 0.5*p_affine_13_2*tmp_34;
      real_t tmp_97 = 0.5*p_affine_13_0*tmp_40 + 0.5*p_affine_13_1*tmp_39 + 0.5*p_affine_13_2*tmp_38;
      real_t a_0_0 = tmp_49*(-tmp_42*tmp_43 + 3.0*tmp_42*tmp_44*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_43*tmp_56 - tmp_45*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_43*tmp_65 - tmp_45*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_43*tmp_74 - tmp_45*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_43*tmp_83 - tmp_45*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_43*tmp_92 - tmp_45*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_1_0 = tmp_49*(3.0*tmp_33*tmp_42*tmp_47 - tmp_33*tmp_45 - tmp_42*tmp_95) + tmp_58*(-tmp_45*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_45*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_45*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_45*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_45*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_2_0 = tmp_49*(3.0*tmp_37*tmp_42*tmp_47 - tmp_37*tmp_45 - tmp_42*tmp_96) + tmp_58*(-tmp_45*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_45*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_45*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_45*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_45*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_3_0 = tmp_49*(3.0*tmp_41*tmp_42*tmp_47 - tmp_41*tmp_45 - tmp_42*tmp_97) + tmp_58*(-tmp_45*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_45*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_45*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_45*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_45*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
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
                                  MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_26 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_27 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_28 = tmp_26*tmp_27;
      real_t tmp_29 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_30 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_31 = tmp_29*tmp_30;
      real_t tmp_32 = tmp_28 - tmp_31;
      real_t tmp_33 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_34 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_35 = tmp_30*tmp_34;
      real_t tmp_36 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_37 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_38 = tmp_27*tmp_37;
      real_t tmp_39 = tmp_26*tmp_34;
      real_t tmp_40 = 1.0 / (tmp_25*tmp_29*tmp_37 - tmp_25*tmp_39 + tmp_28*tmp_36 - tmp_31*tmp_36 + tmp_33*tmp_35 - tmp_33*tmp_38);
      real_t tmp_41 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_42 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_43 = -tmp_42;
      real_t tmp_44 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_45 = 0.091576213509770743*tmp_43 + 0.81684757298045851*tmp_44;
      real_t tmp_46 = tmp_40*(tmp_41 + tmp_45);
      real_t tmp_47 = tmp_29*tmp_37 - tmp_39;
      real_t tmp_48 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_49 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_50 = -tmp_49;
      real_t tmp_51 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_52 = 0.091576213509770743*tmp_50 + 0.81684757298045851*tmp_51;
      real_t tmp_53 = tmp_40*(tmp_48 + tmp_52);
      real_t tmp_54 = tmp_35 - tmp_38;
      real_t tmp_55 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = -tmp_56;
      real_t tmp_58 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_59 = 0.091576213509770743*tmp_57 + 0.81684757298045851*tmp_58;
      real_t tmp_60 = tmp_40*(tmp_55 + tmp_59);
      real_t tmp_61 = tmp_25*tmp_29 - tmp_27*tmp_33;
      real_t tmp_62 = -tmp_29*tmp_36 + tmp_33*tmp_34;
      real_t tmp_63 = -tmp_25*tmp_34 + tmp_27*tmp_36;
      real_t tmp_64 = -tmp_25*tmp_26 + tmp_30*tmp_33;
      real_t tmp_65 = tmp_26*tmp_36 - tmp_33*tmp_37;
      real_t tmp_66 = tmp_25*tmp_37 - tmp_30*tmp_36;
      real_t tmp_67 = tmp_25*(tmp_32*tmp_46 + tmp_47*tmp_53 + tmp_54*tmp_60 - 1.0/4.0) + tmp_27*(tmp_46*tmp_64 + tmp_53*tmp_65 + tmp_60*tmp_66 - 1.0/4.0) + tmp_30*(tmp_46*tmp_61 + tmp_53*tmp_62 + tmp_60*tmp_63 - 1.0/4.0);
      real_t tmp_68 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_69 = tmp_45 + tmp_68;
      real_t tmp_70 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_71 = tmp_52 + tmp_70;
      real_t tmp_72 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_73 = tmp_59 + tmp_72;
      real_t tmp_74 = tmp_17*tmp_73 + tmp_20*tmp_71 + tmp_23*tmp_69;
      real_t tmp_75 = tmp_16*tmp_73 + tmp_19*tmp_71 + tmp_22*tmp_69;
      real_t tmp_76 = tmp_15*tmp_73 + tmp_18*tmp_71 + tmp_21*tmp_69;
      real_t tmp_77 = -tmp_74 - tmp_75 - tmp_76 + 1;
      real_t tmp_78 = tmp_25*tmp_40;
      real_t tmp_79 = tmp_30*tmp_40;
      real_t tmp_80 = tmp_27*tmp_40;
      real_t tmp_81 = 0.5*p_affine_13_0*(tmp_54*tmp_78 + tmp_63*tmp_79 + tmp_66*tmp_80) + 0.5*p_affine_13_1*(tmp_47*tmp_78 + tmp_62*tmp_79 + tmp_65*tmp_80) + 0.5*p_affine_13_2*(tmp_32*tmp_78 + tmp_61*tmp_79 + tmp_64*tmp_80);
      real_t tmp_82 = (std::abs(tmp_42*tmp_51 - tmp_44*tmp_49)*std::abs(tmp_42*tmp_51 - tmp_44*tmp_49)) + (std::abs(tmp_42*tmp_58 - tmp_44*tmp_56)*std::abs(tmp_42*tmp_58 - tmp_44*tmp_56)) + (std::abs(tmp_49*tmp_58 - tmp_51*tmp_56)*std::abs(tmp_49*tmp_58 - tmp_51*tmp_56));
      real_t tmp_83 = 3.0*std::pow(tmp_82, -0.25);
      real_t tmp_84 = tmp_67*tmp_83;
      real_t tmp_85 = 1.0*std::pow(tmp_82, 1.0/2.0);
      real_t tmp_86 = 0.054975871827660928*tmp_85;
      real_t tmp_87 = 0.44594849091596489*tmp_43 + 0.10810301816807022*tmp_44;
      real_t tmp_88 = tmp_40*(tmp_41 + tmp_87);
      real_t tmp_89 = 0.44594849091596489*tmp_50 + 0.10810301816807022*tmp_51;
      real_t tmp_90 = tmp_40*(tmp_48 + tmp_89);
      real_t tmp_91 = 0.44594849091596489*tmp_57 + 0.10810301816807022*tmp_58;
      real_t tmp_92 = tmp_40*(tmp_55 + tmp_91);
      real_t tmp_93 = tmp_25*(tmp_32*tmp_88 + tmp_47*tmp_90 + tmp_54*tmp_92 - 1.0/4.0) + tmp_27*(tmp_64*tmp_88 + tmp_65*tmp_90 + tmp_66*tmp_92 - 1.0/4.0) + tmp_30*(tmp_61*tmp_88 + tmp_62*tmp_90 + tmp_63*tmp_92 - 1.0/4.0);
      real_t tmp_94 = tmp_68 + tmp_87;
      real_t tmp_95 = tmp_70 + tmp_89;
      real_t tmp_96 = tmp_72 + tmp_91;
      real_t tmp_97 = tmp_17*tmp_96 + tmp_20*tmp_95 + tmp_23*tmp_94;
      real_t tmp_98 = tmp_16*tmp_96 + tmp_19*tmp_95 + tmp_22*tmp_94;
      real_t tmp_99 = tmp_15*tmp_96 + tmp_18*tmp_95 + tmp_21*tmp_94;
      real_t tmp_100 = -tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_101 = tmp_83*tmp_93;
      real_t tmp_102 = 0.11169079483900572*tmp_85;
      real_t tmp_103 = 0.81684757298045851*tmp_43 + 0.091576213509770743*tmp_44;
      real_t tmp_104 = tmp_40*(tmp_103 + tmp_41);
      real_t tmp_105 = 0.81684757298045851*tmp_50 + 0.091576213509770743*tmp_51;
      real_t tmp_106 = tmp_40*(tmp_105 + tmp_48);
      real_t tmp_107 = 0.81684757298045851*tmp_57 + 0.091576213509770743*tmp_58;
      real_t tmp_108 = tmp_40*(tmp_107 + tmp_55);
      real_t tmp_109 = tmp_25*(tmp_104*tmp_32 + tmp_106*tmp_47 + tmp_108*tmp_54 - 1.0/4.0) + tmp_27*(tmp_104*tmp_64 + tmp_106*tmp_65 + tmp_108*tmp_66 - 1.0/4.0) + tmp_30*(tmp_104*tmp_61 + tmp_106*tmp_62 + tmp_108*tmp_63 - 1.0/4.0);
      real_t tmp_110 = tmp_103 + tmp_68;
      real_t tmp_111 = tmp_105 + tmp_70;
      real_t tmp_112 = tmp_107 + tmp_72;
      real_t tmp_113 = tmp_110*tmp_23 + tmp_111*tmp_20 + tmp_112*tmp_17;
      real_t tmp_114 = tmp_110*tmp_22 + tmp_111*tmp_19 + tmp_112*tmp_16;
      real_t tmp_115 = tmp_110*tmp_21 + tmp_111*tmp_18 + tmp_112*tmp_15;
      real_t tmp_116 = -tmp_113 - tmp_114 - tmp_115 + 1;
      real_t tmp_117 = tmp_109*tmp_83;
      real_t tmp_118 = 0.054975871827660928*tmp_85;
      real_t tmp_119 = 0.10810301816807022*tmp_43 + 0.44594849091596489*tmp_44;
      real_t tmp_120 = tmp_40*(tmp_119 + tmp_41);
      real_t tmp_121 = 0.10810301816807022*tmp_50 + 0.44594849091596489*tmp_51;
      real_t tmp_122 = tmp_40*(tmp_121 + tmp_48);
      real_t tmp_123 = 0.10810301816807022*tmp_57 + 0.44594849091596489*tmp_58;
      real_t tmp_124 = tmp_40*(tmp_123 + tmp_55);
      real_t tmp_125 = tmp_25*(tmp_120*tmp_32 + tmp_122*tmp_47 + tmp_124*tmp_54 - 1.0/4.0) + tmp_27*(tmp_120*tmp_64 + tmp_122*tmp_65 + tmp_124*tmp_66 - 1.0/4.0) + tmp_30*(tmp_120*tmp_61 + tmp_122*tmp_62 + tmp_124*tmp_63 - 1.0/4.0);
      real_t tmp_126 = tmp_119 + tmp_68;
      real_t tmp_127 = tmp_121 + tmp_70;
      real_t tmp_128 = tmp_123 + tmp_72;
      real_t tmp_129 = tmp_126*tmp_23 + tmp_127*tmp_20 + tmp_128*tmp_17;
      real_t tmp_130 = tmp_126*tmp_22 + tmp_127*tmp_19 + tmp_128*tmp_16;
      real_t tmp_131 = tmp_126*tmp_21 + tmp_127*tmp_18 + tmp_128*tmp_15;
      real_t tmp_132 = -tmp_129 - tmp_130 - tmp_131 + 1;
      real_t tmp_133 = tmp_125*tmp_83;
      real_t tmp_134 = 0.11169079483900572*tmp_85;
      real_t tmp_135 = 0.091576213509770743*tmp_43 + 0.091576213509770743*tmp_44;
      real_t tmp_136 = tmp_40*(tmp_135 + tmp_41);
      real_t tmp_137 = 0.091576213509770743*tmp_50 + 0.091576213509770743*tmp_51;
      real_t tmp_138 = tmp_40*(tmp_137 + tmp_48);
      real_t tmp_139 = 0.091576213509770743*tmp_57 + 0.091576213509770743*tmp_58;
      real_t tmp_140 = tmp_40*(tmp_139 + tmp_55);
      real_t tmp_141 = tmp_25*(tmp_136*tmp_32 + tmp_138*tmp_47 + tmp_140*tmp_54 - 1.0/4.0) + tmp_27*(tmp_136*tmp_64 + tmp_138*tmp_65 + tmp_140*tmp_66 - 1.0/4.0) + tmp_30*(tmp_136*tmp_61 + tmp_138*tmp_62 + tmp_140*tmp_63 - 1.0/4.0);
      real_t tmp_142 = tmp_135 + tmp_68;
      real_t tmp_143 = tmp_137 + tmp_70;
      real_t tmp_144 = tmp_139 + tmp_72;
      real_t tmp_145 = tmp_142*tmp_23 + tmp_143*tmp_20 + tmp_144*tmp_17;
      real_t tmp_146 = tmp_142*tmp_22 + tmp_143*tmp_19 + tmp_144*tmp_16;
      real_t tmp_147 = tmp_142*tmp_21 + tmp_143*tmp_18 + tmp_144*tmp_15;
      real_t tmp_148 = -tmp_145 - tmp_146 - tmp_147 + 1;
      real_t tmp_149 = tmp_141*tmp_83;
      real_t tmp_150 = 0.054975871827660928*tmp_85;
      real_t tmp_151 = 0.44594849091596489*tmp_43 + 0.44594849091596489*tmp_44;
      real_t tmp_152 = tmp_40*(tmp_151 + tmp_41);
      real_t tmp_153 = 0.44594849091596489*tmp_50 + 0.44594849091596489*tmp_51;
      real_t tmp_154 = tmp_40*(tmp_153 + tmp_48);
      real_t tmp_155 = 0.44594849091596489*tmp_57 + 0.44594849091596489*tmp_58;
      real_t tmp_156 = tmp_40*(tmp_155 + tmp_55);
      real_t tmp_157 = tmp_25*(tmp_152*tmp_32 + tmp_154*tmp_47 + tmp_156*tmp_54 - 1.0/4.0) + tmp_27*(tmp_152*tmp_64 + tmp_154*tmp_65 + tmp_156*tmp_66 - 1.0/4.0) + tmp_30*(tmp_152*tmp_61 + tmp_154*tmp_62 + tmp_156*tmp_63 - 1.0/4.0);
      real_t tmp_158 = tmp_151 + tmp_68;
      real_t tmp_159 = tmp_153 + tmp_70;
      real_t tmp_160 = tmp_155 + tmp_72;
      real_t tmp_161 = tmp_158*tmp_23 + tmp_159*tmp_20 + tmp_160*tmp_17;
      real_t tmp_162 = tmp_158*tmp_22 + tmp_159*tmp_19 + tmp_16*tmp_160;
      real_t tmp_163 = tmp_15*tmp_160 + tmp_158*tmp_21 + tmp_159*tmp_18;
      real_t tmp_164 = -tmp_161 - tmp_162 - tmp_163 + 1;
      real_t tmp_165 = tmp_157*tmp_83;
      real_t tmp_166 = 0.11169079483900572*tmp_85;
      real_t tmp_167 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_168 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_169 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t a_0_0 = tmp_102*(-tmp_100*tmp_101 - tmp_100*tmp_81 + 0.5*tmp_24*tmp_93) + tmp_118*(0.5*tmp_109*tmp_24 - tmp_116*tmp_117 - tmp_116*tmp_81) + tmp_134*(0.5*tmp_125*tmp_24 - tmp_132*tmp_133 - tmp_132*tmp_81) + tmp_150*(0.5*tmp_141*tmp_24 - tmp_148*tmp_149 - tmp_148*tmp_81) + tmp_166*(0.5*tmp_157*tmp_24 - tmp_164*tmp_165 - tmp_164*tmp_81) + tmp_86*(0.5*tmp_24*tmp_67 - tmp_77*tmp_81 - tmp_77*tmp_84);
      real_t a_1_0 = tmp_102*(-tmp_101*tmp_97 + 0.5*tmp_167*tmp_93 - tmp_81*tmp_97) + tmp_118*(0.5*tmp_109*tmp_167 - tmp_113*tmp_117 - tmp_113*tmp_81) + tmp_134*(0.5*tmp_125*tmp_167 - tmp_129*tmp_133 - tmp_129*tmp_81) + tmp_150*(0.5*tmp_141*tmp_167 - tmp_145*tmp_149 - tmp_145*tmp_81) + tmp_166*(0.5*tmp_157*tmp_167 - tmp_161*tmp_165 - tmp_161*tmp_81) + tmp_86*(0.5*tmp_167*tmp_67 - tmp_74*tmp_81 - tmp_74*tmp_84);
      real_t a_2_0 = tmp_102*(-tmp_101*tmp_98 + 0.5*tmp_168*tmp_93 - tmp_81*tmp_98) + tmp_118*(0.5*tmp_109*tmp_168 - tmp_114*tmp_117 - tmp_114*tmp_81) + tmp_134*(0.5*tmp_125*tmp_168 - tmp_130*tmp_133 - tmp_130*tmp_81) + tmp_150*(0.5*tmp_141*tmp_168 - tmp_146*tmp_149 - tmp_146*tmp_81) + tmp_166*(0.5*tmp_157*tmp_168 - tmp_162*tmp_165 - tmp_162*tmp_81) + tmp_86*(0.5*tmp_168*tmp_67 - tmp_75*tmp_81 - tmp_75*tmp_84);
      real_t a_3_0 = tmp_102*(-tmp_101*tmp_99 + 0.5*tmp_169*tmp_93 - tmp_81*tmp_99) + tmp_118*(0.5*tmp_109*tmp_169 - tmp_115*tmp_117 - tmp_115*tmp_81) + tmp_134*(0.5*tmp_125*tmp_169 - tmp_131*tmp_133 - tmp_131*tmp_81) + tmp_150*(0.5*tmp_141*tmp_169 - tmp_147*tmp_149 - tmp_147*tmp_81) + tmp_166*(0.5*tmp_157*tmp_169 - tmp_163*tmp_165 - tmp_163*tmp_81) + tmp_86*(0.5*tmp_169*tmp_67 - tmp_76*tmp_81 - tmp_76*tmp_84);
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
                                        MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_26 = -tmp_25;
      real_t tmp_27 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_28 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_29 = 0.091576213509770743*tmp_26 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -tmp_30;
      real_t tmp_32 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_33 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_34 = 0.091576213509770743*tmp_31 + 0.81684757298045851*tmp_32 + tmp_33;
      real_t tmp_35 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_36 = -tmp_35;
      real_t tmp_37 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_38 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_39 = 0.091576213509770743*tmp_36 + 0.81684757298045851*tmp_37 + tmp_38;
      real_t tmp_40 = tmp_17*tmp_39 + tmp_20*tmp_34 + tmp_23*tmp_29;
      real_t tmp_41 = tmp_16*tmp_39 + tmp_19*tmp_34 + tmp_22*tmp_29;
      real_t tmp_42 = tmp_15*tmp_39 + tmp_18*tmp_34 + tmp_21*tmp_29;
      real_t tmp_43 = tmp_0*(tmp_40 - 1.0/4.0) + tmp_10*(tmp_42 - 1.0/4.0) + tmp_3*(tmp_41 - 1.0/4.0);
      real_t tmp_44 = p_affine_13_0*(tmp_0*tmp_17 + tmp_10*tmp_15 + tmp_16*tmp_3) + p_affine_13_1*(tmp_0*tmp_20 + tmp_10*tmp_18 + tmp_19*tmp_3) + p_affine_13_2*(tmp_0*tmp_23 + tmp_10*tmp_21 + tmp_22*tmp_3);
      real_t tmp_45 = -tmp_40 - tmp_41 - tmp_42 + 1;
      real_t tmp_46 = (std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)*std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)) + (std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)*std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)) + (std::abs(tmp_30*tmp_37 - tmp_32*tmp_35)*std::abs(tmp_30*tmp_37 - tmp_32*tmp_35));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_26 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_51 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32 + tmp_33;
      real_t tmp_52 = 0.44594849091596489*tmp_36 + 0.10810301816807022*tmp_37 + tmp_38;
      real_t tmp_53 = tmp_17*tmp_52 + tmp_20*tmp_51 + tmp_23*tmp_50;
      real_t tmp_54 = tmp_16*tmp_52 + tmp_19*tmp_51 + tmp_22*tmp_50;
      real_t tmp_55 = tmp_15*tmp_52 + tmp_18*tmp_51 + tmp_21*tmp_50;
      real_t tmp_56 = tmp_0*(tmp_53 - 1.0/4.0) + tmp_10*(tmp_55 - 1.0/4.0) + tmp_3*(tmp_54 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_60 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_61 = 0.81684757298045851*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_62 = tmp_17*tmp_61 + tmp_20*tmp_60 + tmp_23*tmp_59;
      real_t tmp_63 = tmp_16*tmp_61 + tmp_19*tmp_60 + tmp_22*tmp_59;
      real_t tmp_64 = tmp_15*tmp_61 + tmp_18*tmp_60 + tmp_21*tmp_59;
      real_t tmp_65 = tmp_0*(tmp_62 - 1.0/4.0) + tmp_10*(tmp_64 - 1.0/4.0) + tmp_3*(tmp_63 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_69 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_70 = 0.10810301816807022*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_71 = tmp_17*tmp_70 + tmp_20*tmp_69 + tmp_23*tmp_68;
      real_t tmp_72 = tmp_16*tmp_70 + tmp_19*tmp_69 + tmp_22*tmp_68;
      real_t tmp_73 = tmp_15*tmp_70 + tmp_18*tmp_69 + tmp_21*tmp_68;
      real_t tmp_74 = tmp_0*(tmp_71 - 1.0/4.0) + tmp_10*(tmp_73 - 1.0/4.0) + tmp_3*(tmp_72 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_78 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_79 = 0.091576213509770743*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_80 = tmp_17*tmp_79 + tmp_20*tmp_78 + tmp_23*tmp_77;
      real_t tmp_81 = tmp_16*tmp_79 + tmp_19*tmp_78 + tmp_22*tmp_77;
      real_t tmp_82 = tmp_15*tmp_79 + tmp_18*tmp_78 + tmp_21*tmp_77;
      real_t tmp_83 = tmp_0*(tmp_80 - 1.0/4.0) + tmp_10*(tmp_82 - 1.0/4.0) + tmp_3*(tmp_81 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_87 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_88 = 0.44594849091596489*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_89 = tmp_17*tmp_88 + tmp_20*tmp_87 + tmp_23*tmp_86;
      real_t tmp_90 = tmp_16*tmp_88 + tmp_19*tmp_87 + tmp_22*tmp_86;
      real_t tmp_91 = tmp_15*tmp_88 + tmp_18*tmp_87 + tmp_21*tmp_86;
      real_t tmp_92 = tmp_0*(tmp_89 - 1.0/4.0) + tmp_10*(tmp_91 - 1.0/4.0) + tmp_3*(tmp_90 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_96 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_97 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t a_0_0 = tmp_49*(-tmp_24*tmp_43 + 3.0*tmp_43*tmp_45*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_24*tmp_56 - tmp_44*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_24*tmp_65 - tmp_44*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_24*tmp_74 - tmp_44*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_24*tmp_83 - tmp_44*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_24*tmp_92 - tmp_44*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_1_0 = tmp_49*(3.0*tmp_40*tmp_43*tmp_47 - tmp_40*tmp_44 - tmp_43*tmp_95) + tmp_58*(-tmp_44*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_44*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_44*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_44*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_44*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_2_0 = tmp_49*(3.0*tmp_41*tmp_43*tmp_47 - tmp_41*tmp_44 - tmp_43*tmp_96) + tmp_58*(-tmp_44*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_44*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_44*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_44*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_44*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_3_0 = tmp_49*(3.0*tmp_42*tmp_43*tmp_47 - tmp_42*tmp_44 - tmp_43*tmp_97) + tmp_58*(-tmp_44*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_44*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_44*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_44*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_44*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }

public:



};




class EGVectorLaplaceFormNitscheBC_P1P1_22 : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_P1P1_22()
: callback_Scalar_Variable_Coefficient_3D_g2 ([](const Point3D & p) -> real_t { return 0.; })
    {}



void Scalar_Variable_Coefficient_3D_g2( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g2( Point3D( {in_0, in_1, in_2} ) );
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
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const override
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
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
                                          MatrixXr&                                           elMat ) const override
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
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const override
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
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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

      real_t a_0_0 = 0;
      real_t a_1_0 = 0;
      real_t a_2_0 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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

      real_t Scalar_Variable_Coefficient_3D_g2_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id5 = 0;
      Scalar_Variable_Coefficient_3D_g2( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id0 );
      Scalar_Variable_Coefficient_3D_g2( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id1 );
      Scalar_Variable_Coefficient_3D_g2( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id2 );
      Scalar_Variable_Coefficient_3D_g2( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id3 );
      Scalar_Variable_Coefficient_3D_g2( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id4 );
      Scalar_Variable_Coefficient_3D_g2( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id5 );
      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -tmp_4;
      real_t tmp_9 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_10 = 0.81684757298045851*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_11 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_13 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_14 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_15 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_16 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_19 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_20 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_12*tmp_20;
      real_t tmp_23 = tmp_16*tmp_19;
      real_t tmp_24 = tmp_14*tmp_18;
      real_t tmp_25 = 1.0 / (tmp_11*tmp_12*tmp_18 - tmp_11*tmp_23 + tmp_13*tmp_21 - tmp_13*tmp_24 + tmp_15*tmp_17 - tmp_15*tmp_22);
      real_t tmp_26 = tmp_25*(tmp_11*tmp_12 - tmp_13*tmp_14);
      real_t tmp_27 = -tmp_1;
      real_t tmp_28 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_29 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_30 = tmp_25*(-tmp_11*tmp_16 + tmp_13*tmp_20);
      real_t tmp_31 = -tmp_3;
      real_t tmp_32 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_33 = 0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_34 = tmp_25*(tmp_17 - tmp_22);
      real_t tmp_35 = tmp_10*tmp_26 + tmp_29*tmp_30 + tmp_33*tmp_34;
      real_t tmp_36 = tmp_25*(-tmp_12*tmp_15 + tmp_13*tmp_19);
      real_t tmp_37 = tmp_25*(-tmp_13*tmp_18 + tmp_15*tmp_16);
      real_t tmp_38 = tmp_25*(tmp_12*tmp_18 - tmp_23);
      real_t tmp_39 = tmp_10*tmp_36 + tmp_29*tmp_37 + tmp_33*tmp_38;
      real_t tmp_40 = tmp_25*(-tmp_11*tmp_19 + tmp_14*tmp_15);
      real_t tmp_41 = tmp_25*(tmp_11*tmp_18 - tmp_15*tmp_20);
      real_t tmp_42 = tmp_25*(tmp_21 - tmp_24);
      real_t tmp_43 = tmp_10*tmp_40 + tmp_29*tmp_41 + tmp_33*tmp_42;
      real_t tmp_44 = p_affine_13_0*(-tmp_34 - tmp_38 - tmp_42) + p_affine_13_1*(-tmp_30 - tmp_37 - tmp_41) + p_affine_13_2*(-tmp_26 - tmp_36 - tmp_40);
      real_t tmp_45 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_46 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_g2_out0_id0*tmp_45;
      real_t tmp_47 = 0.10810301816807022*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_48 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_49 = 0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_50 = tmp_26*tmp_47 + tmp_30*tmp_48 + tmp_34*tmp_49;
      real_t tmp_51 = tmp_36*tmp_47 + tmp_37*tmp_48 + tmp_38*tmp_49;
      real_t tmp_52 = tmp_40*tmp_47 + tmp_41*tmp_48 + tmp_42*tmp_49;
      real_t tmp_53 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_g2_out0_id1*tmp_45;
      real_t tmp_54 = 0.091576213509770743*tmp_5 + 0.81684757298045851*tmp_8 + tmp_9;
      real_t tmp_55 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_56 = 0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_31 + tmp_32;
      real_t tmp_57 = tmp_26*tmp_54 + tmp_30*tmp_55 + tmp_34*tmp_56;
      real_t tmp_58 = tmp_36*tmp_54 + tmp_37*tmp_55 + tmp_38*tmp_56;
      real_t tmp_59 = tmp_40*tmp_54 + tmp_41*tmp_55 + tmp_42*tmp_56;
      real_t tmp_60 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_g2_out0_id2*tmp_45;
      real_t tmp_61 = 0.44594849091596489*tmp_5 + 0.10810301816807022*tmp_8 + tmp_9;
      real_t tmp_62 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_63 = 0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_31 + tmp_32;
      real_t tmp_64 = tmp_26*tmp_61 + tmp_30*tmp_62 + tmp_34*tmp_63;
      real_t tmp_65 = tmp_36*tmp_61 + tmp_37*tmp_62 + tmp_38*tmp_63;
      real_t tmp_66 = tmp_40*tmp_61 + tmp_41*tmp_62 + tmp_42*tmp_63;
      real_t tmp_67 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_g2_out0_id3*tmp_45;
      real_t tmp_68 = 0.091576213509770743*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_69 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_70 = 0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_71 = tmp_26*tmp_68 + tmp_30*tmp_69 + tmp_34*tmp_70;
      real_t tmp_72 = tmp_36*tmp_68 + tmp_37*tmp_69 + tmp_38*tmp_70;
      real_t tmp_73 = tmp_40*tmp_68 + tmp_41*tmp_69 + tmp_42*tmp_70;
      real_t tmp_74 = 0.054975871827660928*Scalar_Variable_Coefficient_3D_g2_out0_id4*tmp_45;
      real_t tmp_75 = 0.44594849091596489*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_76 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_77 = 0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_78 = tmp_26*tmp_75 + tmp_30*tmp_76 + tmp_34*tmp_77;
      real_t tmp_79 = tmp_36*tmp_75 + tmp_37*tmp_76 + tmp_38*tmp_77;
      real_t tmp_80 = tmp_40*tmp_75 + tmp_41*tmp_76 + tmp_42*tmp_77;
      real_t tmp_81 = 0.11169079483900572*Scalar_Variable_Coefficient_3D_g2_out0_id5*tmp_45;
      real_t tmp_82 = p_affine_13_0*tmp_34 + p_affine_13_1*tmp_30 + p_affine_13_2*tmp_26;
      real_t tmp_83 = p_affine_13_0*tmp_38 + p_affine_13_1*tmp_37 + p_affine_13_2*tmp_36;
      real_t tmp_84 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_40;
      real_t a_0_0 = tmp_46*(-tmp_44 + 3.0*tmp_7*(-tmp_35 - tmp_39 - tmp_43 + 1)) + tmp_53*(-tmp_44 + 3.0*tmp_7*(-tmp_50 - tmp_51 - tmp_52 + 1)) + tmp_60*(-tmp_44 + 3.0*tmp_7*(-tmp_57 - tmp_58 - tmp_59 + 1)) + tmp_67*(-tmp_44 + 3.0*tmp_7*(-tmp_64 - tmp_65 - tmp_66 + 1)) + tmp_74*(-tmp_44 + 3.0*tmp_7*(-tmp_71 - tmp_72 - tmp_73 + 1)) + tmp_81*(-tmp_44 + 3.0*tmp_7*(-tmp_78 - tmp_79 - tmp_80 + 1));
      real_t a_1_0 = tmp_46*(3.0*tmp_35*tmp_7 - tmp_82) + tmp_53*(3.0*tmp_50*tmp_7 - tmp_82) + tmp_60*(3.0*tmp_57*tmp_7 - tmp_82) + tmp_67*(3.0*tmp_64*tmp_7 - tmp_82) + tmp_74*(3.0*tmp_7*tmp_71 - tmp_82) + tmp_81*(3.0*tmp_7*tmp_78 - tmp_82);
      real_t a_2_0 = tmp_46*(3.0*tmp_39*tmp_7 - tmp_83) + tmp_53*(3.0*tmp_51*tmp_7 - tmp_83) + tmp_60*(3.0*tmp_58*tmp_7 - tmp_83) + tmp_67*(3.0*tmp_65*tmp_7 - tmp_83) + tmp_74*(3.0*tmp_7*tmp_72 - tmp_83) + tmp_81*(3.0*tmp_7*tmp_79 - tmp_83);
      real_t a_3_0 = tmp_46*(3.0*tmp_43*tmp_7 - tmp_84) + tmp_53*(3.0*tmp_52*tmp_7 - tmp_84) + tmp_60*(3.0*tmp_59*tmp_7 - tmp_84) + tmp_67*(3.0*tmp_66*tmp_7 - tmp_84) + tmp_74*(3.0*tmp_7*tmp_73 - tmp_84) + tmp_81*(3.0*tmp_7*tmp_80 - tmp_84);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = tmp_2 - tmp_5;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_10 = tmp_3*tmp_9;
      real_t tmp_11 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_12 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_1*tmp_12;
      real_t tmp_15 = tmp_10*tmp_8 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_13 - tmp_14*tmp_8 + tmp_2*tmp_7 - tmp_5*tmp_7;
      real_t tmp_16 = 1.0 / (tmp_15);
      real_t tmp_17 = tmp_16*tmp_6;
      real_t tmp_18 = tmp_12*tmp_4 - tmp_13;
      real_t tmp_19 = tmp_16*tmp_18;
      real_t tmp_20 = tmp_10 - tmp_14;
      real_t tmp_21 = tmp_16*tmp_20;
      real_t tmp_22 = -tmp_17 - tmp_19 - tmp_21;
      real_t tmp_23 = -tmp_0*tmp_11 + tmp_3*tmp_8;
      real_t tmp_24 = tmp_16*tmp_23;
      real_t tmp_25 = tmp_0*tmp_7 - tmp_12*tmp_8;
      real_t tmp_26 = tmp_16*tmp_25;
      real_t tmp_27 = tmp_11*tmp_12 - tmp_3*tmp_7;
      real_t tmp_28 = tmp_16*tmp_27;
      real_t tmp_29 = -tmp_24 - tmp_26 - tmp_28;
      real_t tmp_30 = -tmp_1*tmp_8 + tmp_11*tmp_4;
      real_t tmp_31 = tmp_16*tmp_30;
      real_t tmp_32 = -tmp_4*tmp_7 + tmp_8*tmp_9;
      real_t tmp_33 = tmp_16*tmp_32;
      real_t tmp_34 = tmp_1*tmp_7 - tmp_11*tmp_9;
      real_t tmp_35 = tmp_16*tmp_34;
      real_t tmp_36 = -tmp_31 - tmp_33 - tmp_35;
      real_t tmp_37 = p_affine_0_0*p_affine_1_1;
      real_t tmp_38 = p_affine_0_0*p_affine_1_2;
      real_t tmp_39 = p_affine_2_1*p_affine_3_2;
      real_t tmp_40 = p_affine_0_1*p_affine_1_0;
      real_t tmp_41 = p_affine_0_1*p_affine_1_2;
      real_t tmp_42 = p_affine_2_2*p_affine_3_0;
      real_t tmp_43 = p_affine_0_2*p_affine_1_0;
      real_t tmp_44 = p_affine_0_2*p_affine_1_1;
      real_t tmp_45 = p_affine_2_0*p_affine_3_1;
      real_t tmp_46 = p_affine_2_2*p_affine_3_1;
      real_t tmp_47 = p_affine_2_0*p_affine_3_2;
      real_t tmp_48 = p_affine_2_1*p_affine_3_0;
      real_t tmp_49 = std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = tmp_49*((tmp_22*tmp_22) + (tmp_29*tmp_29) + (tmp_36*tmp_36));
      real_t tmp_51 = tmp_49*(tmp_21*tmp_22 + tmp_28*tmp_29 + tmp_35*tmp_36);
      real_t tmp_52 = 0.16666666666666666*tmp_51;
      real_t tmp_53 = tmp_49*(tmp_19*tmp_22 + tmp_26*tmp_29 + tmp_33*tmp_36);
      real_t tmp_54 = 0.16666666666666666*tmp_53;
      real_t tmp_55 = tmp_49*(tmp_17*tmp_22 + tmp_24*tmp_29 + tmp_31*tmp_36);
      real_t tmp_56 = 0.16666666666666666*tmp_55;
      real_t tmp_57 = 1.0 / (tmp_15*tmp_15);
      real_t tmp_58 = tmp_49*((tmp_20*tmp_20)*tmp_57 + (tmp_27*tmp_27)*tmp_57 + (tmp_34*tmp_34)*tmp_57);
      real_t tmp_59 = tmp_20*tmp_57;
      real_t tmp_60 = tmp_27*tmp_57;
      real_t tmp_61 = tmp_34*tmp_57;
      real_t tmp_62 = tmp_49*(tmp_18*tmp_59 + tmp_25*tmp_60 + tmp_32*tmp_61);
      real_t tmp_63 = 0.16666666666666666*tmp_62;
      real_t tmp_64 = tmp_49*(tmp_23*tmp_60 + tmp_30*tmp_61 + tmp_59*tmp_6);
      real_t tmp_65 = 0.16666666666666666*tmp_64;
      real_t tmp_66 = tmp_49*((tmp_18*tmp_18)*tmp_57 + (tmp_25*tmp_25)*tmp_57 + (tmp_32*tmp_32)*tmp_57);
      real_t tmp_67 = tmp_49*(tmp_18*tmp_57*tmp_6 + tmp_23*tmp_25*tmp_57 + tmp_30*tmp_32*tmp_57);
      real_t tmp_68 = 0.16666666666666666*tmp_67;
      real_t tmp_69 = tmp_49*((tmp_23*tmp_23)*tmp_57 + (tmp_30*tmp_30)*tmp_57 + tmp_57*(tmp_6*tmp_6));
      real_t a_0_0 = 0.16666666666666666*tmp_50;
      real_t a_0_1 = tmp_52;
      real_t a_0_2 = tmp_54;
      real_t a_0_3 = tmp_56;
      real_t a_1_0 = tmp_52;
      real_t a_1_1 = 0.16666666666666666*tmp_58;
      real_t a_1_2 = tmp_63;
      real_t a_1_3 = tmp_65;
      real_t a_2_0 = tmp_54;
      real_t a_2_1 = tmp_63;
      real_t a_2_2 = 0.16666666666666666*tmp_66;
      real_t a_2_3 = tmp_68;
      real_t a_3_0 = tmp_56;
      real_t a_3_1 = tmp_65;
      real_t a_3_2 = tmp_68;
      real_t a_3_3 = 0.16666666666666666*tmp_69;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const override
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

         real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -tmp_4;
      real_t tmp_9 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_10 = 0.81684757298045851*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_11 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_13 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_14 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_15 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_16 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_19 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_20 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_12*tmp_16;
      real_t tmp_23 = tmp_11*tmp_19;
      real_t tmp_24 = tmp_13*tmp_18;
      real_t tmp_25 = 1.0 / (tmp_11*tmp_12*tmp_18 + tmp_13*tmp_21 - tmp_14*tmp_24 + tmp_15*tmp_17 - tmp_15*tmp_23 - tmp_20*tmp_22);
      real_t tmp_26 = tmp_25*(tmp_11*tmp_12 - tmp_13*tmp_14);
      real_t tmp_27 = -tmp_1;
      real_t tmp_28 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_29 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_30 = tmp_25*(-tmp_11*tmp_15 + tmp_13*tmp_20);
      real_t tmp_31 = -tmp_3;
      real_t tmp_32 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_33 = 0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_34 = tmp_25*(-tmp_12*tmp_20 + tmp_14*tmp_15);
      real_t tmp_35 = tmp_10*tmp_26 + tmp_29*tmp_30 + tmp_33*tmp_34;
      real_t tmp_36 = tmp_25*(tmp_13*tmp_19 - tmp_22);
      real_t tmp_37 = tmp_25*(tmp_15*tmp_16 - tmp_24);
      real_t tmp_38 = tmp_25*(tmp_12*tmp_18 - tmp_15*tmp_19);
      real_t tmp_39 = tmp_10*tmp_36 + tmp_29*tmp_37 + tmp_33*tmp_38;
      real_t tmp_40 = tmp_25*(tmp_17 - tmp_23);
      real_t tmp_41 = tmp_25*(tmp_11*tmp_18 - tmp_16*tmp_20);
      real_t tmp_42 = tmp_25*(-tmp_14*tmp_18 + tmp_21);
      real_t tmp_43 = tmp_10*tmp_40 + tmp_29*tmp_41 + tmp_33*tmp_42;
      real_t tmp_44 = -tmp_35 - tmp_39 - tmp_43 + 1;
      real_t tmp_45 = p_affine_13_0*(-tmp_34 - tmp_38 - tmp_42) + p_affine_13_1*(-tmp_30 - tmp_37 - tmp_41) + p_affine_13_2*(-tmp_26 - tmp_36 - tmp_40);
      real_t tmp_46 = 1.0*tmp_45;
      real_t tmp_47 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_48 = 0.054975871827660928*tmp_47;
      real_t tmp_49 = 0.10810301816807022*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_50 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_51 = 0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_52 = tmp_26*tmp_49 + tmp_30*tmp_50 + tmp_34*tmp_51;
      real_t tmp_53 = tmp_36*tmp_49 + tmp_37*tmp_50 + tmp_38*tmp_51;
      real_t tmp_54 = tmp_40*tmp_49 + tmp_41*tmp_50 + tmp_42*tmp_51;
      real_t tmp_55 = -tmp_52 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11169079483900572*tmp_47;
      real_t tmp_57 = 0.091576213509770743*tmp_5 + 0.81684757298045851*tmp_8 + tmp_9;
      real_t tmp_58 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_59 = 0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_31 + tmp_32;
      real_t tmp_60 = tmp_26*tmp_57 + tmp_30*tmp_58 + tmp_34*tmp_59;
      real_t tmp_61 = tmp_36*tmp_57 + tmp_37*tmp_58 + tmp_38*tmp_59;
      real_t tmp_62 = tmp_40*tmp_57 + tmp_41*tmp_58 + tmp_42*tmp_59;
      real_t tmp_63 = -tmp_60 - tmp_61 - tmp_62 + 1;
      real_t tmp_64 = 0.054975871827660928*tmp_47;
      real_t tmp_65 = 0.44594849091596489*tmp_5 + 0.10810301816807022*tmp_8 + tmp_9;
      real_t tmp_66 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_67 = 0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_31 + tmp_32;
      real_t tmp_68 = tmp_26*tmp_65 + tmp_30*tmp_66 + tmp_34*tmp_67;
      real_t tmp_69 = tmp_36*tmp_65 + tmp_37*tmp_66 + tmp_38*tmp_67;
      real_t tmp_70 = tmp_40*tmp_65 + tmp_41*tmp_66 + tmp_42*tmp_67;
      real_t tmp_71 = -tmp_68 - tmp_69 - tmp_70 + 1;
      real_t tmp_72 = 0.11169079483900572*tmp_47;
      real_t tmp_73 = 0.091576213509770743*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_74 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_75 = 0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_76 = tmp_26*tmp_73 + tmp_30*tmp_74 + tmp_34*tmp_75;
      real_t tmp_77 = tmp_36*tmp_73 + tmp_37*tmp_74 + tmp_38*tmp_75;
      real_t tmp_78 = tmp_40*tmp_73 + tmp_41*tmp_74 + tmp_42*tmp_75;
      real_t tmp_79 = -tmp_76 - tmp_77 - tmp_78 + 1;
      real_t tmp_80 = 0.054975871827660928*tmp_47;
      real_t tmp_81 = 0.44594849091596489*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_82 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_83 = 0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_84 = tmp_26*tmp_81 + tmp_30*tmp_82 + tmp_34*tmp_83;
      real_t tmp_85 = tmp_36*tmp_81 + tmp_37*tmp_82 + tmp_38*tmp_83;
      real_t tmp_86 = tmp_40*tmp_81 + tmp_41*tmp_82 + tmp_42*tmp_83;
      real_t tmp_87 = -tmp_84 - tmp_85 - tmp_86 + 1;
      real_t tmp_88 = 0.11169079483900572*tmp_47;
      real_t tmp_89 = 0.5*tmp_45;
      real_t tmp_90 = p_affine_13_0*tmp_34 + p_affine_13_1*tmp_30 + p_affine_13_2*tmp_26;
      real_t tmp_91 = 0.5*tmp_90;
      real_t tmp_92 = tmp_48*(3.0*tmp_35*tmp_44*tmp_7 - tmp_35*tmp_89 - tmp_44*tmp_91) + tmp_56*(3.0*tmp_52*tmp_55*tmp_7 - tmp_52*tmp_89 - tmp_55*tmp_91) + tmp_64*(3.0*tmp_60*tmp_63*tmp_7 - tmp_60*tmp_89 - tmp_63*tmp_91) + tmp_72*(3.0*tmp_68*tmp_7*tmp_71 - tmp_68*tmp_89 - tmp_71*tmp_91) + tmp_80*(3.0*tmp_7*tmp_76*tmp_79 - tmp_76*tmp_89 - tmp_79*tmp_91) + tmp_88*(3.0*tmp_7*tmp_84*tmp_87 - tmp_84*tmp_89 - tmp_87*tmp_91);
      real_t tmp_93 = p_affine_13_0*tmp_38 + p_affine_13_1*tmp_37 + p_affine_13_2*tmp_36;
      real_t tmp_94 = 0.5*tmp_93;
      real_t tmp_95 = tmp_48*(3.0*tmp_39*tmp_44*tmp_7 - tmp_39*tmp_89 - tmp_44*tmp_94) + tmp_56*(3.0*tmp_53*tmp_55*tmp_7 - tmp_53*tmp_89 - tmp_55*tmp_94) + tmp_64*(3.0*tmp_61*tmp_63*tmp_7 - tmp_61*tmp_89 - tmp_63*tmp_94) + tmp_72*(3.0*tmp_69*tmp_7*tmp_71 - tmp_69*tmp_89 - tmp_71*tmp_94) + tmp_80*(3.0*tmp_7*tmp_77*tmp_79 - tmp_77*tmp_89 - tmp_79*tmp_94) + tmp_88*(3.0*tmp_7*tmp_85*tmp_87 - tmp_85*tmp_89 - tmp_87*tmp_94);
      real_t tmp_96 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_40;
      real_t tmp_97 = 0.5*tmp_96;
      real_t tmp_98 = tmp_48*(3.0*tmp_43*tmp_44*tmp_7 - tmp_43*tmp_89 - tmp_44*tmp_97) + tmp_56*(3.0*tmp_54*tmp_55*tmp_7 - tmp_54*tmp_89 - tmp_55*tmp_97) + tmp_64*(3.0*tmp_62*tmp_63*tmp_7 - tmp_62*tmp_89 - tmp_63*tmp_97) + tmp_72*(3.0*tmp_7*tmp_70*tmp_71 - tmp_70*tmp_89 - tmp_71*tmp_97) + tmp_80*(3.0*tmp_7*tmp_78*tmp_79 - tmp_78*tmp_89 - tmp_79*tmp_97) + tmp_88*(3.0*tmp_7*tmp_86*tmp_87 - tmp_86*tmp_89 - tmp_87*tmp_97);
      real_t tmp_99 = 1.0*tmp_90;
      real_t tmp_100 = tmp_48*(3.0*tmp_35*tmp_39*tmp_7 - tmp_35*tmp_94 - tmp_39*tmp_91) + tmp_56*(3.0*tmp_52*tmp_53*tmp_7 - tmp_52*tmp_94 - tmp_53*tmp_91) + tmp_64*(3.0*tmp_60*tmp_61*tmp_7 - tmp_60*tmp_94 - tmp_61*tmp_91) + tmp_72*(3.0*tmp_68*tmp_69*tmp_7 - tmp_68*tmp_94 - tmp_69*tmp_91) + tmp_80*(3.0*tmp_7*tmp_76*tmp_77 - tmp_76*tmp_94 - tmp_77*tmp_91) + tmp_88*(3.0*tmp_7*tmp_84*tmp_85 - tmp_84*tmp_94 - tmp_85*tmp_91);
      real_t tmp_101 = tmp_48*(3.0*tmp_35*tmp_43*tmp_7 - tmp_35*tmp_97 - tmp_43*tmp_91) + tmp_56*(3.0*tmp_52*tmp_54*tmp_7 - tmp_52*tmp_97 - tmp_54*tmp_91) + tmp_64*(3.0*tmp_60*tmp_62*tmp_7 - tmp_60*tmp_97 - tmp_62*tmp_91) + tmp_72*(3.0*tmp_68*tmp_7*tmp_70 - tmp_68*tmp_97 - tmp_70*tmp_91) + tmp_80*(3.0*tmp_7*tmp_76*tmp_78 - tmp_76*tmp_97 - tmp_78*tmp_91) + tmp_88*(3.0*tmp_7*tmp_84*tmp_86 - tmp_84*tmp_97 - tmp_86*tmp_91);
      real_t tmp_102 = 1.0*tmp_93;
      real_t tmp_103 = tmp_48*(3.0*tmp_39*tmp_43*tmp_7 - tmp_39*tmp_97 - tmp_43*tmp_94) + tmp_56*(3.0*tmp_53*tmp_54*tmp_7 - tmp_53*tmp_97 - tmp_54*tmp_94) + tmp_64*(3.0*tmp_61*tmp_62*tmp_7 - tmp_61*tmp_97 - tmp_62*tmp_94) + tmp_72*(3.0*tmp_69*tmp_7*tmp_70 - tmp_69*tmp_97 - tmp_70*tmp_94) + tmp_80*(3.0*tmp_7*tmp_77*tmp_78 - tmp_77*tmp_97 - tmp_78*tmp_94) + tmp_88*(3.0*tmp_7*tmp_85*tmp_86 - tmp_85*tmp_97 - tmp_86*tmp_94);
      real_t tmp_104 = 1.0*tmp_96;
      real_t a_0_0 = tmp_48*(3.0*(tmp_44*tmp_44)*tmp_7 - tmp_44*tmp_46) + tmp_56*(-tmp_46*tmp_55 + 3.0*(tmp_55*tmp_55)*tmp_7) + tmp_64*(-tmp_46*tmp_63 + 3.0*(tmp_63*tmp_63)*tmp_7) + tmp_72*(-tmp_46*tmp_71 + 3.0*tmp_7*(tmp_71*tmp_71)) + tmp_80*(-tmp_46*tmp_79 + 3.0*tmp_7*(tmp_79*tmp_79)) + tmp_88*(-tmp_46*tmp_87 + 3.0*tmp_7*(tmp_87*tmp_87));
      real_t a_0_1 = tmp_92;
      real_t a_0_2 = tmp_95;
      real_t a_0_3 = tmp_98;
      real_t a_1_0 = tmp_92;
      real_t a_1_1 = tmp_48*(3.0*(tmp_35*tmp_35)*tmp_7 - tmp_35*tmp_99) + tmp_56*(3.0*(tmp_52*tmp_52)*tmp_7 - tmp_52*tmp_99) + tmp_64*(3.0*(tmp_60*tmp_60)*tmp_7 - tmp_60*tmp_99) + tmp_72*(3.0*(tmp_68*tmp_68)*tmp_7 - tmp_68*tmp_99) + tmp_80*(3.0*tmp_7*(tmp_76*tmp_76) - tmp_76*tmp_99) + tmp_88*(3.0*tmp_7*(tmp_84*tmp_84) - tmp_84*tmp_99);
      real_t a_1_2 = tmp_100;
      real_t a_1_3 = tmp_101;
      real_t a_2_0 = tmp_95;
      real_t a_2_1 = tmp_100;
      real_t a_2_2 = tmp_48*(-tmp_102*tmp_39 + 3.0*(tmp_39*tmp_39)*tmp_7) + tmp_56*(-tmp_102*tmp_53 + 3.0*(tmp_53*tmp_53)*tmp_7) + tmp_64*(-tmp_102*tmp_61 + 3.0*(tmp_61*tmp_61)*tmp_7) + tmp_72*(-tmp_102*tmp_69 + 3.0*(tmp_69*tmp_69)*tmp_7) + tmp_80*(-tmp_102*tmp_77 + 3.0*tmp_7*(tmp_77*tmp_77)) + tmp_88*(-tmp_102*tmp_85 + 3.0*tmp_7*(tmp_85*tmp_85));
      real_t a_2_3 = tmp_103;
      real_t a_3_0 = tmp_98;
      real_t a_3_1 = tmp_101;
      real_t a_3_2 = tmp_103;
      real_t a_3_3 = tmp_48*(-tmp_104*tmp_43 + 3.0*(tmp_43*tmp_43)*tmp_7) + tmp_56*(-tmp_104*tmp_54 + 3.0*(tmp_54*tmp_54)*tmp_7) + tmp_64*(-tmp_104*tmp_62 + 3.0*(tmp_62*tmp_62)*tmp_7) + tmp_72*(-tmp_104*tmp_70 + 3.0*tmp_7*(tmp_70*tmp_70)) + tmp_80*(-tmp_104*tmp_78 + 3.0*tmp_7*(tmp_78*tmp_78)) + tmp_88*(-tmp_104*tmp_86 + 3.0*tmp_7*(tmp_86*tmp_86));
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
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
                                  MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_26 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_27 = -tmp_26;
      real_t tmp_28 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_29 = 0.091576213509770743*tmp_27 + 0.81684757298045851*tmp_28;
      real_t tmp_30 = tmp_25 + tmp_29;
      real_t tmp_31 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_32 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_33 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_34 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_35 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_36 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_39 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_40 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_41 = tmp_39*tmp_40;
      real_t tmp_42 = tmp_32*tmp_36;
      real_t tmp_43 = tmp_31*tmp_39;
      real_t tmp_44 = tmp_33*tmp_38;
      real_t tmp_45 = 1.0 / (tmp_31*tmp_32*tmp_38 + tmp_33*tmp_41 - tmp_34*tmp_44 + tmp_35*tmp_37 - tmp_35*tmp_43 - tmp_40*tmp_42);
      real_t tmp_46 = tmp_45*(tmp_31*tmp_32 - tmp_33*tmp_34);
      real_t tmp_47 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_48 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_49 = -tmp_48;
      real_t tmp_50 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_51 = 0.091576213509770743*tmp_49 + 0.81684757298045851*tmp_50;
      real_t tmp_52 = tmp_47 + tmp_51;
      real_t tmp_53 = tmp_45*(-tmp_31*tmp_35 + tmp_33*tmp_40);
      real_t tmp_54 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_55 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_56 = -tmp_55;
      real_t tmp_57 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_58 = 0.091576213509770743*tmp_56 + 0.81684757298045851*tmp_57;
      real_t tmp_59 = tmp_54 + tmp_58;
      real_t tmp_60 = tmp_45*(-tmp_32*tmp_40 + tmp_34*tmp_35);
      real_t tmp_61 = tmp_30*tmp_46 + tmp_52*tmp_53 + tmp_59*tmp_60;
      real_t tmp_62 = tmp_45*(tmp_33*tmp_39 - tmp_42);
      real_t tmp_63 = tmp_45*(tmp_35*tmp_36 - tmp_44);
      real_t tmp_64 = tmp_45*(tmp_32*tmp_38 - tmp_35*tmp_39);
      real_t tmp_65 = tmp_30*tmp_62 + tmp_52*tmp_63 + tmp_59*tmp_64;
      real_t tmp_66 = tmp_45*(tmp_37 - tmp_43);
      real_t tmp_67 = tmp_45*(tmp_31*tmp_38 - tmp_36*tmp_40);
      real_t tmp_68 = tmp_45*(-tmp_34*tmp_38 + tmp_41);
      real_t tmp_69 = tmp_30*tmp_66 + tmp_52*tmp_67 + tmp_59*tmp_68;
      real_t tmp_70 = -tmp_61 - tmp_65 - tmp_69 + 1;
      real_t tmp_71 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_72 = tmp_29 + tmp_71;
      real_t tmp_73 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_74 = tmp_51 + tmp_73;
      real_t tmp_75 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_76 = tmp_58 + tmp_75;
      real_t tmp_77 = tmp_17*tmp_76 + tmp_20*tmp_74 + tmp_23*tmp_72;
      real_t tmp_78 = tmp_16*tmp_76 + tmp_19*tmp_74 + tmp_22*tmp_72;
      real_t tmp_79 = tmp_15*tmp_76 + tmp_18*tmp_74 + tmp_21*tmp_72;
      real_t tmp_80 = -tmp_77 - tmp_78 - tmp_79 + 1;
      real_t tmp_81 = 0.5*p_affine_13_0*(-tmp_60 - tmp_64 - tmp_68) + 0.5*p_affine_13_1*(-tmp_53 - tmp_63 - tmp_67) + 0.5*p_affine_13_2*(-tmp_46 - tmp_62 - tmp_66);
      real_t tmp_82 = (std::abs(tmp_26*tmp_50 - tmp_28*tmp_48)*std::abs(tmp_26*tmp_50 - tmp_28*tmp_48)) + (std::abs(tmp_26*tmp_57 - tmp_28*tmp_55)*std::abs(tmp_26*tmp_57 - tmp_28*tmp_55)) + (std::abs(tmp_48*tmp_57 - tmp_50*tmp_55)*std::abs(tmp_48*tmp_57 - tmp_50*tmp_55));
      real_t tmp_83 = 3.0*std::pow(tmp_82, -0.25);
      real_t tmp_84 = tmp_80*tmp_83;
      real_t tmp_85 = 1.0*std::pow(tmp_82, 1.0/2.0);
      real_t tmp_86 = 0.054975871827660928*tmp_85;
      real_t tmp_87 = 0.44594849091596489*tmp_27 + 0.10810301816807022*tmp_28;
      real_t tmp_88 = tmp_25 + tmp_87;
      real_t tmp_89 = 0.44594849091596489*tmp_49 + 0.10810301816807022*tmp_50;
      real_t tmp_90 = tmp_47 + tmp_89;
      real_t tmp_91 = 0.44594849091596489*tmp_56 + 0.10810301816807022*tmp_57;
      real_t tmp_92 = tmp_54 + tmp_91;
      real_t tmp_93 = tmp_46*tmp_88 + tmp_53*tmp_90 + tmp_60*tmp_92;
      real_t tmp_94 = tmp_62*tmp_88 + tmp_63*tmp_90 + tmp_64*tmp_92;
      real_t tmp_95 = tmp_66*tmp_88 + tmp_67*tmp_90 + tmp_68*tmp_92;
      real_t tmp_96 = -tmp_93 - tmp_94 - tmp_95 + 1;
      real_t tmp_97 = tmp_71 + tmp_87;
      real_t tmp_98 = tmp_73 + tmp_89;
      real_t tmp_99 = tmp_75 + tmp_91;
      real_t tmp_100 = tmp_17*tmp_99 + tmp_20*tmp_98 + tmp_23*tmp_97;
      real_t tmp_101 = tmp_16*tmp_99 + tmp_19*tmp_98 + tmp_22*tmp_97;
      real_t tmp_102 = tmp_15*tmp_99 + tmp_18*tmp_98 + tmp_21*tmp_97;
      real_t tmp_103 = -tmp_100 - tmp_101 - tmp_102 + 1;
      real_t tmp_104 = tmp_103*tmp_83;
      real_t tmp_105 = 0.11169079483900572*tmp_85;
      real_t tmp_106 = 0.81684757298045851*tmp_27 + 0.091576213509770743*tmp_28;
      real_t tmp_107 = tmp_106 + tmp_25;
      real_t tmp_108 = 0.81684757298045851*tmp_49 + 0.091576213509770743*tmp_50;
      real_t tmp_109 = tmp_108 + tmp_47;
      real_t tmp_110 = 0.81684757298045851*tmp_56 + 0.091576213509770743*tmp_57;
      real_t tmp_111 = tmp_110 + tmp_54;
      real_t tmp_112 = tmp_107*tmp_46 + tmp_109*tmp_53 + tmp_111*tmp_60;
      real_t tmp_113 = tmp_107*tmp_62 + tmp_109*tmp_63 + tmp_111*tmp_64;
      real_t tmp_114 = tmp_107*tmp_66 + tmp_109*tmp_67 + tmp_111*tmp_68;
      real_t tmp_115 = -tmp_112 - tmp_113 - tmp_114 + 1;
      real_t tmp_116 = tmp_106 + tmp_71;
      real_t tmp_117 = tmp_108 + tmp_73;
      real_t tmp_118 = tmp_110 + tmp_75;
      real_t tmp_119 = tmp_116*tmp_23 + tmp_117*tmp_20 + tmp_118*tmp_17;
      real_t tmp_120 = tmp_116*tmp_22 + tmp_117*tmp_19 + tmp_118*tmp_16;
      real_t tmp_121 = tmp_116*tmp_21 + tmp_117*tmp_18 + tmp_118*tmp_15;
      real_t tmp_122 = -tmp_119 - tmp_120 - tmp_121 + 1;
      real_t tmp_123 = tmp_122*tmp_83;
      real_t tmp_124 = 0.054975871827660928*tmp_85;
      real_t tmp_125 = 0.10810301816807022*tmp_27 + 0.44594849091596489*tmp_28;
      real_t tmp_126 = tmp_125 + tmp_25;
      real_t tmp_127 = 0.10810301816807022*tmp_49 + 0.44594849091596489*tmp_50;
      real_t tmp_128 = tmp_127 + tmp_47;
      real_t tmp_129 = 0.10810301816807022*tmp_56 + 0.44594849091596489*tmp_57;
      real_t tmp_130 = tmp_129 + tmp_54;
      real_t tmp_131 = tmp_126*tmp_46 + tmp_128*tmp_53 + tmp_130*tmp_60;
      real_t tmp_132 = tmp_126*tmp_62 + tmp_128*tmp_63 + tmp_130*tmp_64;
      real_t tmp_133 = tmp_126*tmp_66 + tmp_128*tmp_67 + tmp_130*tmp_68;
      real_t tmp_134 = -tmp_131 - tmp_132 - tmp_133 + 1;
      real_t tmp_135 = tmp_125 + tmp_71;
      real_t tmp_136 = tmp_127 + tmp_73;
      real_t tmp_137 = tmp_129 + tmp_75;
      real_t tmp_138 = tmp_135*tmp_23 + tmp_136*tmp_20 + tmp_137*tmp_17;
      real_t tmp_139 = tmp_135*tmp_22 + tmp_136*tmp_19 + tmp_137*tmp_16;
      real_t tmp_140 = tmp_135*tmp_21 + tmp_136*tmp_18 + tmp_137*tmp_15;
      real_t tmp_141 = -tmp_138 - tmp_139 - tmp_140 + 1;
      real_t tmp_142 = tmp_141*tmp_83;
      real_t tmp_143 = 0.11169079483900572*tmp_85;
      real_t tmp_144 = 0.091576213509770743*tmp_27 + 0.091576213509770743*tmp_28;
      real_t tmp_145 = tmp_144 + tmp_25;
      real_t tmp_146 = 0.091576213509770743*tmp_49 + 0.091576213509770743*tmp_50;
      real_t tmp_147 = tmp_146 + tmp_47;
      real_t tmp_148 = 0.091576213509770743*tmp_56 + 0.091576213509770743*tmp_57;
      real_t tmp_149 = tmp_148 + tmp_54;
      real_t tmp_150 = tmp_145*tmp_46 + tmp_147*tmp_53 + tmp_149*tmp_60;
      real_t tmp_151 = tmp_145*tmp_62 + tmp_147*tmp_63 + tmp_149*tmp_64;
      real_t tmp_152 = tmp_145*tmp_66 + tmp_147*tmp_67 + tmp_149*tmp_68;
      real_t tmp_153 = -tmp_150 - tmp_151 - tmp_152 + 1;
      real_t tmp_154 = tmp_144 + tmp_71;
      real_t tmp_155 = tmp_146 + tmp_73;
      real_t tmp_156 = tmp_148 + tmp_75;
      real_t tmp_157 = tmp_154*tmp_23 + tmp_155*tmp_20 + tmp_156*tmp_17;
      real_t tmp_158 = tmp_154*tmp_22 + tmp_155*tmp_19 + tmp_156*tmp_16;
      real_t tmp_159 = tmp_15*tmp_156 + tmp_154*tmp_21 + tmp_155*tmp_18;
      real_t tmp_160 = -tmp_157 - tmp_158 - tmp_159 + 1;
      real_t tmp_161 = tmp_160*tmp_83;
      real_t tmp_162 = 0.054975871827660928*tmp_85;
      real_t tmp_163 = 0.44594849091596489*tmp_27 + 0.44594849091596489*tmp_28;
      real_t tmp_164 = tmp_163 + tmp_25;
      real_t tmp_165 = 0.44594849091596489*tmp_49 + 0.44594849091596489*tmp_50;
      real_t tmp_166 = tmp_165 + tmp_47;
      real_t tmp_167 = 0.44594849091596489*tmp_56 + 0.44594849091596489*tmp_57;
      real_t tmp_168 = tmp_167 + tmp_54;
      real_t tmp_169 = tmp_164*tmp_46 + tmp_166*tmp_53 + tmp_168*tmp_60;
      real_t tmp_170 = tmp_164*tmp_62 + tmp_166*tmp_63 + tmp_168*tmp_64;
      real_t tmp_171 = tmp_164*tmp_66 + tmp_166*tmp_67 + tmp_168*tmp_68;
      real_t tmp_172 = -tmp_169 - tmp_170 - tmp_171 + 1;
      real_t tmp_173 = tmp_163 + tmp_71;
      real_t tmp_174 = tmp_165 + tmp_73;
      real_t tmp_175 = tmp_167 + tmp_75;
      real_t tmp_176 = tmp_17*tmp_175 + tmp_173*tmp_23 + tmp_174*tmp_20;
      real_t tmp_177 = tmp_16*tmp_175 + tmp_173*tmp_22 + tmp_174*tmp_19;
      real_t tmp_178 = tmp_15*tmp_175 + tmp_173*tmp_21 + tmp_174*tmp_18;
      real_t tmp_179 = -tmp_176 - tmp_177 - tmp_178 + 1;
      real_t tmp_180 = tmp_179*tmp_83;
      real_t tmp_181 = 0.11169079483900572*tmp_85;
      real_t tmp_182 = 0.5*p_affine_13_0*tmp_60 + 0.5*p_affine_13_1*tmp_53 + 0.5*p_affine_13_2*tmp_46;
      real_t tmp_183 = 0.5*p_affine_13_0*tmp_64 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_62;
      real_t tmp_184 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_67 + 0.5*p_affine_13_2*tmp_66;
      real_t tmp_185 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_186 = tmp_77*tmp_83;
      real_t tmp_187 = tmp_100*tmp_83;
      real_t tmp_188 = tmp_119*tmp_83;
      real_t tmp_189 = tmp_138*tmp_83;
      real_t tmp_190 = tmp_157*tmp_83;
      real_t tmp_191 = tmp_176*tmp_83;
      real_t tmp_192 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_193 = tmp_78*tmp_83;
      real_t tmp_194 = tmp_101*tmp_83;
      real_t tmp_195 = tmp_120*tmp_83;
      real_t tmp_196 = tmp_139*tmp_83;
      real_t tmp_197 = tmp_158*tmp_83;
      real_t tmp_198 = tmp_177*tmp_83;
      real_t tmp_199 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t tmp_200 = tmp_79*tmp_83;
      real_t tmp_201 = tmp_102*tmp_83;
      real_t tmp_202 = tmp_121*tmp_83;
      real_t tmp_203 = tmp_140*tmp_83;
      real_t tmp_204 = tmp_159*tmp_83;
      real_t tmp_205 = tmp_178*tmp_83;
      real_t a_0_0 = tmp_105*(-tmp_103*tmp_81 - tmp_104*tmp_96 + 0.5*tmp_24*tmp_96) + tmp_124*(-tmp_115*tmp_123 + 0.5*tmp_115*tmp_24 - tmp_122*tmp_81) + tmp_143*(-tmp_134*tmp_142 + 0.5*tmp_134*tmp_24 - tmp_141*tmp_81) + tmp_162*(-tmp_153*tmp_161 + 0.5*tmp_153*tmp_24 - tmp_160*tmp_81) + tmp_181*(-tmp_172*tmp_180 + 0.5*tmp_172*tmp_24 - tmp_179*tmp_81) + tmp_86*(0.5*tmp_24*tmp_70 - tmp_70*tmp_84 - tmp_80*tmp_81);
      real_t a_0_1 = tmp_105*(-tmp_103*tmp_182 - tmp_104*tmp_93 + 0.5*tmp_24*tmp_93) + tmp_124*(-tmp_112*tmp_123 + 0.5*tmp_112*tmp_24 - tmp_122*tmp_182) + tmp_143*(-tmp_131*tmp_142 + 0.5*tmp_131*tmp_24 - tmp_141*tmp_182) + tmp_162*(-tmp_150*tmp_161 + 0.5*tmp_150*tmp_24 - tmp_160*tmp_182) + tmp_181*(-tmp_169*tmp_180 + 0.5*tmp_169*tmp_24 - tmp_179*tmp_182) + tmp_86*(-tmp_182*tmp_80 + 0.5*tmp_24*tmp_61 - tmp_61*tmp_84);
      real_t a_0_2 = tmp_105*(-tmp_103*tmp_183 - tmp_104*tmp_94 + 0.5*tmp_24*tmp_94) + tmp_124*(-tmp_113*tmp_123 + 0.5*tmp_113*tmp_24 - tmp_122*tmp_183) + tmp_143*(-tmp_132*tmp_142 + 0.5*tmp_132*tmp_24 - tmp_141*tmp_183) + tmp_162*(-tmp_151*tmp_161 + 0.5*tmp_151*tmp_24 - tmp_160*tmp_183) + tmp_181*(-tmp_170*tmp_180 + 0.5*tmp_170*tmp_24 - tmp_179*tmp_183) + tmp_86*(-tmp_183*tmp_80 + 0.5*tmp_24*tmp_65 - tmp_65*tmp_84);
      real_t a_0_3 = tmp_105*(-tmp_103*tmp_184 - tmp_104*tmp_95 + 0.5*tmp_24*tmp_95) + tmp_124*(-tmp_114*tmp_123 + 0.5*tmp_114*tmp_24 - tmp_122*tmp_184) + tmp_143*(-tmp_133*tmp_142 + 0.5*tmp_133*tmp_24 - tmp_141*tmp_184) + tmp_162*(-tmp_152*tmp_161 + 0.5*tmp_152*tmp_24 - tmp_160*tmp_184) + tmp_181*(-tmp_171*tmp_180 + 0.5*tmp_171*tmp_24 - tmp_179*tmp_184) + tmp_86*(-tmp_184*tmp_80 + 0.5*tmp_24*tmp_69 - tmp_69*tmp_84);
      real_t a_1_0 = tmp_105*(-tmp_100*tmp_81 + 0.5*tmp_185*tmp_96 - tmp_187*tmp_96) + tmp_124*(0.5*tmp_115*tmp_185 - tmp_115*tmp_188 - tmp_119*tmp_81) + tmp_143*(0.5*tmp_134*tmp_185 - tmp_134*tmp_189 - tmp_138*tmp_81) + tmp_162*(0.5*tmp_153*tmp_185 - tmp_153*tmp_190 - tmp_157*tmp_81) + tmp_181*(0.5*tmp_172*tmp_185 - tmp_172*tmp_191 - tmp_176*tmp_81) + tmp_86*(0.5*tmp_185*tmp_70 - tmp_186*tmp_70 - tmp_77*tmp_81);
      real_t a_1_1 = tmp_105*(-tmp_100*tmp_182 + 0.5*tmp_185*tmp_93 - tmp_187*tmp_93) + tmp_124*(0.5*tmp_112*tmp_185 - tmp_112*tmp_188 - tmp_119*tmp_182) + tmp_143*(0.5*tmp_131*tmp_185 - tmp_131*tmp_189 - tmp_138*tmp_182) + tmp_162*(0.5*tmp_150*tmp_185 - tmp_150*tmp_190 - tmp_157*tmp_182) + tmp_181*(0.5*tmp_169*tmp_185 - tmp_169*tmp_191 - tmp_176*tmp_182) + tmp_86*(-tmp_182*tmp_77 + 0.5*tmp_185*tmp_61 - tmp_186*tmp_61);
      real_t a_1_2 = tmp_105*(-tmp_100*tmp_183 + 0.5*tmp_185*tmp_94 - tmp_187*tmp_94) + tmp_124*(0.5*tmp_113*tmp_185 - tmp_113*tmp_188 - tmp_119*tmp_183) + tmp_143*(0.5*tmp_132*tmp_185 - tmp_132*tmp_189 - tmp_138*tmp_183) + tmp_162*(0.5*tmp_151*tmp_185 - tmp_151*tmp_190 - tmp_157*tmp_183) + tmp_181*(0.5*tmp_170*tmp_185 - tmp_170*tmp_191 - tmp_176*tmp_183) + tmp_86*(-tmp_183*tmp_77 + 0.5*tmp_185*tmp_65 - tmp_186*tmp_65);
      real_t a_1_3 = tmp_105*(-tmp_100*tmp_184 + 0.5*tmp_185*tmp_95 - tmp_187*tmp_95) + tmp_124*(0.5*tmp_114*tmp_185 - tmp_114*tmp_188 - tmp_119*tmp_184) + tmp_143*(0.5*tmp_133*tmp_185 - tmp_133*tmp_189 - tmp_138*tmp_184) + tmp_162*(0.5*tmp_152*tmp_185 - tmp_152*tmp_190 - tmp_157*tmp_184) + tmp_181*(0.5*tmp_171*tmp_185 - tmp_171*tmp_191 - tmp_176*tmp_184) + tmp_86*(-tmp_184*tmp_77 + 0.5*tmp_185*tmp_69 - tmp_186*tmp_69);
      real_t a_2_0 = tmp_105*(-tmp_101*tmp_81 + 0.5*tmp_192*tmp_96 - tmp_194*tmp_96) + tmp_124*(0.5*tmp_115*tmp_192 - tmp_115*tmp_195 - tmp_120*tmp_81) + tmp_143*(0.5*tmp_134*tmp_192 - tmp_134*tmp_196 - tmp_139*tmp_81) + tmp_162*(0.5*tmp_153*tmp_192 - tmp_153*tmp_197 - tmp_158*tmp_81) + tmp_181*(0.5*tmp_172*tmp_192 - tmp_172*tmp_198 - tmp_177*tmp_81) + tmp_86*(0.5*tmp_192*tmp_70 - tmp_193*tmp_70 - tmp_78*tmp_81);
      real_t a_2_1 = tmp_105*(-tmp_101*tmp_182 + 0.5*tmp_192*tmp_93 - tmp_194*tmp_93) + tmp_124*(0.5*tmp_112*tmp_192 - tmp_112*tmp_195 - tmp_120*tmp_182) + tmp_143*(0.5*tmp_131*tmp_192 - tmp_131*tmp_196 - tmp_139*tmp_182) + tmp_162*(0.5*tmp_150*tmp_192 - tmp_150*tmp_197 - tmp_158*tmp_182) + tmp_181*(0.5*tmp_169*tmp_192 - tmp_169*tmp_198 - tmp_177*tmp_182) + tmp_86*(-tmp_182*tmp_78 + 0.5*tmp_192*tmp_61 - tmp_193*tmp_61);
      real_t a_2_2 = tmp_105*(-tmp_101*tmp_183 + 0.5*tmp_192*tmp_94 - tmp_194*tmp_94) + tmp_124*(0.5*tmp_113*tmp_192 - tmp_113*tmp_195 - tmp_120*tmp_183) + tmp_143*(0.5*tmp_132*tmp_192 - tmp_132*tmp_196 - tmp_139*tmp_183) + tmp_162*(0.5*tmp_151*tmp_192 - tmp_151*tmp_197 - tmp_158*tmp_183) + tmp_181*(0.5*tmp_170*tmp_192 - tmp_170*tmp_198 - tmp_177*tmp_183) + tmp_86*(-tmp_183*tmp_78 + 0.5*tmp_192*tmp_65 - tmp_193*tmp_65);
      real_t a_2_3 = tmp_105*(-tmp_101*tmp_184 + 0.5*tmp_192*tmp_95 - tmp_194*tmp_95) + tmp_124*(0.5*tmp_114*tmp_192 - tmp_114*tmp_195 - tmp_120*tmp_184) + tmp_143*(0.5*tmp_133*tmp_192 - tmp_133*tmp_196 - tmp_139*tmp_184) + tmp_162*(0.5*tmp_152*tmp_192 - tmp_152*tmp_197 - tmp_158*tmp_184) + tmp_181*(0.5*tmp_171*tmp_192 - tmp_171*tmp_198 - tmp_177*tmp_184) + tmp_86*(-tmp_184*tmp_78 + 0.5*tmp_192*tmp_69 - tmp_193*tmp_69);
      real_t a_3_0 = tmp_105*(-tmp_102*tmp_81 + 0.5*tmp_199*tmp_96 - tmp_201*tmp_96) + tmp_124*(0.5*tmp_115*tmp_199 - tmp_115*tmp_202 - tmp_121*tmp_81) + tmp_143*(0.5*tmp_134*tmp_199 - tmp_134*tmp_203 - tmp_140*tmp_81) + tmp_162*(0.5*tmp_153*tmp_199 - tmp_153*tmp_204 - tmp_159*tmp_81) + tmp_181*(0.5*tmp_172*tmp_199 - tmp_172*tmp_205 - tmp_178*tmp_81) + tmp_86*(0.5*tmp_199*tmp_70 - tmp_200*tmp_70 - tmp_79*tmp_81);
      real_t a_3_1 = tmp_105*(-tmp_102*tmp_182 + 0.5*tmp_199*tmp_93 - tmp_201*tmp_93) + tmp_124*(0.5*tmp_112*tmp_199 - tmp_112*tmp_202 - tmp_121*tmp_182) + tmp_143*(0.5*tmp_131*tmp_199 - tmp_131*tmp_203 - tmp_140*tmp_182) + tmp_162*(0.5*tmp_150*tmp_199 - tmp_150*tmp_204 - tmp_159*tmp_182) + tmp_181*(0.5*tmp_169*tmp_199 - tmp_169*tmp_205 - tmp_178*tmp_182) + tmp_86*(-tmp_182*tmp_79 + 0.5*tmp_199*tmp_61 - tmp_200*tmp_61);
      real_t a_3_2 = tmp_105*(-tmp_102*tmp_183 + 0.5*tmp_199*tmp_94 - tmp_201*tmp_94) + tmp_124*(0.5*tmp_113*tmp_199 - tmp_113*tmp_202 - tmp_121*tmp_183) + tmp_143*(0.5*tmp_132*tmp_199 - tmp_132*tmp_203 - tmp_140*tmp_183) + tmp_162*(0.5*tmp_151*tmp_199 - tmp_151*tmp_204 - tmp_159*tmp_183) + tmp_181*(0.5*tmp_170*tmp_199 - tmp_170*tmp_205 - tmp_178*tmp_183) + tmp_86*(-tmp_183*tmp_79 + 0.5*tmp_199*tmp_65 - tmp_200*tmp_65);
      real_t a_3_3 = tmp_105*(-tmp_102*tmp_184 + 0.5*tmp_199*tmp_95 - tmp_201*tmp_95) + tmp_124*(0.5*tmp_114*tmp_199 - tmp_114*tmp_202 - tmp_121*tmp_184) + tmp_143*(0.5*tmp_133*tmp_199 - tmp_133*tmp_203 - tmp_140*tmp_184) + tmp_162*(0.5*tmp_152*tmp_199 - tmp_152*tmp_204 - tmp_159*tmp_184) + tmp_181*(0.5*tmp_171*tmp_199 - tmp_171*tmp_205 - tmp_178*tmp_184) + tmp_86*(-tmp_184*tmp_79 + 0.5*tmp_199*tmp_69 - tmp_200*tmp_69);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
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
                                        MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = (std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4));
      real_t tmp_7 = std::pow(tmp_6, -0.25);
      real_t tmp_8 = -tmp_4;
      real_t tmp_9 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_10 = 0.81684757298045851*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_11 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_12 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_13 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_14 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_15 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_16 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_19 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_20 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_12*tmp_20;
      real_t tmp_23 = tmp_16*tmp_19;
      real_t tmp_24 = tmp_14*tmp_18;
      real_t tmp_25 = 1.0 / (tmp_11*tmp_12*tmp_18 - tmp_11*tmp_23 + tmp_13*tmp_21 - tmp_13*tmp_24 + tmp_15*tmp_17 - tmp_15*tmp_22);
      real_t tmp_26 = tmp_25*(tmp_11*tmp_12 - tmp_13*tmp_14);
      real_t tmp_27 = -tmp_1;
      real_t tmp_28 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_29 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_30 = tmp_25*(-tmp_11*tmp_16 + tmp_13*tmp_20);
      real_t tmp_31 = -tmp_3;
      real_t tmp_32 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_33 = 0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_34 = tmp_25*(tmp_17 - tmp_22);
      real_t tmp_35 = tmp_10*tmp_26 + tmp_29*tmp_30 + tmp_33*tmp_34;
      real_t tmp_36 = tmp_25*(-tmp_12*tmp_15 + tmp_13*tmp_19);
      real_t tmp_37 = tmp_25*(-tmp_13*tmp_18 + tmp_15*tmp_16);
      real_t tmp_38 = tmp_25*(tmp_12*tmp_18 - tmp_23);
      real_t tmp_39 = tmp_10*tmp_36 + tmp_29*tmp_37 + tmp_33*tmp_38;
      real_t tmp_40 = tmp_25*(-tmp_11*tmp_19 + tmp_14*tmp_15);
      real_t tmp_41 = tmp_25*(tmp_11*tmp_18 - tmp_15*tmp_20);
      real_t tmp_42 = tmp_25*(tmp_21 - tmp_24);
      real_t tmp_43 = tmp_10*tmp_40 + tmp_29*tmp_41 + tmp_33*tmp_42;
      real_t tmp_44 = -tmp_35 - tmp_39 - tmp_43 + 1;
      real_t tmp_45 = p_affine_13_0*(-tmp_34 - tmp_38 - tmp_42) + p_affine_13_1*(-tmp_30 - tmp_37 - tmp_41) + p_affine_13_2*(-tmp_26 - tmp_36 - tmp_40);
      real_t tmp_46 = 2*tmp_45;
      real_t tmp_47 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_48 = 0.054975871827660928*tmp_47;
      real_t tmp_49 = 0.10810301816807022*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_50 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_51 = 0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_52 = tmp_26*tmp_49 + tmp_30*tmp_50 + tmp_34*tmp_51;
      real_t tmp_53 = tmp_36*tmp_49 + tmp_37*tmp_50 + tmp_38*tmp_51;
      real_t tmp_54 = tmp_40*tmp_49 + tmp_41*tmp_50 + tmp_42*tmp_51;
      real_t tmp_55 = -tmp_52 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11169079483900572*tmp_47;
      real_t tmp_57 = 0.091576213509770743*tmp_5 + 0.81684757298045851*tmp_8 + tmp_9;
      real_t tmp_58 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_59 = 0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_31 + tmp_32;
      real_t tmp_60 = tmp_26*tmp_57 + tmp_30*tmp_58 + tmp_34*tmp_59;
      real_t tmp_61 = tmp_36*tmp_57 + tmp_37*tmp_58 + tmp_38*tmp_59;
      real_t tmp_62 = tmp_40*tmp_57 + tmp_41*tmp_58 + tmp_42*tmp_59;
      real_t tmp_63 = -tmp_60 - tmp_61 - tmp_62 + 1;
      real_t tmp_64 = 0.054975871827660928*tmp_47;
      real_t tmp_65 = 0.44594849091596489*tmp_5 + 0.10810301816807022*tmp_8 + tmp_9;
      real_t tmp_66 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_67 = 0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_31 + tmp_32;
      real_t tmp_68 = tmp_26*tmp_65 + tmp_30*tmp_66 + tmp_34*tmp_67;
      real_t tmp_69 = tmp_36*tmp_65 + tmp_37*tmp_66 + tmp_38*tmp_67;
      real_t tmp_70 = tmp_40*tmp_65 + tmp_41*tmp_66 + tmp_42*tmp_67;
      real_t tmp_71 = -tmp_68 - tmp_69 - tmp_70 + 1;
      real_t tmp_72 = 0.11169079483900572*tmp_47;
      real_t tmp_73 = 0.091576213509770743*tmp_5 + 0.091576213509770743*tmp_8 + tmp_9;
      real_t tmp_74 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_75 = 0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_31 + tmp_32;
      real_t tmp_76 = tmp_26*tmp_73 + tmp_30*tmp_74 + tmp_34*tmp_75;
      real_t tmp_77 = tmp_36*tmp_73 + tmp_37*tmp_74 + tmp_38*tmp_75;
      real_t tmp_78 = tmp_40*tmp_73 + tmp_41*tmp_74 + tmp_42*tmp_75;
      real_t tmp_79 = -tmp_76 - tmp_77 - tmp_78 + 1;
      real_t tmp_80 = 0.054975871827660928*tmp_47;
      real_t tmp_81 = 0.44594849091596489*tmp_5 + 0.44594849091596489*tmp_8 + tmp_9;
      real_t tmp_82 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_83 = 0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_31 + tmp_32;
      real_t tmp_84 = tmp_26*tmp_81 + tmp_30*tmp_82 + tmp_34*tmp_83;
      real_t tmp_85 = tmp_36*tmp_81 + tmp_37*tmp_82 + tmp_38*tmp_83;
      real_t tmp_86 = tmp_40*tmp_81 + tmp_41*tmp_82 + tmp_42*tmp_83;
      real_t tmp_87 = -tmp_84 - tmp_85 - tmp_86 + 1;
      real_t tmp_88 = 0.11169079483900572*tmp_47;
      real_t tmp_89 = p_affine_13_0*tmp_34 + p_affine_13_1*tmp_30 + p_affine_13_2*tmp_26;
      real_t tmp_90 = tmp_48*(3.0*tmp_35*tmp_44*tmp_7 - tmp_35*tmp_45 - tmp_44*tmp_89) + tmp_56*(-tmp_45*tmp_52 + 3.0*tmp_52*tmp_55*tmp_7 - tmp_55*tmp_89) + tmp_64*(-tmp_45*tmp_60 + 3.0*tmp_60*tmp_63*tmp_7 - tmp_63*tmp_89) + tmp_72*(-tmp_45*tmp_68 + 3.0*tmp_68*tmp_7*tmp_71 - tmp_71*tmp_89) + tmp_80*(-tmp_45*tmp_76 + 3.0*tmp_7*tmp_76*tmp_79 - tmp_79*tmp_89) + tmp_88*(-tmp_45*tmp_84 + 3.0*tmp_7*tmp_84*tmp_87 - tmp_87*tmp_89);
      real_t tmp_91 = p_affine_13_0*tmp_38 + p_affine_13_1*tmp_37 + p_affine_13_2*tmp_36;
      real_t tmp_92 = tmp_48*(3.0*tmp_39*tmp_44*tmp_7 - tmp_39*tmp_45 - tmp_44*tmp_91) + tmp_56*(-tmp_45*tmp_53 + 3.0*tmp_53*tmp_55*tmp_7 - tmp_55*tmp_91) + tmp_64*(-tmp_45*tmp_61 + 3.0*tmp_61*tmp_63*tmp_7 - tmp_63*tmp_91) + tmp_72*(-tmp_45*tmp_69 + 3.0*tmp_69*tmp_7*tmp_71 - tmp_71*tmp_91) + tmp_80*(-tmp_45*tmp_77 + 3.0*tmp_7*tmp_77*tmp_79 - tmp_79*tmp_91) + tmp_88*(-tmp_45*tmp_85 + 3.0*tmp_7*tmp_85*tmp_87 - tmp_87*tmp_91);
      real_t tmp_93 = p_affine_13_0*tmp_42 + p_affine_13_1*tmp_41 + p_affine_13_2*tmp_40;
      real_t tmp_94 = tmp_48*(3.0*tmp_43*tmp_44*tmp_7 - tmp_43*tmp_45 - tmp_44*tmp_93) + tmp_56*(-tmp_45*tmp_54 + 3.0*tmp_54*tmp_55*tmp_7 - tmp_55*tmp_93) + tmp_64*(-tmp_45*tmp_62 + 3.0*tmp_62*tmp_63*tmp_7 - tmp_63*tmp_93) + tmp_72*(-tmp_45*tmp_70 + 3.0*tmp_7*tmp_70*tmp_71 - tmp_71*tmp_93) + tmp_80*(-tmp_45*tmp_78 + 3.0*tmp_7*tmp_78*tmp_79 - tmp_79*tmp_93) + tmp_88*(-tmp_45*tmp_86 + 3.0*tmp_7*tmp_86*tmp_87 - tmp_87*tmp_93);
      real_t tmp_95 = 2*tmp_89;
      real_t tmp_96 = tmp_48*(3.0*tmp_35*tmp_39*tmp_7 - tmp_35*tmp_91 - tmp_39*tmp_89) + tmp_56*(3.0*tmp_52*tmp_53*tmp_7 - tmp_52*tmp_91 - tmp_53*tmp_89) + tmp_64*(3.0*tmp_60*tmp_61*tmp_7 - tmp_60*tmp_91 - tmp_61*tmp_89) + tmp_72*(3.0*tmp_68*tmp_69*tmp_7 - tmp_68*tmp_91 - tmp_69*tmp_89) + tmp_80*(3.0*tmp_7*tmp_76*tmp_77 - tmp_76*tmp_91 - tmp_77*tmp_89) + tmp_88*(3.0*tmp_7*tmp_84*tmp_85 - tmp_84*tmp_91 - tmp_85*tmp_89);
      real_t tmp_97 = tmp_48*(3.0*tmp_35*tmp_43*tmp_7 - tmp_35*tmp_93 - tmp_43*tmp_89) + tmp_56*(3.0*tmp_52*tmp_54*tmp_7 - tmp_52*tmp_93 - tmp_54*tmp_89) + tmp_64*(3.0*tmp_60*tmp_62*tmp_7 - tmp_60*tmp_93 - tmp_62*tmp_89) + tmp_72*(3.0*tmp_68*tmp_7*tmp_70 - tmp_68*tmp_93 - tmp_70*tmp_89) + tmp_80*(3.0*tmp_7*tmp_76*tmp_78 - tmp_76*tmp_93 - tmp_78*tmp_89) + tmp_88*(3.0*tmp_7*tmp_84*tmp_86 - tmp_84*tmp_93 - tmp_86*tmp_89);
      real_t tmp_98 = 2*tmp_91;
      real_t tmp_99 = tmp_48*(3.0*tmp_39*tmp_43*tmp_7 - tmp_39*tmp_93 - tmp_43*tmp_91) + tmp_56*(3.0*tmp_53*tmp_54*tmp_7 - tmp_53*tmp_93 - tmp_54*tmp_91) + tmp_64*(3.0*tmp_61*tmp_62*tmp_7 - tmp_61*tmp_93 - tmp_62*tmp_91) + tmp_72*(3.0*tmp_69*tmp_7*tmp_70 - tmp_69*tmp_93 - tmp_70*tmp_91) + tmp_80*(3.0*tmp_7*tmp_77*tmp_78 - tmp_77*tmp_93 - tmp_78*tmp_91) + tmp_88*(3.0*tmp_7*tmp_85*tmp_86 - tmp_85*tmp_93 - tmp_86*tmp_91);
      real_t tmp_100 = 2*tmp_93;
      real_t a_0_0 = tmp_48*(3.0*(tmp_44*tmp_44)*tmp_7 - tmp_44*tmp_46) + tmp_56*(-tmp_46*tmp_55 + 3.0*(tmp_55*tmp_55)*tmp_7) + tmp_64*(-tmp_46*tmp_63 + 3.0*(tmp_63*tmp_63)*tmp_7) + tmp_72*(-tmp_46*tmp_71 + 3.0*tmp_7*(tmp_71*tmp_71)) + tmp_80*(-tmp_46*tmp_79 + 3.0*tmp_7*(tmp_79*tmp_79)) + tmp_88*(-tmp_46*tmp_87 + 3.0*tmp_7*(tmp_87*tmp_87));
      real_t a_0_1 = tmp_90;
      real_t a_0_2 = tmp_92;
      real_t a_0_3 = tmp_94;
      real_t a_1_0 = tmp_90;
      real_t a_1_1 = tmp_48*(3.0*(tmp_35*tmp_35)*tmp_7 - tmp_35*tmp_95) + tmp_56*(3.0*(tmp_52*tmp_52)*tmp_7 - tmp_52*tmp_95) + tmp_64*(3.0*(tmp_60*tmp_60)*tmp_7 - tmp_60*tmp_95) + tmp_72*(3.0*(tmp_68*tmp_68)*tmp_7 - tmp_68*tmp_95) + tmp_80*(3.0*tmp_7*(tmp_76*tmp_76) - tmp_76*tmp_95) + tmp_88*(3.0*tmp_7*(tmp_84*tmp_84) - tmp_84*tmp_95);
      real_t a_1_2 = tmp_96;
      real_t a_1_3 = tmp_97;
      real_t a_2_0 = tmp_92;
      real_t a_2_1 = tmp_96;
      real_t a_2_2 = tmp_48*(3.0*(tmp_39*tmp_39)*tmp_7 - tmp_39*tmp_98) + tmp_56*(3.0*(tmp_53*tmp_53)*tmp_7 - tmp_53*tmp_98) + tmp_64*(3.0*(tmp_61*tmp_61)*tmp_7 - tmp_61*tmp_98) + tmp_72*(3.0*(tmp_69*tmp_69)*tmp_7 - tmp_69*tmp_98) + tmp_80*(3.0*tmp_7*(tmp_77*tmp_77) - tmp_77*tmp_98) + tmp_88*(3.0*tmp_7*(tmp_85*tmp_85) - tmp_85*tmp_98);
      real_t a_2_3 = tmp_99;
      real_t a_3_0 = tmp_94;
      real_t a_3_1 = tmp_97;
      real_t a_3_2 = tmp_99;
      real_t a_3_3 = tmp_48*(-tmp_100*tmp_43 + 3.0*(tmp_43*tmp_43)*tmp_7) + tmp_56*(-tmp_100*tmp_54 + 3.0*(tmp_54*tmp_54)*tmp_7) + tmp_64*(-tmp_100*tmp_62 + 3.0*(tmp_62*tmp_62)*tmp_7) + tmp_72*(-tmp_100*tmp_70 + 3.0*tmp_7*(tmp_70*tmp_70)) + tmp_80*(-tmp_100*tmp_78 + 3.0*tmp_7*(tmp_78*tmp_78)) + tmp_88*(-tmp_100*tmp_86 + 3.0*tmp_7*(tmp_86*tmp_86));
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 1, 3) = a_1_3;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
      elMat( 2, 3) = a_2_3;
      elMat( 3, 0) = a_3_0;
      elMat( 3, 1) = a_3_1;
      elMat( 3, 2) = a_3_2;
      elMat( 3, 3) = a_3_3;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g2;

};




class EGVectorLaplaceFormNitscheBC_P1E_2 : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_P1E_2()

    {}





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
                                       MatrixXr&                                           elMat ) const override
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
                                          MatrixXr&                                           elMat ) const override
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
                                                   MatrixXr&                                           elMat ) const override
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
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_8 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_0*tmp_8;
      real_t tmp_13 = tmp_1*tmp_11;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_11*tmp_4 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_2*tmp_6 - tmp_5*tmp_6 + tmp_7*tmp_9);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_11*tmp_4 - tmp_12);
      real_t tmp_17 = tmp_14*(-tmp_13 + tmp_9);
      real_t tmp_18 = tmp_10*tmp_16 + tmp_15*tmp_6 + tmp_17*tmp_7;
      real_t tmp_19 = tmp_14*(-tmp_0*tmp_10 + tmp_3*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_0*tmp_6 - tmp_11*tmp_7);
      real_t tmp_21 = tmp_14*(tmp_10*tmp_11 - tmp_3*tmp_6);
      real_t tmp_22 = tmp_10*tmp_20 + tmp_19*tmp_6 + tmp_21*tmp_7;
      real_t tmp_23 = tmp_14*(-tmp_1*tmp_7 + tmp_10*tmp_4);
      real_t tmp_24 = tmp_14*(-tmp_4*tmp_6 + tmp_7*tmp_8);
      real_t tmp_25 = tmp_14*(tmp_1*tmp_6 - tmp_10*tmp_8);
      real_t tmp_26 = tmp_10*tmp_24 + tmp_23*tmp_6 + tmp_25*tmp_7;
      real_t tmp_27 = p_affine_0_0*p_affine_1_1;
      real_t tmp_28 = p_affine_0_0*p_affine_1_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_2;
      real_t tmp_30 = p_affine_0_1*p_affine_1_0;
      real_t tmp_31 = p_affine_0_1*p_affine_1_2;
      real_t tmp_32 = p_affine_2_2*p_affine_3_0;
      real_t tmp_33 = p_affine_0_2*p_affine_1_0;
      real_t tmp_34 = p_affine_0_2*p_affine_1_1;
      real_t tmp_35 = p_affine_2_0*p_affine_3_1;
      real_t tmp_36 = p_affine_2_2*p_affine_3_1;
      real_t tmp_37 = p_affine_2_0*p_affine_3_2;
      real_t tmp_38 = p_affine_2_1*p_affine_3_0;
      real_t tmp_39 = std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = tmp_39*(tmp_18*(-tmp_15 - tmp_16 - tmp_17) + tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25));
      real_t tmp_41 = tmp_39*(tmp_17*tmp_18 + tmp_21*tmp_22 + tmp_25*tmp_26);
      real_t tmp_42 = tmp_39*(tmp_16*tmp_18 + tmp_20*tmp_22 + tmp_24*tmp_26);
      real_t tmp_43 = tmp_39*(tmp_15*tmp_18 + tmp_19*tmp_22 + tmp_23*tmp_26);
      real_t a_0_0 = 0.16666666666666666*tmp_40;
      real_t a_1_0 = 0.16666666666666666*tmp_41;
      real_t a_2_0 = 0.16666666666666666*tmp_42;
      real_t a_3_0 = 0.16666666666666666*tmp_43;
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
                               MatrixXr&                            elMat ) const override
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

         real_t tmp_0 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_1 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_2 = -tmp_1;
      real_t tmp_3 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_4 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_5 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_3 + tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_7 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_8 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_9 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_10 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_11 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_12 = tmp_11*tmp_9;
      real_t tmp_13 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_14 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_14*tmp_7;
      real_t tmp_17 = tmp_11*tmp_13;
      real_t tmp_18 = tmp_0*tmp_9;
      real_t tmp_19 = 1.0 / (tmp_0*tmp_6*tmp_7 + tmp_10*tmp_12 - tmp_10*tmp_16 + tmp_15*tmp_8 - tmp_17*tmp_6 - tmp_18*tmp_8);
      real_t tmp_20 = tmp_19*(tmp_6*tmp_7 - tmp_8*tmp_9);
      real_t tmp_21 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_22 = -tmp_21;
      real_t tmp_23 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_24 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_25 = 0.091576213509770743*tmp_22 + 0.81684757298045851*tmp_23 + tmp_24;
      real_t tmp_26 = tmp_19*(-tmp_11*tmp_6 + tmp_14*tmp_8);
      real_t tmp_27 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_28 = -tmp_27;
      real_t tmp_29 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_30 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_31 = 0.091576213509770743*tmp_28 + 0.81684757298045851*tmp_29 + tmp_30;
      real_t tmp_32 = tmp_19*(tmp_12 - tmp_16);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_25*tmp_26 + tmp_31*tmp_32;
      real_t tmp_34 = tmp_19*(-tmp_10*tmp_7 + tmp_13*tmp_8);
      real_t tmp_35 = tmp_19*(-tmp_0*tmp_8 + tmp_10*tmp_11);
      real_t tmp_36 = tmp_19*(tmp_0*tmp_7 - tmp_17);
      real_t tmp_37 = tmp_25*tmp_35 + tmp_31*tmp_36 + tmp_34*tmp_5;
      real_t tmp_38 = tmp_19*(tmp_10*tmp_9 - tmp_13*tmp_6);
      real_t tmp_39 = tmp_19*(tmp_0*tmp_6 - tmp_10*tmp_14);
      real_t tmp_40 = tmp_19*(tmp_15 - tmp_18);
      real_t tmp_41 = tmp_25*tmp_39 + tmp_31*tmp_40 + tmp_38*tmp_5;
      real_t tmp_42 = tmp_0*(tmp_33 - 1.0/4.0) + tmp_11*(tmp_41 - 1.0/4.0) + tmp_14*(tmp_37 - 1.0/4.0);
      real_t tmp_43 = 0.5*p_affine_13_0*(-tmp_32 - tmp_36 - tmp_40) + 0.5*p_affine_13_1*(-tmp_26 - tmp_35 - tmp_39) + 0.5*p_affine_13_2*(-tmp_20 - tmp_34 - tmp_38);
      real_t tmp_44 = -tmp_33 - tmp_37 - tmp_41 + 1;
      real_t tmp_45 = 0.5*p_affine_13_0*(tmp_0*tmp_32 + tmp_11*tmp_40 + tmp_14*tmp_36) + 0.5*p_affine_13_1*(tmp_0*tmp_26 + tmp_11*tmp_39 + tmp_14*tmp_35) + 0.5*p_affine_13_2*(tmp_0*tmp_20 + tmp_11*tmp_38 + tmp_14*tmp_34);
      real_t tmp_46 = (std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)*std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)) + (std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)*std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)) + (std::abs(tmp_21*tmp_29 - tmp_23*tmp_27)*std::abs(tmp_21*tmp_29 - tmp_23*tmp_27));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_3 + tmp_4;
      real_t tmp_51 = 0.44594849091596489*tmp_22 + 0.10810301816807022*tmp_23 + tmp_24;
      real_t tmp_52 = 0.44594849091596489*tmp_28 + 0.10810301816807022*tmp_29 + tmp_30;
      real_t tmp_53 = tmp_20*tmp_50 + tmp_26*tmp_51 + tmp_32*tmp_52;
      real_t tmp_54 = tmp_34*tmp_50 + tmp_35*tmp_51 + tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_50 + tmp_39*tmp_51 + tmp_40*tmp_52;
      real_t tmp_56 = tmp_0*(tmp_53 - 1.0/4.0) + tmp_11*(tmp_55 - 1.0/4.0) + tmp_14*(tmp_54 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_60 = 0.81684757298045851*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_61 = 0.81684757298045851*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_62 = tmp_20*tmp_59 + tmp_26*tmp_60 + tmp_32*tmp_61;
      real_t tmp_63 = tmp_34*tmp_59 + tmp_35*tmp_60 + tmp_36*tmp_61;
      real_t tmp_64 = tmp_38*tmp_59 + tmp_39*tmp_60 + tmp_40*tmp_61;
      real_t tmp_65 = tmp_0*(tmp_62 - 1.0/4.0) + tmp_11*(tmp_64 - 1.0/4.0) + tmp_14*(tmp_63 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_69 = 0.10810301816807022*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_70 = 0.10810301816807022*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_71 = tmp_20*tmp_68 + tmp_26*tmp_69 + tmp_32*tmp_70;
      real_t tmp_72 = tmp_34*tmp_68 + tmp_35*tmp_69 + tmp_36*tmp_70;
      real_t tmp_73 = tmp_38*tmp_68 + tmp_39*tmp_69 + tmp_40*tmp_70;
      real_t tmp_74 = tmp_0*(tmp_71 - 1.0/4.0) + tmp_11*(tmp_73 - 1.0/4.0) + tmp_14*(tmp_72 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_78 = 0.091576213509770743*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_79 = 0.091576213509770743*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_80 = tmp_20*tmp_77 + tmp_26*tmp_78 + tmp_32*tmp_79;
      real_t tmp_81 = tmp_34*tmp_77 + tmp_35*tmp_78 + tmp_36*tmp_79;
      real_t tmp_82 = tmp_38*tmp_77 + tmp_39*tmp_78 + tmp_40*tmp_79;
      real_t tmp_83 = tmp_0*(tmp_80 - 1.0/4.0) + tmp_11*(tmp_82 - 1.0/4.0) + tmp_14*(tmp_81 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_87 = 0.44594849091596489*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_88 = 0.44594849091596489*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_89 = tmp_20*tmp_86 + tmp_26*tmp_87 + tmp_32*tmp_88;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = tmp_0*(tmp_89 - 1.0/4.0) + tmp_11*(tmp_91 - 1.0/4.0) + tmp_14*(tmp_90 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = 0.5*p_affine_13_0*tmp_32 + 0.5*p_affine_13_1*tmp_26 + 0.5*p_affine_13_2*tmp_20;
      real_t tmp_96 = 0.5*p_affine_13_0*tmp_36 + 0.5*p_affine_13_1*tmp_35 + 0.5*p_affine_13_2*tmp_34;
      real_t tmp_97 = 0.5*p_affine_13_0*tmp_40 + 0.5*p_affine_13_1*tmp_39 + 0.5*p_affine_13_2*tmp_38;
      real_t a_0_0 = tmp_49*(-tmp_42*tmp_43 + 3.0*tmp_42*tmp_44*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_43*tmp_56 - tmp_45*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_43*tmp_65 - tmp_45*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_43*tmp_74 - tmp_45*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_43*tmp_83 - tmp_45*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_43*tmp_92 - tmp_45*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_1_0 = tmp_49*(3.0*tmp_33*tmp_42*tmp_47 - tmp_33*tmp_45 - tmp_42*tmp_95) + tmp_58*(-tmp_45*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_45*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_45*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_45*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_45*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_2_0 = tmp_49*(3.0*tmp_37*tmp_42*tmp_47 - tmp_37*tmp_45 - tmp_42*tmp_96) + tmp_58*(-tmp_45*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_45*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_45*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_45*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_45*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_3_0 = tmp_49*(3.0*tmp_41*tmp_42*tmp_47 - tmp_41*tmp_45 - tmp_42*tmp_97) + tmp_58*(-tmp_45*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_45*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_45*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_45*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_45*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
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
                                  MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_26 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_27 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_28 = tmp_26*tmp_27;
      real_t tmp_29 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_30 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_31 = tmp_29*tmp_30;
      real_t tmp_32 = tmp_28 - tmp_31;
      real_t tmp_33 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_34 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_35 = tmp_30*tmp_34;
      real_t tmp_36 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_37 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_38 = tmp_27*tmp_37;
      real_t tmp_39 = tmp_26*tmp_34;
      real_t tmp_40 = 1.0 / (tmp_25*tmp_28 - tmp_25*tmp_31 + tmp_29*tmp_36*tmp_37 + tmp_33*tmp_35 - tmp_33*tmp_38 - tmp_36*tmp_39);
      real_t tmp_41 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_42 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_43 = -tmp_42;
      real_t tmp_44 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_45 = 0.091576213509770743*tmp_43 + 0.81684757298045851*tmp_44;
      real_t tmp_46 = tmp_40*(tmp_41 + tmp_45);
      real_t tmp_47 = tmp_29*tmp_37 - tmp_39;
      real_t tmp_48 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_49 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_50 = -tmp_49;
      real_t tmp_51 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_52 = 0.091576213509770743*tmp_50 + 0.81684757298045851*tmp_51;
      real_t tmp_53 = tmp_40*(tmp_48 + tmp_52);
      real_t tmp_54 = tmp_35 - tmp_38;
      real_t tmp_55 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = -tmp_56;
      real_t tmp_58 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_59 = 0.091576213509770743*tmp_57 + 0.81684757298045851*tmp_58;
      real_t tmp_60 = tmp_40*(tmp_55 + tmp_59);
      real_t tmp_61 = -tmp_27*tmp_33 + tmp_29*tmp_36;
      real_t tmp_62 = -tmp_25*tmp_29 + tmp_33*tmp_34;
      real_t tmp_63 = tmp_25*tmp_27 - tmp_34*tmp_36;
      real_t tmp_64 = -tmp_26*tmp_36 + tmp_30*tmp_33;
      real_t tmp_65 = tmp_25*tmp_26 - tmp_33*tmp_37;
      real_t tmp_66 = -tmp_25*tmp_30 + tmp_36*tmp_37;
      real_t tmp_67 = tmp_25*(tmp_32*tmp_46 + tmp_47*tmp_53 + tmp_54*tmp_60 - 1.0/4.0) + tmp_34*(tmp_46*tmp_64 + tmp_53*tmp_65 + tmp_60*tmp_66 - 1.0/4.0) + tmp_37*(tmp_46*tmp_61 + tmp_53*tmp_62 + tmp_60*tmp_63 - 1.0/4.0);
      real_t tmp_68 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_69 = tmp_45 + tmp_68;
      real_t tmp_70 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_71 = tmp_52 + tmp_70;
      real_t tmp_72 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_73 = tmp_59 + tmp_72;
      real_t tmp_74 = tmp_17*tmp_73 + tmp_20*tmp_71 + tmp_23*tmp_69;
      real_t tmp_75 = tmp_16*tmp_73 + tmp_19*tmp_71 + tmp_22*tmp_69;
      real_t tmp_76 = tmp_15*tmp_73 + tmp_18*tmp_71 + tmp_21*tmp_69;
      real_t tmp_77 = -tmp_74 - tmp_75 - tmp_76 + 1;
      real_t tmp_78 = tmp_25*tmp_40;
      real_t tmp_79 = tmp_37*tmp_40;
      real_t tmp_80 = tmp_34*tmp_40;
      real_t tmp_81 = 0.5*p_affine_13_0*(tmp_54*tmp_78 + tmp_63*tmp_79 + tmp_66*tmp_80) + 0.5*p_affine_13_1*(tmp_47*tmp_78 + tmp_62*tmp_79 + tmp_65*tmp_80) + 0.5*p_affine_13_2*(tmp_32*tmp_78 + tmp_61*tmp_79 + tmp_64*tmp_80);
      real_t tmp_82 = (std::abs(tmp_42*tmp_51 - tmp_44*tmp_49)*std::abs(tmp_42*tmp_51 - tmp_44*tmp_49)) + (std::abs(tmp_42*tmp_58 - tmp_44*tmp_56)*std::abs(tmp_42*tmp_58 - tmp_44*tmp_56)) + (std::abs(tmp_49*tmp_58 - tmp_51*tmp_56)*std::abs(tmp_49*tmp_58 - tmp_51*tmp_56));
      real_t tmp_83 = 3.0*std::pow(tmp_82, -0.25);
      real_t tmp_84 = tmp_67*tmp_83;
      real_t tmp_85 = 1.0*std::pow(tmp_82, 1.0/2.0);
      real_t tmp_86 = 0.054975871827660928*tmp_85;
      real_t tmp_87 = 0.44594849091596489*tmp_43 + 0.10810301816807022*tmp_44;
      real_t tmp_88 = tmp_40*(tmp_41 + tmp_87);
      real_t tmp_89 = 0.44594849091596489*tmp_50 + 0.10810301816807022*tmp_51;
      real_t tmp_90 = tmp_40*(tmp_48 + tmp_89);
      real_t tmp_91 = 0.44594849091596489*tmp_57 + 0.10810301816807022*tmp_58;
      real_t tmp_92 = tmp_40*(tmp_55 + tmp_91);
      real_t tmp_93 = tmp_25*(tmp_32*tmp_88 + tmp_47*tmp_90 + tmp_54*tmp_92 - 1.0/4.0) + tmp_34*(tmp_64*tmp_88 + tmp_65*tmp_90 + tmp_66*tmp_92 - 1.0/4.0) + tmp_37*(tmp_61*tmp_88 + tmp_62*tmp_90 + tmp_63*tmp_92 - 1.0/4.0);
      real_t tmp_94 = tmp_68 + tmp_87;
      real_t tmp_95 = tmp_70 + tmp_89;
      real_t tmp_96 = tmp_72 + tmp_91;
      real_t tmp_97 = tmp_17*tmp_96 + tmp_20*tmp_95 + tmp_23*tmp_94;
      real_t tmp_98 = tmp_16*tmp_96 + tmp_19*tmp_95 + tmp_22*tmp_94;
      real_t tmp_99 = tmp_15*tmp_96 + tmp_18*tmp_95 + tmp_21*tmp_94;
      real_t tmp_100 = -tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_101 = tmp_83*tmp_93;
      real_t tmp_102 = 0.11169079483900572*tmp_85;
      real_t tmp_103 = 0.81684757298045851*tmp_43 + 0.091576213509770743*tmp_44;
      real_t tmp_104 = tmp_40*(tmp_103 + tmp_41);
      real_t tmp_105 = 0.81684757298045851*tmp_50 + 0.091576213509770743*tmp_51;
      real_t tmp_106 = tmp_40*(tmp_105 + tmp_48);
      real_t tmp_107 = 0.81684757298045851*tmp_57 + 0.091576213509770743*tmp_58;
      real_t tmp_108 = tmp_40*(tmp_107 + tmp_55);
      real_t tmp_109 = tmp_25*(tmp_104*tmp_32 + tmp_106*tmp_47 + tmp_108*tmp_54 - 1.0/4.0) + tmp_34*(tmp_104*tmp_64 + tmp_106*tmp_65 + tmp_108*tmp_66 - 1.0/4.0) + tmp_37*(tmp_104*tmp_61 + tmp_106*tmp_62 + tmp_108*tmp_63 - 1.0/4.0);
      real_t tmp_110 = tmp_103 + tmp_68;
      real_t tmp_111 = tmp_105 + tmp_70;
      real_t tmp_112 = tmp_107 + tmp_72;
      real_t tmp_113 = tmp_110*tmp_23 + tmp_111*tmp_20 + tmp_112*tmp_17;
      real_t tmp_114 = tmp_110*tmp_22 + tmp_111*tmp_19 + tmp_112*tmp_16;
      real_t tmp_115 = tmp_110*tmp_21 + tmp_111*tmp_18 + tmp_112*tmp_15;
      real_t tmp_116 = -tmp_113 - tmp_114 - tmp_115 + 1;
      real_t tmp_117 = tmp_109*tmp_83;
      real_t tmp_118 = 0.054975871827660928*tmp_85;
      real_t tmp_119 = 0.10810301816807022*tmp_43 + 0.44594849091596489*tmp_44;
      real_t tmp_120 = tmp_40*(tmp_119 + tmp_41);
      real_t tmp_121 = 0.10810301816807022*tmp_50 + 0.44594849091596489*tmp_51;
      real_t tmp_122 = tmp_40*(tmp_121 + tmp_48);
      real_t tmp_123 = 0.10810301816807022*tmp_57 + 0.44594849091596489*tmp_58;
      real_t tmp_124 = tmp_40*(tmp_123 + tmp_55);
      real_t tmp_125 = tmp_25*(tmp_120*tmp_32 + tmp_122*tmp_47 + tmp_124*tmp_54 - 1.0/4.0) + tmp_34*(tmp_120*tmp_64 + tmp_122*tmp_65 + tmp_124*tmp_66 - 1.0/4.0) + tmp_37*(tmp_120*tmp_61 + tmp_122*tmp_62 + tmp_124*tmp_63 - 1.0/4.0);
      real_t tmp_126 = tmp_119 + tmp_68;
      real_t tmp_127 = tmp_121 + tmp_70;
      real_t tmp_128 = tmp_123 + tmp_72;
      real_t tmp_129 = tmp_126*tmp_23 + tmp_127*tmp_20 + tmp_128*tmp_17;
      real_t tmp_130 = tmp_126*tmp_22 + tmp_127*tmp_19 + tmp_128*tmp_16;
      real_t tmp_131 = tmp_126*tmp_21 + tmp_127*tmp_18 + tmp_128*tmp_15;
      real_t tmp_132 = -tmp_129 - tmp_130 - tmp_131 + 1;
      real_t tmp_133 = tmp_125*tmp_83;
      real_t tmp_134 = 0.11169079483900572*tmp_85;
      real_t tmp_135 = 0.091576213509770743*tmp_43 + 0.091576213509770743*tmp_44;
      real_t tmp_136 = tmp_40*(tmp_135 + tmp_41);
      real_t tmp_137 = 0.091576213509770743*tmp_50 + 0.091576213509770743*tmp_51;
      real_t tmp_138 = tmp_40*(tmp_137 + tmp_48);
      real_t tmp_139 = 0.091576213509770743*tmp_57 + 0.091576213509770743*tmp_58;
      real_t tmp_140 = tmp_40*(tmp_139 + tmp_55);
      real_t tmp_141 = tmp_25*(tmp_136*tmp_32 + tmp_138*tmp_47 + tmp_140*tmp_54 - 1.0/4.0) + tmp_34*(tmp_136*tmp_64 + tmp_138*tmp_65 + tmp_140*tmp_66 - 1.0/4.0) + tmp_37*(tmp_136*tmp_61 + tmp_138*tmp_62 + tmp_140*tmp_63 - 1.0/4.0);
      real_t tmp_142 = tmp_135 + tmp_68;
      real_t tmp_143 = tmp_137 + tmp_70;
      real_t tmp_144 = tmp_139 + tmp_72;
      real_t tmp_145 = tmp_142*tmp_23 + tmp_143*tmp_20 + tmp_144*tmp_17;
      real_t tmp_146 = tmp_142*tmp_22 + tmp_143*tmp_19 + tmp_144*tmp_16;
      real_t tmp_147 = tmp_142*tmp_21 + tmp_143*tmp_18 + tmp_144*tmp_15;
      real_t tmp_148 = -tmp_145 - tmp_146 - tmp_147 + 1;
      real_t tmp_149 = tmp_141*tmp_83;
      real_t tmp_150 = 0.054975871827660928*tmp_85;
      real_t tmp_151 = 0.44594849091596489*tmp_43 + 0.44594849091596489*tmp_44;
      real_t tmp_152 = tmp_40*(tmp_151 + tmp_41);
      real_t tmp_153 = 0.44594849091596489*tmp_50 + 0.44594849091596489*tmp_51;
      real_t tmp_154 = tmp_40*(tmp_153 + tmp_48);
      real_t tmp_155 = 0.44594849091596489*tmp_57 + 0.44594849091596489*tmp_58;
      real_t tmp_156 = tmp_40*(tmp_155 + tmp_55);
      real_t tmp_157 = tmp_25*(tmp_152*tmp_32 + tmp_154*tmp_47 + tmp_156*tmp_54 - 1.0/4.0) + tmp_34*(tmp_152*tmp_64 + tmp_154*tmp_65 + tmp_156*tmp_66 - 1.0/4.0) + tmp_37*(tmp_152*tmp_61 + tmp_154*tmp_62 + tmp_156*tmp_63 - 1.0/4.0);
      real_t tmp_158 = tmp_151 + tmp_68;
      real_t tmp_159 = tmp_153 + tmp_70;
      real_t tmp_160 = tmp_155 + tmp_72;
      real_t tmp_161 = tmp_158*tmp_23 + tmp_159*tmp_20 + tmp_160*tmp_17;
      real_t tmp_162 = tmp_158*tmp_22 + tmp_159*tmp_19 + tmp_16*tmp_160;
      real_t tmp_163 = tmp_15*tmp_160 + tmp_158*tmp_21 + tmp_159*tmp_18;
      real_t tmp_164 = -tmp_161 - tmp_162 - tmp_163 + 1;
      real_t tmp_165 = tmp_157*tmp_83;
      real_t tmp_166 = 0.11169079483900572*tmp_85;
      real_t tmp_167 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_168 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_169 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t a_0_0 = tmp_102*(-tmp_100*tmp_101 - tmp_100*tmp_81 + 0.5*tmp_24*tmp_93) + tmp_118*(0.5*tmp_109*tmp_24 - tmp_116*tmp_117 - tmp_116*tmp_81) + tmp_134*(0.5*tmp_125*tmp_24 - tmp_132*tmp_133 - tmp_132*tmp_81) + tmp_150*(0.5*tmp_141*tmp_24 - tmp_148*tmp_149 - tmp_148*tmp_81) + tmp_166*(0.5*tmp_157*tmp_24 - tmp_164*tmp_165 - tmp_164*tmp_81) + tmp_86*(0.5*tmp_24*tmp_67 - tmp_77*tmp_81 - tmp_77*tmp_84);
      real_t a_1_0 = tmp_102*(-tmp_101*tmp_97 + 0.5*tmp_167*tmp_93 - tmp_81*tmp_97) + tmp_118*(0.5*tmp_109*tmp_167 - tmp_113*tmp_117 - tmp_113*tmp_81) + tmp_134*(0.5*tmp_125*tmp_167 - tmp_129*tmp_133 - tmp_129*tmp_81) + tmp_150*(0.5*tmp_141*tmp_167 - tmp_145*tmp_149 - tmp_145*tmp_81) + tmp_166*(0.5*tmp_157*tmp_167 - tmp_161*tmp_165 - tmp_161*tmp_81) + tmp_86*(0.5*tmp_167*tmp_67 - tmp_74*tmp_81 - tmp_74*tmp_84);
      real_t a_2_0 = tmp_102*(-tmp_101*tmp_98 + 0.5*tmp_168*tmp_93 - tmp_81*tmp_98) + tmp_118*(0.5*tmp_109*tmp_168 - tmp_114*tmp_117 - tmp_114*tmp_81) + tmp_134*(0.5*tmp_125*tmp_168 - tmp_130*tmp_133 - tmp_130*tmp_81) + tmp_150*(0.5*tmp_141*tmp_168 - tmp_146*tmp_149 - tmp_146*tmp_81) + tmp_166*(0.5*tmp_157*tmp_168 - tmp_162*tmp_165 - tmp_162*tmp_81) + tmp_86*(0.5*tmp_168*tmp_67 - tmp_75*tmp_81 - tmp_75*tmp_84);
      real_t a_3_0 = tmp_102*(-tmp_101*tmp_99 + 0.5*tmp_169*tmp_93 - tmp_81*tmp_99) + tmp_118*(0.5*tmp_109*tmp_169 - tmp_115*tmp_117 - tmp_115*tmp_81) + tmp_134*(0.5*tmp_125*tmp_169 - tmp_131*tmp_133 - tmp_131*tmp_81) + tmp_150*(0.5*tmp_141*tmp_169 - tmp_147*tmp_149 - tmp_147*tmp_81) + tmp_166*(0.5*tmp_157*tmp_169 - tmp_163*tmp_165 - tmp_163*tmp_81) + tmp_86*(0.5*tmp_169*tmp_67 - tmp_76*tmp_81 - tmp_76*tmp_84);
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
                                        MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_26 = -tmp_25;
      real_t tmp_27 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_28 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_29 = 0.091576213509770743*tmp_26 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -tmp_30;
      real_t tmp_32 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_33 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_34 = 0.091576213509770743*tmp_31 + 0.81684757298045851*tmp_32 + tmp_33;
      real_t tmp_35 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_36 = -tmp_35;
      real_t tmp_37 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_38 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_39 = 0.091576213509770743*tmp_36 + 0.81684757298045851*tmp_37 + tmp_38;
      real_t tmp_40 = tmp_17*tmp_39 + tmp_20*tmp_34 + tmp_23*tmp_29;
      real_t tmp_41 = tmp_16*tmp_39 + tmp_19*tmp_34 + tmp_22*tmp_29;
      real_t tmp_42 = tmp_15*tmp_39 + tmp_18*tmp_34 + tmp_21*tmp_29;
      real_t tmp_43 = tmp_1*(tmp_41 - 1.0/4.0) + tmp_4*(tmp_40 - 1.0/4.0) + tmp_7*(tmp_42 - 1.0/4.0);
      real_t tmp_44 = p_affine_13_0*(tmp_1*tmp_16 + tmp_15*tmp_7 + tmp_17*tmp_4) + p_affine_13_1*(tmp_1*tmp_19 + tmp_18*tmp_7 + tmp_20*tmp_4) + p_affine_13_2*(tmp_1*tmp_22 + tmp_21*tmp_7 + tmp_23*tmp_4);
      real_t tmp_45 = -tmp_40 - tmp_41 - tmp_42 + 1;
      real_t tmp_46 = (std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)*std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)) + (std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)*std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)) + (std::abs(tmp_30*tmp_37 - tmp_32*tmp_35)*std::abs(tmp_30*tmp_37 - tmp_32*tmp_35));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_26 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_51 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32 + tmp_33;
      real_t tmp_52 = 0.44594849091596489*tmp_36 + 0.10810301816807022*tmp_37 + tmp_38;
      real_t tmp_53 = tmp_17*tmp_52 + tmp_20*tmp_51 + tmp_23*tmp_50;
      real_t tmp_54 = tmp_16*tmp_52 + tmp_19*tmp_51 + tmp_22*tmp_50;
      real_t tmp_55 = tmp_15*tmp_52 + tmp_18*tmp_51 + tmp_21*tmp_50;
      real_t tmp_56 = tmp_1*(tmp_54 - 1.0/4.0) + tmp_4*(tmp_53 - 1.0/4.0) + tmp_7*(tmp_55 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_60 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_61 = 0.81684757298045851*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_62 = tmp_17*tmp_61 + tmp_20*tmp_60 + tmp_23*tmp_59;
      real_t tmp_63 = tmp_16*tmp_61 + tmp_19*tmp_60 + tmp_22*tmp_59;
      real_t tmp_64 = tmp_15*tmp_61 + tmp_18*tmp_60 + tmp_21*tmp_59;
      real_t tmp_65 = tmp_1*(tmp_63 - 1.0/4.0) + tmp_4*(tmp_62 - 1.0/4.0) + tmp_7*(tmp_64 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_69 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_70 = 0.10810301816807022*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_71 = tmp_17*tmp_70 + tmp_20*tmp_69 + tmp_23*tmp_68;
      real_t tmp_72 = tmp_16*tmp_70 + tmp_19*tmp_69 + tmp_22*tmp_68;
      real_t tmp_73 = tmp_15*tmp_70 + tmp_18*tmp_69 + tmp_21*tmp_68;
      real_t tmp_74 = tmp_1*(tmp_72 - 1.0/4.0) + tmp_4*(tmp_71 - 1.0/4.0) + tmp_7*(tmp_73 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_78 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_79 = 0.091576213509770743*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_80 = tmp_17*tmp_79 + tmp_20*tmp_78 + tmp_23*tmp_77;
      real_t tmp_81 = tmp_16*tmp_79 + tmp_19*tmp_78 + tmp_22*tmp_77;
      real_t tmp_82 = tmp_15*tmp_79 + tmp_18*tmp_78 + tmp_21*tmp_77;
      real_t tmp_83 = tmp_1*(tmp_81 - 1.0/4.0) + tmp_4*(tmp_80 - 1.0/4.0) + tmp_7*(tmp_82 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_87 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_88 = 0.44594849091596489*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_89 = tmp_17*tmp_88 + tmp_20*tmp_87 + tmp_23*tmp_86;
      real_t tmp_90 = tmp_16*tmp_88 + tmp_19*tmp_87 + tmp_22*tmp_86;
      real_t tmp_91 = tmp_15*tmp_88 + tmp_18*tmp_87 + tmp_21*tmp_86;
      real_t tmp_92 = tmp_1*(tmp_90 - 1.0/4.0) + tmp_4*(tmp_89 - 1.0/4.0) + tmp_7*(tmp_91 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_96 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_97 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t a_0_0 = tmp_49*(-tmp_24*tmp_43 + 3.0*tmp_43*tmp_45*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_24*tmp_56 - tmp_44*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_24*tmp_65 - tmp_44*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_24*tmp_74 - tmp_44*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_24*tmp_83 - tmp_44*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_24*tmp_92 - tmp_44*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_1_0 = tmp_49*(3.0*tmp_40*tmp_43*tmp_47 - tmp_40*tmp_44 - tmp_43*tmp_95) + tmp_58*(-tmp_44*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_44*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_44*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_44*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_44*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_2_0 = tmp_49*(3.0*tmp_41*tmp_43*tmp_47 - tmp_41*tmp_44 - tmp_43*tmp_96) + tmp_58*(-tmp_44*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_44*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_44*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_44*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_44*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_3_0 = tmp_49*(3.0*tmp_42*tmp_43*tmp_47 - tmp_42*tmp_44 - tmp_43*tmp_97) + tmp_58*(-tmp_44*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_44*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_44*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_44*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_44*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }

public:



};




class EGVectorLaplaceFormNitscheBC_EP1_0 : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_EP1_0()

    {}





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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = -tmp_3;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = 1.0 / (tmp_2 + tmp_4*tmp_5);
      real_t tmp_7 = tmp_0*tmp_6;
      real_t tmp_8 = tmp_3*tmp_6;
      real_t tmp_9 = tmp_0*tmp_8 + tmp_4*tmp_7;
      real_t tmp_10 = tmp_1*tmp_6;
      real_t tmp_11 = tmp_5*tmp_6;
      real_t tmp_12 = tmp_11*tmp_4 + tmp_2*tmp_6;
      real_t tmp_13 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = tmp_13*(tmp_12*(-tmp_10 - tmp_11) + tmp_9*(-tmp_7 - tmp_8));
      real_t tmp_15 = tmp_13*(tmp_10*tmp_12 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_13*(tmp_11*tmp_12 + tmp_7*tmp_9);
      real_t a_0_0 = 0.5*tmp_14;
      real_t a_0_1 = 0.5*tmp_15;
      real_t a_0_2 = 0.5*tmp_16;
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
                                       MatrixXr&                                           elMat ) const override
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
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_3 = 0.21132486540518713*tmp_1 + tmp_2;
      real_t tmp_4 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_0*tmp_5;
      real_t tmp_7 = -tmp_4;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = 1.0 / (tmp_6 + tmp_7*tmp_8);
      real_t tmp_10 = tmp_4*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_11 + tmp_12;
      real_t tmp_14 = tmp_5*tmp_9;
      real_t tmp_15 = tmp_10*tmp_3 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_0*tmp_9;
      real_t tmp_17 = tmp_8*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_3;
      real_t tmp_19 = tmp_0*(tmp_15 - 1.0/3.0) + tmp_7*(tmp_18 - 1.0/3.0);
      real_t tmp_20 = 0.5*p_affine_10_0*(-tmp_14 - tmp_17) + 0.5*p_affine_10_1*(-tmp_10 - tmp_16);
      real_t tmp_21 = -tmp_15 - tmp_18 + 1;
      real_t tmp_22 = 0.5*p_affine_10_0*(tmp_17*tmp_7 + tmp_6*tmp_9) + 0.5*p_affine_10_1*(tmp_0*tmp_10 + tmp_16*tmp_7);
      real_t tmp_23 = std::abs(std::pow((tmp_1*tmp_1) + (tmp_11*tmp_11), 1.0/2.0));
      real_t tmp_24 = 1.0 / (tmp_23);
      real_t tmp_25 = 0.5*tmp_23;
      real_t tmp_26 = 0.78867513459481287*tmp_1 + tmp_2;
      real_t tmp_27 = 0.78867513459481287*tmp_11 + tmp_12;
      real_t tmp_28 = tmp_10*tmp_26 + tmp_14*tmp_27;
      real_t tmp_29 = tmp_16*tmp_26 + tmp_17*tmp_27;
      real_t tmp_30 = tmp_0*(tmp_28 - 1.0/3.0) + tmp_7*(tmp_29 - 1.0/3.0);
      real_t tmp_31 = -tmp_28 - tmp_29 + 1;
      real_t tmp_32 = 0.5*tmp_23;
      real_t tmp_33 = 0.5*p_affine_10_0*tmp_14 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_34 = 0.5*p_affine_10_0*tmp_17 + 0.5*p_affine_10_1*tmp_16;
      real_t a_0_0 = tmp_25*(-tmp_19*tmp_20 + 3*tmp_19*tmp_21*tmp_24 - tmp_21*tmp_22) + tmp_32*(-tmp_20*tmp_30 - tmp_22*tmp_31 + 3*tmp_24*tmp_30*tmp_31);
      real_t a_0_1 = tmp_25*(3*tmp_15*tmp_19*tmp_24 - tmp_15*tmp_22 - tmp_19*tmp_33) + tmp_32*(-tmp_22*tmp_28 + 3*tmp_24*tmp_28*tmp_30 - tmp_30*tmp_33);
      real_t a_0_2 = tmp_25*(3*tmp_18*tmp_19*tmp_24 - tmp_18*tmp_22 - tmp_19*tmp_34) + tmp_32*(-tmp_22*tmp_29 + 3*tmp_24*tmp_29*tmp_30 - tmp_30*tmp_34);
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
                                          MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3 = p_affine_6_1 + 0.21132486540518713*tmp_2;
      real_t tmp_4 = tmp_1 + tmp_3;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_0*tmp_5;
      real_t tmp_7 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_8 = -tmp_7;
      real_t tmp_9 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = 1.0 / (tmp_6 + tmp_8*tmp_9);
      real_t tmp_11 = tmp_10*tmp_7;
      real_t tmp_12 = -p_affine_0_0;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + 0.21132486540518713*tmp_13;
      real_t tmp_15 = tmp_10*(tmp_12 + tmp_14);
      real_t tmp_16 = tmp_0*tmp_10;
      real_t tmp_17 = tmp_0*(tmp_11*tmp_4 + tmp_15*tmp_5 - 1.0/3.0) + tmp_8*(tmp_15*tmp_9 + tmp_16*tmp_4 - 1.0/3.0);
      real_t tmp_18 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_19 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_20 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_21 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_22 = 1.0 / (tmp_18*tmp_19 - tmp_20*tmp_21);
      real_t tmp_23 = tmp_18*tmp_22;
      real_t tmp_24 = tmp_21*tmp_22;
      real_t tmp_25 = tmp_19*tmp_22;
      real_t tmp_26 = tmp_20*tmp_22;
      real_t tmp_27 = 0.5*p_affine_10_0*(-tmp_23 - tmp_24) + 0.5*p_affine_10_1*(-tmp_25 - tmp_26);
      real_t tmp_28 = tmp_10*tmp_8;
      real_t tmp_29 = p_affine_10_0*(tmp_10*tmp_6 + tmp_28*tmp_9) + p_affine_10_1*(tmp_0*tmp_11 + tmp_0*tmp_28);
      real_t tmp_30 = -p_affine_3_1;
      real_t tmp_31 = tmp_3 + tmp_30;
      real_t tmp_32 = -p_affine_3_0;
      real_t tmp_33 = tmp_14 + tmp_32;
      real_t tmp_34 = tmp_23*tmp_33 + tmp_26*tmp_31;
      real_t tmp_35 = tmp_24*tmp_33 + tmp_25*tmp_31;
      real_t tmp_36 = -tmp_34 - tmp_35 + 1;
      real_t tmp_37 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_2*tmp_2), 1.0/2.0));
      real_t tmp_38 = 3/tmp_37;
      real_t tmp_39 = tmp_17*tmp_38;
      real_t tmp_40 = 0.5*tmp_37;
      real_t tmp_41 = p_affine_6_1 + 0.78867513459481287*tmp_2;
      real_t tmp_42 = tmp_1 + tmp_41;
      real_t tmp_43 = p_affine_6_0 + 0.78867513459481287*tmp_13;
      real_t tmp_44 = tmp_10*(tmp_12 + tmp_43);
      real_t tmp_45 = tmp_0*(tmp_11*tmp_42 + tmp_44*tmp_5 - 1.0/3.0) + tmp_8*(tmp_16*tmp_42 + tmp_44*tmp_9 - 1.0/3.0);
      real_t tmp_46 = tmp_30 + tmp_41;
      real_t tmp_47 = tmp_32 + tmp_43;
      real_t tmp_48 = tmp_23*tmp_47 + tmp_26*tmp_46;
      real_t tmp_49 = tmp_24*tmp_47 + tmp_25*tmp_46;
      real_t tmp_50 = -tmp_48 - tmp_49 + 1;
      real_t tmp_51 = tmp_38*tmp_45;
      real_t tmp_52 = 0.5*tmp_37;
      real_t tmp_53 = 0.5*p_affine_10_0*tmp_23 + 0.5*p_affine_10_1*tmp_26;
      real_t tmp_54 = 0.5*p_affine_10_0*tmp_24 + 0.5*p_affine_10_1*tmp_25;
      real_t a_0_0 = tmp_40*(-tmp_17*tmp_27 + 0.5*tmp_29*tmp_36 - tmp_36*tmp_39) + tmp_52*(-tmp_27*tmp_45 + 0.5*tmp_29*tmp_50 - tmp_50*tmp_51);
      real_t a_0_1 = tmp_40*(-tmp_17*tmp_53 + 0.5*tmp_29*tmp_34 - tmp_34*tmp_39) + tmp_52*(0.5*tmp_29*tmp_48 - tmp_45*tmp_53 - tmp_48*tmp_51);
      real_t a_0_2 = tmp_40*(-tmp_17*tmp_54 + 0.5*tmp_29*tmp_35 - tmp_35*tmp_39) + tmp_52*(0.5*tmp_29*tmp_49 - tmp_45*tmp_54 - tmp_49*tmp_51);
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
                                                   MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = -tmp_3;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = 1.0 / (tmp_2 + tmp_4*tmp_5);
      real_t tmp_7 = tmp_0*tmp_6;
      real_t tmp_8 = tmp_5*tmp_6;
      real_t tmp_9 = tmp_1*tmp_6;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = p_affine_10_0*(-tmp_7 - tmp_8) + p_affine_10_1*(-tmp_10 - tmp_9);
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_14 = 0.21132486540518713*tmp_12 + tmp_13;
      real_t tmp_15 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_16 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_17 = 0.21132486540518713*tmp_15 + tmp_16;
      real_t tmp_18 = tmp_10*tmp_14 + tmp_17*tmp_7;
      real_t tmp_19 = tmp_14*tmp_9 + tmp_17*tmp_8;
      real_t tmp_20 = tmp_1*(tmp_18 - 1.0/3.0) + tmp_4*(tmp_19 - 1.0/3.0);
      real_t tmp_21 = p_affine_10_0*(tmp_2*tmp_6 + tmp_4*tmp_8) + p_affine_10_1*(tmp_1*tmp_10 + tmp_4*tmp_9);
      real_t tmp_22 = -tmp_18 - tmp_19 + 1;
      real_t tmp_23 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_15*tmp_15), 1.0/2.0));
      real_t tmp_24 = 1.0 / (tmp_23);
      real_t tmp_25 = 0.5*tmp_23;
      real_t tmp_26 = 0.78867513459481287*tmp_12 + tmp_13;
      real_t tmp_27 = 0.78867513459481287*tmp_15 + tmp_16;
      real_t tmp_28 = tmp_10*tmp_26 + tmp_27*tmp_7;
      real_t tmp_29 = tmp_26*tmp_9 + tmp_27*tmp_8;
      real_t tmp_30 = tmp_1*(tmp_28 - 1.0/3.0) + tmp_4*(tmp_29 - 1.0/3.0);
      real_t tmp_31 = -tmp_28 - tmp_29 + 1;
      real_t tmp_32 = 0.5*tmp_23;
      real_t tmp_33 = p_affine_10_0*tmp_7 + p_affine_10_1*tmp_10;
      real_t tmp_34 = p_affine_10_0*tmp_8 + p_affine_10_1*tmp_9;
      real_t a_0_0 = tmp_25*(-tmp_11*tmp_20 + 3*tmp_20*tmp_22*tmp_24 - tmp_21*tmp_22) + tmp_32*(-tmp_11*tmp_30 - tmp_21*tmp_31 + 3*tmp_24*tmp_30*tmp_31);
      real_t a_0_1 = tmp_25*(3*tmp_18*tmp_20*tmp_24 - tmp_18*tmp_21 - tmp_20*tmp_33) + tmp_32*(-tmp_21*tmp_28 + 3*tmp_24*tmp_28*tmp_30 - tmp_30*tmp_33);
      real_t a_0_2 = tmp_25*(3*tmp_19*tmp_20*tmp_24 - tmp_19*tmp_21 - tmp_20*tmp_34) + tmp_32*(-tmp_21*tmp_29 + 3*tmp_24*tmp_29*tmp_30 - tmp_30*tmp_34);
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
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_8 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_0*tmp_8;
      real_t tmp_13 = tmp_1*tmp_11;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_11*tmp_4 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_2*tmp_6 - tmp_5*tmp_6 + tmp_7*tmp_9);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_11*tmp_4 - tmp_12);
      real_t tmp_17 = tmp_14*(-tmp_13 + tmp_9);
      real_t tmp_18 = tmp_0*tmp_17 + tmp_11*tmp_15 + tmp_16*tmp_3;
      real_t tmp_19 = tmp_14*(-tmp_0*tmp_10 + tmp_3*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_0*tmp_6 - tmp_11*tmp_7);
      real_t tmp_21 = tmp_14*(tmp_10*tmp_11 - tmp_3*tmp_6);
      real_t tmp_22 = tmp_0*tmp_21 + tmp_11*tmp_19 + tmp_20*tmp_3;
      real_t tmp_23 = tmp_14*(-tmp_1*tmp_7 + tmp_10*tmp_4);
      real_t tmp_24 = tmp_14*(-tmp_4*tmp_6 + tmp_7*tmp_8);
      real_t tmp_25 = tmp_14*(tmp_1*tmp_6 - tmp_10*tmp_8);
      real_t tmp_26 = tmp_0*tmp_25 + tmp_11*tmp_23 + tmp_24*tmp_3;
      real_t tmp_27 = p_affine_0_0*p_affine_1_1;
      real_t tmp_28 = p_affine_0_0*p_affine_1_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_2;
      real_t tmp_30 = p_affine_0_1*p_affine_1_0;
      real_t tmp_31 = p_affine_0_1*p_affine_1_2;
      real_t tmp_32 = p_affine_2_2*p_affine_3_0;
      real_t tmp_33 = p_affine_0_2*p_affine_1_0;
      real_t tmp_34 = p_affine_0_2*p_affine_1_1;
      real_t tmp_35 = p_affine_2_0*p_affine_3_1;
      real_t tmp_36 = p_affine_2_2*p_affine_3_1;
      real_t tmp_37 = p_affine_2_0*p_affine_3_2;
      real_t tmp_38 = p_affine_2_1*p_affine_3_0;
      real_t tmp_39 = std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = tmp_39*(tmp_18*(-tmp_15 - tmp_16 - tmp_17) + tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25));
      real_t tmp_41 = tmp_39*(tmp_17*tmp_18 + tmp_21*tmp_22 + tmp_25*tmp_26);
      real_t tmp_42 = tmp_39*(tmp_16*tmp_18 + tmp_20*tmp_22 + tmp_24*tmp_26);
      real_t tmp_43 = tmp_39*(tmp_15*tmp_18 + tmp_19*tmp_22 + tmp_23*tmp_26);
      real_t a_0_0 = 0.16666666666666666*tmp_40;
      real_t a_0_1 = 0.16666666666666666*tmp_41;
      real_t a_0_2 = 0.16666666666666666*tmp_42;
      real_t a_0_3 = 0.16666666666666666*tmp_43;
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
                               MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_2 = -tmp_1;
      real_t tmp_3 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_4 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_5 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_3 + tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_7 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_8 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_9 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_13 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_14 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_14*tmp_7;
      real_t tmp_17 = tmp_10*tmp_13;
      real_t tmp_18 = tmp_12*tmp_9;
      real_t tmp_19 = 1.0 / (tmp_0*tmp_11 - tmp_0*tmp_16 + tmp_12*tmp_6*tmp_7 + tmp_15*tmp_8 - tmp_17*tmp_6 - tmp_18*tmp_8);
      real_t tmp_20 = tmp_19*(tmp_6*tmp_7 - tmp_8*tmp_9);
      real_t tmp_21 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_22 = -tmp_21;
      real_t tmp_23 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_24 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_25 = 0.091576213509770743*tmp_22 + 0.81684757298045851*tmp_23 + tmp_24;
      real_t tmp_26 = tmp_19*(-tmp_10*tmp_6 + tmp_14*tmp_8);
      real_t tmp_27 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_28 = -tmp_27;
      real_t tmp_29 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_30 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_31 = 0.091576213509770743*tmp_28 + 0.81684757298045851*tmp_29 + tmp_30;
      real_t tmp_32 = tmp_19*(tmp_11 - tmp_16);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_25*tmp_26 + tmp_31*tmp_32;
      real_t tmp_34 = tmp_19*(-tmp_0*tmp_7 + tmp_13*tmp_8);
      real_t tmp_35 = tmp_19*(tmp_0*tmp_10 - tmp_12*tmp_8);
      real_t tmp_36 = tmp_19*(tmp_12*tmp_7 - tmp_17);
      real_t tmp_37 = tmp_25*tmp_35 + tmp_31*tmp_36 + tmp_34*tmp_5;
      real_t tmp_38 = tmp_19*(tmp_0*tmp_9 - tmp_13*tmp_6);
      real_t tmp_39 = tmp_19*(-tmp_0*tmp_14 + tmp_12*tmp_6);
      real_t tmp_40 = tmp_19*(tmp_15 - tmp_18);
      real_t tmp_41 = tmp_25*tmp_39 + tmp_31*tmp_40 + tmp_38*tmp_5;
      real_t tmp_42 = tmp_0*(tmp_33 - 1.0/4.0) + tmp_6*(tmp_37 - 1.0/4.0) + tmp_8*(tmp_41 - 1.0/4.0);
      real_t tmp_43 = 0.5*p_affine_13_0*(-tmp_32 - tmp_36 - tmp_40) + 0.5*p_affine_13_1*(-tmp_26 - tmp_35 - tmp_39) + 0.5*p_affine_13_2*(-tmp_20 - tmp_34 - tmp_38);
      real_t tmp_44 = -tmp_33 - tmp_37 - tmp_41 + 1;
      real_t tmp_45 = 0.5*p_affine_13_0*(tmp_0*tmp_32 + tmp_36*tmp_6 + tmp_40*tmp_8) + 0.5*p_affine_13_1*(tmp_0*tmp_26 + tmp_35*tmp_6 + tmp_39*tmp_8) + 0.5*p_affine_13_2*(tmp_0*tmp_20 + tmp_34*tmp_6 + tmp_38*tmp_8);
      real_t tmp_46 = (std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)*std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)) + (std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)*std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)) + (std::abs(tmp_21*tmp_29 - tmp_23*tmp_27)*std::abs(tmp_21*tmp_29 - tmp_23*tmp_27));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_3 + tmp_4;
      real_t tmp_51 = 0.44594849091596489*tmp_22 + 0.10810301816807022*tmp_23 + tmp_24;
      real_t tmp_52 = 0.44594849091596489*tmp_28 + 0.10810301816807022*tmp_29 + tmp_30;
      real_t tmp_53 = tmp_20*tmp_50 + tmp_26*tmp_51 + tmp_32*tmp_52;
      real_t tmp_54 = tmp_34*tmp_50 + tmp_35*tmp_51 + tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_50 + tmp_39*tmp_51 + tmp_40*tmp_52;
      real_t tmp_56 = tmp_0*(tmp_53 - 1.0/4.0) + tmp_6*(tmp_54 - 1.0/4.0) + tmp_8*(tmp_55 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_60 = 0.81684757298045851*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_61 = 0.81684757298045851*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_62 = tmp_20*tmp_59 + tmp_26*tmp_60 + tmp_32*tmp_61;
      real_t tmp_63 = tmp_34*tmp_59 + tmp_35*tmp_60 + tmp_36*tmp_61;
      real_t tmp_64 = tmp_38*tmp_59 + tmp_39*tmp_60 + tmp_40*tmp_61;
      real_t tmp_65 = tmp_0*(tmp_62 - 1.0/4.0) + tmp_6*(tmp_63 - 1.0/4.0) + tmp_8*(tmp_64 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_69 = 0.10810301816807022*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_70 = 0.10810301816807022*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_71 = tmp_20*tmp_68 + tmp_26*tmp_69 + tmp_32*tmp_70;
      real_t tmp_72 = tmp_34*tmp_68 + tmp_35*tmp_69 + tmp_36*tmp_70;
      real_t tmp_73 = tmp_38*tmp_68 + tmp_39*tmp_69 + tmp_40*tmp_70;
      real_t tmp_74 = tmp_0*(tmp_71 - 1.0/4.0) + tmp_6*(tmp_72 - 1.0/4.0) + tmp_8*(tmp_73 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_78 = 0.091576213509770743*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_79 = 0.091576213509770743*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_80 = tmp_20*tmp_77 + tmp_26*tmp_78 + tmp_32*tmp_79;
      real_t tmp_81 = tmp_34*tmp_77 + tmp_35*tmp_78 + tmp_36*tmp_79;
      real_t tmp_82 = tmp_38*tmp_77 + tmp_39*tmp_78 + tmp_40*tmp_79;
      real_t tmp_83 = tmp_0*(tmp_80 - 1.0/4.0) + tmp_6*(tmp_81 - 1.0/4.0) + tmp_8*(tmp_82 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_87 = 0.44594849091596489*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_88 = 0.44594849091596489*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_89 = tmp_20*tmp_86 + tmp_26*tmp_87 + tmp_32*tmp_88;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = tmp_0*(tmp_89 - 1.0/4.0) + tmp_6*(tmp_90 - 1.0/4.0) + tmp_8*(tmp_91 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = 0.5*p_affine_13_0*tmp_32 + 0.5*p_affine_13_1*tmp_26 + 0.5*p_affine_13_2*tmp_20;
      real_t tmp_96 = 0.5*p_affine_13_0*tmp_36 + 0.5*p_affine_13_1*tmp_35 + 0.5*p_affine_13_2*tmp_34;
      real_t tmp_97 = 0.5*p_affine_13_0*tmp_40 + 0.5*p_affine_13_1*tmp_39 + 0.5*p_affine_13_2*tmp_38;
      real_t a_0_0 = tmp_49*(-tmp_42*tmp_43 + 3.0*tmp_42*tmp_44*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_43*tmp_56 - tmp_45*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_43*tmp_65 - tmp_45*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_43*tmp_74 - tmp_45*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_43*tmp_83 - tmp_45*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_43*tmp_92 - tmp_45*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_0_1 = tmp_49*(3.0*tmp_33*tmp_42*tmp_47 - tmp_33*tmp_45 - tmp_42*tmp_95) + tmp_58*(-tmp_45*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_45*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_45*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_45*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_45*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_0_2 = tmp_49*(3.0*tmp_37*tmp_42*tmp_47 - tmp_37*tmp_45 - tmp_42*tmp_96) + tmp_58*(-tmp_45*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_45*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_45*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_45*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_45*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_0_3 = tmp_49*(3.0*tmp_41*tmp_42*tmp_47 - tmp_41*tmp_45 - tmp_42*tmp_97) + tmp_58*(-tmp_45*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_45*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_45*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_45*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_45*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
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
                                  MatrixXr&                            elMat ) const override
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
      real_t tmp_19 = 0.091576213509770743*tmp_17 + 0.81684757298045851*tmp_18;
      real_t tmp_20 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_21 = tmp_15*(tmp_19 + tmp_20);
      real_t tmp_22 = -tmp_1*tmp_6 + tmp_10*tmp_3;
      real_t tmp_23 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_24 = -tmp_23;
      real_t tmp_25 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_26 = 0.091576213509770743*tmp_24 + 0.81684757298045851*tmp_25;
      real_t tmp_27 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_28 = tmp_15*(tmp_26 + tmp_27);
      real_t tmp_29 = -tmp_12 + tmp_7;
      real_t tmp_30 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_31 = -tmp_30;
      real_t tmp_32 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_33 = 0.091576213509770743*tmp_31 + 0.81684757298045851*tmp_32;
      real_t tmp_34 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_35 = tmp_15*(tmp_33 + tmp_34);
      real_t tmp_36 = -tmp_0*tmp_2 + tmp_3*tmp_9;
      real_t tmp_37 = tmp_0*tmp_6 - tmp_3*tmp_8;
      real_t tmp_38 = -tmp_13 + tmp_2*tmp_8;
      real_t tmp_39 = tmp_0*tmp_4 - tmp_1*tmp_9;
      real_t tmp_40 = -tmp_0*tmp_10 + tmp_1*tmp_8;
      real_t tmp_41 = tmp_11 - tmp_14;
      real_t tmp_42 = tmp_0*(tmp_21*tmp_5 + tmp_22*tmp_28 + tmp_29*tmp_35 - 1.0/4.0) + tmp_1*(tmp_21*tmp_36 + tmp_28*tmp_37 + tmp_35*tmp_38 - 1.0/4.0) + tmp_3*(tmp_21*tmp_39 + tmp_28*tmp_40 + tmp_35*tmp_41 - 1.0/4.0);
      real_t tmp_43 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_44 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_45 = tmp_43*tmp_44;
      real_t tmp_46 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_47 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_48 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_49 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_50 = tmp_46*tmp_49;
      real_t tmp_51 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_52 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_53 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_54 = tmp_49*tmp_51;
      real_t tmp_55 = tmp_43*tmp_52;
      real_t tmp_56 = tmp_47*tmp_53;
      real_t tmp_57 = 1.0 / (-tmp_44*tmp_54 + tmp_45*tmp_53 - tmp_46*tmp_56 + tmp_47*tmp_51*tmp_52 + tmp_48*tmp_50 - tmp_48*tmp_55);
      real_t tmp_58 = tmp_57*(tmp_45 - tmp_46*tmp_47);
      real_t tmp_59 = tmp_57*(-tmp_43*tmp_48 + tmp_47*tmp_51);
      real_t tmp_60 = tmp_57*(-tmp_44*tmp_51 + tmp_46*tmp_48);
      real_t tmp_61 = tmp_57*(-tmp_44*tmp_49 + tmp_47*tmp_52);
      real_t tmp_62 = tmp_57*(tmp_48*tmp_49 - tmp_56);
      real_t tmp_63 = tmp_57*(tmp_44*tmp_53 - tmp_48*tmp_52);
      real_t tmp_64 = tmp_57*(tmp_50 - tmp_55);
      real_t tmp_65 = tmp_57*(tmp_43*tmp_53 - tmp_54);
      real_t tmp_66 = tmp_57*(-tmp_46*tmp_53 + tmp_51*tmp_52);
      real_t tmp_67 = 0.5*p_affine_13_0*(-tmp_58 - tmp_59 - tmp_60) + 0.5*p_affine_13_1*(-tmp_61 - tmp_62 - tmp_63) + 0.5*p_affine_13_2*(-tmp_64 - tmp_65 - tmp_66);
      real_t tmp_68 = tmp_0*tmp_15;
      real_t tmp_69 = tmp_1*tmp_15;
      real_t tmp_70 = tmp_15*tmp_3;
      real_t tmp_71 = p_affine_13_0*(tmp_29*tmp_68 + tmp_38*tmp_69 + tmp_41*tmp_70) + p_affine_13_1*(tmp_22*tmp_68 + tmp_37*tmp_69 + tmp_40*tmp_70) + p_affine_13_2*(tmp_36*tmp_69 + tmp_39*tmp_70 + tmp_5*tmp_68);
      real_t tmp_72 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_73 = tmp_19 + tmp_72;
      real_t tmp_74 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_75 = tmp_26 + tmp_74;
      real_t tmp_76 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_77 = tmp_33 + tmp_76;
      real_t tmp_78 = tmp_60*tmp_77 + tmp_63*tmp_75 + tmp_66*tmp_73;
      real_t tmp_79 = tmp_59*tmp_77 + tmp_62*tmp_75 + tmp_65*tmp_73;
      real_t tmp_80 = tmp_58*tmp_77 + tmp_61*tmp_75 + tmp_64*tmp_73;
      real_t tmp_81 = -tmp_78 - tmp_79 - tmp_80 + 1;
      real_t tmp_82 = (std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)*std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)) + (std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)*std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)) + (std::abs(tmp_23*tmp_32 - tmp_25*tmp_30)*std::abs(tmp_23*tmp_32 - tmp_25*tmp_30));
      real_t tmp_83 = 3.0*std::pow(tmp_82, -0.25);
      real_t tmp_84 = tmp_42*tmp_83;
      real_t tmp_85 = 1.0*std::pow(tmp_82, 1.0/2.0);
      real_t tmp_86 = 0.054975871827660928*tmp_85;
      real_t tmp_87 = 0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18;
      real_t tmp_88 = tmp_15*(tmp_20 + tmp_87);
      real_t tmp_89 = 0.44594849091596489*tmp_24 + 0.10810301816807022*tmp_25;
      real_t tmp_90 = tmp_15*(tmp_27 + tmp_89);
      real_t tmp_91 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32;
      real_t tmp_92 = tmp_15*(tmp_34 + tmp_91);
      real_t tmp_93 = tmp_0*(tmp_22*tmp_90 + tmp_29*tmp_92 + tmp_5*tmp_88 - 1.0/4.0) + tmp_1*(tmp_36*tmp_88 + tmp_37*tmp_90 + tmp_38*tmp_92 - 1.0/4.0) + tmp_3*(tmp_39*tmp_88 + tmp_40*tmp_90 + tmp_41*tmp_92 - 1.0/4.0);
      real_t tmp_94 = tmp_72 + tmp_87;
      real_t tmp_95 = tmp_74 + tmp_89;
      real_t tmp_96 = tmp_76 + tmp_91;
      real_t tmp_97 = tmp_60*tmp_96 + tmp_63*tmp_95 + tmp_66*tmp_94;
      real_t tmp_98 = tmp_59*tmp_96 + tmp_62*tmp_95 + tmp_65*tmp_94;
      real_t tmp_99 = tmp_58*tmp_96 + tmp_61*tmp_95 + tmp_64*tmp_94;
      real_t tmp_100 = -tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_101 = tmp_83*tmp_93;
      real_t tmp_102 = 0.11169079483900572*tmp_85;
      real_t tmp_103 = 0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_104 = tmp_15*(tmp_103 + tmp_20);
      real_t tmp_105 = 0.81684757298045851*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_106 = tmp_15*(tmp_105 + tmp_27);
      real_t tmp_107 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_108 = tmp_15*(tmp_107 + tmp_34);
      real_t tmp_109 = tmp_0*(tmp_104*tmp_5 + tmp_106*tmp_22 + tmp_108*tmp_29 - 1.0/4.0) + tmp_1*(tmp_104*tmp_36 + tmp_106*tmp_37 + tmp_108*tmp_38 - 1.0/4.0) + tmp_3*(tmp_104*tmp_39 + tmp_106*tmp_40 + tmp_108*tmp_41 - 1.0/4.0);
      real_t tmp_110 = tmp_103 + tmp_72;
      real_t tmp_111 = tmp_105 + tmp_74;
      real_t tmp_112 = tmp_107 + tmp_76;
      real_t tmp_113 = tmp_110*tmp_66 + tmp_111*tmp_63 + tmp_112*tmp_60;
      real_t tmp_114 = tmp_110*tmp_65 + tmp_111*tmp_62 + tmp_112*tmp_59;
      real_t tmp_115 = tmp_110*tmp_64 + tmp_111*tmp_61 + tmp_112*tmp_58;
      real_t tmp_116 = -tmp_113 - tmp_114 - tmp_115 + 1;
      real_t tmp_117 = tmp_109*tmp_83;
      real_t tmp_118 = 0.054975871827660928*tmp_85;
      real_t tmp_119 = 0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_120 = tmp_15*(tmp_119 + tmp_20);
      real_t tmp_121 = 0.10810301816807022*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_122 = tmp_15*(tmp_121 + tmp_27);
      real_t tmp_123 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_124 = tmp_15*(tmp_123 + tmp_34);
      real_t tmp_125 = tmp_0*(tmp_120*tmp_5 + tmp_122*tmp_22 + tmp_124*tmp_29 - 1.0/4.0) + tmp_1*(tmp_120*tmp_36 + tmp_122*tmp_37 + tmp_124*tmp_38 - 1.0/4.0) + tmp_3*(tmp_120*tmp_39 + tmp_122*tmp_40 + tmp_124*tmp_41 - 1.0/4.0);
      real_t tmp_126 = tmp_119 + tmp_72;
      real_t tmp_127 = tmp_121 + tmp_74;
      real_t tmp_128 = tmp_123 + tmp_76;
      real_t tmp_129 = tmp_126*tmp_66 + tmp_127*tmp_63 + tmp_128*tmp_60;
      real_t tmp_130 = tmp_126*tmp_65 + tmp_127*tmp_62 + tmp_128*tmp_59;
      real_t tmp_131 = tmp_126*tmp_64 + tmp_127*tmp_61 + tmp_128*tmp_58;
      real_t tmp_132 = -tmp_129 - tmp_130 - tmp_131 + 1;
      real_t tmp_133 = tmp_125*tmp_83;
      real_t tmp_134 = 0.11169079483900572*tmp_85;
      real_t tmp_135 = 0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_136 = tmp_15*(tmp_135 + tmp_20);
      real_t tmp_137 = 0.091576213509770743*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_138 = tmp_15*(tmp_137 + tmp_27);
      real_t tmp_139 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_140 = tmp_15*(tmp_139 + tmp_34);
      real_t tmp_141 = tmp_0*(tmp_136*tmp_5 + tmp_138*tmp_22 + tmp_140*tmp_29 - 1.0/4.0) + tmp_1*(tmp_136*tmp_36 + tmp_138*tmp_37 + tmp_140*tmp_38 - 1.0/4.0) + tmp_3*(tmp_136*tmp_39 + tmp_138*tmp_40 + tmp_140*tmp_41 - 1.0/4.0);
      real_t tmp_142 = tmp_135 + tmp_72;
      real_t tmp_143 = tmp_137 + tmp_74;
      real_t tmp_144 = tmp_139 + tmp_76;
      real_t tmp_145 = tmp_142*tmp_66 + tmp_143*tmp_63 + tmp_144*tmp_60;
      real_t tmp_146 = tmp_142*tmp_65 + tmp_143*tmp_62 + tmp_144*tmp_59;
      real_t tmp_147 = tmp_142*tmp_64 + tmp_143*tmp_61 + tmp_144*tmp_58;
      real_t tmp_148 = -tmp_145 - tmp_146 - tmp_147 + 1;
      real_t tmp_149 = tmp_141*tmp_83;
      real_t tmp_150 = 0.054975871827660928*tmp_85;
      real_t tmp_151 = 0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_152 = tmp_15*(tmp_151 + tmp_20);
      real_t tmp_153 = 0.44594849091596489*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_154 = tmp_15*(tmp_153 + tmp_27);
      real_t tmp_155 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_156 = tmp_15*(tmp_155 + tmp_34);
      real_t tmp_157 = tmp_0*(tmp_152*tmp_5 + tmp_154*tmp_22 + tmp_156*tmp_29 - 1.0/4.0) + tmp_1*(tmp_152*tmp_36 + tmp_154*tmp_37 + tmp_156*tmp_38 - 1.0/4.0) + tmp_3*(tmp_152*tmp_39 + tmp_154*tmp_40 + tmp_156*tmp_41 - 1.0/4.0);
      real_t tmp_158 = tmp_151 + tmp_72;
      real_t tmp_159 = tmp_153 + tmp_74;
      real_t tmp_160 = tmp_155 + tmp_76;
      real_t tmp_161 = tmp_158*tmp_66 + tmp_159*tmp_63 + tmp_160*tmp_60;
      real_t tmp_162 = tmp_158*tmp_65 + tmp_159*tmp_62 + tmp_160*tmp_59;
      real_t tmp_163 = tmp_158*tmp_64 + tmp_159*tmp_61 + tmp_160*tmp_58;
      real_t tmp_164 = -tmp_161 - tmp_162 - tmp_163 + 1;
      real_t tmp_165 = tmp_157*tmp_83;
      real_t tmp_166 = 0.11169079483900572*tmp_85;
      real_t tmp_167 = 0.5*p_affine_13_0*tmp_60 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_66;
      real_t tmp_168 = 0.5*p_affine_13_0*tmp_59 + 0.5*p_affine_13_1*tmp_62 + 0.5*p_affine_13_2*tmp_65;
      real_t tmp_169 = 0.5*p_affine_13_0*tmp_58 + 0.5*p_affine_13_1*tmp_61 + 0.5*p_affine_13_2*tmp_64;
      real_t a_0_0 = tmp_102*(-tmp_100*tmp_101 + 0.5*tmp_100*tmp_71 - tmp_67*tmp_93) + tmp_118*(-tmp_109*tmp_67 - tmp_116*tmp_117 + 0.5*tmp_116*tmp_71) + tmp_134*(-tmp_125*tmp_67 - tmp_132*tmp_133 + 0.5*tmp_132*tmp_71) + tmp_150*(-tmp_141*tmp_67 - tmp_148*tmp_149 + 0.5*tmp_148*tmp_71) + tmp_166*(-tmp_157*tmp_67 - tmp_164*tmp_165 + 0.5*tmp_164*tmp_71) + tmp_86*(-tmp_42*tmp_67 + 0.5*tmp_71*tmp_81 - tmp_81*tmp_84);
      real_t a_0_1 = tmp_102*(-tmp_101*tmp_97 - tmp_167*tmp_93 + 0.5*tmp_71*tmp_97) + tmp_118*(-tmp_109*tmp_167 - tmp_113*tmp_117 + 0.5*tmp_113*tmp_71) + tmp_134*(-tmp_125*tmp_167 - tmp_129*tmp_133 + 0.5*tmp_129*tmp_71) + tmp_150*(-tmp_141*tmp_167 - tmp_145*tmp_149 + 0.5*tmp_145*tmp_71) + tmp_166*(-tmp_157*tmp_167 - tmp_161*tmp_165 + 0.5*tmp_161*tmp_71) + tmp_86*(-tmp_167*tmp_42 + 0.5*tmp_71*tmp_78 - tmp_78*tmp_84);
      real_t a_0_2 = tmp_102*(-tmp_101*tmp_98 - tmp_168*tmp_93 + 0.5*tmp_71*tmp_98) + tmp_118*(-tmp_109*tmp_168 - tmp_114*tmp_117 + 0.5*tmp_114*tmp_71) + tmp_134*(-tmp_125*tmp_168 - tmp_130*tmp_133 + 0.5*tmp_130*tmp_71) + tmp_150*(-tmp_141*tmp_168 - tmp_146*tmp_149 + 0.5*tmp_146*tmp_71) + tmp_166*(-tmp_157*tmp_168 - tmp_162*tmp_165 + 0.5*tmp_162*tmp_71) + tmp_86*(-tmp_168*tmp_42 + 0.5*tmp_71*tmp_79 - tmp_79*tmp_84);
      real_t a_0_3 = tmp_102*(-tmp_101*tmp_99 - tmp_169*tmp_93 + 0.5*tmp_71*tmp_99) + tmp_118*(-tmp_109*tmp_169 - tmp_115*tmp_117 + 0.5*tmp_115*tmp_71) + tmp_134*(-tmp_125*tmp_169 - tmp_131*tmp_133 + 0.5*tmp_131*tmp_71) + tmp_150*(-tmp_141*tmp_169 - tmp_147*tmp_149 + 0.5*tmp_147*tmp_71) + tmp_166*(-tmp_157*tmp_169 - tmp_163*tmp_165 + 0.5*tmp_163*tmp_71) + tmp_86*(-tmp_169*tmp_42 + 0.5*tmp_71*tmp_80 - tmp_80*tmp_84);
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
                                        MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_26 = -tmp_25;
      real_t tmp_27 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_28 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_29 = 0.091576213509770743*tmp_26 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -tmp_30;
      real_t tmp_32 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_33 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_34 = 0.091576213509770743*tmp_31 + 0.81684757298045851*tmp_32 + tmp_33;
      real_t tmp_35 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_36 = -tmp_35;
      real_t tmp_37 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_38 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_39 = 0.091576213509770743*tmp_36 + 0.81684757298045851*tmp_37 + tmp_38;
      real_t tmp_40 = tmp_17*tmp_39 + tmp_20*tmp_34 + tmp_23*tmp_29;
      real_t tmp_41 = tmp_16*tmp_39 + tmp_19*tmp_34 + tmp_22*tmp_29;
      real_t tmp_42 = tmp_15*tmp_39 + tmp_18*tmp_34 + tmp_21*tmp_29;
      real_t tmp_43 = tmp_11*(tmp_42 - 1.0/4.0) + tmp_6*(tmp_40 - 1.0/4.0) + tmp_9*(tmp_41 - 1.0/4.0);
      real_t tmp_44 = p_affine_13_0*(tmp_11*tmp_15 + tmp_16*tmp_9 + tmp_17*tmp_6) + p_affine_13_1*(tmp_11*tmp_18 + tmp_19*tmp_9 + tmp_20*tmp_6) + p_affine_13_2*(tmp_11*tmp_21 + tmp_22*tmp_9 + tmp_23*tmp_6);
      real_t tmp_45 = -tmp_40 - tmp_41 - tmp_42 + 1;
      real_t tmp_46 = (std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)*std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)) + (std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)*std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)) + (std::abs(tmp_30*tmp_37 - tmp_32*tmp_35)*std::abs(tmp_30*tmp_37 - tmp_32*tmp_35));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_26 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_51 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32 + tmp_33;
      real_t tmp_52 = 0.44594849091596489*tmp_36 + 0.10810301816807022*tmp_37 + tmp_38;
      real_t tmp_53 = tmp_17*tmp_52 + tmp_20*tmp_51 + tmp_23*tmp_50;
      real_t tmp_54 = tmp_16*tmp_52 + tmp_19*tmp_51 + tmp_22*tmp_50;
      real_t tmp_55 = tmp_15*tmp_52 + tmp_18*tmp_51 + tmp_21*tmp_50;
      real_t tmp_56 = tmp_11*(tmp_55 - 1.0/4.0) + tmp_6*(tmp_53 - 1.0/4.0) + tmp_9*(tmp_54 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_60 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_61 = 0.81684757298045851*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_62 = tmp_17*tmp_61 + tmp_20*tmp_60 + tmp_23*tmp_59;
      real_t tmp_63 = tmp_16*tmp_61 + tmp_19*tmp_60 + tmp_22*tmp_59;
      real_t tmp_64 = tmp_15*tmp_61 + tmp_18*tmp_60 + tmp_21*tmp_59;
      real_t tmp_65 = tmp_11*(tmp_64 - 1.0/4.0) + tmp_6*(tmp_62 - 1.0/4.0) + tmp_9*(tmp_63 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_69 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_70 = 0.10810301816807022*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_71 = tmp_17*tmp_70 + tmp_20*tmp_69 + tmp_23*tmp_68;
      real_t tmp_72 = tmp_16*tmp_70 + tmp_19*tmp_69 + tmp_22*tmp_68;
      real_t tmp_73 = tmp_15*tmp_70 + tmp_18*tmp_69 + tmp_21*tmp_68;
      real_t tmp_74 = tmp_11*(tmp_73 - 1.0/4.0) + tmp_6*(tmp_71 - 1.0/4.0) + tmp_9*(tmp_72 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_78 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_79 = 0.091576213509770743*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_80 = tmp_17*tmp_79 + tmp_20*tmp_78 + tmp_23*tmp_77;
      real_t tmp_81 = tmp_16*tmp_79 + tmp_19*tmp_78 + tmp_22*tmp_77;
      real_t tmp_82 = tmp_15*tmp_79 + tmp_18*tmp_78 + tmp_21*tmp_77;
      real_t tmp_83 = tmp_11*(tmp_82 - 1.0/4.0) + tmp_6*(tmp_80 - 1.0/4.0) + tmp_9*(tmp_81 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_87 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_88 = 0.44594849091596489*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_89 = tmp_17*tmp_88 + tmp_20*tmp_87 + tmp_23*tmp_86;
      real_t tmp_90 = tmp_16*tmp_88 + tmp_19*tmp_87 + tmp_22*tmp_86;
      real_t tmp_91 = tmp_15*tmp_88 + tmp_18*tmp_87 + tmp_21*tmp_86;
      real_t tmp_92 = tmp_11*(tmp_91 - 1.0/4.0) + tmp_6*(tmp_89 - 1.0/4.0) + tmp_9*(tmp_90 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_96 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_97 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t a_0_0 = tmp_49*(-tmp_24*tmp_43 + 3.0*tmp_43*tmp_45*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_24*tmp_56 - tmp_44*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_24*tmp_65 - tmp_44*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_24*tmp_74 - tmp_44*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_24*tmp_83 - tmp_44*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_24*tmp_92 - tmp_44*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_0_1 = tmp_49*(3.0*tmp_40*tmp_43*tmp_47 - tmp_40*tmp_44 - tmp_43*tmp_95) + tmp_58*(-tmp_44*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_44*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_44*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_44*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_44*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_0_2 = tmp_49*(3.0*tmp_41*tmp_43*tmp_47 - tmp_41*tmp_44 - tmp_43*tmp_96) + tmp_58*(-tmp_44*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_44*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_44*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_44*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_44*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_0_3 = tmp_49*(3.0*tmp_42*tmp_43*tmp_47 - tmp_42*tmp_44 - tmp_43*tmp_97) + tmp_58*(-tmp_44*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_44*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_44*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_44*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_44*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }

public:



};




class EGVectorLaplaceFormNitscheBC_EP1_1 : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_EP1_1()

    {}





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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_5 = -tmp_4;
      real_t tmp_6 = 1.0 / (tmp_2 + tmp_3*tmp_5);
      real_t tmp_7 = tmp_0*tmp_6;
      real_t tmp_8 = tmp_3*tmp_6;
      real_t tmp_9 = tmp_2*tmp_6 + tmp_5*tmp_8;
      real_t tmp_10 = tmp_1*tmp_6;
      real_t tmp_11 = tmp_4*tmp_6;
      real_t tmp_12 = tmp_1*tmp_11 + tmp_10*tmp_5;
      real_t tmp_13 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = tmp_13*(tmp_12*(-tmp_10 - tmp_11) + tmp_9*(-tmp_7 - tmp_8));
      real_t tmp_15 = tmp_13*(tmp_10*tmp_12 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_13*(tmp_11*tmp_12 + tmp_7*tmp_9);
      real_t a_0_0 = 0.5*tmp_14;
      real_t a_0_1 = 0.5*tmp_15;
      real_t a_0_2 = 0.5*tmp_16;
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
                                       MatrixXr&                                           elMat ) const override
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
      real_t tmp_2 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_4 = 0.21132486540518713*tmp_2 + tmp_3;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = 1.0 / (tmp_1*tmp_5 + tmp_8);
      real_t tmp_10 = tmp_5*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = 0.21132486540518713*tmp_11 + tmp_12;
      real_t tmp_14 = tmp_7*tmp_9;
      real_t tmp_15 = tmp_10*tmp_4 + tmp_13*tmp_14;
      real_t tmp_16 = tmp_6*tmp_9;
      real_t tmp_17 = tmp_0*tmp_9;
      real_t tmp_18 = tmp_13*tmp_17 + tmp_16*tmp_4;
      real_t tmp_19 = tmp_1*(tmp_15 - 1.0/3.0) + tmp_7*(tmp_18 - 1.0/3.0);
      real_t tmp_20 = 0.5*p_affine_10_0*(-tmp_14 - tmp_17) + 0.5*p_affine_10_1*(-tmp_10 - tmp_16);
      real_t tmp_21 = -tmp_15 - tmp_18 + 1;
      real_t tmp_22 = 0.5*p_affine_10_0*(tmp_1*tmp_14 + tmp_17*tmp_7) + 0.5*p_affine_10_1*(tmp_1*tmp_10 + tmp_8*tmp_9);
      real_t tmp_23 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_2*tmp_2), 1.0/2.0));
      real_t tmp_24 = 1.0 / (tmp_23);
      real_t tmp_25 = 0.5*tmp_23;
      real_t tmp_26 = 0.78867513459481287*tmp_2 + tmp_3;
      real_t tmp_27 = 0.78867513459481287*tmp_11 + tmp_12;
      real_t tmp_28 = tmp_10*tmp_26 + tmp_14*tmp_27;
      real_t tmp_29 = tmp_16*tmp_26 + tmp_17*tmp_27;
      real_t tmp_30 = tmp_1*(tmp_28 - 1.0/3.0) + tmp_7*(tmp_29 - 1.0/3.0);
      real_t tmp_31 = -tmp_28 - tmp_29 + 1;
      real_t tmp_32 = 0.5*tmp_23;
      real_t tmp_33 = 0.5*p_affine_10_0*tmp_14 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_34 = 0.5*p_affine_10_0*tmp_17 + 0.5*p_affine_10_1*tmp_16;
      real_t a_0_0 = tmp_25*(-tmp_19*tmp_20 + 3*tmp_19*tmp_21*tmp_24 - tmp_21*tmp_22) + tmp_32*(-tmp_20*tmp_30 - tmp_22*tmp_31 + 3*tmp_24*tmp_30*tmp_31);
      real_t a_0_1 = tmp_25*(3*tmp_15*tmp_19*tmp_24 - tmp_15*tmp_22 - tmp_19*tmp_33) + tmp_32*(-tmp_22*tmp_28 + 3*tmp_24*tmp_28*tmp_30 - tmp_30*tmp_33);
      real_t a_0_2 = tmp_25*(3*tmp_18*tmp_19*tmp_24 - tmp_18*tmp_22 - tmp_19*tmp_34) + tmp_32*(-tmp_22*tmp_29 + 3*tmp_24*tmp_29*tmp_30 - tmp_30*tmp_34);
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
                                          MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_4 = p_affine_6_1 + 0.21132486540518713*tmp_3;
      real_t tmp_5 = tmp_2 + tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = 1.0 / (tmp_1*tmp_9 + tmp_8);
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = -p_affine_0_0;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + 0.21132486540518713*tmp_13;
      real_t tmp_15 = tmp_12 + tmp_14;
      real_t tmp_16 = tmp_10*tmp_7;
      real_t tmp_17 = tmp_10*tmp_6;
      real_t tmp_18 = tmp_0*tmp_10;
      real_t tmp_19 = tmp_1*(tmp_11*tmp_5 + tmp_15*tmp_16 - 1.0/3.0) + tmp_7*(tmp_15*tmp_18 + tmp_17*tmp_5 - 1.0/3.0);
      real_t tmp_20 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_21 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_22 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_23 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_24 = 1.0 / (tmp_20*tmp_21 - tmp_22*tmp_23);
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = tmp_23*tmp_24;
      real_t tmp_27 = tmp_21*tmp_24;
      real_t tmp_28 = tmp_22*tmp_24;
      real_t tmp_29 = 0.5*p_affine_10_0*(-tmp_25 - tmp_26) + 0.5*p_affine_10_1*(-tmp_27 - tmp_28);
      real_t tmp_30 = p_affine_10_0*(tmp_0*tmp_16 + tmp_1*tmp_16) + p_affine_10_1*(tmp_1*tmp_11 + tmp_10*tmp_8);
      real_t tmp_31 = -p_affine_3_1;
      real_t tmp_32 = tmp_31 + tmp_4;
      real_t tmp_33 = -p_affine_3_0;
      real_t tmp_34 = tmp_14 + tmp_33;
      real_t tmp_35 = tmp_25*tmp_34 + tmp_28*tmp_32;
      real_t tmp_36 = tmp_26*tmp_34 + tmp_27*tmp_32;
      real_t tmp_37 = -tmp_35 - tmp_36 + 1;
      real_t tmp_38 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_3*tmp_3), 1.0/2.0));
      real_t tmp_39 = 3/tmp_38;
      real_t tmp_40 = tmp_19*tmp_39;
      real_t tmp_41 = 0.5*tmp_38;
      real_t tmp_42 = p_affine_6_1 + 0.78867513459481287*tmp_3;
      real_t tmp_43 = tmp_2 + tmp_42;
      real_t tmp_44 = p_affine_6_0 + 0.78867513459481287*tmp_13;
      real_t tmp_45 = tmp_12 + tmp_44;
      real_t tmp_46 = tmp_1*(tmp_11*tmp_43 + tmp_16*tmp_45 - 1.0/3.0) + tmp_7*(tmp_17*tmp_43 + tmp_18*tmp_45 - 1.0/3.0);
      real_t tmp_47 = tmp_31 + tmp_42;
      real_t tmp_48 = tmp_33 + tmp_44;
      real_t tmp_49 = tmp_25*tmp_48 + tmp_28*tmp_47;
      real_t tmp_50 = tmp_26*tmp_48 + tmp_27*tmp_47;
      real_t tmp_51 = -tmp_49 - tmp_50 + 1;
      real_t tmp_52 = tmp_39*tmp_46;
      real_t tmp_53 = 0.5*tmp_38;
      real_t tmp_54 = 0.5*p_affine_10_0*tmp_25 + 0.5*p_affine_10_1*tmp_28;
      real_t tmp_55 = 0.5*p_affine_10_0*tmp_26 + 0.5*p_affine_10_1*tmp_27;
      real_t a_0_0 = tmp_41*(-tmp_19*tmp_29 + 0.5*tmp_30*tmp_37 - tmp_37*tmp_40) + tmp_53*(-tmp_29*tmp_46 + 0.5*tmp_30*tmp_51 - tmp_51*tmp_52);
      real_t a_0_1 = tmp_41*(-tmp_19*tmp_54 + 0.5*tmp_30*tmp_35 - tmp_35*tmp_40) + tmp_53*(0.5*tmp_30*tmp_49 - tmp_46*tmp_54 - tmp_49*tmp_52);
      real_t a_0_2 = tmp_41*(-tmp_19*tmp_55 + 0.5*tmp_30*tmp_36 - tmp_36*tmp_40) + tmp_53*(0.5*tmp_30*tmp_50 - tmp_46*tmp_55 - tmp_50*tmp_52);
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
                                                   MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_4 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_5 = -tmp_4;
      real_t tmp_6 = 1.0 / (tmp_2 + tmp_3*tmp_5);
      real_t tmp_7 = tmp_0*tmp_6;
      real_t tmp_8 = tmp_4*tmp_6;
      real_t tmp_9 = tmp_1*tmp_6;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = p_affine_10_0*(-tmp_7 - tmp_8) + p_affine_10_1*(-tmp_10 - tmp_9);
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_14 = 0.21132486540518713*tmp_12 + tmp_13;
      real_t tmp_15 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_16 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_17 = 0.21132486540518713*tmp_15 + tmp_16;
      real_t tmp_18 = tmp_10*tmp_14 + tmp_17*tmp_7;
      real_t tmp_19 = tmp_14*tmp_9 + tmp_17*tmp_8;
      real_t tmp_20 = tmp_0*(tmp_19 - 1.0/3.0) + tmp_5*(tmp_18 - 1.0/3.0);
      real_t tmp_21 = p_affine_10_0*(tmp_0*tmp_8 + tmp_5*tmp_7) + p_affine_10_1*(tmp_10*tmp_5 + tmp_2*tmp_6);
      real_t tmp_22 = -tmp_18 - tmp_19 + 1;
      real_t tmp_23 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_15*tmp_15), 1.0/2.0));
      real_t tmp_24 = 1.0 / (tmp_23);
      real_t tmp_25 = 0.5*tmp_23;
      real_t tmp_26 = 0.78867513459481287*tmp_12 + tmp_13;
      real_t tmp_27 = 0.78867513459481287*tmp_15 + tmp_16;
      real_t tmp_28 = tmp_10*tmp_26 + tmp_27*tmp_7;
      real_t tmp_29 = tmp_26*tmp_9 + tmp_27*tmp_8;
      real_t tmp_30 = tmp_0*(tmp_29 - 1.0/3.0) + tmp_5*(tmp_28 - 1.0/3.0);
      real_t tmp_31 = -tmp_28 - tmp_29 + 1;
      real_t tmp_32 = 0.5*tmp_23;
      real_t tmp_33 = p_affine_10_0*tmp_7 + p_affine_10_1*tmp_10;
      real_t tmp_34 = p_affine_10_0*tmp_8 + p_affine_10_1*tmp_9;
      real_t a_0_0 = tmp_25*(-tmp_11*tmp_20 + 3*tmp_20*tmp_22*tmp_24 - tmp_21*tmp_22) + tmp_32*(-tmp_11*tmp_30 - tmp_21*tmp_31 + 3*tmp_24*tmp_30*tmp_31);
      real_t a_0_1 = tmp_25*(3*tmp_18*tmp_20*tmp_24 - tmp_18*tmp_21 - tmp_20*tmp_33) + tmp_32*(-tmp_21*tmp_28 + 3*tmp_24*tmp_28*tmp_30 - tmp_30*tmp_33);
      real_t a_0_2 = tmp_25*(3*tmp_19*tmp_20*tmp_24 - tmp_19*tmp_21 - tmp_20*tmp_34) + tmp_32*(-tmp_21*tmp_29 + 3*tmp_24*tmp_29*tmp_30 - tmp_30*tmp_34);
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
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_8 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_0*tmp_8;
      real_t tmp_13 = tmp_1*tmp_11;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_11*tmp_4 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_2*tmp_6 - tmp_5*tmp_6 + tmp_7*tmp_9);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_11*tmp_4 - tmp_12);
      real_t tmp_17 = tmp_14*(-tmp_13 + tmp_9);
      real_t tmp_18 = tmp_1*tmp_16 + tmp_15*tmp_8 + tmp_17*tmp_4;
      real_t tmp_19 = tmp_14*(-tmp_0*tmp_10 + tmp_3*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_0*tmp_6 - tmp_11*tmp_7);
      real_t tmp_21 = tmp_14*(tmp_10*tmp_11 - tmp_3*tmp_6);
      real_t tmp_22 = tmp_1*tmp_20 + tmp_19*tmp_8 + tmp_21*tmp_4;
      real_t tmp_23 = tmp_14*(-tmp_1*tmp_7 + tmp_10*tmp_4);
      real_t tmp_24 = tmp_14*(-tmp_4*tmp_6 + tmp_7*tmp_8);
      real_t tmp_25 = tmp_14*(tmp_1*tmp_6 - tmp_10*tmp_8);
      real_t tmp_26 = tmp_1*tmp_24 + tmp_23*tmp_8 + tmp_25*tmp_4;
      real_t tmp_27 = p_affine_0_0*p_affine_1_1;
      real_t tmp_28 = p_affine_0_0*p_affine_1_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_2;
      real_t tmp_30 = p_affine_0_1*p_affine_1_0;
      real_t tmp_31 = p_affine_0_1*p_affine_1_2;
      real_t tmp_32 = p_affine_2_2*p_affine_3_0;
      real_t tmp_33 = p_affine_0_2*p_affine_1_0;
      real_t tmp_34 = p_affine_0_2*p_affine_1_1;
      real_t tmp_35 = p_affine_2_0*p_affine_3_1;
      real_t tmp_36 = p_affine_2_2*p_affine_3_1;
      real_t tmp_37 = p_affine_2_0*p_affine_3_2;
      real_t tmp_38 = p_affine_2_1*p_affine_3_0;
      real_t tmp_39 = std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = tmp_39*(tmp_18*(-tmp_15 - tmp_16 - tmp_17) + tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25));
      real_t tmp_41 = tmp_39*(tmp_17*tmp_18 + tmp_21*tmp_22 + tmp_25*tmp_26);
      real_t tmp_42 = tmp_39*(tmp_16*tmp_18 + tmp_20*tmp_22 + tmp_24*tmp_26);
      real_t tmp_43 = tmp_39*(tmp_15*tmp_18 + tmp_19*tmp_22 + tmp_23*tmp_26);
      real_t a_0_0 = 0.16666666666666666*tmp_40;
      real_t a_0_1 = 0.16666666666666666*tmp_41;
      real_t a_0_2 = 0.16666666666666666*tmp_42;
      real_t a_0_3 = 0.16666666666666666*tmp_43;
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
                               MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_2 = -tmp_1;
      real_t tmp_3 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_4 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_5 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_3 + tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_7 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_8 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_9 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_10 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_11 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_12 = tmp_11*tmp_9;
      real_t tmp_13 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_14 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_15 = tmp_0*tmp_14;
      real_t tmp_16 = tmp_14*tmp_7;
      real_t tmp_17 = tmp_0*tmp_11;
      real_t tmp_18 = tmp_13*tmp_9;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_16 + tmp_13*tmp_6*tmp_7 + tmp_15*tmp_8 - tmp_17*tmp_6 - tmp_18*tmp_8);
      real_t tmp_20 = tmp_19*(tmp_6*tmp_7 - tmp_8*tmp_9);
      real_t tmp_21 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_22 = -tmp_21;
      real_t tmp_23 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_24 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_25 = 0.091576213509770743*tmp_22 + 0.81684757298045851*tmp_23 + tmp_24;
      real_t tmp_26 = tmp_19*(-tmp_11*tmp_6 + tmp_14*tmp_8);
      real_t tmp_27 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_28 = -tmp_27;
      real_t tmp_29 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_30 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_31 = 0.091576213509770743*tmp_28 + 0.81684757298045851*tmp_29 + tmp_30;
      real_t tmp_32 = tmp_19*(tmp_12 - tmp_16);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_25*tmp_26 + tmp_31*tmp_32;
      real_t tmp_34 = tmp_19*(tmp_0*tmp_8 - tmp_10*tmp_7);
      real_t tmp_35 = tmp_19*(tmp_10*tmp_11 - tmp_13*tmp_8);
      real_t tmp_36 = tmp_19*(tmp_13*tmp_7 - tmp_17);
      real_t tmp_37 = tmp_25*tmp_35 + tmp_31*tmp_36 + tmp_34*tmp_5;
      real_t tmp_38 = tmp_19*(-tmp_0*tmp_6 + tmp_10*tmp_9);
      real_t tmp_39 = tmp_19*(-tmp_10*tmp_14 + tmp_13*tmp_6);
      real_t tmp_40 = tmp_19*(tmp_15 - tmp_18);
      real_t tmp_41 = tmp_25*tmp_39 + tmp_31*tmp_40 + tmp_38*tmp_5;
      real_t tmp_42 = tmp_0*(tmp_33 - 1.0/4.0) + tmp_7*(tmp_41 - 1.0/4.0) + tmp_9*(tmp_37 - 1.0/4.0);
      real_t tmp_43 = 0.5*p_affine_13_0*(-tmp_32 - tmp_36 - tmp_40) + 0.5*p_affine_13_1*(-tmp_26 - tmp_35 - tmp_39) + 0.5*p_affine_13_2*(-tmp_20 - tmp_34 - tmp_38);
      real_t tmp_44 = -tmp_33 - tmp_37 - tmp_41 + 1;
      real_t tmp_45 = 0.5*p_affine_13_0*(tmp_0*tmp_32 + tmp_36*tmp_9 + tmp_40*tmp_7) + 0.5*p_affine_13_1*(tmp_0*tmp_26 + tmp_35*tmp_9 + tmp_39*tmp_7) + 0.5*p_affine_13_2*(tmp_0*tmp_20 + tmp_34*tmp_9 + tmp_38*tmp_7);
      real_t tmp_46 = (std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)*std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)) + (std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)*std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)) + (std::abs(tmp_21*tmp_29 - tmp_23*tmp_27)*std::abs(tmp_21*tmp_29 - tmp_23*tmp_27));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_3 + tmp_4;
      real_t tmp_51 = 0.44594849091596489*tmp_22 + 0.10810301816807022*tmp_23 + tmp_24;
      real_t tmp_52 = 0.44594849091596489*tmp_28 + 0.10810301816807022*tmp_29 + tmp_30;
      real_t tmp_53 = tmp_20*tmp_50 + tmp_26*tmp_51 + tmp_32*tmp_52;
      real_t tmp_54 = tmp_34*tmp_50 + tmp_35*tmp_51 + tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_50 + tmp_39*tmp_51 + tmp_40*tmp_52;
      real_t tmp_56 = tmp_0*(tmp_53 - 1.0/4.0) + tmp_7*(tmp_55 - 1.0/4.0) + tmp_9*(tmp_54 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_60 = 0.81684757298045851*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_61 = 0.81684757298045851*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_62 = tmp_20*tmp_59 + tmp_26*tmp_60 + tmp_32*tmp_61;
      real_t tmp_63 = tmp_34*tmp_59 + tmp_35*tmp_60 + tmp_36*tmp_61;
      real_t tmp_64 = tmp_38*tmp_59 + tmp_39*tmp_60 + tmp_40*tmp_61;
      real_t tmp_65 = tmp_0*(tmp_62 - 1.0/4.0) + tmp_7*(tmp_64 - 1.0/4.0) + tmp_9*(tmp_63 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_69 = 0.10810301816807022*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_70 = 0.10810301816807022*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_71 = tmp_20*tmp_68 + tmp_26*tmp_69 + tmp_32*tmp_70;
      real_t tmp_72 = tmp_34*tmp_68 + tmp_35*tmp_69 + tmp_36*tmp_70;
      real_t tmp_73 = tmp_38*tmp_68 + tmp_39*tmp_69 + tmp_40*tmp_70;
      real_t tmp_74 = tmp_0*(tmp_71 - 1.0/4.0) + tmp_7*(tmp_73 - 1.0/4.0) + tmp_9*(tmp_72 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_78 = 0.091576213509770743*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_79 = 0.091576213509770743*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_80 = tmp_20*tmp_77 + tmp_26*tmp_78 + tmp_32*tmp_79;
      real_t tmp_81 = tmp_34*tmp_77 + tmp_35*tmp_78 + tmp_36*tmp_79;
      real_t tmp_82 = tmp_38*tmp_77 + tmp_39*tmp_78 + tmp_40*tmp_79;
      real_t tmp_83 = tmp_0*(tmp_80 - 1.0/4.0) + tmp_7*(tmp_82 - 1.0/4.0) + tmp_9*(tmp_81 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_87 = 0.44594849091596489*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_88 = 0.44594849091596489*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_89 = tmp_20*tmp_86 + tmp_26*tmp_87 + tmp_32*tmp_88;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = tmp_0*(tmp_89 - 1.0/4.0) + tmp_7*(tmp_91 - 1.0/4.0) + tmp_9*(tmp_90 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = 0.5*p_affine_13_0*tmp_32 + 0.5*p_affine_13_1*tmp_26 + 0.5*p_affine_13_2*tmp_20;
      real_t tmp_96 = 0.5*p_affine_13_0*tmp_36 + 0.5*p_affine_13_1*tmp_35 + 0.5*p_affine_13_2*tmp_34;
      real_t tmp_97 = 0.5*p_affine_13_0*tmp_40 + 0.5*p_affine_13_1*tmp_39 + 0.5*p_affine_13_2*tmp_38;
      real_t a_0_0 = tmp_49*(-tmp_42*tmp_43 + 3.0*tmp_42*tmp_44*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_43*tmp_56 - tmp_45*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_43*tmp_65 - tmp_45*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_43*tmp_74 - tmp_45*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_43*tmp_83 - tmp_45*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_43*tmp_92 - tmp_45*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_0_1 = tmp_49*(3.0*tmp_33*tmp_42*tmp_47 - tmp_33*tmp_45 - tmp_42*tmp_95) + tmp_58*(-tmp_45*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_45*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_45*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_45*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_45*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_0_2 = tmp_49*(3.0*tmp_37*tmp_42*tmp_47 - tmp_37*tmp_45 - tmp_42*tmp_96) + tmp_58*(-tmp_45*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_45*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_45*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_45*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_45*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_0_3 = tmp_49*(3.0*tmp_41*tmp_42*tmp_47 - tmp_41*tmp_45 - tmp_42*tmp_97) + tmp_58*(-tmp_45*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_45*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_45*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_45*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_45*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
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
                                  MatrixXr&                            elMat ) const override
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
      real_t tmp_36 = tmp_0*tmp_3 - tmp_2*tmp_6;
      real_t tmp_37 = -tmp_3*tmp_9 + tmp_6*tmp_7;
      real_t tmp_38 = -tmp_13 + tmp_2*tmp_9;
      real_t tmp_39 = -tmp_0*tmp_1 + tmp_4*tmp_6;
      real_t tmp_40 = tmp_1*tmp_9 - tmp_10*tmp_6;
      real_t tmp_41 = tmp_11 - tmp_14;
      real_t tmp_42 = tmp_0*(tmp_21*tmp_5 + tmp_22*tmp_28 + tmp_29*tmp_35 - 1.0/4.0) + tmp_2*(tmp_21*tmp_39 + tmp_28*tmp_40 + tmp_35*tmp_41 - 1.0/4.0) + tmp_4*(tmp_21*tmp_36 + tmp_28*tmp_37 + tmp_35*tmp_38 - 1.0/4.0);
      real_t tmp_43 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_44 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_45 = tmp_43*tmp_44;
      real_t tmp_46 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_47 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_48 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_49 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_50 = tmp_46*tmp_49;
      real_t tmp_51 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_52 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_53 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_54 = tmp_49*tmp_51;
      real_t tmp_55 = tmp_43*tmp_52;
      real_t tmp_56 = tmp_47*tmp_53;
      real_t tmp_57 = 1.0 / (-tmp_44*tmp_54 + tmp_45*tmp_53 - tmp_46*tmp_56 + tmp_47*tmp_51*tmp_52 + tmp_48*tmp_50 - tmp_48*tmp_55);
      real_t tmp_58 = tmp_57*(tmp_45 - tmp_46*tmp_47);
      real_t tmp_59 = tmp_57*(-tmp_43*tmp_48 + tmp_47*tmp_51);
      real_t tmp_60 = tmp_57*(-tmp_44*tmp_51 + tmp_46*tmp_48);
      real_t tmp_61 = tmp_57*(-tmp_44*tmp_49 + tmp_47*tmp_52);
      real_t tmp_62 = tmp_57*(tmp_48*tmp_49 - tmp_56);
      real_t tmp_63 = tmp_57*(tmp_44*tmp_53 - tmp_48*tmp_52);
      real_t tmp_64 = tmp_57*(tmp_50 - tmp_55);
      real_t tmp_65 = tmp_57*(tmp_43*tmp_53 - tmp_54);
      real_t tmp_66 = tmp_57*(-tmp_46*tmp_53 + tmp_51*tmp_52);
      real_t tmp_67 = 0.5*p_affine_13_0*(-tmp_58 - tmp_59 - tmp_60) + 0.5*p_affine_13_1*(-tmp_61 - tmp_62 - tmp_63) + 0.5*p_affine_13_2*(-tmp_64 - tmp_65 - tmp_66);
      real_t tmp_68 = tmp_0*tmp_15;
      real_t tmp_69 = tmp_15*tmp_4;
      real_t tmp_70 = tmp_15*tmp_2;
      real_t tmp_71 = p_affine_13_0*(tmp_29*tmp_68 + tmp_38*tmp_69 + tmp_41*tmp_70) + p_affine_13_1*(tmp_22*tmp_68 + tmp_37*tmp_69 + tmp_40*tmp_70) + p_affine_13_2*(tmp_36*tmp_69 + tmp_39*tmp_70 + tmp_5*tmp_68);
      real_t tmp_72 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_73 = tmp_19 + tmp_72;
      real_t tmp_74 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_75 = tmp_26 + tmp_74;
      real_t tmp_76 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_77 = tmp_33 + tmp_76;
      real_t tmp_78 = tmp_60*tmp_77 + tmp_63*tmp_75 + tmp_66*tmp_73;
      real_t tmp_79 = tmp_59*tmp_77 + tmp_62*tmp_75 + tmp_65*tmp_73;
      real_t tmp_80 = tmp_58*tmp_77 + tmp_61*tmp_75 + tmp_64*tmp_73;
      real_t tmp_81 = -tmp_78 - tmp_79 - tmp_80 + 1;
      real_t tmp_82 = (std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)*std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)) + (std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)*std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)) + (std::abs(tmp_23*tmp_32 - tmp_25*tmp_30)*std::abs(tmp_23*tmp_32 - tmp_25*tmp_30));
      real_t tmp_83 = 3.0*std::pow(tmp_82, -0.25);
      real_t tmp_84 = tmp_42*tmp_83;
      real_t tmp_85 = 1.0*std::pow(tmp_82, 1.0/2.0);
      real_t tmp_86 = 0.054975871827660928*tmp_85;
      real_t tmp_87 = 0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18;
      real_t tmp_88 = tmp_15*(tmp_20 + tmp_87);
      real_t tmp_89 = 0.44594849091596489*tmp_24 + 0.10810301816807022*tmp_25;
      real_t tmp_90 = tmp_15*(tmp_27 + tmp_89);
      real_t tmp_91 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32;
      real_t tmp_92 = tmp_15*(tmp_34 + tmp_91);
      real_t tmp_93 = tmp_0*(tmp_22*tmp_90 + tmp_29*tmp_92 + tmp_5*tmp_88 - 1.0/4.0) + tmp_2*(tmp_39*tmp_88 + tmp_40*tmp_90 + tmp_41*tmp_92 - 1.0/4.0) + tmp_4*(tmp_36*tmp_88 + tmp_37*tmp_90 + tmp_38*tmp_92 - 1.0/4.0);
      real_t tmp_94 = tmp_72 + tmp_87;
      real_t tmp_95 = tmp_74 + tmp_89;
      real_t tmp_96 = tmp_76 + tmp_91;
      real_t tmp_97 = tmp_60*tmp_96 + tmp_63*tmp_95 + tmp_66*tmp_94;
      real_t tmp_98 = tmp_59*tmp_96 + tmp_62*tmp_95 + tmp_65*tmp_94;
      real_t tmp_99 = tmp_58*tmp_96 + tmp_61*tmp_95 + tmp_64*tmp_94;
      real_t tmp_100 = -tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_101 = tmp_83*tmp_93;
      real_t tmp_102 = 0.11169079483900572*tmp_85;
      real_t tmp_103 = 0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_104 = tmp_15*(tmp_103 + tmp_20);
      real_t tmp_105 = 0.81684757298045851*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_106 = tmp_15*(tmp_105 + tmp_27);
      real_t tmp_107 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_108 = tmp_15*(tmp_107 + tmp_34);
      real_t tmp_109 = tmp_0*(tmp_104*tmp_5 + tmp_106*tmp_22 + tmp_108*tmp_29 - 1.0/4.0) + tmp_2*(tmp_104*tmp_39 + tmp_106*tmp_40 + tmp_108*tmp_41 - 1.0/4.0) + tmp_4*(tmp_104*tmp_36 + tmp_106*tmp_37 + tmp_108*tmp_38 - 1.0/4.0);
      real_t tmp_110 = tmp_103 + tmp_72;
      real_t tmp_111 = tmp_105 + tmp_74;
      real_t tmp_112 = tmp_107 + tmp_76;
      real_t tmp_113 = tmp_110*tmp_66 + tmp_111*tmp_63 + tmp_112*tmp_60;
      real_t tmp_114 = tmp_110*tmp_65 + tmp_111*tmp_62 + tmp_112*tmp_59;
      real_t tmp_115 = tmp_110*tmp_64 + tmp_111*tmp_61 + tmp_112*tmp_58;
      real_t tmp_116 = -tmp_113 - tmp_114 - tmp_115 + 1;
      real_t tmp_117 = tmp_109*tmp_83;
      real_t tmp_118 = 0.054975871827660928*tmp_85;
      real_t tmp_119 = 0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_120 = tmp_15*(tmp_119 + tmp_20);
      real_t tmp_121 = 0.10810301816807022*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_122 = tmp_15*(tmp_121 + tmp_27);
      real_t tmp_123 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_124 = tmp_15*(tmp_123 + tmp_34);
      real_t tmp_125 = tmp_0*(tmp_120*tmp_5 + tmp_122*tmp_22 + tmp_124*tmp_29 - 1.0/4.0) + tmp_2*(tmp_120*tmp_39 + tmp_122*tmp_40 + tmp_124*tmp_41 - 1.0/4.0) + tmp_4*(tmp_120*tmp_36 + tmp_122*tmp_37 + tmp_124*tmp_38 - 1.0/4.0);
      real_t tmp_126 = tmp_119 + tmp_72;
      real_t tmp_127 = tmp_121 + tmp_74;
      real_t tmp_128 = tmp_123 + tmp_76;
      real_t tmp_129 = tmp_126*tmp_66 + tmp_127*tmp_63 + tmp_128*tmp_60;
      real_t tmp_130 = tmp_126*tmp_65 + tmp_127*tmp_62 + tmp_128*tmp_59;
      real_t tmp_131 = tmp_126*tmp_64 + tmp_127*tmp_61 + tmp_128*tmp_58;
      real_t tmp_132 = -tmp_129 - tmp_130 - tmp_131 + 1;
      real_t tmp_133 = tmp_125*tmp_83;
      real_t tmp_134 = 0.11169079483900572*tmp_85;
      real_t tmp_135 = 0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_136 = tmp_15*(tmp_135 + tmp_20);
      real_t tmp_137 = 0.091576213509770743*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_138 = tmp_15*(tmp_137 + tmp_27);
      real_t tmp_139 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_140 = tmp_15*(tmp_139 + tmp_34);
      real_t tmp_141 = tmp_0*(tmp_136*tmp_5 + tmp_138*tmp_22 + tmp_140*tmp_29 - 1.0/4.0) + tmp_2*(tmp_136*tmp_39 + tmp_138*tmp_40 + tmp_140*tmp_41 - 1.0/4.0) + tmp_4*(tmp_136*tmp_36 + tmp_138*tmp_37 + tmp_140*tmp_38 - 1.0/4.0);
      real_t tmp_142 = tmp_135 + tmp_72;
      real_t tmp_143 = tmp_137 + tmp_74;
      real_t tmp_144 = tmp_139 + tmp_76;
      real_t tmp_145 = tmp_142*tmp_66 + tmp_143*tmp_63 + tmp_144*tmp_60;
      real_t tmp_146 = tmp_142*tmp_65 + tmp_143*tmp_62 + tmp_144*tmp_59;
      real_t tmp_147 = tmp_142*tmp_64 + tmp_143*tmp_61 + tmp_144*tmp_58;
      real_t tmp_148 = -tmp_145 - tmp_146 - tmp_147 + 1;
      real_t tmp_149 = tmp_141*tmp_83;
      real_t tmp_150 = 0.054975871827660928*tmp_85;
      real_t tmp_151 = 0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_152 = tmp_15*(tmp_151 + tmp_20);
      real_t tmp_153 = 0.44594849091596489*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_154 = tmp_15*(tmp_153 + tmp_27);
      real_t tmp_155 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_156 = tmp_15*(tmp_155 + tmp_34);
      real_t tmp_157 = tmp_0*(tmp_152*tmp_5 + tmp_154*tmp_22 + tmp_156*tmp_29 - 1.0/4.0) + tmp_2*(tmp_152*tmp_39 + tmp_154*tmp_40 + tmp_156*tmp_41 - 1.0/4.0) + tmp_4*(tmp_152*tmp_36 + tmp_154*tmp_37 + tmp_156*tmp_38 - 1.0/4.0);
      real_t tmp_158 = tmp_151 + tmp_72;
      real_t tmp_159 = tmp_153 + tmp_74;
      real_t tmp_160 = tmp_155 + tmp_76;
      real_t tmp_161 = tmp_158*tmp_66 + tmp_159*tmp_63 + tmp_160*tmp_60;
      real_t tmp_162 = tmp_158*tmp_65 + tmp_159*tmp_62 + tmp_160*tmp_59;
      real_t tmp_163 = tmp_158*tmp_64 + tmp_159*tmp_61 + tmp_160*tmp_58;
      real_t tmp_164 = -tmp_161 - tmp_162 - tmp_163 + 1;
      real_t tmp_165 = tmp_157*tmp_83;
      real_t tmp_166 = 0.11169079483900572*tmp_85;
      real_t tmp_167 = 0.5*p_affine_13_0*tmp_60 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_66;
      real_t tmp_168 = 0.5*p_affine_13_0*tmp_59 + 0.5*p_affine_13_1*tmp_62 + 0.5*p_affine_13_2*tmp_65;
      real_t tmp_169 = 0.5*p_affine_13_0*tmp_58 + 0.5*p_affine_13_1*tmp_61 + 0.5*p_affine_13_2*tmp_64;
      real_t a_0_0 = tmp_102*(-tmp_100*tmp_101 + 0.5*tmp_100*tmp_71 - tmp_67*tmp_93) + tmp_118*(-tmp_109*tmp_67 - tmp_116*tmp_117 + 0.5*tmp_116*tmp_71) + tmp_134*(-tmp_125*tmp_67 - tmp_132*tmp_133 + 0.5*tmp_132*tmp_71) + tmp_150*(-tmp_141*tmp_67 - tmp_148*tmp_149 + 0.5*tmp_148*tmp_71) + tmp_166*(-tmp_157*tmp_67 - tmp_164*tmp_165 + 0.5*tmp_164*tmp_71) + tmp_86*(-tmp_42*tmp_67 + 0.5*tmp_71*tmp_81 - tmp_81*tmp_84);
      real_t a_0_1 = tmp_102*(-tmp_101*tmp_97 - tmp_167*tmp_93 + 0.5*tmp_71*tmp_97) + tmp_118*(-tmp_109*tmp_167 - tmp_113*tmp_117 + 0.5*tmp_113*tmp_71) + tmp_134*(-tmp_125*tmp_167 - tmp_129*tmp_133 + 0.5*tmp_129*tmp_71) + tmp_150*(-tmp_141*tmp_167 - tmp_145*tmp_149 + 0.5*tmp_145*tmp_71) + tmp_166*(-tmp_157*tmp_167 - tmp_161*tmp_165 + 0.5*tmp_161*tmp_71) + tmp_86*(-tmp_167*tmp_42 + 0.5*tmp_71*tmp_78 - tmp_78*tmp_84);
      real_t a_0_2 = tmp_102*(-tmp_101*tmp_98 - tmp_168*tmp_93 + 0.5*tmp_71*tmp_98) + tmp_118*(-tmp_109*tmp_168 - tmp_114*tmp_117 + 0.5*tmp_114*tmp_71) + tmp_134*(-tmp_125*tmp_168 - tmp_130*tmp_133 + 0.5*tmp_130*tmp_71) + tmp_150*(-tmp_141*tmp_168 - tmp_146*tmp_149 + 0.5*tmp_146*tmp_71) + tmp_166*(-tmp_157*tmp_168 - tmp_162*tmp_165 + 0.5*tmp_162*tmp_71) + tmp_86*(-tmp_168*tmp_42 + 0.5*tmp_71*tmp_79 - tmp_79*tmp_84);
      real_t a_0_3 = tmp_102*(-tmp_101*tmp_99 - tmp_169*tmp_93 + 0.5*tmp_71*tmp_99) + tmp_118*(-tmp_109*tmp_169 - tmp_115*tmp_117 + 0.5*tmp_115*tmp_71) + tmp_134*(-tmp_125*tmp_169 - tmp_131*tmp_133 + 0.5*tmp_131*tmp_71) + tmp_150*(-tmp_141*tmp_169 - tmp_147*tmp_149 + 0.5*tmp_147*tmp_71) + tmp_166*(-tmp_157*tmp_169 - tmp_163*tmp_165 + 0.5*tmp_163*tmp_71) + tmp_86*(-tmp_169*tmp_42 + 0.5*tmp_71*tmp_80 - tmp_80*tmp_84);
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
                                        MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_26 = -tmp_25;
      real_t tmp_27 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_28 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_29 = 0.091576213509770743*tmp_26 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -tmp_30;
      real_t tmp_32 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_33 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_34 = 0.091576213509770743*tmp_31 + 0.81684757298045851*tmp_32 + tmp_33;
      real_t tmp_35 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_36 = -tmp_35;
      real_t tmp_37 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_38 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_39 = 0.091576213509770743*tmp_36 + 0.81684757298045851*tmp_37 + tmp_38;
      real_t tmp_40 = tmp_17*tmp_39 + tmp_20*tmp_34 + tmp_23*tmp_29;
      real_t tmp_41 = tmp_16*tmp_39 + tmp_19*tmp_34 + tmp_22*tmp_29;
      real_t tmp_42 = tmp_15*tmp_39 + tmp_18*tmp_34 + tmp_21*tmp_29;
      real_t tmp_43 = tmp_0*(tmp_40 - 1.0/4.0) + tmp_10*(tmp_42 - 1.0/4.0) + tmp_3*(tmp_41 - 1.0/4.0);
      real_t tmp_44 = p_affine_13_0*(tmp_0*tmp_17 + tmp_10*tmp_15 + tmp_16*tmp_3) + p_affine_13_1*(tmp_0*tmp_20 + tmp_10*tmp_18 + tmp_19*tmp_3) + p_affine_13_2*(tmp_0*tmp_23 + tmp_10*tmp_21 + tmp_22*tmp_3);
      real_t tmp_45 = -tmp_40 - tmp_41 - tmp_42 + 1;
      real_t tmp_46 = (std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)*std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)) + (std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)*std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)) + (std::abs(tmp_30*tmp_37 - tmp_32*tmp_35)*std::abs(tmp_30*tmp_37 - tmp_32*tmp_35));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_26 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_51 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32 + tmp_33;
      real_t tmp_52 = 0.44594849091596489*tmp_36 + 0.10810301816807022*tmp_37 + tmp_38;
      real_t tmp_53 = tmp_17*tmp_52 + tmp_20*tmp_51 + tmp_23*tmp_50;
      real_t tmp_54 = tmp_16*tmp_52 + tmp_19*tmp_51 + tmp_22*tmp_50;
      real_t tmp_55 = tmp_15*tmp_52 + tmp_18*tmp_51 + tmp_21*tmp_50;
      real_t tmp_56 = tmp_0*(tmp_53 - 1.0/4.0) + tmp_10*(tmp_55 - 1.0/4.0) + tmp_3*(tmp_54 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_60 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_61 = 0.81684757298045851*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_62 = tmp_17*tmp_61 + tmp_20*tmp_60 + tmp_23*tmp_59;
      real_t tmp_63 = tmp_16*tmp_61 + tmp_19*tmp_60 + tmp_22*tmp_59;
      real_t tmp_64 = tmp_15*tmp_61 + tmp_18*tmp_60 + tmp_21*tmp_59;
      real_t tmp_65 = tmp_0*(tmp_62 - 1.0/4.0) + tmp_10*(tmp_64 - 1.0/4.0) + tmp_3*(tmp_63 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_69 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_70 = 0.10810301816807022*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_71 = tmp_17*tmp_70 + tmp_20*tmp_69 + tmp_23*tmp_68;
      real_t tmp_72 = tmp_16*tmp_70 + tmp_19*tmp_69 + tmp_22*tmp_68;
      real_t tmp_73 = tmp_15*tmp_70 + tmp_18*tmp_69 + tmp_21*tmp_68;
      real_t tmp_74 = tmp_0*(tmp_71 - 1.0/4.0) + tmp_10*(tmp_73 - 1.0/4.0) + tmp_3*(tmp_72 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_78 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_79 = 0.091576213509770743*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_80 = tmp_17*tmp_79 + tmp_20*tmp_78 + tmp_23*tmp_77;
      real_t tmp_81 = tmp_16*tmp_79 + tmp_19*tmp_78 + tmp_22*tmp_77;
      real_t tmp_82 = tmp_15*tmp_79 + tmp_18*tmp_78 + tmp_21*tmp_77;
      real_t tmp_83 = tmp_0*(tmp_80 - 1.0/4.0) + tmp_10*(tmp_82 - 1.0/4.0) + tmp_3*(tmp_81 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_87 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_88 = 0.44594849091596489*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_89 = tmp_17*tmp_88 + tmp_20*tmp_87 + tmp_23*tmp_86;
      real_t tmp_90 = tmp_16*tmp_88 + tmp_19*tmp_87 + tmp_22*tmp_86;
      real_t tmp_91 = tmp_15*tmp_88 + tmp_18*tmp_87 + tmp_21*tmp_86;
      real_t tmp_92 = tmp_0*(tmp_89 - 1.0/4.0) + tmp_10*(tmp_91 - 1.0/4.0) + tmp_3*(tmp_90 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_96 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_97 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t a_0_0 = tmp_49*(-tmp_24*tmp_43 + 3.0*tmp_43*tmp_45*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_24*tmp_56 - tmp_44*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_24*tmp_65 - tmp_44*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_24*tmp_74 - tmp_44*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_24*tmp_83 - tmp_44*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_24*tmp_92 - tmp_44*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_0_1 = tmp_49*(3.0*tmp_40*tmp_43*tmp_47 - tmp_40*tmp_44 - tmp_43*tmp_95) + tmp_58*(-tmp_44*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_44*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_44*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_44*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_44*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_0_2 = tmp_49*(3.0*tmp_41*tmp_43*tmp_47 - tmp_41*tmp_44 - tmp_43*tmp_96) + tmp_58*(-tmp_44*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_44*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_44*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_44*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_44*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_0_3 = tmp_49*(3.0*tmp_42*tmp_43*tmp_47 - tmp_42*tmp_44 - tmp_43*tmp_97) + tmp_58*(-tmp_44*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_44*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_44*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_44*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_44*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }

public:



};




class EGVectorLaplaceFormNitscheBC_EP1_2 : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_EP1_2()

    {}





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
                                       MatrixXr&                                           elMat ) const override
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
                                          MatrixXr&                                           elMat ) const override
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
                                                   MatrixXr&                                           elMat ) const override
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
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_8 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_0*tmp_8;
      real_t tmp_13 = tmp_1*tmp_11;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_11*tmp_4 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_2*tmp_6 - tmp_5*tmp_6 + tmp_7*tmp_9);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_11*tmp_4 - tmp_12);
      real_t tmp_17 = tmp_14*(-tmp_13 + tmp_9);
      real_t tmp_18 = tmp_10*tmp_16 + tmp_15*tmp_6 + tmp_17*tmp_7;
      real_t tmp_19 = tmp_14*(-tmp_0*tmp_10 + tmp_3*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_0*tmp_6 - tmp_11*tmp_7);
      real_t tmp_21 = tmp_14*(tmp_10*tmp_11 - tmp_3*tmp_6);
      real_t tmp_22 = tmp_10*tmp_20 + tmp_19*tmp_6 + tmp_21*tmp_7;
      real_t tmp_23 = tmp_14*(-tmp_1*tmp_7 + tmp_10*tmp_4);
      real_t tmp_24 = tmp_14*(-tmp_4*tmp_6 + tmp_7*tmp_8);
      real_t tmp_25 = tmp_14*(tmp_1*tmp_6 - tmp_10*tmp_8);
      real_t tmp_26 = tmp_10*tmp_24 + tmp_23*tmp_6 + tmp_25*tmp_7;
      real_t tmp_27 = p_affine_0_0*p_affine_1_1;
      real_t tmp_28 = p_affine_0_0*p_affine_1_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_2;
      real_t tmp_30 = p_affine_0_1*p_affine_1_0;
      real_t tmp_31 = p_affine_0_1*p_affine_1_2;
      real_t tmp_32 = p_affine_2_2*p_affine_3_0;
      real_t tmp_33 = p_affine_0_2*p_affine_1_0;
      real_t tmp_34 = p_affine_0_2*p_affine_1_1;
      real_t tmp_35 = p_affine_2_0*p_affine_3_1;
      real_t tmp_36 = p_affine_2_2*p_affine_3_1;
      real_t tmp_37 = p_affine_2_0*p_affine_3_2;
      real_t tmp_38 = p_affine_2_1*p_affine_3_0;
      real_t tmp_39 = std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = tmp_39*(tmp_18*(-tmp_15 - tmp_16 - tmp_17) + tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25));
      real_t tmp_41 = tmp_39*(tmp_17*tmp_18 + tmp_21*tmp_22 + tmp_25*tmp_26);
      real_t tmp_42 = tmp_39*(tmp_16*tmp_18 + tmp_20*tmp_22 + tmp_24*tmp_26);
      real_t tmp_43 = tmp_39*(tmp_15*tmp_18 + tmp_19*tmp_22 + tmp_23*tmp_26);
      real_t a_0_0 = 0.16666666666666666*tmp_40;
      real_t a_0_1 = 0.16666666666666666*tmp_41;
      real_t a_0_2 = 0.16666666666666666*tmp_42;
      real_t a_0_3 = 0.16666666666666666*tmp_43;
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
                               MatrixXr&                            elMat ) const override
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

         real_t tmp_0 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_1 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_2 = -tmp_1;
      real_t tmp_3 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_4 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_5 = 0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_3 + tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_7 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_8 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_9 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_10 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_11 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_12 = tmp_11*tmp_9;
      real_t tmp_13 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_14 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_14*tmp_7;
      real_t tmp_17 = tmp_11*tmp_13;
      real_t tmp_18 = tmp_0*tmp_9;
      real_t tmp_19 = 1.0 / (tmp_0*tmp_6*tmp_7 + tmp_10*tmp_12 - tmp_10*tmp_16 + tmp_15*tmp_8 - tmp_17*tmp_6 - tmp_18*tmp_8);
      real_t tmp_20 = tmp_19*(tmp_6*tmp_7 - tmp_8*tmp_9);
      real_t tmp_21 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_22 = -tmp_21;
      real_t tmp_23 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_24 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_25 = 0.091576213509770743*tmp_22 + 0.81684757298045851*tmp_23 + tmp_24;
      real_t tmp_26 = tmp_19*(-tmp_11*tmp_6 + tmp_14*tmp_8);
      real_t tmp_27 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_28 = -tmp_27;
      real_t tmp_29 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_30 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_31 = 0.091576213509770743*tmp_28 + 0.81684757298045851*tmp_29 + tmp_30;
      real_t tmp_32 = tmp_19*(tmp_12 - tmp_16);
      real_t tmp_33 = tmp_20*tmp_5 + tmp_25*tmp_26 + tmp_31*tmp_32;
      real_t tmp_34 = tmp_19*(-tmp_10*tmp_7 + tmp_13*tmp_8);
      real_t tmp_35 = tmp_19*(-tmp_0*tmp_8 + tmp_10*tmp_11);
      real_t tmp_36 = tmp_19*(tmp_0*tmp_7 - tmp_17);
      real_t tmp_37 = tmp_25*tmp_35 + tmp_31*tmp_36 + tmp_34*tmp_5;
      real_t tmp_38 = tmp_19*(tmp_10*tmp_9 - tmp_13*tmp_6);
      real_t tmp_39 = tmp_19*(tmp_0*tmp_6 - tmp_10*tmp_14);
      real_t tmp_40 = tmp_19*(tmp_15 - tmp_18);
      real_t tmp_41 = tmp_25*tmp_39 + tmp_31*tmp_40 + tmp_38*tmp_5;
      real_t tmp_42 = tmp_0*(tmp_33 - 1.0/4.0) + tmp_11*(tmp_41 - 1.0/4.0) + tmp_14*(tmp_37 - 1.0/4.0);
      real_t tmp_43 = 0.5*p_affine_13_0*(-tmp_32 - tmp_36 - tmp_40) + 0.5*p_affine_13_1*(-tmp_26 - tmp_35 - tmp_39) + 0.5*p_affine_13_2*(-tmp_20 - tmp_34 - tmp_38);
      real_t tmp_44 = -tmp_33 - tmp_37 - tmp_41 + 1;
      real_t tmp_45 = 0.5*p_affine_13_0*(tmp_0*tmp_32 + tmp_11*tmp_40 + tmp_14*tmp_36) + 0.5*p_affine_13_1*(tmp_0*tmp_26 + tmp_11*tmp_39 + tmp_14*tmp_35) + 0.5*p_affine_13_2*(tmp_0*tmp_20 + tmp_11*tmp_38 + tmp_14*tmp_34);
      real_t tmp_46 = (std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)*std::abs(tmp_1*tmp_23 - tmp_21*tmp_3)) + (std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)*std::abs(tmp_1*tmp_29 - tmp_27*tmp_3)) + (std::abs(tmp_21*tmp_29 - tmp_23*tmp_27)*std::abs(tmp_21*tmp_29 - tmp_23*tmp_27));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_3 + tmp_4;
      real_t tmp_51 = 0.44594849091596489*tmp_22 + 0.10810301816807022*tmp_23 + tmp_24;
      real_t tmp_52 = 0.44594849091596489*tmp_28 + 0.10810301816807022*tmp_29 + tmp_30;
      real_t tmp_53 = tmp_20*tmp_50 + tmp_26*tmp_51 + tmp_32*tmp_52;
      real_t tmp_54 = tmp_34*tmp_50 + tmp_35*tmp_51 + tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_50 + tmp_39*tmp_51 + tmp_40*tmp_52;
      real_t tmp_56 = tmp_0*(tmp_53 - 1.0/4.0) + tmp_11*(tmp_55 - 1.0/4.0) + tmp_14*(tmp_54 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_60 = 0.81684757298045851*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_61 = 0.81684757298045851*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_62 = tmp_20*tmp_59 + tmp_26*tmp_60 + tmp_32*tmp_61;
      real_t tmp_63 = tmp_34*tmp_59 + tmp_35*tmp_60 + tmp_36*tmp_61;
      real_t tmp_64 = tmp_38*tmp_59 + tmp_39*tmp_60 + tmp_40*tmp_61;
      real_t tmp_65 = tmp_0*(tmp_62 - 1.0/4.0) + tmp_11*(tmp_64 - 1.0/4.0) + tmp_14*(tmp_63 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_69 = 0.10810301816807022*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_70 = 0.10810301816807022*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_71 = tmp_20*tmp_68 + tmp_26*tmp_69 + tmp_32*tmp_70;
      real_t tmp_72 = tmp_34*tmp_68 + tmp_35*tmp_69 + tmp_36*tmp_70;
      real_t tmp_73 = tmp_38*tmp_68 + tmp_39*tmp_69 + tmp_40*tmp_70;
      real_t tmp_74 = tmp_0*(tmp_71 - 1.0/4.0) + tmp_11*(tmp_73 - 1.0/4.0) + tmp_14*(tmp_72 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_3 + tmp_4;
      real_t tmp_78 = 0.091576213509770743*tmp_22 + 0.091576213509770743*tmp_23 + tmp_24;
      real_t tmp_79 = 0.091576213509770743*tmp_28 + 0.091576213509770743*tmp_29 + tmp_30;
      real_t tmp_80 = tmp_20*tmp_77 + tmp_26*tmp_78 + tmp_32*tmp_79;
      real_t tmp_81 = tmp_34*tmp_77 + tmp_35*tmp_78 + tmp_36*tmp_79;
      real_t tmp_82 = tmp_38*tmp_77 + tmp_39*tmp_78 + tmp_40*tmp_79;
      real_t tmp_83 = tmp_0*(tmp_80 - 1.0/4.0) + tmp_11*(tmp_82 - 1.0/4.0) + tmp_14*(tmp_81 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_3 + tmp_4;
      real_t tmp_87 = 0.44594849091596489*tmp_22 + 0.44594849091596489*tmp_23 + tmp_24;
      real_t tmp_88 = 0.44594849091596489*tmp_28 + 0.44594849091596489*tmp_29 + tmp_30;
      real_t tmp_89 = tmp_20*tmp_86 + tmp_26*tmp_87 + tmp_32*tmp_88;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = tmp_0*(tmp_89 - 1.0/4.0) + tmp_11*(tmp_91 - 1.0/4.0) + tmp_14*(tmp_90 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = 0.5*p_affine_13_0*tmp_32 + 0.5*p_affine_13_1*tmp_26 + 0.5*p_affine_13_2*tmp_20;
      real_t tmp_96 = 0.5*p_affine_13_0*tmp_36 + 0.5*p_affine_13_1*tmp_35 + 0.5*p_affine_13_2*tmp_34;
      real_t tmp_97 = 0.5*p_affine_13_0*tmp_40 + 0.5*p_affine_13_1*tmp_39 + 0.5*p_affine_13_2*tmp_38;
      real_t a_0_0 = tmp_49*(-tmp_42*tmp_43 + 3.0*tmp_42*tmp_44*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_43*tmp_56 - tmp_45*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_43*tmp_65 - tmp_45*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_43*tmp_74 - tmp_45*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_43*tmp_83 - tmp_45*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_43*tmp_92 - tmp_45*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_0_1 = tmp_49*(3.0*tmp_33*tmp_42*tmp_47 - tmp_33*tmp_45 - tmp_42*tmp_95) + tmp_58*(-tmp_45*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_45*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_45*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_45*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_45*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_0_2 = tmp_49*(3.0*tmp_37*tmp_42*tmp_47 - tmp_37*tmp_45 - tmp_42*tmp_96) + tmp_58*(-tmp_45*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_45*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_45*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_45*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_45*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_0_3 = tmp_49*(3.0*tmp_41*tmp_42*tmp_47 - tmp_41*tmp_45 - tmp_42*tmp_97) + tmp_58*(-tmp_45*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_45*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_45*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_45*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_45*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
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
                                  MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_1*tmp_2 - tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_10*tmp_2;
      real_t tmp_13 = tmp_7*tmp_9;
      real_t tmp_14 = tmp_0*tmp_4;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_2 - tmp_1*tmp_13 + tmp_11*tmp_3 - tmp_12*tmp_6 - tmp_14*tmp_3 + tmp_6*tmp_8);
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
      real_t tmp_36 = -tmp_2*tmp_6 + tmp_3*tmp_9;
      real_t tmp_37 = -tmp_0*tmp_3 + tmp_6*tmp_7;
      real_t tmp_38 = tmp_0*tmp_2 - tmp_13;
      real_t tmp_39 = -tmp_1*tmp_9 + tmp_4*tmp_6;
      real_t tmp_40 = tmp_0*tmp_1 - tmp_10*tmp_6;
      real_t tmp_41 = tmp_11 - tmp_14;
      real_t tmp_42 = tmp_0*(tmp_21*tmp_5 + tmp_22*tmp_28 + tmp_29*tmp_35 - 1.0/4.0) + tmp_10*(tmp_21*tmp_36 + tmp_28*tmp_37 + tmp_35*tmp_38 - 1.0/4.0) + tmp_7*(tmp_21*tmp_39 + tmp_28*tmp_40 + tmp_35*tmp_41 - 1.0/4.0);
      real_t tmp_43 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_44 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_45 = tmp_43*tmp_44;
      real_t tmp_46 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_47 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_48 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_49 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_50 = tmp_46*tmp_49;
      real_t tmp_51 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_52 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_53 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_54 = tmp_49*tmp_51;
      real_t tmp_55 = tmp_43*tmp_52;
      real_t tmp_56 = tmp_47*tmp_53;
      real_t tmp_57 = 1.0 / (-tmp_44*tmp_54 + tmp_45*tmp_53 - tmp_46*tmp_56 + tmp_47*tmp_51*tmp_52 + tmp_48*tmp_50 - tmp_48*tmp_55);
      real_t tmp_58 = tmp_57*(tmp_45 - tmp_46*tmp_47);
      real_t tmp_59 = tmp_57*(-tmp_43*tmp_48 + tmp_47*tmp_51);
      real_t tmp_60 = tmp_57*(-tmp_44*tmp_51 + tmp_46*tmp_48);
      real_t tmp_61 = tmp_57*(-tmp_44*tmp_49 + tmp_47*tmp_52);
      real_t tmp_62 = tmp_57*(tmp_48*tmp_49 - tmp_56);
      real_t tmp_63 = tmp_57*(tmp_44*tmp_53 - tmp_48*tmp_52);
      real_t tmp_64 = tmp_57*(tmp_50 - tmp_55);
      real_t tmp_65 = tmp_57*(tmp_43*tmp_53 - tmp_54);
      real_t tmp_66 = tmp_57*(-tmp_46*tmp_53 + tmp_51*tmp_52);
      real_t tmp_67 = 0.5*p_affine_13_0*(-tmp_58 - tmp_59 - tmp_60) + 0.5*p_affine_13_1*(-tmp_61 - tmp_62 - tmp_63) + 0.5*p_affine_13_2*(-tmp_64 - tmp_65 - tmp_66);
      real_t tmp_68 = tmp_0*tmp_15;
      real_t tmp_69 = tmp_10*tmp_15;
      real_t tmp_70 = tmp_15*tmp_7;
      real_t tmp_71 = p_affine_13_0*(tmp_29*tmp_68 + tmp_38*tmp_69 + tmp_41*tmp_70) + p_affine_13_1*(tmp_22*tmp_68 + tmp_37*tmp_69 + tmp_40*tmp_70) + p_affine_13_2*(tmp_36*tmp_69 + tmp_39*tmp_70 + tmp_5*tmp_68);
      real_t tmp_72 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_73 = tmp_19 + tmp_72;
      real_t tmp_74 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_75 = tmp_26 + tmp_74;
      real_t tmp_76 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_77 = tmp_33 + tmp_76;
      real_t tmp_78 = tmp_60*tmp_77 + tmp_63*tmp_75 + tmp_66*tmp_73;
      real_t tmp_79 = tmp_59*tmp_77 + tmp_62*tmp_75 + tmp_65*tmp_73;
      real_t tmp_80 = tmp_58*tmp_77 + tmp_61*tmp_75 + tmp_64*tmp_73;
      real_t tmp_81 = -tmp_78 - tmp_79 - tmp_80 + 1;
      real_t tmp_82 = (std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)*std::abs(tmp_16*tmp_25 - tmp_18*tmp_23)) + (std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)*std::abs(tmp_16*tmp_32 - tmp_18*tmp_30)) + (std::abs(tmp_23*tmp_32 - tmp_25*tmp_30)*std::abs(tmp_23*tmp_32 - tmp_25*tmp_30));
      real_t tmp_83 = 3.0*std::pow(tmp_82, -0.25);
      real_t tmp_84 = tmp_42*tmp_83;
      real_t tmp_85 = 1.0*std::pow(tmp_82, 1.0/2.0);
      real_t tmp_86 = 0.054975871827660928*tmp_85;
      real_t tmp_87 = 0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18;
      real_t tmp_88 = tmp_15*(tmp_20 + tmp_87);
      real_t tmp_89 = 0.44594849091596489*tmp_24 + 0.10810301816807022*tmp_25;
      real_t tmp_90 = tmp_15*(tmp_27 + tmp_89);
      real_t tmp_91 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32;
      real_t tmp_92 = tmp_15*(tmp_34 + tmp_91);
      real_t tmp_93 = tmp_0*(tmp_22*tmp_90 + tmp_29*tmp_92 + tmp_5*tmp_88 - 1.0/4.0) + tmp_10*(tmp_36*tmp_88 + tmp_37*tmp_90 + tmp_38*tmp_92 - 1.0/4.0) + tmp_7*(tmp_39*tmp_88 + tmp_40*tmp_90 + tmp_41*tmp_92 - 1.0/4.0);
      real_t tmp_94 = tmp_72 + tmp_87;
      real_t tmp_95 = tmp_74 + tmp_89;
      real_t tmp_96 = tmp_76 + tmp_91;
      real_t tmp_97 = tmp_60*tmp_96 + tmp_63*tmp_95 + tmp_66*tmp_94;
      real_t tmp_98 = tmp_59*tmp_96 + tmp_62*tmp_95 + tmp_65*tmp_94;
      real_t tmp_99 = tmp_58*tmp_96 + tmp_61*tmp_95 + tmp_64*tmp_94;
      real_t tmp_100 = -tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_101 = tmp_83*tmp_93;
      real_t tmp_102 = 0.11169079483900572*tmp_85;
      real_t tmp_103 = 0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_104 = tmp_15*(tmp_103 + tmp_20);
      real_t tmp_105 = 0.81684757298045851*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_106 = tmp_15*(tmp_105 + tmp_27);
      real_t tmp_107 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_108 = tmp_15*(tmp_107 + tmp_34);
      real_t tmp_109 = tmp_0*(tmp_104*tmp_5 + tmp_106*tmp_22 + tmp_108*tmp_29 - 1.0/4.0) + tmp_10*(tmp_104*tmp_36 + tmp_106*tmp_37 + tmp_108*tmp_38 - 1.0/4.0) + tmp_7*(tmp_104*tmp_39 + tmp_106*tmp_40 + tmp_108*tmp_41 - 1.0/4.0);
      real_t tmp_110 = tmp_103 + tmp_72;
      real_t tmp_111 = tmp_105 + tmp_74;
      real_t tmp_112 = tmp_107 + tmp_76;
      real_t tmp_113 = tmp_110*tmp_66 + tmp_111*tmp_63 + tmp_112*tmp_60;
      real_t tmp_114 = tmp_110*tmp_65 + tmp_111*tmp_62 + tmp_112*tmp_59;
      real_t tmp_115 = tmp_110*tmp_64 + tmp_111*tmp_61 + tmp_112*tmp_58;
      real_t tmp_116 = -tmp_113 - tmp_114 - tmp_115 + 1;
      real_t tmp_117 = tmp_109*tmp_83;
      real_t tmp_118 = 0.054975871827660928*tmp_85;
      real_t tmp_119 = 0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_120 = tmp_15*(tmp_119 + tmp_20);
      real_t tmp_121 = 0.10810301816807022*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_122 = tmp_15*(tmp_121 + tmp_27);
      real_t tmp_123 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_124 = tmp_15*(tmp_123 + tmp_34);
      real_t tmp_125 = tmp_0*(tmp_120*tmp_5 + tmp_122*tmp_22 + tmp_124*tmp_29 - 1.0/4.0) + tmp_10*(tmp_120*tmp_36 + tmp_122*tmp_37 + tmp_124*tmp_38 - 1.0/4.0) + tmp_7*(tmp_120*tmp_39 + tmp_122*tmp_40 + tmp_124*tmp_41 - 1.0/4.0);
      real_t tmp_126 = tmp_119 + tmp_72;
      real_t tmp_127 = tmp_121 + tmp_74;
      real_t tmp_128 = tmp_123 + tmp_76;
      real_t tmp_129 = tmp_126*tmp_66 + tmp_127*tmp_63 + tmp_128*tmp_60;
      real_t tmp_130 = tmp_126*tmp_65 + tmp_127*tmp_62 + tmp_128*tmp_59;
      real_t tmp_131 = tmp_126*tmp_64 + tmp_127*tmp_61 + tmp_128*tmp_58;
      real_t tmp_132 = -tmp_129 - tmp_130 - tmp_131 + 1;
      real_t tmp_133 = tmp_125*tmp_83;
      real_t tmp_134 = 0.11169079483900572*tmp_85;
      real_t tmp_135 = 0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18;
      real_t tmp_136 = tmp_15*(tmp_135 + tmp_20);
      real_t tmp_137 = 0.091576213509770743*tmp_24 + 0.091576213509770743*tmp_25;
      real_t tmp_138 = tmp_15*(tmp_137 + tmp_27);
      real_t tmp_139 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32;
      real_t tmp_140 = tmp_15*(tmp_139 + tmp_34);
      real_t tmp_141 = tmp_0*(tmp_136*tmp_5 + tmp_138*tmp_22 + tmp_140*tmp_29 - 1.0/4.0) + tmp_10*(tmp_136*tmp_36 + tmp_138*tmp_37 + tmp_140*tmp_38 - 1.0/4.0) + tmp_7*(tmp_136*tmp_39 + tmp_138*tmp_40 + tmp_140*tmp_41 - 1.0/4.0);
      real_t tmp_142 = tmp_135 + tmp_72;
      real_t tmp_143 = tmp_137 + tmp_74;
      real_t tmp_144 = tmp_139 + tmp_76;
      real_t tmp_145 = tmp_142*tmp_66 + tmp_143*tmp_63 + tmp_144*tmp_60;
      real_t tmp_146 = tmp_142*tmp_65 + tmp_143*tmp_62 + tmp_144*tmp_59;
      real_t tmp_147 = tmp_142*tmp_64 + tmp_143*tmp_61 + tmp_144*tmp_58;
      real_t tmp_148 = -tmp_145 - tmp_146 - tmp_147 + 1;
      real_t tmp_149 = tmp_141*tmp_83;
      real_t tmp_150 = 0.054975871827660928*tmp_85;
      real_t tmp_151 = 0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_152 = tmp_15*(tmp_151 + tmp_20);
      real_t tmp_153 = 0.44594849091596489*tmp_24 + 0.44594849091596489*tmp_25;
      real_t tmp_154 = tmp_15*(tmp_153 + tmp_27);
      real_t tmp_155 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32;
      real_t tmp_156 = tmp_15*(tmp_155 + tmp_34);
      real_t tmp_157 = tmp_0*(tmp_152*tmp_5 + tmp_154*tmp_22 + tmp_156*tmp_29 - 1.0/4.0) + tmp_10*(tmp_152*tmp_36 + tmp_154*tmp_37 + tmp_156*tmp_38 - 1.0/4.0) + tmp_7*(tmp_152*tmp_39 + tmp_154*tmp_40 + tmp_156*tmp_41 - 1.0/4.0);
      real_t tmp_158 = tmp_151 + tmp_72;
      real_t tmp_159 = tmp_153 + tmp_74;
      real_t tmp_160 = tmp_155 + tmp_76;
      real_t tmp_161 = tmp_158*tmp_66 + tmp_159*tmp_63 + tmp_160*tmp_60;
      real_t tmp_162 = tmp_158*tmp_65 + tmp_159*tmp_62 + tmp_160*tmp_59;
      real_t tmp_163 = tmp_158*tmp_64 + tmp_159*tmp_61 + tmp_160*tmp_58;
      real_t tmp_164 = -tmp_161 - tmp_162 - tmp_163 + 1;
      real_t tmp_165 = tmp_157*tmp_83;
      real_t tmp_166 = 0.11169079483900572*tmp_85;
      real_t tmp_167 = 0.5*p_affine_13_0*tmp_60 + 0.5*p_affine_13_1*tmp_63 + 0.5*p_affine_13_2*tmp_66;
      real_t tmp_168 = 0.5*p_affine_13_0*tmp_59 + 0.5*p_affine_13_1*tmp_62 + 0.5*p_affine_13_2*tmp_65;
      real_t tmp_169 = 0.5*p_affine_13_0*tmp_58 + 0.5*p_affine_13_1*tmp_61 + 0.5*p_affine_13_2*tmp_64;
      real_t a_0_0 = tmp_102*(-tmp_100*tmp_101 + 0.5*tmp_100*tmp_71 - tmp_67*tmp_93) + tmp_118*(-tmp_109*tmp_67 - tmp_116*tmp_117 + 0.5*tmp_116*tmp_71) + tmp_134*(-tmp_125*tmp_67 - tmp_132*tmp_133 + 0.5*tmp_132*tmp_71) + tmp_150*(-tmp_141*tmp_67 - tmp_148*tmp_149 + 0.5*tmp_148*tmp_71) + tmp_166*(-tmp_157*tmp_67 - tmp_164*tmp_165 + 0.5*tmp_164*tmp_71) + tmp_86*(-tmp_42*tmp_67 + 0.5*tmp_71*tmp_81 - tmp_81*tmp_84);
      real_t a_0_1 = tmp_102*(-tmp_101*tmp_97 - tmp_167*tmp_93 + 0.5*tmp_71*tmp_97) + tmp_118*(-tmp_109*tmp_167 - tmp_113*tmp_117 + 0.5*tmp_113*tmp_71) + tmp_134*(-tmp_125*tmp_167 - tmp_129*tmp_133 + 0.5*tmp_129*tmp_71) + tmp_150*(-tmp_141*tmp_167 - tmp_145*tmp_149 + 0.5*tmp_145*tmp_71) + tmp_166*(-tmp_157*tmp_167 - tmp_161*tmp_165 + 0.5*tmp_161*tmp_71) + tmp_86*(-tmp_167*tmp_42 + 0.5*tmp_71*tmp_78 - tmp_78*tmp_84);
      real_t a_0_2 = tmp_102*(-tmp_101*tmp_98 - tmp_168*tmp_93 + 0.5*tmp_71*tmp_98) + tmp_118*(-tmp_109*tmp_168 - tmp_114*tmp_117 + 0.5*tmp_114*tmp_71) + tmp_134*(-tmp_125*tmp_168 - tmp_130*tmp_133 + 0.5*tmp_130*tmp_71) + tmp_150*(-tmp_141*tmp_168 - tmp_146*tmp_149 + 0.5*tmp_146*tmp_71) + tmp_166*(-tmp_157*tmp_168 - tmp_162*tmp_165 + 0.5*tmp_162*tmp_71) + tmp_86*(-tmp_168*tmp_42 + 0.5*tmp_71*tmp_79 - tmp_79*tmp_84);
      real_t a_0_3 = tmp_102*(-tmp_101*tmp_99 - tmp_169*tmp_93 + 0.5*tmp_71*tmp_99) + tmp_118*(-tmp_109*tmp_169 - tmp_115*tmp_117 + 0.5*tmp_115*tmp_71) + tmp_134*(-tmp_125*tmp_169 - tmp_131*tmp_133 + 0.5*tmp_131*tmp_71) + tmp_150*(-tmp_141*tmp_169 - tmp_147*tmp_149 + 0.5*tmp_147*tmp_71) + tmp_166*(-tmp_157*tmp_169 - tmp_163*tmp_165 + 0.5*tmp_163*tmp_71) + tmp_86*(-tmp_169*tmp_42 + 0.5*tmp_71*tmp_80 - tmp_80*tmp_84);
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
                                        MatrixXr&                            elMat ) const override
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
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = tmp_14*(-tmp_1*tmp_6 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_14*(-tmp_11*tmp_4 + tmp_6*tmp_7);
      real_t tmp_20 = tmp_14*(tmp_1*tmp_11 - tmp_7*tmp_9);
      real_t tmp_21 = tmp_14*(-tmp_0*tmp_9 + tmp_3*tmp_6);
      real_t tmp_22 = tmp_14*(tmp_0*tmp_11 - tmp_10*tmp_6);
      real_t tmp_23 = tmp_14*(tmp_10*tmp_9 - tmp_11*tmp_3);
      real_t tmp_24 = p_affine_13_0*(-tmp_15 - tmp_16 - tmp_17) + p_affine_13_1*(-tmp_18 - tmp_19 - tmp_20) + p_affine_13_2*(-tmp_21 - tmp_22 - tmp_23);
      real_t tmp_25 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_26 = -tmp_25;
      real_t tmp_27 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_28 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_29 = 0.091576213509770743*tmp_26 + 0.81684757298045851*tmp_27 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -tmp_30;
      real_t tmp_32 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_33 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_34 = 0.091576213509770743*tmp_31 + 0.81684757298045851*tmp_32 + tmp_33;
      real_t tmp_35 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_36 = -tmp_35;
      real_t tmp_37 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_38 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_39 = 0.091576213509770743*tmp_36 + 0.81684757298045851*tmp_37 + tmp_38;
      real_t tmp_40 = tmp_17*tmp_39 + tmp_20*tmp_34 + tmp_23*tmp_29;
      real_t tmp_41 = tmp_16*tmp_39 + tmp_19*tmp_34 + tmp_22*tmp_29;
      real_t tmp_42 = tmp_15*tmp_39 + tmp_18*tmp_34 + tmp_21*tmp_29;
      real_t tmp_43 = tmp_1*(tmp_41 - 1.0/4.0) + tmp_4*(tmp_40 - 1.0/4.0) + tmp_7*(tmp_42 - 1.0/4.0);
      real_t tmp_44 = p_affine_13_0*(tmp_1*tmp_16 + tmp_15*tmp_7 + tmp_17*tmp_4) + p_affine_13_1*(tmp_1*tmp_19 + tmp_18*tmp_7 + tmp_20*tmp_4) + p_affine_13_2*(tmp_1*tmp_22 + tmp_21*tmp_7 + tmp_23*tmp_4);
      real_t tmp_45 = -tmp_40 - tmp_41 - tmp_42 + 1;
      real_t tmp_46 = (std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)*std::abs(tmp_25*tmp_32 - tmp_27*tmp_30)) + (std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)*std::abs(tmp_25*tmp_37 - tmp_27*tmp_35)) + (std::abs(tmp_30*tmp_37 - tmp_32*tmp_35)*std::abs(tmp_30*tmp_37 - tmp_32*tmp_35));
      real_t tmp_47 = std::pow(tmp_46, -0.25);
      real_t tmp_48 = 1.0*std::pow(tmp_46, 1.0/2.0);
      real_t tmp_49 = 0.054975871827660928*tmp_48;
      real_t tmp_50 = 0.44594849091596489*tmp_26 + 0.10810301816807022*tmp_27 + tmp_28;
      real_t tmp_51 = 0.44594849091596489*tmp_31 + 0.10810301816807022*tmp_32 + tmp_33;
      real_t tmp_52 = 0.44594849091596489*tmp_36 + 0.10810301816807022*tmp_37 + tmp_38;
      real_t tmp_53 = tmp_17*tmp_52 + tmp_20*tmp_51 + tmp_23*tmp_50;
      real_t tmp_54 = tmp_16*tmp_52 + tmp_19*tmp_51 + tmp_22*tmp_50;
      real_t tmp_55 = tmp_15*tmp_52 + tmp_18*tmp_51 + tmp_21*tmp_50;
      real_t tmp_56 = tmp_1*(tmp_54 - 1.0/4.0) + tmp_4*(tmp_53 - 1.0/4.0) + tmp_7*(tmp_55 - 1.0/4.0);
      real_t tmp_57 = -tmp_53 - tmp_54 - tmp_55 + 1;
      real_t tmp_58 = 0.11169079483900572*tmp_48;
      real_t tmp_59 = 0.81684757298045851*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_60 = 0.81684757298045851*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_61 = 0.81684757298045851*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_62 = tmp_17*tmp_61 + tmp_20*tmp_60 + tmp_23*tmp_59;
      real_t tmp_63 = tmp_16*tmp_61 + tmp_19*tmp_60 + tmp_22*tmp_59;
      real_t tmp_64 = tmp_15*tmp_61 + tmp_18*tmp_60 + tmp_21*tmp_59;
      real_t tmp_65 = tmp_1*(tmp_63 - 1.0/4.0) + tmp_4*(tmp_62 - 1.0/4.0) + tmp_7*(tmp_64 - 1.0/4.0);
      real_t tmp_66 = -tmp_62 - tmp_63 - tmp_64 + 1;
      real_t tmp_67 = 0.054975871827660928*tmp_48;
      real_t tmp_68 = 0.10810301816807022*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_69 = 0.10810301816807022*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_70 = 0.10810301816807022*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_71 = tmp_17*tmp_70 + tmp_20*tmp_69 + tmp_23*tmp_68;
      real_t tmp_72 = tmp_16*tmp_70 + tmp_19*tmp_69 + tmp_22*tmp_68;
      real_t tmp_73 = tmp_15*tmp_70 + tmp_18*tmp_69 + tmp_21*tmp_68;
      real_t tmp_74 = tmp_1*(tmp_72 - 1.0/4.0) + tmp_4*(tmp_71 - 1.0/4.0) + tmp_7*(tmp_73 - 1.0/4.0);
      real_t tmp_75 = -tmp_71 - tmp_72 - tmp_73 + 1;
      real_t tmp_76 = 0.11169079483900572*tmp_48;
      real_t tmp_77 = 0.091576213509770743*tmp_26 + 0.091576213509770743*tmp_27 + tmp_28;
      real_t tmp_78 = 0.091576213509770743*tmp_31 + 0.091576213509770743*tmp_32 + tmp_33;
      real_t tmp_79 = 0.091576213509770743*tmp_36 + 0.091576213509770743*tmp_37 + tmp_38;
      real_t tmp_80 = tmp_17*tmp_79 + tmp_20*tmp_78 + tmp_23*tmp_77;
      real_t tmp_81 = tmp_16*tmp_79 + tmp_19*tmp_78 + tmp_22*tmp_77;
      real_t tmp_82 = tmp_15*tmp_79 + tmp_18*tmp_78 + tmp_21*tmp_77;
      real_t tmp_83 = tmp_1*(tmp_81 - 1.0/4.0) + tmp_4*(tmp_80 - 1.0/4.0) + tmp_7*(tmp_82 - 1.0/4.0);
      real_t tmp_84 = -tmp_80 - tmp_81 - tmp_82 + 1;
      real_t tmp_85 = 0.054975871827660928*tmp_48;
      real_t tmp_86 = 0.44594849091596489*tmp_26 + 0.44594849091596489*tmp_27 + tmp_28;
      real_t tmp_87 = 0.44594849091596489*tmp_31 + 0.44594849091596489*tmp_32 + tmp_33;
      real_t tmp_88 = 0.44594849091596489*tmp_36 + 0.44594849091596489*tmp_37 + tmp_38;
      real_t tmp_89 = tmp_17*tmp_88 + tmp_20*tmp_87 + tmp_23*tmp_86;
      real_t tmp_90 = tmp_16*tmp_88 + tmp_19*tmp_87 + tmp_22*tmp_86;
      real_t tmp_91 = tmp_15*tmp_88 + tmp_18*tmp_87 + tmp_21*tmp_86;
      real_t tmp_92 = tmp_1*(tmp_90 - 1.0/4.0) + tmp_4*(tmp_89 - 1.0/4.0) + tmp_7*(tmp_91 - 1.0/4.0);
      real_t tmp_93 = -tmp_89 - tmp_90 - tmp_91 + 1;
      real_t tmp_94 = 0.11169079483900572*tmp_48;
      real_t tmp_95 = p_affine_13_0*tmp_17 + p_affine_13_1*tmp_20 + p_affine_13_2*tmp_23;
      real_t tmp_96 = p_affine_13_0*tmp_16 + p_affine_13_1*tmp_19 + p_affine_13_2*tmp_22;
      real_t tmp_97 = p_affine_13_0*tmp_15 + p_affine_13_1*tmp_18 + p_affine_13_2*tmp_21;
      real_t a_0_0 = tmp_49*(-tmp_24*tmp_43 + 3.0*tmp_43*tmp_45*tmp_47 - tmp_44*tmp_45) + tmp_58*(-tmp_24*tmp_56 - tmp_44*tmp_57 + 3.0*tmp_47*tmp_56*tmp_57) + tmp_67*(-tmp_24*tmp_65 - tmp_44*tmp_66 + 3.0*tmp_47*tmp_65*tmp_66) + tmp_76*(-tmp_24*tmp_74 - tmp_44*tmp_75 + 3.0*tmp_47*tmp_74*tmp_75) + tmp_85*(-tmp_24*tmp_83 - tmp_44*tmp_84 + 3.0*tmp_47*tmp_83*tmp_84) + tmp_94*(-tmp_24*tmp_92 - tmp_44*tmp_93 + 3.0*tmp_47*tmp_92*tmp_93);
      real_t a_0_1 = tmp_49*(3.0*tmp_40*tmp_43*tmp_47 - tmp_40*tmp_44 - tmp_43*tmp_95) + tmp_58*(-tmp_44*tmp_53 + 3.0*tmp_47*tmp_53*tmp_56 - tmp_56*tmp_95) + tmp_67*(-tmp_44*tmp_62 + 3.0*tmp_47*tmp_62*tmp_65 - tmp_65*tmp_95) + tmp_76*(-tmp_44*tmp_71 + 3.0*tmp_47*tmp_71*tmp_74 - tmp_74*tmp_95) + tmp_85*(-tmp_44*tmp_80 + 3.0*tmp_47*tmp_80*tmp_83 - tmp_83*tmp_95) + tmp_94*(-tmp_44*tmp_89 + 3.0*tmp_47*tmp_89*tmp_92 - tmp_92*tmp_95);
      real_t a_0_2 = tmp_49*(3.0*tmp_41*tmp_43*tmp_47 - tmp_41*tmp_44 - tmp_43*tmp_96) + tmp_58*(-tmp_44*tmp_54 + 3.0*tmp_47*tmp_54*tmp_56 - tmp_56*tmp_96) + tmp_67*(-tmp_44*tmp_63 + 3.0*tmp_47*tmp_63*tmp_65 - tmp_65*tmp_96) + tmp_76*(-tmp_44*tmp_72 + 3.0*tmp_47*tmp_72*tmp_74 - tmp_74*tmp_96) + tmp_85*(-tmp_44*tmp_81 + 3.0*tmp_47*tmp_81*tmp_83 - tmp_83*tmp_96) + tmp_94*(-tmp_44*tmp_90 + 3.0*tmp_47*tmp_90*tmp_92 - tmp_92*tmp_96);
      real_t a_0_3 = tmp_49*(3.0*tmp_42*tmp_43*tmp_47 - tmp_42*tmp_44 - tmp_43*tmp_97) + tmp_58*(-tmp_44*tmp_55 + 3.0*tmp_47*tmp_55*tmp_56 - tmp_56*tmp_97) + tmp_67*(-tmp_44*tmp_64 + 3.0*tmp_47*tmp_64*tmp_65 - tmp_65*tmp_97) + tmp_76*(-tmp_44*tmp_73 + 3.0*tmp_47*tmp_73*tmp_74 - tmp_74*tmp_97) + tmp_85*(-tmp_44*tmp_82 + 3.0*tmp_47*tmp_82*tmp_83 - tmp_83*tmp_97) + tmp_94*(-tmp_44*tmp_91 + 3.0*tmp_47*tmp_91*tmp_92 - tmp_92*tmp_97);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }

public:



};




class EGVectorLaplaceFormNitscheBC_EE : public hyteg::dg::DGForm
{

 public:
    EGVectorLaplaceFormNitscheBC_EE()
: callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g2 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_g0( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g0( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g1( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g1( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g2( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g2( Point3D( {in_0, in_1, in_2} ) );
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

      real_t tmp_0 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = tmp_2*tmp_3;
      real_t tmp_5 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_6 = -tmp_5;
      real_t tmp_7 = 1.0 / (-tmp_1*tmp_6 + tmp_4);
      real_t tmp_8 = tmp_2*tmp_7;
      real_t tmp_9 = tmp_4*tmp_7;
      real_t tmp_10 = tmp_5*tmp_7;
      real_t tmp_11 = tmp_6*tmp_7;
      real_t tmp_12 = (((tmp_0*tmp_11 + tmp_9)*(tmp_0*tmp_11 + tmp_9)) + ((tmp_0*tmp_8 + tmp_1*tmp_8)*(tmp_0*tmp_8 + tmp_1*tmp_8)) + ((tmp_1*tmp_10 + tmp_9)*(tmp_1*tmp_10 + tmp_9)) + ((tmp_10*tmp_3 + tmp_11*tmp_3)*(tmp_10*tmp_3 + tmp_11*tmp_3)))*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.5*tmp_12;
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
                                       MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_3*tmp_5;
      real_t tmp_7 = -tmp_4;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = -tmp_8;
      real_t tmp_10 = 1.0 / (tmp_6 - tmp_7*tmp_9);
      real_t tmp_11 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_12 = tmp_10*(0.21132486540518713*tmp_1 + tmp_11);
      real_t tmp_13 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_14 = tmp_10*(0.21132486540518713*tmp_0 + tmp_13);
      real_t tmp_15 = tmp_12*tmp_4 + tmp_14*tmp_5 - 1.0/3.0;
      real_t tmp_16 = tmp_12*tmp_3 + tmp_14*tmp_8 - 1.0/3.0;
      real_t tmp_17 = tmp_15*tmp_3 + tmp_16*tmp_7;
      real_t tmp_18 = tmp_10*tmp_6;
      real_t tmp_19 = tmp_10*tmp_7;
      real_t tmp_20 = tmp_10*tmp_4;
      real_t tmp_21 = 1.0*p_affine_10_0*(tmp_18 + tmp_19*tmp_8) + 1.0*p_affine_10_1*(tmp_19*tmp_3 + tmp_20*tmp_3);
      real_t tmp_22 = tmp_15*tmp_9 + tmp_16*tmp_5;
      real_t tmp_23 = tmp_10*tmp_5;
      real_t tmp_24 = 1.0*p_affine_10_0*(tmp_23*tmp_8 + tmp_23*tmp_9) + 1.0*p_affine_10_1*(tmp_18 + tmp_20*tmp_9);
      real_t tmp_25 = 1.0 / (tmp_2);
      real_t tmp_26 = 0.78867513459481287*tmp_1 + tmp_11;
      real_t tmp_27 = 0.78867513459481287*tmp_0 + tmp_13;
      real_t tmp_28 = tmp_20*tmp_26 + tmp_23*tmp_27 - 1.0/3.0;
      real_t tmp_29 = tmp_10*tmp_26*tmp_3 + tmp_10*tmp_27*tmp_8 - 1.0/3.0;
      real_t tmp_30 = tmp_28*tmp_3 + tmp_29*tmp_7;
      real_t tmp_31 = tmp_28*tmp_9 + tmp_29*tmp_5;
      real_t a_0_0 = 0.5*tmp_2*(-tmp_17*tmp_21 - tmp_22*tmp_24 + 3*tmp_25*((tmp_17*tmp_17) + (tmp_22*tmp_22))) + 0.5*tmp_2*(-tmp_21*tmp_30 - tmp_24*tmp_31 + 3*tmp_25*((tmp_30*tmp_30) + (tmp_31*tmp_31)));
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
                                          MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_7 = -tmp_6;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = -tmp_8;
      real_t tmp_10 = 1.0 / (tmp_5 - tmp_7*tmp_9);
      real_t tmp_11 = tmp_10*tmp_5;
      real_t tmp_12 = tmp_10*tmp_7;
      real_t tmp_13 = tmp_10*tmp_6;
      real_t tmp_14 = p_affine_10_0*(tmp_11 + tmp_12*tmp_8) + p_affine_10_1*(tmp_12*tmp_3 + tmp_13*tmp_3);
      real_t tmp_15 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_16 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_17 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_18 = tmp_15*tmp_17;
      real_t tmp_19 = -tmp_16;
      real_t tmp_20 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_21 = -tmp_20;
      real_t tmp_22 = 1.0 / (tmp_18 - tmp_19*tmp_21);
      real_t tmp_23 = -p_affine_3_1;
      real_t tmp_24 = p_affine_6_1 + 0.21132486540518713*tmp_1;
      real_t tmp_25 = tmp_22*(tmp_23 + tmp_24);
      real_t tmp_26 = -p_affine_3_0;
      real_t tmp_27 = p_affine_6_0 + 0.21132486540518713*tmp_0;
      real_t tmp_28 = tmp_22*(tmp_26 + tmp_27);
      real_t tmp_29 = tmp_16*tmp_25 + tmp_17*tmp_28 - 1.0/3.0;
      real_t tmp_30 = tmp_15*tmp_25 + tmp_20*tmp_28 - 1.0/3.0;
      real_t tmp_31 = tmp_15*tmp_29 + tmp_19*tmp_30;
      real_t tmp_32 = tmp_10*tmp_4;
      real_t tmp_33 = p_affine_10_0*(tmp_32*tmp_8 + tmp_32*tmp_9) + p_affine_10_1*(tmp_11 + tmp_13*tmp_9);
      real_t tmp_34 = tmp_17*tmp_30 + tmp_21*tmp_29;
      real_t tmp_35 = -p_affine_0_1;
      real_t tmp_36 = tmp_10*(tmp_24 + tmp_35);
      real_t tmp_37 = -p_affine_0_0;
      real_t tmp_38 = tmp_10*(tmp_27 + tmp_37);
      real_t tmp_39 = tmp_36*tmp_6 + tmp_38*tmp_4 - 1.0/3.0;
      real_t tmp_40 = tmp_3*tmp_36 + tmp_38*tmp_8 - 1.0/3.0;
      real_t tmp_41 = tmp_3*tmp_39 + tmp_40*tmp_7;
      real_t tmp_42 = tmp_18*tmp_22;
      real_t tmp_43 = tmp_19*tmp_22;
      real_t tmp_44 = tmp_16*tmp_22;
      real_t tmp_45 = 0.5*p_affine_10_0*(tmp_20*tmp_43 + tmp_42) + 0.5*p_affine_10_1*(tmp_15*tmp_43 + tmp_15*tmp_44);
      real_t tmp_46 = tmp_39*tmp_9 + tmp_4*tmp_40;
      real_t tmp_47 = tmp_17*tmp_22;
      real_t tmp_48 = 0.5*p_affine_10_0*(tmp_20*tmp_47 + tmp_21*tmp_47) + 0.5*p_affine_10_1*(tmp_21*tmp_44 + tmp_42);
      real_t tmp_49 = 3/tmp_2;
      real_t tmp_50 = p_affine_6_1 + 0.78867513459481287*tmp_1;
      real_t tmp_51 = tmp_23 + tmp_50;
      real_t tmp_52 = p_affine_6_0 + 0.78867513459481287*tmp_0;
      real_t tmp_53 = tmp_26 + tmp_52;
      real_t tmp_54 = tmp_44*tmp_51 + tmp_47*tmp_53 - 1.0/3.0;
      real_t tmp_55 = tmp_15*tmp_22*tmp_51 + tmp_20*tmp_22*tmp_53 - 1.0/3.0;
      real_t tmp_56 = tmp_15*tmp_54 + tmp_19*tmp_55;
      real_t tmp_57 = tmp_17*tmp_55 + tmp_21*tmp_54;
      real_t tmp_58 = tmp_35 + tmp_50;
      real_t tmp_59 = tmp_37 + tmp_52;
      real_t tmp_60 = tmp_13*tmp_58 + tmp_32*tmp_59 - 1.0/3.0;
      real_t tmp_61 = tmp_10*tmp_3*tmp_58 + tmp_10*tmp_59*tmp_8 - 1.0/3.0;
      real_t tmp_62 = tmp_3*tmp_60 + tmp_61*tmp_7;
      real_t tmp_63 = tmp_4*tmp_61 + tmp_60*tmp_9;
      real_t a_0_0 = 0.5*tmp_2*(0.5*tmp_14*tmp_31 + 0.5*tmp_33*tmp_34 - tmp_41*tmp_45 - tmp_46*tmp_48 - tmp_49*(tmp_31*tmp_41 + tmp_34*tmp_46)) + 0.5*tmp_2*(0.5*tmp_14*tmp_56 + 0.5*tmp_33*tmp_57 - tmp_45*tmp_62 - tmp_48*tmp_63 - tmp_49*(tmp_56*tmp_62 + tmp_57*tmp_63));
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
                                                   MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_5 = 0.21132486540518713*tmp_1 + tmp_4;
      real_t tmp_6 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_9 = -tmp_8;
      real_t tmp_10 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_11 = -tmp_10;
      real_t tmp_12 = 1.0 / (-tmp_11*tmp_9 + tmp_7);
      real_t tmp_13 = tmp_12*tmp_8;
      real_t tmp_14 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_15 = tmp_12*(0.21132486540518713*tmp_0 + tmp_14);
      real_t tmp_16 = tmp_13*tmp_5 + tmp_15*tmp_6 - 1.0/3.0;
      real_t tmp_17 = tmp_12*tmp_3;
      real_t tmp_18 = tmp_10*tmp_15 + tmp_17*tmp_5 - 1.0/3.0;
      real_t tmp_19 = tmp_16*tmp_3 + tmp_18*tmp_9;
      real_t tmp_20 = tmp_12*tmp_7;
      real_t tmp_21 = tmp_12*tmp_9;
      real_t tmp_22 = 2*p_affine_10_0*(tmp_10*tmp_21 + tmp_20) + 2*p_affine_10_1*(tmp_13*tmp_3 + tmp_21*tmp_3);
      real_t tmp_23 = tmp_11*tmp_16 + tmp_18*tmp_6;
      real_t tmp_24 = tmp_12*tmp_6;
      real_t tmp_25 = 2*p_affine_10_0*(tmp_10*tmp_24 + tmp_11*tmp_24) + 2*p_affine_10_1*(tmp_11*tmp_13 + tmp_20);
      real_t tmp_26 = 1.0 / (tmp_2);
      real_t tmp_27 = 0.78867513459481287*tmp_1 + tmp_4;
      real_t tmp_28 = 0.78867513459481287*tmp_0 + tmp_14;
      real_t tmp_29 = tmp_13*tmp_27 + tmp_24*tmp_28 - 1.0/3.0;
      real_t tmp_30 = tmp_10*tmp_12*tmp_28 + tmp_17*tmp_27 - 1.0/3.0;
      real_t tmp_31 = tmp_29*tmp_3 + tmp_30*tmp_9;
      real_t tmp_32 = tmp_11*tmp_29 + tmp_30*tmp_6;
      real_t a_0_0 = 0.5*tmp_2*(-tmp_19*tmp_22 - tmp_23*tmp_25 + 3*tmp_26*((tmp_19*tmp_19) + (tmp_23*tmp_23))) + 0.5*tmp_2*(-tmp_22*tmp_31 - tmp_25*tmp_32 + 3*tmp_26*((tmp_31*tmp_31) + (tmp_32*tmp_32)));
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
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id3 = 0;
      Scalar_Variable_Coefficient_2D_g0( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id0 );
      Scalar_Variable_Coefficient_2D_g1( 0.78867513459481287*p_affine_6_0 + 0.21132486540518713*p_affine_7_0, 0.78867513459481287*p_affine_6_1 + 0.21132486540518713*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id1 );
      Scalar_Variable_Coefficient_2D_g0( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id2 );
      Scalar_Variable_Coefficient_2D_g1( 0.21132486540518713*p_affine_6_0 + 0.78867513459481287*p_affine_7_0, 0.21132486540518713*p_affine_6_1 + 0.78867513459481287*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id3 );
      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = 1.0 / (tmp_2);
      real_t tmp_4 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_6 = 0.21132486540518713*tmp_1 + tmp_5;
      real_t tmp_7 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = -tmp_9;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = -tmp_11;
      real_t tmp_13 = 1.0 / (-tmp_10*tmp_12 + tmp_8);
      real_t tmp_14 = tmp_13*tmp_9;
      real_t tmp_15 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_16 = tmp_13*(0.21132486540518713*tmp_0 + tmp_15);
      real_t tmp_17 = tmp_14*tmp_6 + tmp_16*tmp_7 - 1.0/3.0;
      real_t tmp_18 = tmp_13*tmp_4;
      real_t tmp_19 = tmp_11*tmp_16 + tmp_18*tmp_6 - 1.0/3.0;
      real_t tmp_20 = tmp_13*tmp_8;
      real_t tmp_21 = tmp_10*tmp_13;
      real_t tmp_22 = p_affine_10_0*(tmp_11*tmp_21 + tmp_20) + p_affine_10_1*(tmp_14*tmp_4 + tmp_21*tmp_4);
      real_t tmp_23 = tmp_13*tmp_7;
      real_t tmp_24 = p_affine_10_0*(tmp_11*tmp_23 + tmp_12*tmp_23) + p_affine_10_1*(tmp_12*tmp_14 + tmp_20);
      real_t tmp_25 = 0.78867513459481287*tmp_1 + tmp_5;
      real_t tmp_26 = 0.78867513459481287*tmp_0 + tmp_15;
      real_t tmp_27 = tmp_14*tmp_25 + tmp_23*tmp_26 - 1.0/3.0;
      real_t tmp_28 = tmp_11*tmp_13*tmp_26 + tmp_18*tmp_25 - 1.0/3.0;
      real_t a_0_0 = 0.5*tmp_2*(Scalar_Variable_Coefficient_2D_g0_out0_id0*(-tmp_22 + 3*tmp_3*(tmp_10*tmp_19 + tmp_17*tmp_4)) + Scalar_Variable_Coefficient_2D_g1_out0_id1*(-tmp_24 + 3*tmp_3*(tmp_12*tmp_17 + tmp_19*tmp_7))) + 0.5*tmp_2*(Scalar_Variable_Coefficient_2D_g0_out0_id2*(-tmp_22 + 3*tmp_3*(tmp_10*tmp_28 + tmp_27*tmp_4)) + Scalar_Variable_Coefficient_2D_g1_out0_id3*(-tmp_24 + 3*tmp_3*(tmp_12*tmp_27 + tmp_28*tmp_7)));
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
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id13 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id14 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id15 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id16 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id17 = 0;
      Scalar_Variable_Coefficient_3D_g0( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id0 );
      Scalar_Variable_Coefficient_3D_g1( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id1 );
      Scalar_Variable_Coefficient_3D_g2( 0.81684757298045851*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.81684757298045851*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.81684757298045851*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id2 );
      Scalar_Variable_Coefficient_3D_g0( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id3 );
      Scalar_Variable_Coefficient_3D_g1( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id4 );
      Scalar_Variable_Coefficient_3D_g2( 0.10810301816807022*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.10810301816807022*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.10810301816807022*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id5 );
      Scalar_Variable_Coefficient_3D_g0( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id6 );
      Scalar_Variable_Coefficient_3D_g1( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id7 );
      Scalar_Variable_Coefficient_3D_g2( 0.091576213509770743*p_affine_10_0 + 0.091576213509770743*p_affine_8_0 + 0.81684757298045851*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.091576213509770743*p_affine_8_1 + 0.81684757298045851*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.091576213509770743*p_affine_8_2 + 0.81684757298045851*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id8 );
      Scalar_Variable_Coefficient_3D_g0( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id9 );
      Scalar_Variable_Coefficient_3D_g1( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id10 );
      Scalar_Variable_Coefficient_3D_g2( 0.44594849091596489*p_affine_10_0 + 0.44594849091596489*p_affine_8_0 + 0.10810301816807022*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.44594849091596489*p_affine_8_1 + 0.10810301816807022*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.44594849091596489*p_affine_8_2 + 0.10810301816807022*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id11 );
      Scalar_Variable_Coefficient_3D_g0( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id12 );
      Scalar_Variable_Coefficient_3D_g1( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id13 );
      Scalar_Variable_Coefficient_3D_g2( 0.091576213509770743*p_affine_10_0 + 0.81684757298045851*p_affine_8_0 + 0.091576213509770743*p_affine_9_0, 0.091576213509770743*p_affine_10_1 + 0.81684757298045851*p_affine_8_1 + 0.091576213509770743*p_affine_9_1, 0.091576213509770743*p_affine_10_2 + 0.81684757298045851*p_affine_8_2 + 0.091576213509770743*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id14 );
      Scalar_Variable_Coefficient_3D_g0( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id15 );
      Scalar_Variable_Coefficient_3D_g1( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id16 );
      Scalar_Variable_Coefficient_3D_g2( 0.44594849091596489*p_affine_10_0 + 0.10810301816807022*p_affine_8_0 + 0.44594849091596489*p_affine_9_0, 0.44594849091596489*p_affine_10_1 + 0.10810301816807022*p_affine_8_1 + 0.44594849091596489*p_affine_9_1, 0.44594849091596489*p_affine_10_2 + 0.10810301816807022*p_affine_8_2 + 0.44594849091596489*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id17 );
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
      real_t tmp_26 = tmp_23*(0.091576213509770743*tmp_24 + tmp_25 + 0.81684757298045851*tmp_5);
      real_t tmp_27 = tmp_11*tmp_18 - tmp_14*tmp_9;
      real_t tmp_28 = -tmp_1;
      real_t tmp_29 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_30 = tmp_23*(0.81684757298045851*tmp_2 + 0.091576213509770743*tmp_28 + tmp_29);
      real_t tmp_31 = tmp_15 - tmp_20;
      real_t tmp_32 = -tmp_3;
      real_t tmp_33 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_34 = tmp_23*(0.81684757298045851*tmp_0 + 0.091576213509770743*tmp_32 + tmp_33);
      real_t tmp_35 = tmp_13*tmp_26 + tmp_27*tmp_30 + tmp_31*tmp_34 - 1.0/4.0;
      real_t tmp_36 = -tmp_10*tmp_8 + tmp_11*tmp_17;
      real_t tmp_37 = -tmp_11*tmp_16 + tmp_14*tmp_8;
      real_t tmp_38 = tmp_10*tmp_16 - tmp_21;
      real_t tmp_39 = tmp_26*tmp_36 + tmp_30*tmp_37 + tmp_34*tmp_38 - 1.0/4.0;
      real_t tmp_40 = tmp_12*tmp_8 - tmp_17*tmp_9;
      real_t tmp_41 = tmp_16*tmp_9 - tmp_18*tmp_8;
      real_t tmp_42 = tmp_19 - tmp_22;
      real_t tmp_43 = tmp_26*tmp_40 + tmp_30*tmp_41 + tmp_34*tmp_42 - 1.0/4.0;
      real_t tmp_44 = tmp_23*tmp_8;
      real_t tmp_45 = tmp_23*tmp_9;
      real_t tmp_46 = tmp_11*tmp_23;
      real_t tmp_47 = p_affine_13_0*(tmp_31*tmp_44 + tmp_38*tmp_45 + tmp_42*tmp_46) + p_affine_13_1*(tmp_27*tmp_44 + tmp_37*tmp_45 + tmp_41*tmp_46) + p_affine_13_2*(tmp_13*tmp_44 + tmp_36*tmp_45 + tmp_40*tmp_46);
      real_t tmp_48 = tmp_17*tmp_23;
      real_t tmp_49 = tmp_12*tmp_23;
      real_t tmp_50 = tmp_10*tmp_23;
      real_t tmp_51 = p_affine_13_0*(tmp_31*tmp_48 + tmp_38*tmp_49 + tmp_42*tmp_50) + p_affine_13_1*(tmp_27*tmp_48 + tmp_37*tmp_49 + tmp_41*tmp_50) + p_affine_13_2*(tmp_13*tmp_48 + tmp_36*tmp_49 + tmp_40*tmp_50);
      real_t tmp_52 = tmp_16*tmp_23;
      real_t tmp_53 = tmp_18*tmp_23;
      real_t tmp_54 = tmp_14*tmp_23;
      real_t tmp_55 = p_affine_13_0*(tmp_31*tmp_52 + tmp_38*tmp_53 + tmp_42*tmp_54) + p_affine_13_1*(tmp_27*tmp_52 + tmp_37*tmp_53 + tmp_41*tmp_54) + p_affine_13_2*(tmp_13*tmp_52 + tmp_36*tmp_53 + tmp_40*tmp_54);
      real_t tmp_56 = 1.0*std::pow(tmp_6, 1.0/2.0);
      real_t tmp_57 = tmp_23*(0.44594849091596489*tmp_24 + tmp_25 + 0.10810301816807022*tmp_5);
      real_t tmp_58 = tmp_23*(0.10810301816807022*tmp_2 + 0.44594849091596489*tmp_28 + tmp_29);
      real_t tmp_59 = tmp_23*(0.10810301816807022*tmp_0 + 0.44594849091596489*tmp_32 + tmp_33);
      real_t tmp_60 = tmp_13*tmp_57 + tmp_27*tmp_58 + tmp_31*tmp_59 - 1.0/4.0;
      real_t tmp_61 = tmp_36*tmp_57 + tmp_37*tmp_58 + tmp_38*tmp_59 - 1.0/4.0;
      real_t tmp_62 = tmp_40*tmp_57 + tmp_41*tmp_58 + tmp_42*tmp_59 - 1.0/4.0;
      real_t tmp_63 = tmp_23*(0.81684757298045851*tmp_24 + tmp_25 + 0.091576213509770743*tmp_5);
      real_t tmp_64 = tmp_23*(0.091576213509770743*tmp_2 + 0.81684757298045851*tmp_28 + tmp_29);
      real_t tmp_65 = tmp_23*(0.091576213509770743*tmp_0 + 0.81684757298045851*tmp_32 + tmp_33);
      real_t tmp_66 = tmp_13*tmp_63 + tmp_27*tmp_64 + tmp_31*tmp_65 - 1.0/4.0;
      real_t tmp_67 = tmp_36*tmp_63 + tmp_37*tmp_64 + tmp_38*tmp_65 - 1.0/4.0;
      real_t tmp_68 = tmp_40*tmp_63 + tmp_41*tmp_64 + tmp_42*tmp_65 - 1.0/4.0;
      real_t tmp_69 = tmp_23*(0.10810301816807022*tmp_24 + tmp_25 + 0.44594849091596489*tmp_5);
      real_t tmp_70 = tmp_23*(0.44594849091596489*tmp_2 + 0.10810301816807022*tmp_28 + tmp_29);
      real_t tmp_71 = tmp_23*(0.44594849091596489*tmp_0 + 0.10810301816807022*tmp_32 + tmp_33);
      real_t tmp_72 = tmp_13*tmp_69 + tmp_27*tmp_70 + tmp_31*tmp_71 - 1.0/4.0;
      real_t tmp_73 = tmp_36*tmp_69 + tmp_37*tmp_70 + tmp_38*tmp_71 - 1.0/4.0;
      real_t tmp_74 = tmp_40*tmp_69 + tmp_41*tmp_70 + tmp_42*tmp_71 - 1.0/4.0;
      real_t tmp_75 = tmp_23*(0.091576213509770743*tmp_24 + tmp_25 + 0.091576213509770743*tmp_5);
      real_t tmp_76 = tmp_23*(0.091576213509770743*tmp_2 + 0.091576213509770743*tmp_28 + tmp_29);
      real_t tmp_77 = tmp_23*(0.091576213509770743*tmp_0 + 0.091576213509770743*tmp_32 + tmp_33);
      real_t tmp_78 = tmp_13*tmp_75 + tmp_27*tmp_76 + tmp_31*tmp_77 - 1.0/4.0;
      real_t tmp_79 = tmp_36*tmp_75 + tmp_37*tmp_76 + tmp_38*tmp_77 - 1.0/4.0;
      real_t tmp_80 = tmp_40*tmp_75 + tmp_41*tmp_76 + tmp_42*tmp_77 - 1.0/4.0;
      real_t tmp_81 = tmp_23*(0.44594849091596489*tmp_24 + tmp_25 + 0.44594849091596489*tmp_5);
      real_t tmp_82 = tmp_23*(0.44594849091596489*tmp_2 + 0.44594849091596489*tmp_28 + tmp_29);
      real_t tmp_83 = tmp_23*(0.44594849091596489*tmp_0 + 0.44594849091596489*tmp_32 + tmp_33);
      real_t tmp_84 = tmp_13*tmp_81 + tmp_27*tmp_82 + tmp_31*tmp_83 - 1.0/4.0;
      real_t tmp_85 = tmp_36*tmp_81 + tmp_37*tmp_82 + tmp_38*tmp_83 - 1.0/4.0;
      real_t tmp_86 = tmp_40*tmp_81 + tmp_41*tmp_82 + tmp_42*tmp_83 - 1.0/4.0;
      real_t a_0_0 = 0.054975871827660928*tmp_56*(Scalar_Variable_Coefficient_3D_g0_out0_id0*(-tmp_47 + 3.0*tmp_7*(tmp_11*tmp_43 + tmp_35*tmp_8 + tmp_39*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id1*(-tmp_51 + 3.0*tmp_7*(tmp_10*tmp_43 + tmp_12*tmp_39 + tmp_17*tmp_35)) + Scalar_Variable_Coefficient_3D_g2_out0_id2*(-tmp_55 + 3.0*tmp_7*(tmp_14*tmp_43 + tmp_16*tmp_35 + tmp_18*tmp_39))) + 0.054975871827660928*tmp_56*(Scalar_Variable_Coefficient_3D_g0_out0_id12*(-tmp_47 + 3.0*tmp_7*(tmp_11*tmp_80 + tmp_78*tmp_8 + tmp_79*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id13*(-tmp_51 + 3.0*tmp_7*(tmp_10*tmp_80 + tmp_12*tmp_79 + tmp_17*tmp_78)) + Scalar_Variable_Coefficient_3D_g2_out0_id14*(-tmp_55 + 3.0*tmp_7*(tmp_14*tmp_80 + tmp_16*tmp_78 + tmp_18*tmp_79))) + 0.11169079483900572*tmp_56*(Scalar_Variable_Coefficient_3D_g0_out0_id15*(-tmp_47 + 3.0*tmp_7*(tmp_11*tmp_86 + tmp_8*tmp_84 + tmp_85*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id16*(-tmp_51 + 3.0*tmp_7*(tmp_10*tmp_86 + tmp_12*tmp_85 + tmp_17*tmp_84)) + Scalar_Variable_Coefficient_3D_g2_out0_id17*(-tmp_55 + 3.0*tmp_7*(tmp_14*tmp_86 + tmp_16*tmp_84 + tmp_18*tmp_85))) + 0.11169079483900572*tmp_56*(Scalar_Variable_Coefficient_3D_g0_out0_id3*(-tmp_47 + 3.0*tmp_7*(tmp_11*tmp_62 + tmp_60*tmp_8 + tmp_61*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id4*(-tmp_51 + 3.0*tmp_7*(tmp_10*tmp_62 + tmp_12*tmp_61 + tmp_17*tmp_60)) + Scalar_Variable_Coefficient_3D_g2_out0_id5*(-tmp_55 + 3.0*tmp_7*(tmp_14*tmp_62 + tmp_16*tmp_60 + tmp_18*tmp_61))) + 0.054975871827660928*tmp_56*(Scalar_Variable_Coefficient_3D_g0_out0_id6*(-tmp_47 + 3.0*tmp_7*(tmp_11*tmp_68 + tmp_66*tmp_8 + tmp_67*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id7*(-tmp_51 + 3.0*tmp_7*(tmp_10*tmp_68 + tmp_12*tmp_67 + tmp_17*tmp_66)) + Scalar_Variable_Coefficient_3D_g2_out0_id8*(-tmp_55 + 3.0*tmp_7*(tmp_14*tmp_68 + tmp_16*tmp_66 + tmp_18*tmp_67))) + 0.11169079483900572*tmp_56*(Scalar_Variable_Coefficient_3D_g0_out0_id9*(-tmp_47 + 3.0*tmp_7*(tmp_11*tmp_74 + tmp_72*tmp_8 + tmp_73*tmp_9)) + Scalar_Variable_Coefficient_3D_g1_out0_id10*(-tmp_51 + 3.0*tmp_7*(tmp_10*tmp_74 + tmp_12*tmp_73 + tmp_17*tmp_72)) + Scalar_Variable_Coefficient_3D_g2_out0_id11*(-tmp_55 + 3.0*tmp_7*(tmp_14*tmp_74 + tmp_16*tmp_72 + tmp_18*tmp_73)));
      elMat( 0, 0) = a_0_0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
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

      real_t tmp_0 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = tmp_2 - tmp_5;
      real_t tmp_7 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_8 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_9 = tmp_4*tmp_7;
      real_t tmp_10 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_11 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_12 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_13 = tmp_1*tmp_7;
      real_t tmp_14 = tmp_0*tmp_12;
      real_t tmp_15 = 1.0 / (tmp_10*tmp_2 - tmp_10*tmp_5 + tmp_11*tmp_12*tmp_3 - tmp_11*tmp_13 - tmp_14*tmp_8 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_15*tmp_7;
      real_t tmp_17 = tmp_12*tmp_3 - tmp_13;
      real_t tmp_18 = tmp_0*tmp_15;
      real_t tmp_19 = -tmp_14 + tmp_9;
      real_t tmp_20 = tmp_15*tmp_3;
      real_t tmp_21 = -tmp_0*tmp_8 + tmp_11*tmp_3;
      real_t tmp_22 = -tmp_10*tmp_3 + tmp_7*tmp_8;
      real_t tmp_23 = tmp_0*tmp_10 - tmp_11*tmp_7;
      real_t tmp_24 = -tmp_1*tmp_11 + tmp_4*tmp_8;
      real_t tmp_25 = tmp_1*tmp_10 - tmp_12*tmp_8;
      real_t tmp_26 = -tmp_10*tmp_4 + tmp_11*tmp_12;
      real_t tmp_27 = tmp_12*tmp_15;
      real_t tmp_28 = tmp_15*tmp_4;
      real_t tmp_29 = tmp_1*tmp_15;
      real_t tmp_30 = tmp_10*tmp_15;
      real_t tmp_31 = tmp_11*tmp_15;
      real_t tmp_32 = tmp_15*tmp_8;
      real_t tmp_33 = p_affine_0_0*p_affine_1_1;
      real_t tmp_34 = p_affine_0_0*p_affine_1_2;
      real_t tmp_35 = p_affine_2_1*p_affine_3_2;
      real_t tmp_36 = p_affine_0_1*p_affine_1_0;
      real_t tmp_37 = p_affine_0_1*p_affine_1_2;
      real_t tmp_38 = p_affine_2_2*p_affine_3_0;
      real_t tmp_39 = p_affine_0_2*p_affine_1_0;
      real_t tmp_40 = p_affine_0_2*p_affine_1_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_1;
      real_t tmp_42 = p_affine_2_2*p_affine_3_1;
      real_t tmp_43 = p_affine_2_0*p_affine_3_2;
      real_t tmp_44 = p_affine_2_1*p_affine_3_0;
      real_t tmp_45 = (((tmp_16*tmp_21 + tmp_18*tmp_22 + tmp_20*tmp_23)*(tmp_16*tmp_21 + tmp_18*tmp_22 + tmp_20*tmp_23)) + ((tmp_16*tmp_24 + tmp_18*tmp_25 + tmp_20*tmp_26)*(tmp_16*tmp_24 + tmp_18*tmp_25 + tmp_20*tmp_26)) + ((tmp_16*tmp_6 + tmp_17*tmp_18 + tmp_19*tmp_20)*(tmp_16*tmp_6 + tmp_17*tmp_18 + tmp_19*tmp_20)) + ((tmp_17*tmp_28 + tmp_19*tmp_29 + tmp_27*tmp_6)*(tmp_17*tmp_28 + tmp_19*tmp_29 + tmp_27*tmp_6)) + ((tmp_17*tmp_31 + tmp_19*tmp_32 + tmp_30*tmp_6)*(tmp_17*tmp_31 + tmp_19*tmp_32 + tmp_30*tmp_6)) + ((tmp_21*tmp_27 + tmp_22*tmp_28 + tmp_23*tmp_29)*(tmp_21*tmp_27 + tmp_22*tmp_28 + tmp_23*tmp_29)) + ((tmp_21*tmp_30 + tmp_22*tmp_31 + tmp_23*tmp_32)*(tmp_21*tmp_30 + tmp_22*tmp_31 + tmp_23*tmp_32)) + ((tmp_24*tmp_27 + tmp_25*tmp_28 + tmp_26*tmp_29)*(tmp_24*tmp_27 + tmp_25*tmp_28 + tmp_26*tmp_29)) + ((tmp_24*tmp_30 + tmp_25*tmp_31 + tmp_26*tmp_32)*(tmp_24*tmp_30 + tmp_25*tmp_31 + tmp_26*tmp_32)))*std::abs(p_affine_0_0*tmp_35 - p_affine_0_0*tmp_42 + p_affine_0_1*tmp_38 - p_affine_0_1*tmp_43 + p_affine_0_2*tmp_41 - p_affine_0_2*tmp_44 - p_affine_1_0*tmp_35 + p_affine_1_0*tmp_42 - p_affine_1_1*tmp_38 + p_affine_1_1*tmp_43 - p_affine_1_2*tmp_41 + p_affine_1_2*tmp_44 + p_affine_2_0*tmp_37 - p_affine_2_0*tmp_40 - p_affine_2_1*tmp_34 + p_affine_2_1*tmp_39 + p_affine_2_2*tmp_33 - p_affine_2_2*tmp_36 - p_affine_3_0*tmp_37 + p_affine_3_0*tmp_40 + p_affine_3_1*tmp_34 - p_affine_3_1*tmp_39 - p_affine_3_2*tmp_33 + p_affine_3_2*tmp_36);
      real_t a_0_0 = 0.16666666666666666*tmp_45;
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
                               MatrixXr&                            elMat ) const override
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
      real_t tmp_43 = tmp_0*tmp_15;
      real_t tmp_44 = tmp_1*tmp_15;
      real_t tmp_45 = tmp_15*tmp_4;
      real_t tmp_46 = 1.0*p_affine_13_0*(tmp_27*tmp_43 + tmp_36*tmp_44 + tmp_40*tmp_45) + 1.0*p_affine_13_1*(tmp_21*tmp_43 + tmp_35*tmp_44 + tmp_39*tmp_45) + 1.0*p_affine_13_2*(tmp_34*tmp_44 + tmp_38*tmp_45 + tmp_43*tmp_7);
      real_t tmp_47 = tmp_11*tmp_33 + tmp_2*tmp_41 + tmp_37*tmp_5;
      real_t tmp_48 = tmp_11*tmp_15;
      real_t tmp_49 = tmp_15*tmp_5;
      real_t tmp_50 = tmp_15*tmp_2;
      real_t tmp_51 = 1.0*p_affine_13_0*(tmp_27*tmp_48 + tmp_36*tmp_49 + tmp_40*tmp_50) + 1.0*p_affine_13_1*(tmp_21*tmp_48 + tmp_35*tmp_49 + tmp_39*tmp_50) + 1.0*p_affine_13_2*(tmp_34*tmp_49 + tmp_38*tmp_50 + tmp_48*tmp_7);
      real_t tmp_52 = tmp_10*tmp_33 + tmp_12*tmp_37 + tmp_41*tmp_8;
      real_t tmp_53 = tmp_10*tmp_15;
      real_t tmp_54 = tmp_12*tmp_15;
      real_t tmp_55 = tmp_15*tmp_8;
      real_t tmp_56 = 1.0*p_affine_13_0*(tmp_27*tmp_53 + tmp_36*tmp_54 + tmp_40*tmp_55) + 1.0*p_affine_13_1*(tmp_21*tmp_53 + tmp_35*tmp_54 + tmp_39*tmp_55) + 1.0*p_affine_13_2*(tmp_34*tmp_54 + tmp_38*tmp_55 + tmp_53*tmp_7);
      real_t tmp_57 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_58 = std::pow(tmp_57, -0.25);
      real_t tmp_59 = 1.0*std::pow(tmp_57, 1.0/2.0);
      real_t tmp_60 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_61 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_62 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_63 = tmp_21*tmp_61 + tmp_27*tmp_62 + tmp_60*tmp_7 - 1.0/4.0;
      real_t tmp_64 = tmp_34*tmp_60 + tmp_35*tmp_61 + tmp_36*tmp_62 - 1.0/4.0;
      real_t tmp_65 = tmp_38*tmp_60 + tmp_39*tmp_61 + tmp_40*tmp_62 - 1.0/4.0;
      real_t tmp_66 = tmp_0*tmp_63 + tmp_1*tmp_64 + tmp_4*tmp_65;
      real_t tmp_67 = tmp_11*tmp_63 + tmp_2*tmp_65 + tmp_5*tmp_64;
      real_t tmp_68 = tmp_10*tmp_63 + tmp_12*tmp_64 + tmp_65*tmp_8;
      real_t tmp_69 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_70 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_71 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_72 = tmp_21*tmp_70 + tmp_27*tmp_71 + tmp_69*tmp_7 - 1.0/4.0;
      real_t tmp_73 = tmp_34*tmp_69 + tmp_35*tmp_70 + tmp_36*tmp_71 - 1.0/4.0;
      real_t tmp_74 = tmp_38*tmp_69 + tmp_39*tmp_70 + tmp_40*tmp_71 - 1.0/4.0;
      real_t tmp_75 = tmp_0*tmp_72 + tmp_1*tmp_73 + tmp_4*tmp_74;
      real_t tmp_76 = tmp_11*tmp_72 + tmp_2*tmp_74 + tmp_5*tmp_73;
      real_t tmp_77 = tmp_10*tmp_72 + tmp_12*tmp_73 + tmp_74*tmp_8;
      real_t tmp_78 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_79 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_80 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_81 = tmp_21*tmp_79 + tmp_27*tmp_80 + tmp_7*tmp_78 - 1.0/4.0;
      real_t tmp_82 = tmp_34*tmp_78 + tmp_35*tmp_79 + tmp_36*tmp_80 - 1.0/4.0;
      real_t tmp_83 = tmp_38*tmp_78 + tmp_39*tmp_79 + tmp_40*tmp_80 - 1.0/4.0;
      real_t tmp_84 = tmp_0*tmp_81 + tmp_1*tmp_82 + tmp_4*tmp_83;
      real_t tmp_85 = tmp_11*tmp_81 + tmp_2*tmp_83 + tmp_5*tmp_82;
      real_t tmp_86 = tmp_10*tmp_81 + tmp_12*tmp_82 + tmp_8*tmp_83;
      real_t tmp_87 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_88 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_89 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_90 = tmp_21*tmp_88 + tmp_27*tmp_89 + tmp_7*tmp_87 - 1.0/4.0;
      real_t tmp_91 = tmp_34*tmp_87 + tmp_35*tmp_88 + tmp_36*tmp_89 - 1.0/4.0;
      real_t tmp_92 = tmp_38*tmp_87 + tmp_39*tmp_88 + tmp_40*tmp_89 - 1.0/4.0;
      real_t tmp_93 = tmp_0*tmp_90 + tmp_1*tmp_91 + tmp_4*tmp_92;
      real_t tmp_94 = tmp_11*tmp_90 + tmp_2*tmp_92 + tmp_5*tmp_91;
      real_t tmp_95 = tmp_10*tmp_90 + tmp_12*tmp_91 + tmp_8*tmp_92;
      real_t tmp_96 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_97 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_98 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_99 = tmp_21*tmp_97 + tmp_27*tmp_98 + tmp_7*tmp_96 - 1.0/4.0;
      real_t tmp_100 = tmp_34*tmp_96 + tmp_35*tmp_97 + tmp_36*tmp_98 - 1.0/4.0;
      real_t tmp_101 = tmp_38*tmp_96 + tmp_39*tmp_97 + tmp_40*tmp_98 - 1.0/4.0;
      real_t tmp_102 = tmp_0*tmp_99 + tmp_1*tmp_100 + tmp_101*tmp_4;
      real_t tmp_103 = tmp_100*tmp_5 + tmp_101*tmp_2 + tmp_11*tmp_99;
      real_t tmp_104 = tmp_10*tmp_99 + tmp_100*tmp_12 + tmp_101*tmp_8;
      real_t a_0_0 = 0.11169079483900572*tmp_59*(-tmp_102*tmp_46 - tmp_103*tmp_51 - tmp_104*tmp_56 + 3.0*tmp_58*((tmp_102*tmp_102) + (tmp_103*tmp_103) + (tmp_104*tmp_104))) + 0.054975871827660928*tmp_59*(-tmp_42*tmp_46 - tmp_47*tmp_51 - tmp_52*tmp_56 + 3.0*tmp_58*((tmp_42*tmp_42) + (tmp_47*tmp_47) + (tmp_52*tmp_52))) + 0.11169079483900572*tmp_59*(-tmp_46*tmp_66 - tmp_51*tmp_67 - tmp_56*tmp_68 + 3.0*tmp_58*((tmp_66*tmp_66) + (tmp_67*tmp_67) + (tmp_68*tmp_68))) + 0.054975871827660928*tmp_59*(-tmp_46*tmp_75 - tmp_51*tmp_76 - tmp_56*tmp_77 + 3.0*tmp_58*((tmp_75*tmp_75) + (tmp_76*tmp_76) + (tmp_77*tmp_77))) + 0.11169079483900572*tmp_59*(-tmp_46*tmp_84 - tmp_51*tmp_85 - tmp_56*tmp_86 + 3.0*tmp_58*((tmp_84*tmp_84) + (tmp_85*tmp_85) + (tmp_86*tmp_86))) + 0.054975871827660928*tmp_59*(-tmp_46*tmp_93 - tmp_51*tmp_94 - tmp_56*tmp_95 + 3.0*tmp_58*((tmp_93*tmp_93) + (tmp_94*tmp_94) + (tmp_95*tmp_95)));
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
                                  MatrixXr&                            elMat ) const override
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


      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = tmp_2 - tmp_5;
      real_t tmp_7 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_8 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_9 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_10 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_11 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = tmp_1*tmp_11;
      real_t tmp_14 = tmp_0*tmp_9;
      real_t tmp_15 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_14 - tmp_13*tmp_8 + tmp_2*tmp_7 + tmp_3*tmp_8*tmp_9 - tmp_5*tmp_7);
      real_t tmp_16 = tmp_15*tmp_7;
      real_t tmp_17 = -tmp_13 + tmp_3*tmp_9;
      real_t tmp_18 = tmp_15*tmp_8;
      real_t tmp_19 = tmp_12 - tmp_14;
      real_t tmp_20 = tmp_10*tmp_15;
      real_t tmp_21 = -tmp_1*tmp_8 + tmp_10*tmp_4;
      real_t tmp_22 = tmp_1*tmp_7 - tmp_10*tmp_9;
      real_t tmp_23 = -tmp_4*tmp_7 + tmp_8*tmp_9;
      real_t tmp_24 = -tmp_0*tmp_10 + tmp_3*tmp_8;
      real_t tmp_25 = tmp_10*tmp_11 - tmp_3*tmp_7;
      real_t tmp_26 = tmp_0*tmp_7 - tmp_11*tmp_8;
      real_t tmp_27 = p_affine_13_0*(tmp_16*tmp_6 + tmp_17*tmp_18 + tmp_19*tmp_20) + p_affine_13_1*(tmp_16*tmp_21 + tmp_18*tmp_22 + tmp_20*tmp_23) + p_affine_13_2*(tmp_16*tmp_24 + tmp_18*tmp_25 + tmp_20*tmp_26);
      real_t tmp_28 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_29 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_30 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_31 = tmp_29*tmp_30;
      real_t tmp_32 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_33 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_34 = tmp_32*tmp_33;
      real_t tmp_35 = tmp_31 - tmp_34;
      real_t tmp_36 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_37 = tmp_33*tmp_36;
      real_t tmp_38 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_39 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_40 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_41 = tmp_30*tmp_40;
      real_t tmp_42 = tmp_29*tmp_36;
      real_t tmp_43 = 1.0 / (tmp_28*tmp_37 - tmp_28*tmp_41 + tmp_31*tmp_38 + tmp_32*tmp_39*tmp_40 - tmp_34*tmp_38 - tmp_39*tmp_42);
      real_t tmp_44 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_45 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_46 = -tmp_45;
      real_t tmp_47 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_48 = 0.091576213509770743*tmp_46 + 0.81684757298045851*tmp_47;
      real_t tmp_49 = tmp_43*(tmp_44 + tmp_48);
      real_t tmp_50 = tmp_32*tmp_40 - tmp_42;
      real_t tmp_51 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = -tmp_52;
      real_t tmp_54 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_55 = 0.091576213509770743*tmp_53 + 0.81684757298045851*tmp_54;
      real_t tmp_56 = tmp_43*(tmp_51 + tmp_55);
      real_t tmp_57 = tmp_37 - tmp_41;
      real_t tmp_58 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_59 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_60 = -tmp_59;
      real_t tmp_61 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_62 = 0.091576213509770743*tmp_60 + 0.81684757298045851*tmp_61;
      real_t tmp_63 = tmp_43*(tmp_58 + tmp_62);
      real_t tmp_64 = tmp_35*tmp_49 + tmp_50*tmp_56 + tmp_57*tmp_63 - 1.0/4.0;
      real_t tmp_65 = -tmp_28*tmp_30 + tmp_32*tmp_39;
      real_t tmp_66 = tmp_28*tmp_36 - tmp_32*tmp_38;
      real_t tmp_67 = tmp_30*tmp_38 - tmp_36*tmp_39;
      real_t tmp_68 = tmp_49*tmp_65 + tmp_56*tmp_66 + tmp_63*tmp_67 - 1.0/4.0;
      real_t tmp_69 = tmp_28*tmp_33 - tmp_29*tmp_39;
      real_t tmp_70 = -tmp_28*tmp_40 + tmp_29*tmp_38;
      real_t tmp_71 = -tmp_33*tmp_38 + tmp_39*tmp_40;
      real_t tmp_72 = tmp_49*tmp_69 + tmp_56*tmp_70 + tmp_63*tmp_71 - 1.0/4.0;
      real_t tmp_73 = tmp_28*tmp_64 + tmp_29*tmp_68 + tmp_32*tmp_72;
      real_t tmp_74 = tmp_11*tmp_15;
      real_t tmp_75 = tmp_0*tmp_15;
      real_t tmp_76 = tmp_15*tmp_3;
      real_t tmp_77 = p_affine_13_0*(tmp_17*tmp_75 + tmp_19*tmp_76 + tmp_6*tmp_74) + p_affine_13_1*(tmp_21*tmp_74 + tmp_22*tmp_75 + tmp_23*tmp_76) + p_affine_13_2*(tmp_24*tmp_74 + tmp_25*tmp_75 + tmp_26*tmp_76);
      real_t tmp_78 = tmp_30*tmp_72 + tmp_33*tmp_68 + tmp_39*tmp_64;
      real_t tmp_79 = tmp_15*tmp_9;
      real_t tmp_80 = tmp_15*tmp_4;
      real_t tmp_81 = tmp_1*tmp_15;
      real_t tmp_82 = p_affine_13_0*(tmp_17*tmp_80 + tmp_19*tmp_81 + tmp_6*tmp_79) + p_affine_13_1*(tmp_21*tmp_79 + tmp_22*tmp_80 + tmp_23*tmp_81) + p_affine_13_2*(tmp_24*tmp_79 + tmp_25*tmp_80 + tmp_26*tmp_81);
      real_t tmp_83 = tmp_36*tmp_72 + tmp_38*tmp_64 + tmp_40*tmp_68;
      real_t tmp_84 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_85 = tmp_15*(tmp_48 + tmp_84);
      real_t tmp_86 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_87 = tmp_15*(tmp_55 + tmp_86);
      real_t tmp_88 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_89 = tmp_15*(tmp_62 + tmp_88);
      real_t tmp_90 = tmp_21*tmp_87 + tmp_24*tmp_85 + tmp_6*tmp_89 - 1.0/4.0;
      real_t tmp_91 = tmp_17*tmp_89 + tmp_22*tmp_87 + tmp_25*tmp_85 - 1.0/4.0;
      real_t tmp_92 = tmp_19*tmp_89 + tmp_23*tmp_87 + tmp_26*tmp_85 - 1.0/4.0;
      real_t tmp_93 = tmp_10*tmp_92 + tmp_7*tmp_90 + tmp_8*tmp_91;
      real_t tmp_94 = tmp_28*tmp_43;
      real_t tmp_95 = tmp_29*tmp_43;
      real_t tmp_96 = tmp_32*tmp_43;
      real_t tmp_97 = 0.5*p_affine_13_0*(tmp_57*tmp_94 + tmp_67*tmp_95 + tmp_71*tmp_96) + 0.5*p_affine_13_1*(tmp_50*tmp_94 + tmp_66*tmp_95 + tmp_70*tmp_96) + 0.5*p_affine_13_2*(tmp_35*tmp_94 + tmp_65*tmp_95 + tmp_69*tmp_96);
      real_t tmp_98 = tmp_0*tmp_91 + tmp_11*tmp_90 + tmp_3*tmp_92;
      real_t tmp_99 = tmp_39*tmp_43;
      real_t tmp_100 = tmp_33*tmp_43;
      real_t tmp_101 = tmp_30*tmp_43;
      real_t tmp_102 = 0.5*p_affine_13_0*(tmp_100*tmp_67 + tmp_101*tmp_71 + tmp_57*tmp_99) + 0.5*p_affine_13_1*(tmp_100*tmp_66 + tmp_101*tmp_70 + tmp_50*tmp_99) + 0.5*p_affine_13_2*(tmp_100*tmp_65 + tmp_101*tmp_69 + tmp_35*tmp_99);
      real_t tmp_103 = tmp_1*tmp_92 + tmp_4*tmp_91 + tmp_9*tmp_90;
      real_t tmp_104 = tmp_38*tmp_43;
      real_t tmp_105 = tmp_40*tmp_43;
      real_t tmp_106 = tmp_36*tmp_43;
      real_t tmp_107 = 0.5*p_affine_13_0*(tmp_104*tmp_57 + tmp_105*tmp_67 + tmp_106*tmp_71) + 0.5*p_affine_13_1*(tmp_104*tmp_50 + tmp_105*tmp_66 + tmp_106*tmp_70) + 0.5*p_affine_13_2*(tmp_104*tmp_35 + tmp_105*tmp_65 + tmp_106*tmp_69);
      real_t tmp_108 = (std::abs(tmp_45*tmp_54 - tmp_47*tmp_52)*std::abs(tmp_45*tmp_54 - tmp_47*tmp_52)) + (std::abs(tmp_45*tmp_61 - tmp_47*tmp_59)*std::abs(tmp_45*tmp_61 - tmp_47*tmp_59)) + (std::abs(tmp_52*tmp_61 - tmp_54*tmp_59)*std::abs(tmp_52*tmp_61 - tmp_54*tmp_59));
      real_t tmp_109 = 3.0*std::pow(tmp_108, -0.25);
      real_t tmp_110 = 1.0*std::pow(tmp_108, 1.0/2.0);
      real_t tmp_111 = 0.44594849091596489*tmp_46 + 0.10810301816807022*tmp_47;
      real_t tmp_112 = tmp_43*(tmp_111 + tmp_44);
      real_t tmp_113 = 0.44594849091596489*tmp_53 + 0.10810301816807022*tmp_54;
      real_t tmp_114 = tmp_43*(tmp_113 + tmp_51);
      real_t tmp_115 = 0.44594849091596489*tmp_60 + 0.10810301816807022*tmp_61;
      real_t tmp_116 = tmp_43*(tmp_115 + tmp_58);
      real_t tmp_117 = tmp_112*tmp_35 + tmp_114*tmp_50 + tmp_116*tmp_57 - 1.0/4.0;
      real_t tmp_118 = tmp_112*tmp_65 + tmp_114*tmp_66 + tmp_116*tmp_67 - 1.0/4.0;
      real_t tmp_119 = tmp_112*tmp_69 + tmp_114*tmp_70 + tmp_116*tmp_71 - 1.0/4.0;
      real_t tmp_120 = tmp_117*tmp_28 + tmp_118*tmp_29 + tmp_119*tmp_32;
      real_t tmp_121 = tmp_117*tmp_39 + tmp_118*tmp_33 + tmp_119*tmp_30;
      real_t tmp_122 = tmp_117*tmp_38 + tmp_118*tmp_40 + tmp_119*tmp_36;
      real_t tmp_123 = tmp_15*(tmp_111 + tmp_84);
      real_t tmp_124 = tmp_15*(tmp_113 + tmp_86);
      real_t tmp_125 = tmp_15*(tmp_115 + tmp_88);
      real_t tmp_126 = tmp_123*tmp_24 + tmp_124*tmp_21 + tmp_125*tmp_6 - 1.0/4.0;
      real_t tmp_127 = tmp_123*tmp_25 + tmp_124*tmp_22 + tmp_125*tmp_17 - 1.0/4.0;
      real_t tmp_128 = tmp_123*tmp_26 + tmp_124*tmp_23 + tmp_125*tmp_19 - 1.0/4.0;
      real_t tmp_129 = tmp_10*tmp_128 + tmp_126*tmp_7 + tmp_127*tmp_8;
      real_t tmp_130 = tmp_0*tmp_127 + tmp_11*tmp_126 + tmp_128*tmp_3;
      real_t tmp_131 = tmp_1*tmp_128 + tmp_126*tmp_9 + tmp_127*tmp_4;
      real_t tmp_132 = 0.81684757298045851*tmp_46 + 0.091576213509770743*tmp_47;
      real_t tmp_133 = tmp_43*(tmp_132 + tmp_44);
      real_t tmp_134 = 0.81684757298045851*tmp_53 + 0.091576213509770743*tmp_54;
      real_t tmp_135 = tmp_43*(tmp_134 + tmp_51);
      real_t tmp_136 = 0.81684757298045851*tmp_60 + 0.091576213509770743*tmp_61;
      real_t tmp_137 = tmp_43*(tmp_136 + tmp_58);
      real_t tmp_138 = tmp_133*tmp_35 + tmp_135*tmp_50 + tmp_137*tmp_57 - 1.0/4.0;
      real_t tmp_139 = tmp_133*tmp_65 + tmp_135*tmp_66 + tmp_137*tmp_67 - 1.0/4.0;
      real_t tmp_140 = tmp_133*tmp_69 + tmp_135*tmp_70 + tmp_137*tmp_71 - 1.0/4.0;
      real_t tmp_141 = tmp_138*tmp_28 + tmp_139*tmp_29 + tmp_140*tmp_32;
      real_t tmp_142 = tmp_138*tmp_39 + tmp_139*tmp_33 + tmp_140*tmp_30;
      real_t tmp_143 = tmp_138*tmp_38 + tmp_139*tmp_40 + tmp_140*tmp_36;
      real_t tmp_144 = tmp_15*(tmp_132 + tmp_84);
      real_t tmp_145 = tmp_15*(tmp_134 + tmp_86);
      real_t tmp_146 = tmp_15*(tmp_136 + tmp_88);
      real_t tmp_147 = tmp_144*tmp_24 + tmp_145*tmp_21 + tmp_146*tmp_6 - 1.0/4.0;
      real_t tmp_148 = tmp_144*tmp_25 + tmp_145*tmp_22 + tmp_146*tmp_17 - 1.0/4.0;
      real_t tmp_149 = tmp_144*tmp_26 + tmp_145*tmp_23 + tmp_146*tmp_19 - 1.0/4.0;
      real_t tmp_150 = tmp_10*tmp_149 + tmp_147*tmp_7 + tmp_148*tmp_8;
      real_t tmp_151 = tmp_0*tmp_148 + tmp_11*tmp_147 + tmp_149*tmp_3;
      real_t tmp_152 = tmp_1*tmp_149 + tmp_147*tmp_9 + tmp_148*tmp_4;
      real_t tmp_153 = 0.10810301816807022*tmp_46 + 0.44594849091596489*tmp_47;
      real_t tmp_154 = tmp_43*(tmp_153 + tmp_44);
      real_t tmp_155 = 0.10810301816807022*tmp_53 + 0.44594849091596489*tmp_54;
      real_t tmp_156 = tmp_43*(tmp_155 + tmp_51);
      real_t tmp_157 = 0.10810301816807022*tmp_60 + 0.44594849091596489*tmp_61;
      real_t tmp_158 = tmp_43*(tmp_157 + tmp_58);
      real_t tmp_159 = tmp_154*tmp_35 + tmp_156*tmp_50 + tmp_158*tmp_57 - 1.0/4.0;
      real_t tmp_160 = tmp_154*tmp_65 + tmp_156*tmp_66 + tmp_158*tmp_67 - 1.0/4.0;
      real_t tmp_161 = tmp_154*tmp_69 + tmp_156*tmp_70 + tmp_158*tmp_71 - 1.0/4.0;
      real_t tmp_162 = tmp_159*tmp_28 + tmp_160*tmp_29 + tmp_161*tmp_32;
      real_t tmp_163 = tmp_159*tmp_39 + tmp_160*tmp_33 + tmp_161*tmp_30;
      real_t tmp_164 = tmp_159*tmp_38 + tmp_160*tmp_40 + tmp_161*tmp_36;
      real_t tmp_165 = tmp_15*(tmp_153 + tmp_84);
      real_t tmp_166 = tmp_15*(tmp_155 + tmp_86);
      real_t tmp_167 = tmp_15*(tmp_157 + tmp_88);
      real_t tmp_168 = tmp_165*tmp_24 + tmp_166*tmp_21 + tmp_167*tmp_6 - 1.0/4.0;
      real_t tmp_169 = tmp_165*tmp_25 + tmp_166*tmp_22 + tmp_167*tmp_17 - 1.0/4.0;
      real_t tmp_170 = tmp_165*tmp_26 + tmp_166*tmp_23 + tmp_167*tmp_19 - 1.0/4.0;
      real_t tmp_171 = tmp_10*tmp_170 + tmp_168*tmp_7 + tmp_169*tmp_8;
      real_t tmp_172 = tmp_0*tmp_169 + tmp_11*tmp_168 + tmp_170*tmp_3;
      real_t tmp_173 = tmp_1*tmp_170 + tmp_168*tmp_9 + tmp_169*tmp_4;
      real_t tmp_174 = 0.091576213509770743*tmp_46 + 0.091576213509770743*tmp_47;
      real_t tmp_175 = tmp_43*(tmp_174 + tmp_44);
      real_t tmp_176 = 0.091576213509770743*tmp_53 + 0.091576213509770743*tmp_54;
      real_t tmp_177 = tmp_43*(tmp_176 + tmp_51);
      real_t tmp_178 = 0.091576213509770743*tmp_60 + 0.091576213509770743*tmp_61;
      real_t tmp_179 = tmp_43*(tmp_178 + tmp_58);
      real_t tmp_180 = tmp_175*tmp_35 + tmp_177*tmp_50 + tmp_179*tmp_57 - 1.0/4.0;
      real_t tmp_181 = tmp_175*tmp_65 + tmp_177*tmp_66 + tmp_179*tmp_67 - 1.0/4.0;
      real_t tmp_182 = tmp_175*tmp_69 + tmp_177*tmp_70 + tmp_179*tmp_71 - 1.0/4.0;
      real_t tmp_183 = tmp_180*tmp_28 + tmp_181*tmp_29 + tmp_182*tmp_32;
      real_t tmp_184 = tmp_180*tmp_39 + tmp_181*tmp_33 + tmp_182*tmp_30;
      real_t tmp_185 = tmp_180*tmp_38 + tmp_181*tmp_40 + tmp_182*tmp_36;
      real_t tmp_186 = tmp_15*(tmp_174 + tmp_84);
      real_t tmp_187 = tmp_15*(tmp_176 + tmp_86);
      real_t tmp_188 = tmp_15*(tmp_178 + tmp_88);
      real_t tmp_189 = tmp_186*tmp_24 + tmp_187*tmp_21 + tmp_188*tmp_6 - 1.0/4.0;
      real_t tmp_190 = tmp_17*tmp_188 + tmp_186*tmp_25 + tmp_187*tmp_22 - 1.0/4.0;
      real_t tmp_191 = tmp_186*tmp_26 + tmp_187*tmp_23 + tmp_188*tmp_19 - 1.0/4.0;
      real_t tmp_192 = tmp_10*tmp_191 + tmp_189*tmp_7 + tmp_190*tmp_8;
      real_t tmp_193 = tmp_0*tmp_190 + tmp_11*tmp_189 + tmp_191*tmp_3;
      real_t tmp_194 = tmp_1*tmp_191 + tmp_189*tmp_9 + tmp_190*tmp_4;
      real_t tmp_195 = 0.44594849091596489*tmp_46 + 0.44594849091596489*tmp_47;
      real_t tmp_196 = tmp_43*(tmp_195 + tmp_44);
      real_t tmp_197 = 0.44594849091596489*tmp_53 + 0.44594849091596489*tmp_54;
      real_t tmp_198 = tmp_43*(tmp_197 + tmp_51);
      real_t tmp_199 = 0.44594849091596489*tmp_60 + 0.44594849091596489*tmp_61;
      real_t tmp_200 = tmp_43*(tmp_199 + tmp_58);
      real_t tmp_201 = tmp_196*tmp_35 + tmp_198*tmp_50 + tmp_200*tmp_57 - 1.0/4.0;
      real_t tmp_202 = tmp_196*tmp_65 + tmp_198*tmp_66 + tmp_200*tmp_67 - 1.0/4.0;
      real_t tmp_203 = tmp_196*tmp_69 + tmp_198*tmp_70 + tmp_200*tmp_71 - 1.0/4.0;
      real_t tmp_204 = tmp_201*tmp_28 + tmp_202*tmp_29 + tmp_203*tmp_32;
      real_t tmp_205 = tmp_201*tmp_39 + tmp_202*tmp_33 + tmp_203*tmp_30;
      real_t tmp_206 = tmp_201*tmp_38 + tmp_202*tmp_40 + tmp_203*tmp_36;
      real_t tmp_207 = tmp_15*(tmp_195 + tmp_84);
      real_t tmp_208 = tmp_15*(tmp_197 + tmp_86);
      real_t tmp_209 = tmp_15*(tmp_199 + tmp_88);
      real_t tmp_210 = tmp_207*tmp_24 + tmp_208*tmp_21 + tmp_209*tmp_6 - 1.0/4.0;
      real_t tmp_211 = tmp_17*tmp_209 + tmp_207*tmp_25 + tmp_208*tmp_22 - 1.0/4.0;
      real_t tmp_212 = tmp_19*tmp_209 + tmp_207*tmp_26 + tmp_208*tmp_23 - 1.0/4.0;
      real_t tmp_213 = tmp_10*tmp_212 + tmp_210*tmp_7 + tmp_211*tmp_8;
      real_t tmp_214 = tmp_0*tmp_211 + tmp_11*tmp_210 + tmp_212*tmp_3;
      real_t tmp_215 = tmp_1*tmp_212 + tmp_210*tmp_9 + tmp_211*tmp_4;
      real_t a_0_0 = 0.11169079483900572*tmp_110*(-tmp_102*tmp_130 - tmp_107*tmp_131 - tmp_109*(tmp_120*tmp_129 + tmp_121*tmp_130 + tmp_122*tmp_131) + 0.5*tmp_120*tmp_27 + 0.5*tmp_121*tmp_77 + 0.5*tmp_122*tmp_82 - tmp_129*tmp_97) + 0.054975871827660928*tmp_110*(-tmp_102*tmp_151 - tmp_107*tmp_152 - tmp_109*(tmp_141*tmp_150 + tmp_142*tmp_151 + tmp_143*tmp_152) + 0.5*tmp_141*tmp_27 + 0.5*tmp_142*tmp_77 + 0.5*tmp_143*tmp_82 - tmp_150*tmp_97) + 0.11169079483900572*tmp_110*(-tmp_102*tmp_172 - tmp_107*tmp_173 - tmp_109*(tmp_162*tmp_171 + tmp_163*tmp_172 + tmp_164*tmp_173) + 0.5*tmp_162*tmp_27 + 0.5*tmp_163*tmp_77 + 0.5*tmp_164*tmp_82 - tmp_171*tmp_97) + 0.054975871827660928*tmp_110*(-tmp_102*tmp_193 - tmp_107*tmp_194 - tmp_109*(tmp_183*tmp_192 + tmp_184*tmp_193 + tmp_185*tmp_194) + 0.5*tmp_183*tmp_27 + 0.5*tmp_184*tmp_77 + 0.5*tmp_185*tmp_82 - tmp_192*tmp_97) + 0.11169079483900572*tmp_110*(-tmp_102*tmp_214 - tmp_107*tmp_215 - tmp_109*(tmp_204*tmp_213 + tmp_205*tmp_214 + tmp_206*tmp_215) + 0.5*tmp_204*tmp_27 + 0.5*tmp_205*tmp_77 + 0.5*tmp_206*tmp_82 - tmp_213*tmp_97) + 0.054975871827660928*tmp_110*(-tmp_102*tmp_98 - tmp_103*tmp_107 - tmp_109*(tmp_103*tmp_83 + tmp_73*tmp_93 + tmp_78*tmp_98) + 0.5*tmp_27*tmp_73 + 0.5*tmp_77*tmp_78 + 0.5*tmp_82*tmp_83 - tmp_93*tmp_97);
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
                                        MatrixXr&                            elMat ) const override
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
      real_t tmp_43 = tmp_0*tmp_15;
      real_t tmp_44 = tmp_1*tmp_15;
      real_t tmp_45 = tmp_15*tmp_3;
      real_t tmp_46 = 2*p_affine_13_0*(tmp_27*tmp_43 + tmp_36*tmp_44 + tmp_40*tmp_45) + 2*p_affine_13_1*(tmp_21*tmp_43 + tmp_35*tmp_44 + tmp_39*tmp_45) + 2*p_affine_13_2*(tmp_34*tmp_44 + tmp_38*tmp_45 + tmp_43*tmp_5);
      real_t tmp_47 = tmp_2*tmp_41 + tmp_33*tmp_9 + tmp_37*tmp_4;
      real_t tmp_48 = tmp_15*tmp_9;
      real_t tmp_49 = tmp_15*tmp_4;
      real_t tmp_50 = tmp_15*tmp_2;
      real_t tmp_51 = 2*p_affine_13_0*(tmp_27*tmp_48 + tmp_36*tmp_49 + tmp_40*tmp_50) + 2*p_affine_13_1*(tmp_21*tmp_48 + tmp_35*tmp_49 + tmp_39*tmp_50) + 2*p_affine_13_2*(tmp_34*tmp_49 + tmp_38*tmp_50 + tmp_48*tmp_5);
      real_t tmp_52 = tmp_10*tmp_37 + tmp_33*tmp_8 + tmp_41*tmp_6;
      real_t tmp_53 = tmp_15*tmp_8;
      real_t tmp_54 = tmp_10*tmp_15;
      real_t tmp_55 = tmp_15*tmp_6;
      real_t tmp_56 = 2*p_affine_13_0*(tmp_27*tmp_53 + tmp_36*tmp_54 + tmp_40*tmp_55) + 2*p_affine_13_1*(tmp_21*tmp_53 + tmp_35*tmp_54 + tmp_39*tmp_55) + 2*p_affine_13_2*(tmp_34*tmp_54 + tmp_38*tmp_55 + tmp_5*tmp_53);
      real_t tmp_57 = (std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28));
      real_t tmp_58 = std::pow(tmp_57, -0.25);
      real_t tmp_59 = 1.0*std::pow(tmp_57, 1.0/2.0);
      real_t tmp_60 = tmp_15*(0.44594849091596489*tmp_17 + 0.10810301816807022*tmp_18 + tmp_19);
      real_t tmp_61 = tmp_15*(0.44594849091596489*tmp_23 + 0.10810301816807022*tmp_24 + tmp_25);
      real_t tmp_62 = tmp_15*(0.44594849091596489*tmp_29 + 0.10810301816807022*tmp_30 + tmp_31);
      real_t tmp_63 = tmp_21*tmp_61 + tmp_27*tmp_62 + tmp_5*tmp_60 - 1.0/4.0;
      real_t tmp_64 = tmp_34*tmp_60 + tmp_35*tmp_61 + tmp_36*tmp_62 - 1.0/4.0;
      real_t tmp_65 = tmp_38*tmp_60 + tmp_39*tmp_61 + tmp_40*tmp_62 - 1.0/4.0;
      real_t tmp_66 = tmp_0*tmp_63 + tmp_1*tmp_64 + tmp_3*tmp_65;
      real_t tmp_67 = tmp_2*tmp_65 + tmp_4*tmp_64 + tmp_63*tmp_9;
      real_t tmp_68 = tmp_10*tmp_64 + tmp_6*tmp_65 + tmp_63*tmp_8;
      real_t tmp_69 = tmp_15*(0.81684757298045851*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_70 = tmp_15*(0.81684757298045851*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_71 = tmp_15*(0.81684757298045851*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_72 = tmp_21*tmp_70 + tmp_27*tmp_71 + tmp_5*tmp_69 - 1.0/4.0;
      real_t tmp_73 = tmp_34*tmp_69 + tmp_35*tmp_70 + tmp_36*tmp_71 - 1.0/4.0;
      real_t tmp_74 = tmp_38*tmp_69 + tmp_39*tmp_70 + tmp_40*tmp_71 - 1.0/4.0;
      real_t tmp_75 = tmp_0*tmp_72 + tmp_1*tmp_73 + tmp_3*tmp_74;
      real_t tmp_76 = tmp_2*tmp_74 + tmp_4*tmp_73 + tmp_72*tmp_9;
      real_t tmp_77 = tmp_10*tmp_73 + tmp_6*tmp_74 + tmp_72*tmp_8;
      real_t tmp_78 = tmp_15*(0.10810301816807022*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_79 = tmp_15*(0.10810301816807022*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_80 = tmp_15*(0.10810301816807022*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_81 = tmp_21*tmp_79 + tmp_27*tmp_80 + tmp_5*tmp_78 - 1.0/4.0;
      real_t tmp_82 = tmp_34*tmp_78 + tmp_35*tmp_79 + tmp_36*tmp_80 - 1.0/4.0;
      real_t tmp_83 = tmp_38*tmp_78 + tmp_39*tmp_79 + tmp_40*tmp_80 - 1.0/4.0;
      real_t tmp_84 = tmp_0*tmp_81 + tmp_1*tmp_82 + tmp_3*tmp_83;
      real_t tmp_85 = tmp_2*tmp_83 + tmp_4*tmp_82 + tmp_81*tmp_9;
      real_t tmp_86 = tmp_10*tmp_82 + tmp_6*tmp_83 + tmp_8*tmp_81;
      real_t tmp_87 = tmp_15*(0.091576213509770743*tmp_17 + 0.091576213509770743*tmp_18 + tmp_19);
      real_t tmp_88 = tmp_15*(0.091576213509770743*tmp_23 + 0.091576213509770743*tmp_24 + tmp_25);
      real_t tmp_89 = tmp_15*(0.091576213509770743*tmp_29 + 0.091576213509770743*tmp_30 + tmp_31);
      real_t tmp_90 = tmp_21*tmp_88 + tmp_27*tmp_89 + tmp_5*tmp_87 - 1.0/4.0;
      real_t tmp_91 = tmp_34*tmp_87 + tmp_35*tmp_88 + tmp_36*tmp_89 - 1.0/4.0;
      real_t tmp_92 = tmp_38*tmp_87 + tmp_39*tmp_88 + tmp_40*tmp_89 - 1.0/4.0;
      real_t tmp_93 = tmp_0*tmp_90 + tmp_1*tmp_91 + tmp_3*tmp_92;
      real_t tmp_94 = tmp_2*tmp_92 + tmp_4*tmp_91 + tmp_9*tmp_90;
      real_t tmp_95 = tmp_10*tmp_91 + tmp_6*tmp_92 + tmp_8*tmp_90;
      real_t tmp_96 = tmp_15*(0.44594849091596489*tmp_17 + 0.44594849091596489*tmp_18 + tmp_19);
      real_t tmp_97 = tmp_15*(0.44594849091596489*tmp_23 + 0.44594849091596489*tmp_24 + tmp_25);
      real_t tmp_98 = tmp_15*(0.44594849091596489*tmp_29 + 0.44594849091596489*tmp_30 + tmp_31);
      real_t tmp_99 = tmp_21*tmp_97 + tmp_27*tmp_98 + tmp_5*tmp_96 - 1.0/4.0;
      real_t tmp_100 = tmp_34*tmp_96 + tmp_35*tmp_97 + tmp_36*tmp_98 - 1.0/4.0;
      real_t tmp_101 = tmp_38*tmp_96 + tmp_39*tmp_97 + tmp_40*tmp_98 - 1.0/4.0;
      real_t tmp_102 = tmp_0*tmp_99 + tmp_1*tmp_100 + tmp_101*tmp_3;
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_2 + tmp_9*tmp_99;
      real_t tmp_104 = tmp_10*tmp_100 + tmp_101*tmp_6 + tmp_8*tmp_99;
      real_t a_0_0 = 0.11169079483900572*tmp_59*(-tmp_102*tmp_46 - tmp_103*tmp_51 - tmp_104*tmp_56 + 3.0*tmp_58*((tmp_102*tmp_102) + (tmp_103*tmp_103) + (tmp_104*tmp_104))) + 0.054975871827660928*tmp_59*(-tmp_42*tmp_46 - tmp_47*tmp_51 - tmp_52*tmp_56 + 3.0*tmp_58*((tmp_42*tmp_42) + (tmp_47*tmp_47) + (tmp_52*tmp_52))) + 0.11169079483900572*tmp_59*(-tmp_46*tmp_66 - tmp_51*tmp_67 - tmp_56*tmp_68 + 3.0*tmp_58*((tmp_66*tmp_66) + (tmp_67*tmp_67) + (tmp_68*tmp_68))) + 0.054975871827660928*tmp_59*(-tmp_46*tmp_75 - tmp_51*tmp_76 - tmp_56*tmp_77 + 3.0*tmp_58*((tmp_75*tmp_75) + (tmp_76*tmp_76) + (tmp_77*tmp_77))) + 0.11169079483900572*tmp_59*(-tmp_46*tmp_84 - tmp_51*tmp_85 - tmp_56*tmp_86 + 3.0*tmp_58*((tmp_84*tmp_84) + (tmp_85*tmp_85) + (tmp_86*tmp_86))) + 0.054975871827660928*tmp_59*(-tmp_46*tmp_93 - tmp_51*tmp_94 - tmp_56*tmp_95 + 3.0*tmp_58*((tmp_93*tmp_93) + (tmp_94*tmp_94) + (tmp_95*tmp_95)));
      elMat( 0, 0) = a_0_0;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g2;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;

};


} //eg
} // dg
} // hyteg
