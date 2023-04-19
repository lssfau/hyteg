
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

class EGNIPGVectorLaplaceFormP1E_0 : public hyteg::dg::DGForm
{
 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_1_1 + tmp_2));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = tmp_1*tmp_8 + tmp_5*tmp_7;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_11*tmp_5 + tmp_4*tmp_6;
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

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_1_1 + tmp_2));
      real_t tmp_6 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7 = p_affine_6_1 + tmp_2;
      real_t tmp_8 = tmp_5*(0.046910077030668018*tmp_6 + tmp_7);
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5*(0.046910077030668018*tmp_12 + tmp_13);
      real_t tmp_15 = tmp_14*tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 0.35539032758428413*tmp_1*(tmp_18 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_19 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 0.71794300574904923*tmp_1*(tmp_27 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_28 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 0.8533333333333335*tmp_1*(tmp_36 - 1.0/3.0) + 0.8533333333333335*tmp_4*(tmp_37 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 0.71794300574904923*tmp_1*(tmp_45 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_46 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.35539032758428413*tmp_1*(tmp_54 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_55 - 1.0/3.0);
      real_t a_0_0 = tmp_20*(-tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1) + tmp_29*(-tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1) + tmp_38*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_47*(-tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1) + tmp_56*(-tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1);
      real_t a_1_0 = tmp_18*tmp_20 + tmp_27*tmp_29 + tmp_36*tmp_38 + tmp_45*tmp_47 + tmp_54*tmp_56;
      real_t a_2_0 = tmp_19*tmp_20 + tmp_28*tmp_29 + tmp_37*tmp_38 + tmp_46*tmp_47 + tmp_55*tmp_56;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + 0.046910077030668018*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -p_affine_3_0;
      real_t tmp_18 = p_affine_4_0 + tmp_17;
      real_t tmp_19 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_20 = -p_affine_3_1;
      real_t tmp_21 = p_affine_5_1 + tmp_20;
      real_t tmp_22 = p_affine_5_0 + tmp_17;
      real_t tmp_23 = 1.0 / (tmp_18*tmp_21 - tmp_22*(p_affine_4_1 + tmp_20));
      real_t tmp_24 = tmp_23*(tmp_20 + tmp_6);
      real_t tmp_25 = tmp_23*(tmp_12 + tmp_17);
      real_t tmp_26 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_27 = 0.35539032758428413*tmp_18*(tmp_19*tmp_24 + tmp_21*tmp_25 - 1.0/3.0) + 0.35539032758428413*tmp_22*(tmp_18*tmp_24 + tmp_25*tmp_26 - 1.0/3.0);
      real_t tmp_28 = p_affine_6_1 + 0.23076534494715845*tmp_5;
      real_t tmp_29 = tmp_4*(tmp_2 + tmp_28);
      real_t tmp_30 = tmp_1*tmp_29;
      real_t tmp_31 = tmp_29*tmp_9;
      real_t tmp_32 = p_affine_6_0 + 0.23076534494715845*tmp_11;
      real_t tmp_33 = tmp_4*(tmp_0 + tmp_32);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_15*tmp_33;
      real_t tmp_36 = tmp_23*(tmp_20 + tmp_28);
      real_t tmp_37 = tmp_23*(tmp_17 + tmp_32);
      real_t tmp_38 = 0.71794300574904923*tmp_18*(tmp_19*tmp_36 + tmp_21*tmp_37 - 1.0/3.0) + 0.71794300574904923*tmp_22*(tmp_18*tmp_36 + tmp_26*tmp_37 - 1.0/3.0);
      real_t tmp_39 = p_affine_6_1 + 0.5*tmp_5;
      real_t tmp_40 = tmp_4*(tmp_2 + tmp_39);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = p_affine_6_0 + 0.5*tmp_11;
      real_t tmp_44 = tmp_4*(tmp_0 + tmp_43);
      real_t tmp_45 = tmp_3*tmp_44;
      real_t tmp_46 = tmp_15*tmp_44;
      real_t tmp_47 = tmp_23*(tmp_20 + tmp_39);
      real_t tmp_48 = tmp_23*(tmp_17 + tmp_43);
      real_t tmp_49 = 0.8533333333333335*tmp_18*(tmp_19*tmp_47 + tmp_21*tmp_48 - 1.0/3.0) + 0.8533333333333335*tmp_22*(tmp_18*tmp_47 + tmp_26*tmp_48 - 1.0/3.0);
      real_t tmp_50 = p_affine_6_1 + 0.7692346550528415*tmp_5;
      real_t tmp_51 = tmp_4*(tmp_2 + tmp_50);
      real_t tmp_52 = tmp_1*tmp_51;
      real_t tmp_53 = tmp_51*tmp_9;
      real_t tmp_54 = p_affine_6_0 + 0.7692346550528415*tmp_11;
      real_t tmp_55 = tmp_4*(tmp_0 + tmp_54);
      real_t tmp_56 = tmp_3*tmp_55;
      real_t tmp_57 = tmp_15*tmp_55;
      real_t tmp_58 = tmp_23*(tmp_20 + tmp_50);
      real_t tmp_59 = tmp_23*(tmp_17 + tmp_54);
      real_t tmp_60 = 0.71794300574904923*tmp_18*(tmp_19*tmp_58 + tmp_21*tmp_59 - 1.0/3.0) + 0.71794300574904923*tmp_22*(tmp_18*tmp_58 + tmp_26*tmp_59 - 1.0/3.0);
      real_t tmp_61 = p_affine_6_1 + 0.95308992296933193*tmp_5;
      real_t tmp_62 = tmp_4*(tmp_2 + tmp_61);
      real_t tmp_63 = tmp_1*tmp_62;
      real_t tmp_64 = tmp_62*tmp_9;
      real_t tmp_65 = p_affine_6_0 + 0.95308992296933193*tmp_11;
      real_t tmp_66 = tmp_4*(tmp_0 + tmp_65);
      real_t tmp_67 = tmp_3*tmp_66;
      real_t tmp_68 = tmp_15*tmp_66;
      real_t tmp_69 = tmp_23*(tmp_20 + tmp_61);
      real_t tmp_70 = tmp_23*(tmp_17 + tmp_65);
      real_t tmp_71 = 0.35539032758428413*tmp_18*(tmp_19*tmp_69 + tmp_21*tmp_70 - 1.0/3.0) + 0.35539032758428413*tmp_22*(tmp_18*tmp_69 + tmp_26*tmp_70 - 1.0/3.0);
      real_t a_0_0 = -tmp_27*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) - tmp_38*(-tmp_30 - tmp_31 - tmp_34 - tmp_35 + 1) - tmp_49*(-tmp_41 - tmp_42 - tmp_45 - tmp_46 + 1) - tmp_60*(-tmp_52 - tmp_53 - tmp_56 - tmp_57 + 1) - tmp_71*(-tmp_63 - tmp_64 - tmp_67 - tmp_68 + 1);
      real_t a_1_0 = -tmp_27*(tmp_10 + tmp_14) - tmp_38*(tmp_31 + tmp_34) - tmp_49*(tmp_42 + tmp_45) - tmp_60*(tmp_53 + tmp_56) - tmp_71*(tmp_64 + tmp_67);
      real_t a_2_0 = -tmp_27*(tmp_16 + tmp_8) - tmp_38*(tmp_30 + tmp_35) - tmp_49*(tmp_41 + tmp_46) - tmp_60*(tmp_52 + tmp_57) - tmp_71*(tmp_63 + tmp_68);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_1_1 + tmp_2));
      real_t tmp_6 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7 = p_affine_6_1 + tmp_2;
      real_t tmp_8 = tmp_5*(0.046910077030668018*tmp_6 + tmp_7);
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5*(0.046910077030668018*tmp_12 + tmp_13);
      real_t tmp_15 = tmp_14*tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 0.35539032758428413*tmp_1*(tmp_18 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_19 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 0.71794300574904923*tmp_1*(tmp_27 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_28 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 0.8533333333333335*tmp_1*(tmp_36 - 1.0/3.0) + 0.8533333333333335*tmp_4*(tmp_37 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 0.71794300574904923*tmp_1*(tmp_45 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_46 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.35539032758428413*tmp_1*(tmp_54 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_55 - 1.0/3.0);
      real_t a_0_0 = tmp_20*(-tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1) + tmp_29*(-tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1) + tmp_38*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_47*(-tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1) + tmp_56*(-tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1);
      real_t a_1_0 = tmp_18*tmp_20 + tmp_27*tmp_29 + tmp_36*tmp_38 + tmp_45*tmp_47 + tmp_54*tmp_56;
      real_t a_2_0 = tmp_19*tmp_20 + tmp_28*tmp_29 + tmp_37*tmp_38 + tmp_46*tmp_47 + tmp_55*tmp_56;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_2;
      real_t tmp_9 = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_11;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_16 + tmp_4*tmp_9 - tmp_7*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_21 = tmp_18*(tmp_12 - tmp_17);
      real_t tmp_22 = tmp_1*tmp_21 + tmp_14*tmp_19 + tmp_20*tmp_5;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_1*tmp_25 + tmp_14*tmp_23 + tmp_24*tmp_5;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_1*tmp_29 + tmp_14*tmp_27 + tmp_28*tmp_5;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_1_0 = 0.1666666666666668*tmp_45;
      real_t a_2_0 = 0.1666666666666668*tmp_46;
      real_t a_3_0 = 0.1666666666666668*tmp_47;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 3.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_52 - 1.0/4.0) + tmp_13*(tmp_54 - 1.0/4.0) + tmp_4*(tmp_53 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_72 - 1.0/4.0) + tmp_13*(tmp_74 - 1.0/4.0) + tmp_4*(tmp_73 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_88 - 1.0/4.0) + tmp_13*(tmp_90 - 1.0/4.0) + tmp_4*(tmp_89 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_104 - 1.0/4.0) + tmp_13*(tmp_106 - 1.0/4.0) + tmp_4*(tmp_105 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_120 - 1.0/4.0) + tmp_13*(tmp_122 - 1.0/4.0) + tmp_4*(tmp_121 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_136 - 1.0/4.0) + tmp_13*(tmp_138 - 1.0/4.0) + tmp_4*(tmp_137 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_152 - 1.0/4.0) + tmp_13*(tmp_154 - 1.0/4.0) + tmp_4*(tmp_153 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_168 - 1.0/4.0) + tmp_13*(tmp_170 - 1.0/4.0) + tmp_4*(tmp_169 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_184 - 1.0/4.0) + tmp_13*(tmp_186 - 1.0/4.0) + tmp_4*(tmp_185 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_200 - 1.0/4.0) + tmp_13*(tmp_202 - 1.0/4.0) + tmp_4*(tmp_201 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_216 - 1.0/4.0) + tmp_13*(tmp_218 - 1.0/4.0) + tmp_4*(tmp_217 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_232 - 1.0/4.0) + tmp_13*(tmp_234 - 1.0/4.0) + tmp_4*(tmp_233 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_248 - 1.0/4.0) + tmp_13*(tmp_250 - 1.0/4.0) + tmp_4*(tmp_249 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_264 - 1.0/4.0) + tmp_13*(tmp_266 - 1.0/4.0) + tmp_4*(tmp_265 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_280 - 1.0/4.0) + tmp_13*(tmp_282 - 1.0/4.0) + tmp_4*(tmp_281 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_296 - 1.0/4.0) + tmp_13*(tmp_298 - 1.0/4.0) + tmp_4*(tmp_297 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_312 - 1.0/4.0) + tmp_13*(tmp_314 - 1.0/4.0) + tmp_4*(tmp_313 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_328 - 1.0/4.0) + tmp_13*(tmp_330 - 1.0/4.0) + tmp_4*(tmp_329 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_344 - 1.0/4.0) + tmp_13*(tmp_346 - 1.0/4.0) + tmp_4*(tmp_345 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_360 - 1.0/4.0) + tmp_13*(tmp_362 - 1.0/4.0) + tmp_4*(tmp_361 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_376 - 1.0/4.0) + tmp_13*(tmp_378 - 1.0/4.0) + tmp_4*(tmp_377 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_1_0 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_2_0 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_3_0 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = 0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22;
      real_t tmp_24 = p_affine_8_2 + tmp_9;
      real_t tmp_25 = tmp_19*(tmp_23 + tmp_24);
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_28 = tmp_25*tmp_27;
      real_t tmp_29 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_30 = -p_affine_8_1;
      real_t tmp_31 = p_affine_9_1 + tmp_30;
      real_t tmp_32 = p_affine_10_1 + tmp_30;
      real_t tmp_33 = 0.031405749086161582*tmp_31 + 0.93718850182767688*tmp_32;
      real_t tmp_34 = p_affine_8_1 + tmp_2;
      real_t tmp_35 = tmp_19*(tmp_33 + tmp_34);
      real_t tmp_36 = tmp_29*tmp_35;
      real_t tmp_37 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_40 = tmp_25*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_42 = tmp_35*tmp_41;
      real_t tmp_43 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_44 = -p_affine_8_0;
      real_t tmp_45 = p_affine_9_0 + tmp_44;
      real_t tmp_46 = p_affine_10_0 + tmp_44;
      real_t tmp_47 = 0.031405749086161582*tmp_45 + 0.93718850182767688*tmp_46;
      real_t tmp_48 = p_affine_8_0 + tmp_0;
      real_t tmp_49 = tmp_19*(tmp_47 + tmp_48);
      real_t tmp_50 = tmp_43*tmp_49;
      real_t tmp_51 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_52 = tmp_49*tmp_51;
      real_t tmp_53 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_54 = tmp_49*tmp_53;
      real_t tmp_55 = -p_affine_4_0;
      real_t tmp_56 = p_affine_5_0 + tmp_55;
      real_t tmp_57 = p_affine_6_0 + tmp_55;
      real_t tmp_58 = -p_affine_4_1;
      real_t tmp_59 = p_affine_7_1 + tmp_58;
      real_t tmp_60 = tmp_57*tmp_59;
      real_t tmp_61 = p_affine_7_0 + tmp_55;
      real_t tmp_62 = p_affine_6_1 + tmp_58;
      real_t tmp_63 = tmp_61*tmp_62;
      real_t tmp_64 = tmp_60 - tmp_63;
      real_t tmp_65 = -p_affine_4_2;
      real_t tmp_66 = p_affine_7_2 + tmp_65;
      real_t tmp_67 = tmp_62*tmp_66;
      real_t tmp_68 = p_affine_5_2 + tmp_65;
      real_t tmp_69 = p_affine_5_1 + tmp_58;
      real_t tmp_70 = p_affine_6_2 + tmp_65;
      real_t tmp_71 = tmp_61*tmp_70;
      real_t tmp_72 = tmp_59*tmp_70;
      real_t tmp_73 = tmp_57*tmp_66;
      real_t tmp_74 = 1.0 / (tmp_56*tmp_67 - tmp_56*tmp_72 + tmp_60*tmp_68 - tmp_63*tmp_68 + tmp_69*tmp_71 - tmp_69*tmp_73);
      real_t tmp_75 = p_affine_8_2 + tmp_65;
      real_t tmp_76 = tmp_74*(tmp_23 + tmp_75);
      real_t tmp_77 = tmp_71 - tmp_73;
      real_t tmp_78 = p_affine_8_1 + tmp_58;
      real_t tmp_79 = tmp_74*(tmp_33 + tmp_78);
      real_t tmp_80 = tmp_67 - tmp_72;
      real_t tmp_81 = p_affine_8_0 + tmp_55;
      real_t tmp_82 = tmp_74*(tmp_47 + tmp_81);
      real_t tmp_83 = -tmp_56*tmp_59 + tmp_61*tmp_69;
      real_t tmp_84 = tmp_56*tmp_66 - tmp_61*tmp_68;
      real_t tmp_85 = tmp_59*tmp_68 - tmp_66*tmp_69;
      real_t tmp_86 = tmp_56*tmp_62 - tmp_57*tmp_69;
      real_t tmp_87 = -tmp_56*tmp_70 + tmp_57*tmp_68;
      real_t tmp_88 = -tmp_62*tmp_68 + tmp_69*tmp_70;
      real_t tmp_89 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_90 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_91 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_92 = 3.0*std::pow((std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)*std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)) + (std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)*std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)) + (std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)*std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)), 0.25);
      real_t tmp_93 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_64*tmp_76 + tmp_77*tmp_79 + tmp_80*tmp_82 - 1.0/4.0) + tmp_57*(tmp_76*tmp_83 + tmp_79*tmp_84 + tmp_82*tmp_85 - 1.0/4.0) + tmp_61*(tmp_76*tmp_86 + tmp_79*tmp_87 + tmp_82*tmp_88 - 1.0/4.0));
      real_t tmp_94 = 0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22;
      real_t tmp_95 = tmp_19*(tmp_24 + tmp_94);
      real_t tmp_96 = tmp_8*tmp_95;
      real_t tmp_97 = tmp_27*tmp_95;
      real_t tmp_98 = 0.19601935860219369*tmp_31 + 0.60796128279561268*tmp_32;
      real_t tmp_99 = tmp_19*(tmp_34 + tmp_98);
      real_t tmp_100 = tmp_29*tmp_99;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = tmp_39*tmp_95;
      real_t tmp_103 = tmp_41*tmp_99;
      real_t tmp_104 = 0.19601935860219369*tmp_45 + 0.60796128279561268*tmp_46;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_48);
      real_t tmp_106 = tmp_105*tmp_43;
      real_t tmp_107 = tmp_105*tmp_51;
      real_t tmp_108 = tmp_105*tmp_53;
      real_t tmp_109 = tmp_74*(tmp_75 + tmp_94);
      real_t tmp_110 = tmp_74*(tmp_78 + tmp_98);
      real_t tmp_111 = tmp_74*(tmp_104 + tmp_81);
      real_t tmp_112 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_109*tmp_64 + tmp_110*tmp_77 + tmp_111*tmp_80 - 1.0/4.0) + tmp_57*(tmp_109*tmp_83 + tmp_110*tmp_84 + tmp_111*tmp_85 - 1.0/4.0) + tmp_61*(tmp_109*tmp_86 + tmp_110*tmp_87 + tmp_111*tmp_88 - 1.0/4.0));
      real_t tmp_113 = 0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_114 = tmp_19*(tmp_113 + tmp_24);
      real_t tmp_115 = tmp_114*tmp_8;
      real_t tmp_116 = tmp_114*tmp_27;
      real_t tmp_117 = 0.37605877282253791*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_118 = tmp_19*(tmp_117 + tmp_34);
      real_t tmp_119 = tmp_118*tmp_29;
      real_t tmp_120 = tmp_118*tmp_37;
      real_t tmp_121 = tmp_114*tmp_39;
      real_t tmp_122 = tmp_118*tmp_41;
      real_t tmp_123 = 0.37605877282253791*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_124 = tmp_19*(tmp_123 + tmp_48);
      real_t tmp_125 = tmp_124*tmp_43;
      real_t tmp_126 = tmp_124*tmp_51;
      real_t tmp_127 = tmp_124*tmp_53;
      real_t tmp_128 = tmp_74*(tmp_113 + tmp_75);
      real_t tmp_129 = tmp_74*(tmp_117 + tmp_78);
      real_t tmp_130 = tmp_74*(tmp_123 + tmp_81);
      real_t tmp_131 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_128*tmp_64 + tmp_129*tmp_77 + tmp_130*tmp_80 - 1.0/4.0) + tmp_57*(tmp_128*tmp_83 + tmp_129*tmp_84 + tmp_130*tmp_85 - 1.0/4.0) + tmp_61*(tmp_128*tmp_86 + tmp_129*tmp_87 + tmp_130*tmp_88 - 1.0/4.0));
      real_t tmp_132 = 0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_133 = tmp_19*(tmp_132 + tmp_24);
      real_t tmp_134 = tmp_133*tmp_8;
      real_t tmp_135 = tmp_133*tmp_27;
      real_t tmp_136 = 0.78764240869137092*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_137 = tmp_19*(tmp_136 + tmp_34);
      real_t tmp_138 = tmp_137*tmp_29;
      real_t tmp_139 = tmp_137*tmp_37;
      real_t tmp_140 = tmp_133*tmp_39;
      real_t tmp_141 = tmp_137*tmp_41;
      real_t tmp_142 = 0.78764240869137092*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_143 = tmp_19*(tmp_142 + tmp_48);
      real_t tmp_144 = tmp_143*tmp_43;
      real_t tmp_145 = tmp_143*tmp_51;
      real_t tmp_146 = tmp_143*tmp_53;
      real_t tmp_147 = tmp_74*(tmp_132 + tmp_75);
      real_t tmp_148 = tmp_74*(tmp_136 + tmp_78);
      real_t tmp_149 = tmp_74*(tmp_142 + tmp_81);
      real_t tmp_150 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_147*tmp_64 + tmp_148*tmp_77 + tmp_149*tmp_80 - 1.0/4.0) + tmp_57*(tmp_147*tmp_83 + tmp_148*tmp_84 + tmp_149*tmp_85 - 1.0/4.0) + tmp_61*(tmp_147*tmp_86 + tmp_148*tmp_87 + tmp_149*tmp_88 - 1.0/4.0));
      real_t tmp_151 = 0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_152 = tmp_19*(tmp_151 + tmp_24);
      real_t tmp_153 = tmp_152*tmp_8;
      real_t tmp_154 = tmp_152*tmp_27;
      real_t tmp_155 = 0.58463275527740355*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_156 = tmp_19*(tmp_155 + tmp_34);
      real_t tmp_157 = tmp_156*tmp_29;
      real_t tmp_158 = tmp_156*tmp_37;
      real_t tmp_159 = tmp_152*tmp_39;
      real_t tmp_160 = tmp_156*tmp_41;
      real_t tmp_161 = 0.58463275527740355*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_162 = tmp_19*(tmp_161 + tmp_48);
      real_t tmp_163 = tmp_162*tmp_43;
      real_t tmp_164 = tmp_162*tmp_51;
      real_t tmp_165 = tmp_162*tmp_53;
      real_t tmp_166 = tmp_74*(tmp_151 + tmp_75);
      real_t tmp_167 = tmp_74*(tmp_155 + tmp_78);
      real_t tmp_168 = tmp_74*(tmp_161 + tmp_81);
      real_t tmp_169 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_166*tmp_64 + tmp_167*tmp_77 + tmp_168*tmp_80 - 1.0/4.0) + tmp_57*(tmp_166*tmp_83 + tmp_167*tmp_84 + tmp_168*tmp_85 - 1.0/4.0) + tmp_61*(tmp_166*tmp_86 + tmp_167*tmp_87 + tmp_168*tmp_88 - 1.0/4.0));
      real_t tmp_170 = 0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_24);
      real_t tmp_172 = tmp_171*tmp_8;
      real_t tmp_173 = tmp_171*tmp_27;
      real_t tmp_174 = 0.041227165399737475*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_34);
      real_t tmp_176 = tmp_175*tmp_29;
      real_t tmp_177 = tmp_175*tmp_37;
      real_t tmp_178 = tmp_171*tmp_39;
      real_t tmp_179 = tmp_175*tmp_41;
      real_t tmp_180 = 0.041227165399737475*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_181 = tmp_19*(tmp_180 + tmp_48);
      real_t tmp_182 = tmp_181*tmp_43;
      real_t tmp_183 = tmp_181*tmp_51;
      real_t tmp_184 = tmp_181*tmp_53;
      real_t tmp_185 = tmp_74*(tmp_170 + tmp_75);
      real_t tmp_186 = tmp_74*(tmp_174 + tmp_78);
      real_t tmp_187 = tmp_74*(tmp_180 + tmp_81);
      real_t tmp_188 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_185*tmp_64 + tmp_186*tmp_77 + tmp_187*tmp_80 - 1.0/4.0) + tmp_57*(tmp_185*tmp_83 + tmp_186*tmp_84 + tmp_187*tmp_85 - 1.0/4.0) + tmp_61*(tmp_185*tmp_86 + tmp_186*tmp_87 + tmp_187*tmp_88 - 1.0/4.0));
      real_t tmp_189 = 0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_190 = tmp_19*(tmp_189 + tmp_24);
      real_t tmp_191 = tmp_190*tmp_8;
      real_t tmp_192 = tmp_190*tmp_27;
      real_t tmp_193 = 0.039308471900058539*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_194 = tmp_19*(tmp_193 + tmp_34);
      real_t tmp_195 = tmp_194*tmp_29;
      real_t tmp_196 = tmp_194*tmp_37;
      real_t tmp_197 = tmp_190*tmp_39;
      real_t tmp_198 = tmp_194*tmp_41;
      real_t tmp_199 = 0.039308471900058539*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_200 = tmp_19*(tmp_199 + tmp_48);
      real_t tmp_201 = tmp_200*tmp_43;
      real_t tmp_202 = tmp_200*tmp_51;
      real_t tmp_203 = tmp_200*tmp_53;
      real_t tmp_204 = tmp_74*(tmp_189 + tmp_75);
      real_t tmp_205 = tmp_74*(tmp_193 + tmp_78);
      real_t tmp_206 = tmp_74*(tmp_199 + tmp_81);
      real_t tmp_207 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_204*tmp_64 + tmp_205*tmp_77 + tmp_206*tmp_80 - 1.0/4.0) + tmp_57*(tmp_204*tmp_83 + tmp_205*tmp_84 + tmp_206*tmp_85 - 1.0/4.0) + tmp_61*(tmp_204*tmp_86 + tmp_205*tmp_87 + tmp_206*tmp_88 - 1.0/4.0));
      real_t tmp_208 = 0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_209 = tmp_19*(tmp_208 + tmp_24);
      real_t tmp_210 = tmp_209*tmp_8;
      real_t tmp_211 = tmp_209*tmp_27;
      real_t tmp_212 = 0.78764240869137092*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_213 = tmp_19*(tmp_212 + tmp_34);
      real_t tmp_214 = tmp_213*tmp_29;
      real_t tmp_215 = tmp_213*tmp_37;
      real_t tmp_216 = tmp_209*tmp_39;
      real_t tmp_217 = tmp_213*tmp_41;
      real_t tmp_218 = 0.78764240869137092*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_48);
      real_t tmp_220 = tmp_219*tmp_43;
      real_t tmp_221 = tmp_219*tmp_51;
      real_t tmp_222 = tmp_219*tmp_53;
      real_t tmp_223 = tmp_74*(tmp_208 + tmp_75);
      real_t tmp_224 = tmp_74*(tmp_212 + tmp_78);
      real_t tmp_225 = tmp_74*(tmp_218 + tmp_81);
      real_t tmp_226 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_223*tmp_64 + tmp_224*tmp_77 + tmp_225*tmp_80 - 1.0/4.0) + tmp_57*(tmp_223*tmp_83 + tmp_224*tmp_84 + tmp_225*tmp_85 - 1.0/4.0) + tmp_61*(tmp_223*tmp_86 + tmp_224*tmp_87 + tmp_225*tmp_88 - 1.0/4.0));
      real_t tmp_227 = 0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_228 = tmp_19*(tmp_227 + tmp_24);
      real_t tmp_229 = tmp_228*tmp_8;
      real_t tmp_230 = tmp_228*tmp_27;
      real_t tmp_231 = 0.58463275527740355*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_232 = tmp_19*(tmp_231 + tmp_34);
      real_t tmp_233 = tmp_232*tmp_29;
      real_t tmp_234 = tmp_232*tmp_37;
      real_t tmp_235 = tmp_228*tmp_39;
      real_t tmp_236 = tmp_232*tmp_41;
      real_t tmp_237 = 0.58463275527740355*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_238 = tmp_19*(tmp_237 + tmp_48);
      real_t tmp_239 = tmp_238*tmp_43;
      real_t tmp_240 = tmp_238*tmp_51;
      real_t tmp_241 = tmp_238*tmp_53;
      real_t tmp_242 = tmp_74*(tmp_227 + tmp_75);
      real_t tmp_243 = tmp_74*(tmp_231 + tmp_78);
      real_t tmp_244 = tmp_74*(tmp_237 + tmp_81);
      real_t tmp_245 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_242*tmp_64 + tmp_243*tmp_77 + tmp_244*tmp_80 - 1.0/4.0) + tmp_57*(tmp_242*tmp_83 + tmp_243*tmp_84 + tmp_244*tmp_85 - 1.0/4.0) + tmp_61*(tmp_242*tmp_86 + tmp_243*tmp_87 + tmp_244*tmp_88 - 1.0/4.0));
      real_t tmp_246 = 0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_247 = tmp_19*(tmp_24 + tmp_246);
      real_t tmp_248 = tmp_247*tmp_8;
      real_t tmp_249 = tmp_247*tmp_27;
      real_t tmp_250 = 0.1711304259088916*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_251 = tmp_19*(tmp_250 + tmp_34);
      real_t tmp_252 = tmp_251*tmp_29;
      real_t tmp_253 = tmp_251*tmp_37;
      real_t tmp_254 = tmp_247*tmp_39;
      real_t tmp_255 = tmp_251*tmp_41;
      real_t tmp_256 = 0.1711304259088916*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_257 = tmp_19*(tmp_256 + tmp_48);
      real_t tmp_258 = tmp_257*tmp_43;
      real_t tmp_259 = tmp_257*tmp_51;
      real_t tmp_260 = tmp_257*tmp_53;
      real_t tmp_261 = tmp_74*(tmp_246 + tmp_75);
      real_t tmp_262 = tmp_74*(tmp_250 + tmp_78);
      real_t tmp_263 = tmp_74*(tmp_256 + tmp_81);
      real_t tmp_264 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_261*tmp_64 + tmp_262*tmp_77 + tmp_263*tmp_80 - 1.0/4.0) + tmp_57*(tmp_261*tmp_83 + tmp_262*tmp_84 + tmp_263*tmp_85 - 1.0/4.0) + tmp_61*(tmp_261*tmp_86 + tmp_262*tmp_87 + tmp_263*tmp_88 - 1.0/4.0));
      real_t tmp_265 = 0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_266 = tmp_19*(tmp_24 + tmp_265);
      real_t tmp_267 = tmp_266*tmp_8;
      real_t tmp_268 = tmp_266*tmp_27;
      real_t tmp_269 = 0.37605877282253791*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_270 = tmp_19*(tmp_269 + tmp_34);
      real_t tmp_271 = tmp_270*tmp_29;
      real_t tmp_272 = tmp_270*tmp_37;
      real_t tmp_273 = tmp_266*tmp_39;
      real_t tmp_274 = tmp_270*tmp_41;
      real_t tmp_275 = 0.37605877282253791*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_276 = tmp_19*(tmp_275 + tmp_48);
      real_t tmp_277 = tmp_276*tmp_43;
      real_t tmp_278 = tmp_276*tmp_51;
      real_t tmp_279 = tmp_276*tmp_53;
      real_t tmp_280 = tmp_74*(tmp_265 + tmp_75);
      real_t tmp_281 = tmp_74*(tmp_269 + tmp_78);
      real_t tmp_282 = tmp_74*(tmp_275 + tmp_81);
      real_t tmp_283 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_280*tmp_64 + tmp_281*tmp_77 + tmp_282*tmp_80 - 1.0/4.0) + tmp_57*(tmp_280*tmp_83 + tmp_281*tmp_84 + tmp_282*tmp_85 - 1.0/4.0) + tmp_61*(tmp_280*tmp_86 + tmp_281*tmp_87 + tmp_282*tmp_88 - 1.0/4.0));
      real_t tmp_284 = 0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_285 = tmp_19*(tmp_24 + tmp_284);
      real_t tmp_286 = tmp_285*tmp_8;
      real_t tmp_287 = tmp_27*tmp_285;
      real_t tmp_288 = 0.041227165399737475*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_289 = tmp_19*(tmp_288 + tmp_34);
      real_t tmp_290 = tmp_289*tmp_29;
      real_t tmp_291 = tmp_289*tmp_37;
      real_t tmp_292 = tmp_285*tmp_39;
      real_t tmp_293 = tmp_289*tmp_41;
      real_t tmp_294 = 0.041227165399737475*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_295 = tmp_19*(tmp_294 + tmp_48);
      real_t tmp_296 = tmp_295*tmp_43;
      real_t tmp_297 = tmp_295*tmp_51;
      real_t tmp_298 = tmp_295*tmp_53;
      real_t tmp_299 = tmp_74*(tmp_284 + tmp_75);
      real_t tmp_300 = tmp_74*(tmp_288 + tmp_78);
      real_t tmp_301 = tmp_74*(tmp_294 + tmp_81);
      real_t tmp_302 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_299*tmp_64 + tmp_300*tmp_77 + tmp_301*tmp_80 - 1.0/4.0) + tmp_57*(tmp_299*tmp_83 + tmp_300*tmp_84 + tmp_301*tmp_85 - 1.0/4.0) + tmp_61*(tmp_299*tmp_86 + tmp_300*tmp_87 + tmp_301*tmp_88 - 1.0/4.0));
      real_t tmp_303 = 0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22;
      real_t tmp_304 = tmp_19*(tmp_24 + tmp_303);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_27*tmp_304;
      real_t tmp_307 = 0.40446199974765351*tmp_31 + 0.19107600050469298*tmp_32;
      real_t tmp_308 = tmp_19*(tmp_307 + tmp_34);
      real_t tmp_309 = tmp_29*tmp_308;
      real_t tmp_310 = tmp_308*tmp_37;
      real_t tmp_311 = tmp_304*tmp_39;
      real_t tmp_312 = tmp_308*tmp_41;
      real_t tmp_313 = 0.40446199974765351*tmp_45 + 0.19107600050469298*tmp_46;
      real_t tmp_314 = tmp_19*(tmp_313 + tmp_48);
      real_t tmp_315 = tmp_314*tmp_43;
      real_t tmp_316 = tmp_314*tmp_51;
      real_t tmp_317 = tmp_314*tmp_53;
      real_t tmp_318 = tmp_74*(tmp_303 + tmp_75);
      real_t tmp_319 = tmp_74*(tmp_307 + tmp_78);
      real_t tmp_320 = tmp_74*(tmp_313 + tmp_81);
      real_t tmp_321 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_318*tmp_64 + tmp_319*tmp_77 + tmp_320*tmp_80 - 1.0/4.0) + tmp_57*(tmp_318*tmp_83 + tmp_319*tmp_84 + tmp_320*tmp_85 - 1.0/4.0) + tmp_61*(tmp_318*tmp_86 + tmp_319*tmp_87 + tmp_320*tmp_88 - 1.0/4.0));
      real_t tmp_322 = 0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_323 = tmp_19*(tmp_24 + tmp_322);
      real_t tmp_324 = tmp_323*tmp_8;
      real_t tmp_325 = tmp_27*tmp_323;
      real_t tmp_326 = 0.039308471900058539*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_327 = tmp_19*(tmp_326 + tmp_34);
      real_t tmp_328 = tmp_29*tmp_327;
      real_t tmp_329 = tmp_327*tmp_37;
      real_t tmp_330 = tmp_323*tmp_39;
      real_t tmp_331 = tmp_327*tmp_41;
      real_t tmp_332 = 0.039308471900058539*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_333 = tmp_19*(tmp_332 + tmp_48);
      real_t tmp_334 = tmp_333*tmp_43;
      real_t tmp_335 = tmp_333*tmp_51;
      real_t tmp_336 = tmp_333*tmp_53;
      real_t tmp_337 = tmp_74*(tmp_322 + tmp_75);
      real_t tmp_338 = tmp_74*(tmp_326 + tmp_78);
      real_t tmp_339 = tmp_74*(tmp_332 + tmp_81);
      real_t tmp_340 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_337*tmp_64 + tmp_338*tmp_77 + tmp_339*tmp_80 - 1.0/4.0) + tmp_57*(tmp_337*tmp_83 + tmp_338*tmp_84 + tmp_339*tmp_85 - 1.0/4.0) + tmp_61*(tmp_337*tmp_86 + tmp_338*tmp_87 + tmp_339*tmp_88 - 1.0/4.0));
      real_t tmp_341 = 0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_342 = tmp_19*(tmp_24 + tmp_341);
      real_t tmp_343 = tmp_342*tmp_8;
      real_t tmp_344 = tmp_27*tmp_342;
      real_t tmp_345 = 0.93718850182767688*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_346 = tmp_19*(tmp_34 + tmp_345);
      real_t tmp_347 = tmp_29*tmp_346;
      real_t tmp_348 = tmp_346*tmp_37;
      real_t tmp_349 = tmp_342*tmp_39;
      real_t tmp_350 = tmp_346*tmp_41;
      real_t tmp_351 = 0.93718850182767688*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_352 = tmp_19*(tmp_351 + tmp_48);
      real_t tmp_353 = tmp_352*tmp_43;
      real_t tmp_354 = tmp_352*tmp_51;
      real_t tmp_355 = tmp_352*tmp_53;
      real_t tmp_356 = tmp_74*(tmp_341 + tmp_75);
      real_t tmp_357 = tmp_74*(tmp_345 + tmp_78);
      real_t tmp_358 = tmp_74*(tmp_351 + tmp_81);
      real_t tmp_359 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_356*tmp_64 + tmp_357*tmp_77 + tmp_358*tmp_80 - 1.0/4.0) + tmp_57*(tmp_356*tmp_83 + tmp_357*tmp_84 + tmp_358*tmp_85 - 1.0/4.0) + tmp_61*(tmp_356*tmp_86 + tmp_357*tmp_87 + tmp_358*tmp_88 - 1.0/4.0));
      real_t tmp_360 = 0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_361 = tmp_19*(tmp_24 + tmp_360);
      real_t tmp_362 = tmp_361*tmp_8;
      real_t tmp_363 = tmp_27*tmp_361;
      real_t tmp_364 = 0.60796128279561268*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_365 = tmp_19*(tmp_34 + tmp_364);
      real_t tmp_366 = tmp_29*tmp_365;
      real_t tmp_367 = tmp_365*tmp_37;
      real_t tmp_368 = tmp_361*tmp_39;
      real_t tmp_369 = tmp_365*tmp_41;
      real_t tmp_370 = 0.60796128279561268*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_371 = tmp_19*(tmp_370 + tmp_48);
      real_t tmp_372 = tmp_371*tmp_43;
      real_t tmp_373 = tmp_371*tmp_51;
      real_t tmp_374 = tmp_371*tmp_53;
      real_t tmp_375 = tmp_74*(tmp_360 + tmp_75);
      real_t tmp_376 = tmp_74*(tmp_364 + tmp_78);
      real_t tmp_377 = tmp_74*(tmp_370 + tmp_81);
      real_t tmp_378 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_375*tmp_64 + tmp_376*tmp_77 + tmp_377*tmp_80 - 1.0/4.0) + tmp_57*(tmp_375*tmp_83 + tmp_376*tmp_84 + tmp_377*tmp_85 - 1.0/4.0) + tmp_61*(tmp_375*tmp_86 + tmp_376*tmp_87 + tmp_377*tmp_88 - 1.0/4.0));
      real_t tmp_379 = 0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_380 = tmp_19*(tmp_24 + tmp_379);
      real_t tmp_381 = tmp_380*tmp_8;
      real_t tmp_382 = tmp_27*tmp_380;
      real_t tmp_383 = 0.19107600050469298*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_384 = tmp_19*(tmp_34 + tmp_383);
      real_t tmp_385 = tmp_29*tmp_384;
      real_t tmp_386 = tmp_37*tmp_384;
      real_t tmp_387 = tmp_380*tmp_39;
      real_t tmp_388 = tmp_384*tmp_41;
      real_t tmp_389 = 0.19107600050469298*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_390 = tmp_19*(tmp_389 + tmp_48);
      real_t tmp_391 = tmp_390*tmp_43;
      real_t tmp_392 = tmp_390*tmp_51;
      real_t tmp_393 = tmp_390*tmp_53;
      real_t tmp_394 = tmp_74*(tmp_379 + tmp_75);
      real_t tmp_395 = tmp_74*(tmp_383 + tmp_78);
      real_t tmp_396 = tmp_74*(tmp_389 + tmp_81);
      real_t tmp_397 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_394*tmp_64 + tmp_395*tmp_77 + tmp_396*tmp_80 - 1.0/4.0) + tmp_57*(tmp_394*tmp_83 + tmp_395*tmp_84 + tmp_396*tmp_85 - 1.0/4.0) + tmp_61*(tmp_394*tmp_86 + tmp_395*tmp_87 + tmp_396*tmp_88 - 1.0/4.0));
      real_t tmp_398 = 0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_399 = tmp_19*(tmp_24 + tmp_398);
      real_t tmp_400 = tmp_399*tmp_8;
      real_t tmp_401 = tmp_27*tmp_399;
      real_t tmp_402 = 0.031405749086161582*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_403 = tmp_19*(tmp_34 + tmp_402);
      real_t tmp_404 = tmp_29*tmp_403;
      real_t tmp_405 = tmp_37*tmp_403;
      real_t tmp_406 = tmp_39*tmp_399;
      real_t tmp_407 = tmp_403*tmp_41;
      real_t tmp_408 = 0.031405749086161582*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_409 = tmp_19*(tmp_408 + tmp_48);
      real_t tmp_410 = tmp_409*tmp_43;
      real_t tmp_411 = tmp_409*tmp_51;
      real_t tmp_412 = tmp_409*tmp_53;
      real_t tmp_413 = tmp_74*(tmp_398 + tmp_75);
      real_t tmp_414 = tmp_74*(tmp_402 + tmp_78);
      real_t tmp_415 = tmp_74*(tmp_408 + tmp_81);
      real_t tmp_416 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_413*tmp_64 + tmp_414*tmp_77 + tmp_415*tmp_80 - 1.0/4.0) + tmp_57*(tmp_413*tmp_83 + tmp_414*tmp_84 + tmp_415*tmp_85 - 1.0/4.0) + tmp_61*(tmp_413*tmp_86 + tmp_414*tmp_87 + tmp_415*tmp_88 - 1.0/4.0));
      real_t tmp_417 = 0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_418 = tmp_19*(tmp_24 + tmp_417);
      real_t tmp_419 = tmp_418*tmp_8;
      real_t tmp_420 = tmp_27*tmp_418;
      real_t tmp_421 = 0.19601935860219369*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_422 = tmp_19*(tmp_34 + tmp_421);
      real_t tmp_423 = tmp_29*tmp_422;
      real_t tmp_424 = tmp_37*tmp_422;
      real_t tmp_425 = tmp_39*tmp_418;
      real_t tmp_426 = tmp_41*tmp_422;
      real_t tmp_427 = 0.19601935860219369*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_428 = tmp_19*(tmp_427 + tmp_48);
      real_t tmp_429 = tmp_428*tmp_43;
      real_t tmp_430 = tmp_428*tmp_51;
      real_t tmp_431 = tmp_428*tmp_53;
      real_t tmp_432 = tmp_74*(tmp_417 + tmp_75);
      real_t tmp_433 = tmp_74*(tmp_421 + tmp_78);
      real_t tmp_434 = tmp_74*(tmp_427 + tmp_81);
      real_t tmp_435 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_432*tmp_64 + tmp_433*tmp_77 + tmp_434*tmp_80 - 1.0/4.0) + tmp_57*(tmp_432*tmp_83 + tmp_433*tmp_84 + tmp_434*tmp_85 - 1.0/4.0) + tmp_61*(tmp_432*tmp_86 + tmp_433*tmp_87 + tmp_434*tmp_88 - 1.0/4.0));
      real_t tmp_436 = 0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_437 = tmp_19*(tmp_24 + tmp_436);
      real_t tmp_438 = tmp_437*tmp_8;
      real_t tmp_439 = tmp_27*tmp_437;
      real_t tmp_440 = 0.40446199974765351*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_441 = tmp_19*(tmp_34 + tmp_440);
      real_t tmp_442 = tmp_29*tmp_441;
      real_t tmp_443 = tmp_37*tmp_441;
      real_t tmp_444 = tmp_39*tmp_437;
      real_t tmp_445 = tmp_41*tmp_441;
      real_t tmp_446 = 0.40446199974765351*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_447 = tmp_19*(tmp_446 + tmp_48);
      real_t tmp_448 = tmp_43*tmp_447;
      real_t tmp_449 = tmp_447*tmp_51;
      real_t tmp_450 = tmp_447*tmp_53;
      real_t tmp_451 = tmp_74*(tmp_436 + tmp_75);
      real_t tmp_452 = tmp_74*(tmp_440 + tmp_78);
      real_t tmp_453 = tmp_74*(tmp_446 + tmp_81);
      real_t tmp_454 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_451*tmp_64 + tmp_452*tmp_77 + tmp_453*tmp_80 - 1.0/4.0) + tmp_57*(tmp_451*tmp_83 + tmp_452*tmp_84 + tmp_453*tmp_85 - 1.0/4.0) + tmp_61*(tmp_451*tmp_86 + tmp_452*tmp_87 + tmp_453*tmp_88 - 1.0/4.0));
      real_t tmp_455 = 0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_456 = tmp_19*(tmp_24 + tmp_455);
      real_t tmp_457 = tmp_456*tmp_8;
      real_t tmp_458 = tmp_27*tmp_456;
      real_t tmp_459 = 0.1711304259088916*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_460 = tmp_19*(tmp_34 + tmp_459);
      real_t tmp_461 = tmp_29*tmp_460;
      real_t tmp_462 = tmp_37*tmp_460;
      real_t tmp_463 = tmp_39*tmp_456;
      real_t tmp_464 = tmp_41*tmp_460;
      real_t tmp_465 = 0.1711304259088916*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_466 = tmp_19*(tmp_465 + tmp_48);
      real_t tmp_467 = tmp_43*tmp_466;
      real_t tmp_468 = tmp_466*tmp_51;
      real_t tmp_469 = tmp_466*tmp_53;
      real_t tmp_470 = tmp_74*(tmp_455 + tmp_75);
      real_t tmp_471 = tmp_74*(tmp_459 + tmp_78);
      real_t tmp_472 = tmp_74*(tmp_465 + tmp_81);
      real_t tmp_473 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_470*tmp_64 + tmp_471*tmp_77 + tmp_472*tmp_80 - 1.0/4.0) + tmp_57*(tmp_470*tmp_83 + tmp_471*tmp_84 + tmp_472*tmp_85 - 1.0/4.0) + tmp_61*(tmp_470*tmp_86 + tmp_471*tmp_87 + tmp_472*tmp_88 - 1.0/4.0));
      real_t a_0_0 = -tmp_112*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_106 - tmp_107 - tmp_108 - tmp_96 - tmp_97 + 1) - tmp_131*(-tmp_115 - tmp_116 - tmp_119 - tmp_120 - tmp_121 - tmp_122 - tmp_125 - tmp_126 - tmp_127 + 1) - tmp_150*(-tmp_134 - tmp_135 - tmp_138 - tmp_139 - tmp_140 - tmp_141 - tmp_144 - tmp_145 - tmp_146 + 1) - tmp_169*(-tmp_153 - tmp_154 - tmp_157 - tmp_158 - tmp_159 - tmp_160 - tmp_163 - tmp_164 - tmp_165 + 1) - tmp_188*(-tmp_172 - tmp_173 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_182 - tmp_183 - tmp_184 + 1) - tmp_207*(-tmp_191 - tmp_192 - tmp_195 - tmp_196 - tmp_197 - tmp_198 - tmp_201 - tmp_202 - tmp_203 + 1) - tmp_226*(-tmp_210 - tmp_211 - tmp_214 - tmp_215 - tmp_216 - tmp_217 - tmp_220 - tmp_221 - tmp_222 + 1) - tmp_245*(-tmp_229 - tmp_230 - tmp_233 - tmp_234 - tmp_235 - tmp_236 - tmp_239 - tmp_240 - tmp_241 + 1) - tmp_264*(-tmp_248 - tmp_249 - tmp_252 - tmp_253 - tmp_254 - tmp_255 - tmp_258 - tmp_259 - tmp_260 + 1) - tmp_283*(-tmp_267 - tmp_268 - tmp_271 - tmp_272 - tmp_273 - tmp_274 - tmp_277 - tmp_278 - tmp_279 + 1) - tmp_302*(-tmp_286 - tmp_287 - tmp_290 - tmp_291 - tmp_292 - tmp_293 - tmp_296 - tmp_297 - tmp_298 + 1) - tmp_321*(-tmp_305 - tmp_306 - tmp_309 - tmp_310 - tmp_311 - tmp_312 - tmp_315 - tmp_316 - tmp_317 + 1) - tmp_340*(-tmp_324 - tmp_325 - tmp_328 - tmp_329 - tmp_330 - tmp_331 - tmp_334 - tmp_335 - tmp_336 + 1) - tmp_359*(-tmp_343 - tmp_344 - tmp_347 - tmp_348 - tmp_349 - tmp_350 - tmp_353 - tmp_354 - tmp_355 + 1) - tmp_378*(-tmp_362 - tmp_363 - tmp_366 - tmp_367 - tmp_368 - tmp_369 - tmp_372 - tmp_373 - tmp_374 + 1) - tmp_397*(-tmp_381 - tmp_382 - tmp_385 - tmp_386 - tmp_387 - tmp_388 - tmp_391 - tmp_392 - tmp_393 + 1) - tmp_416*(-tmp_400 - tmp_401 - tmp_404 - tmp_405 - tmp_406 - tmp_407 - tmp_410 - tmp_411 - tmp_412 + 1) - tmp_435*(-tmp_419 - tmp_420 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_429 - tmp_430 - tmp_431 + 1) - tmp_454*(-tmp_438 - tmp_439 - tmp_442 - tmp_443 - tmp_444 - tmp_445 - tmp_448 - tmp_449 - tmp_450 + 1) - tmp_473*(-tmp_457 - tmp_458 - tmp_461 - tmp_462 - tmp_463 - tmp_464 - tmp_467 - tmp_468 - tmp_469 + 1) - tmp_93*(-tmp_26 - tmp_28 - tmp_36 - tmp_38 - tmp_40 - tmp_42 - tmp_50 - tmp_52 - tmp_54 + 1);
      real_t a_1_0 = -tmp_112*(tmp_102 + tmp_103 + tmp_108) - tmp_131*(tmp_121 + tmp_122 + tmp_127) - tmp_150*(tmp_140 + tmp_141 + tmp_146) - tmp_169*(tmp_159 + tmp_160 + tmp_165) - tmp_188*(tmp_178 + tmp_179 + tmp_184) - tmp_207*(tmp_197 + tmp_198 + tmp_203) - tmp_226*(tmp_216 + tmp_217 + tmp_222) - tmp_245*(tmp_235 + tmp_236 + tmp_241) - tmp_264*(tmp_254 + tmp_255 + tmp_260) - tmp_283*(tmp_273 + tmp_274 + tmp_279) - tmp_302*(tmp_292 + tmp_293 + tmp_298) - tmp_321*(tmp_311 + tmp_312 + tmp_317) - tmp_340*(tmp_330 + tmp_331 + tmp_336) - tmp_359*(tmp_349 + tmp_350 + tmp_355) - tmp_378*(tmp_368 + tmp_369 + tmp_374) - tmp_397*(tmp_387 + tmp_388 + tmp_393) - tmp_416*(tmp_406 + tmp_407 + tmp_412) - tmp_435*(tmp_425 + tmp_426 + tmp_431) - tmp_454*(tmp_444 + tmp_445 + tmp_450) - tmp_473*(tmp_463 + tmp_464 + tmp_469) - tmp_93*(tmp_40 + tmp_42 + tmp_54);
      real_t a_2_0 = -tmp_112*(tmp_101 + tmp_107 + tmp_97) - tmp_131*(tmp_116 + tmp_120 + tmp_126) - tmp_150*(tmp_135 + tmp_139 + tmp_145) - tmp_169*(tmp_154 + tmp_158 + tmp_164) - tmp_188*(tmp_173 + tmp_177 + tmp_183) - tmp_207*(tmp_192 + tmp_196 + tmp_202) - tmp_226*(tmp_211 + tmp_215 + tmp_221) - tmp_245*(tmp_230 + tmp_234 + tmp_240) - tmp_264*(tmp_249 + tmp_253 + tmp_259) - tmp_283*(tmp_268 + tmp_272 + tmp_278) - tmp_302*(tmp_287 + tmp_291 + tmp_297) - tmp_321*(tmp_306 + tmp_310 + tmp_316) - tmp_340*(tmp_325 + tmp_329 + tmp_335) - tmp_359*(tmp_344 + tmp_348 + tmp_354) - tmp_378*(tmp_363 + tmp_367 + tmp_373) - tmp_397*(tmp_382 + tmp_386 + tmp_392) - tmp_416*(tmp_401 + tmp_405 + tmp_411) - tmp_435*(tmp_420 + tmp_424 + tmp_430) - tmp_454*(tmp_439 + tmp_443 + tmp_449) - tmp_473*(tmp_458 + tmp_462 + tmp_468) - tmp_93*(tmp_28 + tmp_38 + tmp_52);
      real_t a_3_0 = -tmp_112*(tmp_100 + tmp_106 + tmp_96) - tmp_131*(tmp_115 + tmp_119 + tmp_125) - tmp_150*(tmp_134 + tmp_138 + tmp_144) - tmp_169*(tmp_153 + tmp_157 + tmp_163) - tmp_188*(tmp_172 + tmp_176 + tmp_182) - tmp_207*(tmp_191 + tmp_195 + tmp_201) - tmp_226*(tmp_210 + tmp_214 + tmp_220) - tmp_245*(tmp_229 + tmp_233 + tmp_239) - tmp_264*(tmp_248 + tmp_252 + tmp_258) - tmp_283*(tmp_267 + tmp_271 + tmp_277) - tmp_302*(tmp_286 + tmp_290 + tmp_296) - tmp_321*(tmp_305 + tmp_309 + tmp_315) - tmp_340*(tmp_324 + tmp_328 + tmp_334) - tmp_359*(tmp_343 + tmp_347 + tmp_353) - tmp_378*(tmp_362 + tmp_366 + tmp_372) - tmp_397*(tmp_381 + tmp_385 + tmp_391) - tmp_416*(tmp_400 + tmp_404 + tmp_410) - tmp_435*(tmp_419 + tmp_423 + tmp_429) - tmp_454*(tmp_438 + tmp_442 + tmp_448) - tmp_473*(tmp_457 + tmp_461 + tmp_467) - tmp_93*(tmp_26 + tmp_36 + tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 3.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_52 - 1.0/4.0) + tmp_13*(tmp_54 - 1.0/4.0) + tmp_4*(tmp_53 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_72 - 1.0/4.0) + tmp_13*(tmp_74 - 1.0/4.0) + tmp_4*(tmp_73 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_88 - 1.0/4.0) + tmp_13*(tmp_90 - 1.0/4.0) + tmp_4*(tmp_89 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_104 - 1.0/4.0) + tmp_13*(tmp_106 - 1.0/4.0) + tmp_4*(tmp_105 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_120 - 1.0/4.0) + tmp_13*(tmp_122 - 1.0/4.0) + tmp_4*(tmp_121 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_136 - 1.0/4.0) + tmp_13*(tmp_138 - 1.0/4.0) + tmp_4*(tmp_137 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_152 - 1.0/4.0) + tmp_13*(tmp_154 - 1.0/4.0) + tmp_4*(tmp_153 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_168 - 1.0/4.0) + tmp_13*(tmp_170 - 1.0/4.0) + tmp_4*(tmp_169 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_184 - 1.0/4.0) + tmp_13*(tmp_186 - 1.0/4.0) + tmp_4*(tmp_185 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_200 - 1.0/4.0) + tmp_13*(tmp_202 - 1.0/4.0) + tmp_4*(tmp_201 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_216 - 1.0/4.0) + tmp_13*(tmp_218 - 1.0/4.0) + tmp_4*(tmp_217 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_232 - 1.0/4.0) + tmp_13*(tmp_234 - 1.0/4.0) + tmp_4*(tmp_233 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_248 - 1.0/4.0) + tmp_13*(tmp_250 - 1.0/4.0) + tmp_4*(tmp_249 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_264 - 1.0/4.0) + tmp_13*(tmp_266 - 1.0/4.0) + tmp_4*(tmp_265 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_280 - 1.0/4.0) + tmp_13*(tmp_282 - 1.0/4.0) + tmp_4*(tmp_281 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_296 - 1.0/4.0) + tmp_13*(tmp_298 - 1.0/4.0) + tmp_4*(tmp_297 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_312 - 1.0/4.0) + tmp_13*(tmp_314 - 1.0/4.0) + tmp_4*(tmp_313 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_328 - 1.0/4.0) + tmp_13*(tmp_330 - 1.0/4.0) + tmp_4*(tmp_329 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_344 - 1.0/4.0) + tmp_13*(tmp_346 - 1.0/4.0) + tmp_4*(tmp_345 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_360 - 1.0/4.0) + tmp_13*(tmp_362 - 1.0/4.0) + tmp_4*(tmp_361 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_376 - 1.0/4.0) + tmp_13*(tmp_378 - 1.0/4.0) + tmp_4*(tmp_377 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_1_0 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_2_0 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_3_0 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }


};




class EGNIPGVectorLaplaceFormP1E_1 : public hyteg::dg::DGForm
{
 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_2_0 + tmp_0));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = tmp_4*tmp_6 + tmp_5*tmp_8;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_10*tmp_5 + tmp_11*tmp_3;
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

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_1_1 + tmp_2;
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_2_0 + tmp_0));
      real_t tmp_6 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7 = p_affine_6_1 + tmp_2;
      real_t tmp_8 = tmp_5*(0.046910077030668018*tmp_6 + tmp_7);
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5*(0.046910077030668018*tmp_12 + tmp_13);
      real_t tmp_15 = tmp_14*tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 0.35539032758428413*tmp_3*(tmp_19 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_18 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 0.71794300574904923*tmp_3*(tmp_28 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_27 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 0.8533333333333335*tmp_3*(tmp_37 - 1.0/3.0) + 0.8533333333333335*tmp_4*(tmp_36 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 0.71794300574904923*tmp_3*(tmp_46 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_45 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.35539032758428413*tmp_3*(tmp_55 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_54 - 1.0/3.0);
      real_t a_0_0 = tmp_20*(-tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1) + tmp_29*(-tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1) + tmp_38*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_47*(-tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1) + tmp_56*(-tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1);
      real_t a_1_0 = tmp_18*tmp_20 + tmp_27*tmp_29 + tmp_36*tmp_38 + tmp_45*tmp_47 + tmp_54*tmp_56;
      real_t a_2_0 = tmp_19*tmp_20 + tmp_28*tmp_29 + tmp_37*tmp_38 + tmp_46*tmp_47 + tmp_55*tmp_56;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + 0.046910077030668018*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -p_affine_3_1;
      real_t tmp_18 = p_affine_4_1 + tmp_17;
      real_t tmp_19 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_20 = -p_affine_3_0;
      real_t tmp_21 = p_affine_4_0 + tmp_20;
      real_t tmp_22 = p_affine_5_1 + tmp_17;
      real_t tmp_23 = 1.0 / (-tmp_18*(p_affine_5_0 + tmp_20) + tmp_21*tmp_22);
      real_t tmp_24 = tmp_23*(tmp_17 + tmp_6);
      real_t tmp_25 = tmp_23*(tmp_12 + tmp_20);
      real_t tmp_26 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_27 = 0.35539032758428413*tmp_18*(tmp_19*tmp_24 + tmp_22*tmp_25 - 1.0/3.0) + 0.35539032758428413*tmp_22*(tmp_21*tmp_24 + tmp_25*tmp_26 - 1.0/3.0);
      real_t tmp_28 = p_affine_6_1 + 0.23076534494715845*tmp_5;
      real_t tmp_29 = tmp_4*(tmp_2 + tmp_28);
      real_t tmp_30 = tmp_1*tmp_29;
      real_t tmp_31 = tmp_29*tmp_9;
      real_t tmp_32 = p_affine_6_0 + 0.23076534494715845*tmp_11;
      real_t tmp_33 = tmp_4*(tmp_0 + tmp_32);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_15*tmp_33;
      real_t tmp_36 = tmp_23*(tmp_17 + tmp_28);
      real_t tmp_37 = tmp_23*(tmp_20 + tmp_32);
      real_t tmp_38 = 0.71794300574904923*tmp_18*(tmp_19*tmp_36 + tmp_22*tmp_37 - 1.0/3.0) + 0.71794300574904923*tmp_22*(tmp_21*tmp_36 + tmp_26*tmp_37 - 1.0/3.0);
      real_t tmp_39 = p_affine_6_1 + 0.5*tmp_5;
      real_t tmp_40 = tmp_4*(tmp_2 + tmp_39);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = p_affine_6_0 + 0.5*tmp_11;
      real_t tmp_44 = tmp_4*(tmp_0 + tmp_43);
      real_t tmp_45 = tmp_3*tmp_44;
      real_t tmp_46 = tmp_15*tmp_44;
      real_t tmp_47 = tmp_23*(tmp_17 + tmp_39);
      real_t tmp_48 = tmp_23*(tmp_20 + tmp_43);
      real_t tmp_49 = 0.8533333333333335*tmp_18*(tmp_19*tmp_47 + tmp_22*tmp_48 - 1.0/3.0) + 0.8533333333333335*tmp_22*(tmp_21*tmp_47 + tmp_26*tmp_48 - 1.0/3.0);
      real_t tmp_50 = p_affine_6_1 + 0.7692346550528415*tmp_5;
      real_t tmp_51 = tmp_4*(tmp_2 + tmp_50);
      real_t tmp_52 = tmp_1*tmp_51;
      real_t tmp_53 = tmp_51*tmp_9;
      real_t tmp_54 = p_affine_6_0 + 0.7692346550528415*tmp_11;
      real_t tmp_55 = tmp_4*(tmp_0 + tmp_54);
      real_t tmp_56 = tmp_3*tmp_55;
      real_t tmp_57 = tmp_15*tmp_55;
      real_t tmp_58 = tmp_23*(tmp_17 + tmp_50);
      real_t tmp_59 = tmp_23*(tmp_20 + tmp_54);
      real_t tmp_60 = 0.71794300574904923*tmp_18*(tmp_19*tmp_58 + tmp_22*tmp_59 - 1.0/3.0) + 0.71794300574904923*tmp_22*(tmp_21*tmp_58 + tmp_26*tmp_59 - 1.0/3.0);
      real_t tmp_61 = p_affine_6_1 + 0.95308992296933193*tmp_5;
      real_t tmp_62 = tmp_4*(tmp_2 + tmp_61);
      real_t tmp_63 = tmp_1*tmp_62;
      real_t tmp_64 = tmp_62*tmp_9;
      real_t tmp_65 = p_affine_6_0 + 0.95308992296933193*tmp_11;
      real_t tmp_66 = tmp_4*(tmp_0 + tmp_65);
      real_t tmp_67 = tmp_3*tmp_66;
      real_t tmp_68 = tmp_15*tmp_66;
      real_t tmp_69 = tmp_23*(tmp_17 + tmp_61);
      real_t tmp_70 = tmp_23*(tmp_20 + tmp_65);
      real_t tmp_71 = 0.35539032758428413*tmp_18*(tmp_19*tmp_69 + tmp_22*tmp_70 - 1.0/3.0) + 0.35539032758428413*tmp_22*(tmp_21*tmp_69 + tmp_26*tmp_70 - 1.0/3.0);
      real_t a_0_0 = -tmp_27*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) - tmp_38*(-tmp_30 - tmp_31 - tmp_34 - tmp_35 + 1) - tmp_49*(-tmp_41 - tmp_42 - tmp_45 - tmp_46 + 1) - tmp_60*(-tmp_52 - tmp_53 - tmp_56 - tmp_57 + 1) - tmp_71*(-tmp_63 - tmp_64 - tmp_67 - tmp_68 + 1);
      real_t a_1_0 = -tmp_27*(tmp_10 + tmp_14) - tmp_38*(tmp_31 + tmp_34) - tmp_49*(tmp_42 + tmp_45) - tmp_60*(tmp_53 + tmp_56) - tmp_71*(tmp_64 + tmp_67);
      real_t a_2_0 = -tmp_27*(tmp_16 + tmp_8) - tmp_38*(tmp_30 + tmp_35) - tmp_49*(tmp_41 + tmp_46) - tmp_60*(tmp_52 + tmp_57) - tmp_71*(tmp_63 + tmp_68);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_1_1 + tmp_2;
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_2_0 + tmp_0));
      real_t tmp_6 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7 = p_affine_6_1 + tmp_2;
      real_t tmp_8 = tmp_5*(0.046910077030668018*tmp_6 + tmp_7);
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5*(0.046910077030668018*tmp_12 + tmp_13);
      real_t tmp_15 = tmp_14*tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 0.35539032758428413*tmp_3*(tmp_19 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_18 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 0.71794300574904923*tmp_3*(tmp_28 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_27 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 0.8533333333333335*tmp_3*(tmp_37 - 1.0/3.0) + 0.8533333333333335*tmp_4*(tmp_36 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 0.71794300574904923*tmp_3*(tmp_46 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_45 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.35539032758428413*tmp_3*(tmp_55 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_54 - 1.0/3.0);
      real_t a_0_0 = tmp_20*(-tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1) + tmp_29*(-tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1) + tmp_38*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_47*(-tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1) + tmp_56*(-tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1);
      real_t a_1_0 = tmp_18*tmp_20 + tmp_27*tmp_29 + tmp_36*tmp_38 + tmp_45*tmp_47 + tmp_54*tmp_56;
      real_t a_2_0 = tmp_19*tmp_20 + tmp_28*tmp_29 + tmp_37*tmp_38 + tmp_46*tmp_47 + tmp_55*tmp_56;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_2;
      real_t tmp_9 = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_11;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_16 + tmp_4*tmp_9 - tmp_7*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_21 = tmp_18*(tmp_12 - tmp_17);
      real_t tmp_22 = tmp_11*tmp_19 + tmp_20*tmp_3 + tmp_21*tmp_6;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_11*tmp_23 + tmp_24*tmp_3 + tmp_25*tmp_6;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_11*tmp_27 + tmp_28*tmp_3 + tmp_29*tmp_6;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_1_0 = 0.1666666666666668*tmp_45;
      real_t a_2_0 = 0.1666666666666668*tmp_46;
      real_t a_3_0 = 0.1666666666666668*tmp_47;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 3.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_54 - 1.0/4.0) + tmp_3*(tmp_53 - 1.0/4.0) + tmp_5*(tmp_52 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_74 - 1.0/4.0) + tmp_3*(tmp_73 - 1.0/4.0) + tmp_5*(tmp_72 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_90 - 1.0/4.0) + tmp_3*(tmp_89 - 1.0/4.0) + tmp_5*(tmp_88 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_106 - 1.0/4.0) + tmp_3*(tmp_105 - 1.0/4.0) + tmp_5*(tmp_104 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_122 - 1.0/4.0) + tmp_3*(tmp_121 - 1.0/4.0) + tmp_5*(tmp_120 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_138 - 1.0/4.0) + tmp_3*(tmp_137 - 1.0/4.0) + tmp_5*(tmp_136 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_154 - 1.0/4.0) + tmp_3*(tmp_153 - 1.0/4.0) + tmp_5*(tmp_152 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_170 - 1.0/4.0) + tmp_3*(tmp_169 - 1.0/4.0) + tmp_5*(tmp_168 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_186 - 1.0/4.0) + tmp_3*(tmp_185 - 1.0/4.0) + tmp_5*(tmp_184 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_202 - 1.0/4.0) + tmp_3*(tmp_201 - 1.0/4.0) + tmp_5*(tmp_200 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_218 - 1.0/4.0) + tmp_3*(tmp_217 - 1.0/4.0) + tmp_5*(tmp_216 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_234 - 1.0/4.0) + tmp_3*(tmp_233 - 1.0/4.0) + tmp_5*(tmp_232 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_250 - 1.0/4.0) + tmp_3*(tmp_249 - 1.0/4.0) + tmp_5*(tmp_248 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_266 - 1.0/4.0) + tmp_3*(tmp_265 - 1.0/4.0) + tmp_5*(tmp_264 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_282 - 1.0/4.0) + tmp_3*(tmp_281 - 1.0/4.0) + tmp_5*(tmp_280 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_298 - 1.0/4.0) + tmp_3*(tmp_297 - 1.0/4.0) + tmp_5*(tmp_296 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_314 - 1.0/4.0) + tmp_3*(tmp_313 - 1.0/4.0) + tmp_5*(tmp_312 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_330 - 1.0/4.0) + tmp_3*(tmp_329 - 1.0/4.0) + tmp_5*(tmp_328 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_346 - 1.0/4.0) + tmp_3*(tmp_345 - 1.0/4.0) + tmp_5*(tmp_344 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_362 - 1.0/4.0) + tmp_3*(tmp_361 - 1.0/4.0) + tmp_5*(tmp_360 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_378 - 1.0/4.0) + tmp_3*(tmp_377 - 1.0/4.0) + tmp_5*(tmp_376 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_1_0 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_2_0 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_3_0 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = 0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22;
      real_t tmp_24 = p_affine_8_2 + tmp_9;
      real_t tmp_25 = tmp_19*(tmp_23 + tmp_24);
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_28 = tmp_25*tmp_27;
      real_t tmp_29 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_30 = -p_affine_8_1;
      real_t tmp_31 = p_affine_9_1 + tmp_30;
      real_t tmp_32 = p_affine_10_1 + tmp_30;
      real_t tmp_33 = 0.031405749086161582*tmp_31 + 0.93718850182767688*tmp_32;
      real_t tmp_34 = p_affine_8_1 + tmp_2;
      real_t tmp_35 = tmp_19*(tmp_33 + tmp_34);
      real_t tmp_36 = tmp_29*tmp_35;
      real_t tmp_37 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_40 = tmp_25*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_42 = tmp_35*tmp_41;
      real_t tmp_43 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_44 = -p_affine_8_0;
      real_t tmp_45 = p_affine_9_0 + tmp_44;
      real_t tmp_46 = p_affine_10_0 + tmp_44;
      real_t tmp_47 = 0.031405749086161582*tmp_45 + 0.93718850182767688*tmp_46;
      real_t tmp_48 = p_affine_8_0 + tmp_0;
      real_t tmp_49 = tmp_19*(tmp_47 + tmp_48);
      real_t tmp_50 = tmp_43*tmp_49;
      real_t tmp_51 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_52 = tmp_49*tmp_51;
      real_t tmp_53 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_54 = tmp_49*tmp_53;
      real_t tmp_55 = -p_affine_4_1;
      real_t tmp_56 = p_affine_5_1 + tmp_55;
      real_t tmp_57 = -p_affine_4_0;
      real_t tmp_58 = p_affine_6_0 + tmp_57;
      real_t tmp_59 = p_affine_7_1 + tmp_55;
      real_t tmp_60 = tmp_58*tmp_59;
      real_t tmp_61 = p_affine_7_0 + tmp_57;
      real_t tmp_62 = p_affine_6_1 + tmp_55;
      real_t tmp_63 = tmp_61*tmp_62;
      real_t tmp_64 = tmp_60 - tmp_63;
      real_t tmp_65 = p_affine_5_0 + tmp_57;
      real_t tmp_66 = -p_affine_4_2;
      real_t tmp_67 = p_affine_7_2 + tmp_66;
      real_t tmp_68 = tmp_62*tmp_67;
      real_t tmp_69 = p_affine_5_2 + tmp_66;
      real_t tmp_70 = p_affine_6_2 + tmp_66;
      real_t tmp_71 = tmp_61*tmp_70;
      real_t tmp_72 = tmp_59*tmp_70;
      real_t tmp_73 = tmp_58*tmp_67;
      real_t tmp_74 = 1.0 / (tmp_56*tmp_71 - tmp_56*tmp_73 + tmp_60*tmp_69 - tmp_63*tmp_69 + tmp_65*tmp_68 - tmp_65*tmp_72);
      real_t tmp_75 = p_affine_8_2 + tmp_66;
      real_t tmp_76 = tmp_74*(tmp_23 + tmp_75);
      real_t tmp_77 = tmp_71 - tmp_73;
      real_t tmp_78 = p_affine_8_1 + tmp_55;
      real_t tmp_79 = tmp_74*(tmp_33 + tmp_78);
      real_t tmp_80 = tmp_68 - tmp_72;
      real_t tmp_81 = p_affine_8_0 + tmp_57;
      real_t tmp_82 = tmp_74*(tmp_47 + tmp_81);
      real_t tmp_83 = tmp_56*tmp_61 - tmp_59*tmp_65;
      real_t tmp_84 = -tmp_61*tmp_69 + tmp_65*tmp_67;
      real_t tmp_85 = -tmp_56*tmp_67 + tmp_59*tmp_69;
      real_t tmp_86 = -tmp_56*tmp_58 + tmp_62*tmp_65;
      real_t tmp_87 = tmp_58*tmp_69 - tmp_65*tmp_70;
      real_t tmp_88 = tmp_56*tmp_70 - tmp_62*tmp_69;
      real_t tmp_89 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_90 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_91 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_92 = 3.0*std::pow((std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)*std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)) + (std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)*std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)) + (std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)*std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)), 0.25);
      real_t tmp_93 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_64*tmp_76 + tmp_77*tmp_79 + tmp_80*tmp_82 - 1.0/4.0) + tmp_59*(tmp_76*tmp_86 + tmp_79*tmp_87 + tmp_82*tmp_88 - 1.0/4.0) + tmp_62*(tmp_76*tmp_83 + tmp_79*tmp_84 + tmp_82*tmp_85 - 1.0/4.0));
      real_t tmp_94 = 0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22;
      real_t tmp_95 = tmp_19*(tmp_24 + tmp_94);
      real_t tmp_96 = tmp_8*tmp_95;
      real_t tmp_97 = tmp_27*tmp_95;
      real_t tmp_98 = 0.19601935860219369*tmp_31 + 0.60796128279561268*tmp_32;
      real_t tmp_99 = tmp_19*(tmp_34 + tmp_98);
      real_t tmp_100 = tmp_29*tmp_99;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = tmp_39*tmp_95;
      real_t tmp_103 = tmp_41*tmp_99;
      real_t tmp_104 = 0.19601935860219369*tmp_45 + 0.60796128279561268*tmp_46;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_48);
      real_t tmp_106 = tmp_105*tmp_43;
      real_t tmp_107 = tmp_105*tmp_51;
      real_t tmp_108 = tmp_105*tmp_53;
      real_t tmp_109 = tmp_74*(tmp_75 + tmp_94);
      real_t tmp_110 = tmp_74*(tmp_78 + tmp_98);
      real_t tmp_111 = tmp_74*(tmp_104 + tmp_81);
      real_t tmp_112 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_109*tmp_64 + tmp_110*tmp_77 + tmp_111*tmp_80 - 1.0/4.0) + tmp_59*(tmp_109*tmp_86 + tmp_110*tmp_87 + tmp_111*tmp_88 - 1.0/4.0) + tmp_62*(tmp_109*tmp_83 + tmp_110*tmp_84 + tmp_111*tmp_85 - 1.0/4.0));
      real_t tmp_113 = 0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_114 = tmp_19*(tmp_113 + tmp_24);
      real_t tmp_115 = tmp_114*tmp_8;
      real_t tmp_116 = tmp_114*tmp_27;
      real_t tmp_117 = 0.37605877282253791*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_118 = tmp_19*(tmp_117 + tmp_34);
      real_t tmp_119 = tmp_118*tmp_29;
      real_t tmp_120 = tmp_118*tmp_37;
      real_t tmp_121 = tmp_114*tmp_39;
      real_t tmp_122 = tmp_118*tmp_41;
      real_t tmp_123 = 0.37605877282253791*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_124 = tmp_19*(tmp_123 + tmp_48);
      real_t tmp_125 = tmp_124*tmp_43;
      real_t tmp_126 = tmp_124*tmp_51;
      real_t tmp_127 = tmp_124*tmp_53;
      real_t tmp_128 = tmp_74*(tmp_113 + tmp_75);
      real_t tmp_129 = tmp_74*(tmp_117 + tmp_78);
      real_t tmp_130 = tmp_74*(tmp_123 + tmp_81);
      real_t tmp_131 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_128*tmp_64 + tmp_129*tmp_77 + tmp_130*tmp_80 - 1.0/4.0) + tmp_59*(tmp_128*tmp_86 + tmp_129*tmp_87 + tmp_130*tmp_88 - 1.0/4.0) + tmp_62*(tmp_128*tmp_83 + tmp_129*tmp_84 + tmp_130*tmp_85 - 1.0/4.0));
      real_t tmp_132 = 0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_133 = tmp_19*(tmp_132 + tmp_24);
      real_t tmp_134 = tmp_133*tmp_8;
      real_t tmp_135 = tmp_133*tmp_27;
      real_t tmp_136 = 0.78764240869137092*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_137 = tmp_19*(tmp_136 + tmp_34);
      real_t tmp_138 = tmp_137*tmp_29;
      real_t tmp_139 = tmp_137*tmp_37;
      real_t tmp_140 = tmp_133*tmp_39;
      real_t tmp_141 = tmp_137*tmp_41;
      real_t tmp_142 = 0.78764240869137092*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_143 = tmp_19*(tmp_142 + tmp_48);
      real_t tmp_144 = tmp_143*tmp_43;
      real_t tmp_145 = tmp_143*tmp_51;
      real_t tmp_146 = tmp_143*tmp_53;
      real_t tmp_147 = tmp_74*(tmp_132 + tmp_75);
      real_t tmp_148 = tmp_74*(tmp_136 + tmp_78);
      real_t tmp_149 = tmp_74*(tmp_142 + tmp_81);
      real_t tmp_150 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_147*tmp_64 + tmp_148*tmp_77 + tmp_149*tmp_80 - 1.0/4.0) + tmp_59*(tmp_147*tmp_86 + tmp_148*tmp_87 + tmp_149*tmp_88 - 1.0/4.0) + tmp_62*(tmp_147*tmp_83 + tmp_148*tmp_84 + tmp_149*tmp_85 - 1.0/4.0));
      real_t tmp_151 = 0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_152 = tmp_19*(tmp_151 + tmp_24);
      real_t tmp_153 = tmp_152*tmp_8;
      real_t tmp_154 = tmp_152*tmp_27;
      real_t tmp_155 = 0.58463275527740355*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_156 = tmp_19*(tmp_155 + tmp_34);
      real_t tmp_157 = tmp_156*tmp_29;
      real_t tmp_158 = tmp_156*tmp_37;
      real_t tmp_159 = tmp_152*tmp_39;
      real_t tmp_160 = tmp_156*tmp_41;
      real_t tmp_161 = 0.58463275527740355*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_162 = tmp_19*(tmp_161 + tmp_48);
      real_t tmp_163 = tmp_162*tmp_43;
      real_t tmp_164 = tmp_162*tmp_51;
      real_t tmp_165 = tmp_162*tmp_53;
      real_t tmp_166 = tmp_74*(tmp_151 + tmp_75);
      real_t tmp_167 = tmp_74*(tmp_155 + tmp_78);
      real_t tmp_168 = tmp_74*(tmp_161 + tmp_81);
      real_t tmp_169 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_166*tmp_64 + tmp_167*tmp_77 + tmp_168*tmp_80 - 1.0/4.0) + tmp_59*(tmp_166*tmp_86 + tmp_167*tmp_87 + tmp_168*tmp_88 - 1.0/4.0) + tmp_62*(tmp_166*tmp_83 + tmp_167*tmp_84 + tmp_168*tmp_85 - 1.0/4.0));
      real_t tmp_170 = 0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_24);
      real_t tmp_172 = tmp_171*tmp_8;
      real_t tmp_173 = tmp_171*tmp_27;
      real_t tmp_174 = 0.041227165399737475*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_34);
      real_t tmp_176 = tmp_175*tmp_29;
      real_t tmp_177 = tmp_175*tmp_37;
      real_t tmp_178 = tmp_171*tmp_39;
      real_t tmp_179 = tmp_175*tmp_41;
      real_t tmp_180 = 0.041227165399737475*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_181 = tmp_19*(tmp_180 + tmp_48);
      real_t tmp_182 = tmp_181*tmp_43;
      real_t tmp_183 = tmp_181*tmp_51;
      real_t tmp_184 = tmp_181*tmp_53;
      real_t tmp_185 = tmp_74*(tmp_170 + tmp_75);
      real_t tmp_186 = tmp_74*(tmp_174 + tmp_78);
      real_t tmp_187 = tmp_74*(tmp_180 + tmp_81);
      real_t tmp_188 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_185*tmp_64 + tmp_186*tmp_77 + tmp_187*tmp_80 - 1.0/4.0) + tmp_59*(tmp_185*tmp_86 + tmp_186*tmp_87 + tmp_187*tmp_88 - 1.0/4.0) + tmp_62*(tmp_185*tmp_83 + tmp_186*tmp_84 + tmp_187*tmp_85 - 1.0/4.0));
      real_t tmp_189 = 0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_190 = tmp_19*(tmp_189 + tmp_24);
      real_t tmp_191 = tmp_190*tmp_8;
      real_t tmp_192 = tmp_190*tmp_27;
      real_t tmp_193 = 0.039308471900058539*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_194 = tmp_19*(tmp_193 + tmp_34);
      real_t tmp_195 = tmp_194*tmp_29;
      real_t tmp_196 = tmp_194*tmp_37;
      real_t tmp_197 = tmp_190*tmp_39;
      real_t tmp_198 = tmp_194*tmp_41;
      real_t tmp_199 = 0.039308471900058539*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_200 = tmp_19*(tmp_199 + tmp_48);
      real_t tmp_201 = tmp_200*tmp_43;
      real_t tmp_202 = tmp_200*tmp_51;
      real_t tmp_203 = tmp_200*tmp_53;
      real_t tmp_204 = tmp_74*(tmp_189 + tmp_75);
      real_t tmp_205 = tmp_74*(tmp_193 + tmp_78);
      real_t tmp_206 = tmp_74*(tmp_199 + tmp_81);
      real_t tmp_207 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_204*tmp_64 + tmp_205*tmp_77 + tmp_206*tmp_80 - 1.0/4.0) + tmp_59*(tmp_204*tmp_86 + tmp_205*tmp_87 + tmp_206*tmp_88 - 1.0/4.0) + tmp_62*(tmp_204*tmp_83 + tmp_205*tmp_84 + tmp_206*tmp_85 - 1.0/4.0));
      real_t tmp_208 = 0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_209 = tmp_19*(tmp_208 + tmp_24);
      real_t tmp_210 = tmp_209*tmp_8;
      real_t tmp_211 = tmp_209*tmp_27;
      real_t tmp_212 = 0.78764240869137092*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_213 = tmp_19*(tmp_212 + tmp_34);
      real_t tmp_214 = tmp_213*tmp_29;
      real_t tmp_215 = tmp_213*tmp_37;
      real_t tmp_216 = tmp_209*tmp_39;
      real_t tmp_217 = tmp_213*tmp_41;
      real_t tmp_218 = 0.78764240869137092*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_48);
      real_t tmp_220 = tmp_219*tmp_43;
      real_t tmp_221 = tmp_219*tmp_51;
      real_t tmp_222 = tmp_219*tmp_53;
      real_t tmp_223 = tmp_74*(tmp_208 + tmp_75);
      real_t tmp_224 = tmp_74*(tmp_212 + tmp_78);
      real_t tmp_225 = tmp_74*(tmp_218 + tmp_81);
      real_t tmp_226 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_223*tmp_64 + tmp_224*tmp_77 + tmp_225*tmp_80 - 1.0/4.0) + tmp_59*(tmp_223*tmp_86 + tmp_224*tmp_87 + tmp_225*tmp_88 - 1.0/4.0) + tmp_62*(tmp_223*tmp_83 + tmp_224*tmp_84 + tmp_225*tmp_85 - 1.0/4.0));
      real_t tmp_227 = 0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_228 = tmp_19*(tmp_227 + tmp_24);
      real_t tmp_229 = tmp_228*tmp_8;
      real_t tmp_230 = tmp_228*tmp_27;
      real_t tmp_231 = 0.58463275527740355*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_232 = tmp_19*(tmp_231 + tmp_34);
      real_t tmp_233 = tmp_232*tmp_29;
      real_t tmp_234 = tmp_232*tmp_37;
      real_t tmp_235 = tmp_228*tmp_39;
      real_t tmp_236 = tmp_232*tmp_41;
      real_t tmp_237 = 0.58463275527740355*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_238 = tmp_19*(tmp_237 + tmp_48);
      real_t tmp_239 = tmp_238*tmp_43;
      real_t tmp_240 = tmp_238*tmp_51;
      real_t tmp_241 = tmp_238*tmp_53;
      real_t tmp_242 = tmp_74*(tmp_227 + tmp_75);
      real_t tmp_243 = tmp_74*(tmp_231 + tmp_78);
      real_t tmp_244 = tmp_74*(tmp_237 + tmp_81);
      real_t tmp_245 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_242*tmp_64 + tmp_243*tmp_77 + tmp_244*tmp_80 - 1.0/4.0) + tmp_59*(tmp_242*tmp_86 + tmp_243*tmp_87 + tmp_244*tmp_88 - 1.0/4.0) + tmp_62*(tmp_242*tmp_83 + tmp_243*tmp_84 + tmp_244*tmp_85 - 1.0/4.0));
      real_t tmp_246 = 0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_247 = tmp_19*(tmp_24 + tmp_246);
      real_t tmp_248 = tmp_247*tmp_8;
      real_t tmp_249 = tmp_247*tmp_27;
      real_t tmp_250 = 0.1711304259088916*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_251 = tmp_19*(tmp_250 + tmp_34);
      real_t tmp_252 = tmp_251*tmp_29;
      real_t tmp_253 = tmp_251*tmp_37;
      real_t tmp_254 = tmp_247*tmp_39;
      real_t tmp_255 = tmp_251*tmp_41;
      real_t tmp_256 = 0.1711304259088916*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_257 = tmp_19*(tmp_256 + tmp_48);
      real_t tmp_258 = tmp_257*tmp_43;
      real_t tmp_259 = tmp_257*tmp_51;
      real_t tmp_260 = tmp_257*tmp_53;
      real_t tmp_261 = tmp_74*(tmp_246 + tmp_75);
      real_t tmp_262 = tmp_74*(tmp_250 + tmp_78);
      real_t tmp_263 = tmp_74*(tmp_256 + tmp_81);
      real_t tmp_264 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_261*tmp_64 + tmp_262*tmp_77 + tmp_263*tmp_80 - 1.0/4.0) + tmp_59*(tmp_261*tmp_86 + tmp_262*tmp_87 + tmp_263*tmp_88 - 1.0/4.0) + tmp_62*(tmp_261*tmp_83 + tmp_262*tmp_84 + tmp_263*tmp_85 - 1.0/4.0));
      real_t tmp_265 = 0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_266 = tmp_19*(tmp_24 + tmp_265);
      real_t tmp_267 = tmp_266*tmp_8;
      real_t tmp_268 = tmp_266*tmp_27;
      real_t tmp_269 = 0.37605877282253791*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_270 = tmp_19*(tmp_269 + tmp_34);
      real_t tmp_271 = tmp_270*tmp_29;
      real_t tmp_272 = tmp_270*tmp_37;
      real_t tmp_273 = tmp_266*tmp_39;
      real_t tmp_274 = tmp_270*tmp_41;
      real_t tmp_275 = 0.37605877282253791*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_276 = tmp_19*(tmp_275 + tmp_48);
      real_t tmp_277 = tmp_276*tmp_43;
      real_t tmp_278 = tmp_276*tmp_51;
      real_t tmp_279 = tmp_276*tmp_53;
      real_t tmp_280 = tmp_74*(tmp_265 + tmp_75);
      real_t tmp_281 = tmp_74*(tmp_269 + tmp_78);
      real_t tmp_282 = tmp_74*(tmp_275 + tmp_81);
      real_t tmp_283 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_280*tmp_64 + tmp_281*tmp_77 + tmp_282*tmp_80 - 1.0/4.0) + tmp_59*(tmp_280*tmp_86 + tmp_281*tmp_87 + tmp_282*tmp_88 - 1.0/4.0) + tmp_62*(tmp_280*tmp_83 + tmp_281*tmp_84 + tmp_282*tmp_85 - 1.0/4.0));
      real_t tmp_284 = 0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_285 = tmp_19*(tmp_24 + tmp_284);
      real_t tmp_286 = tmp_285*tmp_8;
      real_t tmp_287 = tmp_27*tmp_285;
      real_t tmp_288 = 0.041227165399737475*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_289 = tmp_19*(tmp_288 + tmp_34);
      real_t tmp_290 = tmp_289*tmp_29;
      real_t tmp_291 = tmp_289*tmp_37;
      real_t tmp_292 = tmp_285*tmp_39;
      real_t tmp_293 = tmp_289*tmp_41;
      real_t tmp_294 = 0.041227165399737475*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_295 = tmp_19*(tmp_294 + tmp_48);
      real_t tmp_296 = tmp_295*tmp_43;
      real_t tmp_297 = tmp_295*tmp_51;
      real_t tmp_298 = tmp_295*tmp_53;
      real_t tmp_299 = tmp_74*(tmp_284 + tmp_75);
      real_t tmp_300 = tmp_74*(tmp_288 + tmp_78);
      real_t tmp_301 = tmp_74*(tmp_294 + tmp_81);
      real_t tmp_302 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_299*tmp_64 + tmp_300*tmp_77 + tmp_301*tmp_80 - 1.0/4.0) + tmp_59*(tmp_299*tmp_86 + tmp_300*tmp_87 + tmp_301*tmp_88 - 1.0/4.0) + tmp_62*(tmp_299*tmp_83 + tmp_300*tmp_84 + tmp_301*tmp_85 - 1.0/4.0));
      real_t tmp_303 = 0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22;
      real_t tmp_304 = tmp_19*(tmp_24 + tmp_303);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_27*tmp_304;
      real_t tmp_307 = 0.40446199974765351*tmp_31 + 0.19107600050469298*tmp_32;
      real_t tmp_308 = tmp_19*(tmp_307 + tmp_34);
      real_t tmp_309 = tmp_29*tmp_308;
      real_t tmp_310 = tmp_308*tmp_37;
      real_t tmp_311 = tmp_304*tmp_39;
      real_t tmp_312 = tmp_308*tmp_41;
      real_t tmp_313 = 0.40446199974765351*tmp_45 + 0.19107600050469298*tmp_46;
      real_t tmp_314 = tmp_19*(tmp_313 + tmp_48);
      real_t tmp_315 = tmp_314*tmp_43;
      real_t tmp_316 = tmp_314*tmp_51;
      real_t tmp_317 = tmp_314*tmp_53;
      real_t tmp_318 = tmp_74*(tmp_303 + tmp_75);
      real_t tmp_319 = tmp_74*(tmp_307 + tmp_78);
      real_t tmp_320 = tmp_74*(tmp_313 + tmp_81);
      real_t tmp_321 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_318*tmp_64 + tmp_319*tmp_77 + tmp_320*tmp_80 - 1.0/4.0) + tmp_59*(tmp_318*tmp_86 + tmp_319*tmp_87 + tmp_320*tmp_88 - 1.0/4.0) + tmp_62*(tmp_318*tmp_83 + tmp_319*tmp_84 + tmp_320*tmp_85 - 1.0/4.0));
      real_t tmp_322 = 0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_323 = tmp_19*(tmp_24 + tmp_322);
      real_t tmp_324 = tmp_323*tmp_8;
      real_t tmp_325 = tmp_27*tmp_323;
      real_t tmp_326 = 0.039308471900058539*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_327 = tmp_19*(tmp_326 + tmp_34);
      real_t tmp_328 = tmp_29*tmp_327;
      real_t tmp_329 = tmp_327*tmp_37;
      real_t tmp_330 = tmp_323*tmp_39;
      real_t tmp_331 = tmp_327*tmp_41;
      real_t tmp_332 = 0.039308471900058539*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_333 = tmp_19*(tmp_332 + tmp_48);
      real_t tmp_334 = tmp_333*tmp_43;
      real_t tmp_335 = tmp_333*tmp_51;
      real_t tmp_336 = tmp_333*tmp_53;
      real_t tmp_337 = tmp_74*(tmp_322 + tmp_75);
      real_t tmp_338 = tmp_74*(tmp_326 + tmp_78);
      real_t tmp_339 = tmp_74*(tmp_332 + tmp_81);
      real_t tmp_340 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_337*tmp_64 + tmp_338*tmp_77 + tmp_339*tmp_80 - 1.0/4.0) + tmp_59*(tmp_337*tmp_86 + tmp_338*tmp_87 + tmp_339*tmp_88 - 1.0/4.0) + tmp_62*(tmp_337*tmp_83 + tmp_338*tmp_84 + tmp_339*tmp_85 - 1.0/4.0));
      real_t tmp_341 = 0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_342 = tmp_19*(tmp_24 + tmp_341);
      real_t tmp_343 = tmp_342*tmp_8;
      real_t tmp_344 = tmp_27*tmp_342;
      real_t tmp_345 = 0.93718850182767688*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_346 = tmp_19*(tmp_34 + tmp_345);
      real_t tmp_347 = tmp_29*tmp_346;
      real_t tmp_348 = tmp_346*tmp_37;
      real_t tmp_349 = tmp_342*tmp_39;
      real_t tmp_350 = tmp_346*tmp_41;
      real_t tmp_351 = 0.93718850182767688*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_352 = tmp_19*(tmp_351 + tmp_48);
      real_t tmp_353 = tmp_352*tmp_43;
      real_t tmp_354 = tmp_352*tmp_51;
      real_t tmp_355 = tmp_352*tmp_53;
      real_t tmp_356 = tmp_74*(tmp_341 + tmp_75);
      real_t tmp_357 = tmp_74*(tmp_345 + tmp_78);
      real_t tmp_358 = tmp_74*(tmp_351 + tmp_81);
      real_t tmp_359 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_356*tmp_64 + tmp_357*tmp_77 + tmp_358*tmp_80 - 1.0/4.0) + tmp_59*(tmp_356*tmp_86 + tmp_357*tmp_87 + tmp_358*tmp_88 - 1.0/4.0) + tmp_62*(tmp_356*tmp_83 + tmp_357*tmp_84 + tmp_358*tmp_85 - 1.0/4.0));
      real_t tmp_360 = 0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_361 = tmp_19*(tmp_24 + tmp_360);
      real_t tmp_362 = tmp_361*tmp_8;
      real_t tmp_363 = tmp_27*tmp_361;
      real_t tmp_364 = 0.60796128279561268*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_365 = tmp_19*(tmp_34 + tmp_364);
      real_t tmp_366 = tmp_29*tmp_365;
      real_t tmp_367 = tmp_365*tmp_37;
      real_t tmp_368 = tmp_361*tmp_39;
      real_t tmp_369 = tmp_365*tmp_41;
      real_t tmp_370 = 0.60796128279561268*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_371 = tmp_19*(tmp_370 + tmp_48);
      real_t tmp_372 = tmp_371*tmp_43;
      real_t tmp_373 = tmp_371*tmp_51;
      real_t tmp_374 = tmp_371*tmp_53;
      real_t tmp_375 = tmp_74*(tmp_360 + tmp_75);
      real_t tmp_376 = tmp_74*(tmp_364 + tmp_78);
      real_t tmp_377 = tmp_74*(tmp_370 + tmp_81);
      real_t tmp_378 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_375*tmp_64 + tmp_376*tmp_77 + tmp_377*tmp_80 - 1.0/4.0) + tmp_59*(tmp_375*tmp_86 + tmp_376*tmp_87 + tmp_377*tmp_88 - 1.0/4.0) + tmp_62*(tmp_375*tmp_83 + tmp_376*tmp_84 + tmp_377*tmp_85 - 1.0/4.0));
      real_t tmp_379 = 0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_380 = tmp_19*(tmp_24 + tmp_379);
      real_t tmp_381 = tmp_380*tmp_8;
      real_t tmp_382 = tmp_27*tmp_380;
      real_t tmp_383 = 0.19107600050469298*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_384 = tmp_19*(tmp_34 + tmp_383);
      real_t tmp_385 = tmp_29*tmp_384;
      real_t tmp_386 = tmp_37*tmp_384;
      real_t tmp_387 = tmp_380*tmp_39;
      real_t tmp_388 = tmp_384*tmp_41;
      real_t tmp_389 = 0.19107600050469298*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_390 = tmp_19*(tmp_389 + tmp_48);
      real_t tmp_391 = tmp_390*tmp_43;
      real_t tmp_392 = tmp_390*tmp_51;
      real_t tmp_393 = tmp_390*tmp_53;
      real_t tmp_394 = tmp_74*(tmp_379 + tmp_75);
      real_t tmp_395 = tmp_74*(tmp_383 + tmp_78);
      real_t tmp_396 = tmp_74*(tmp_389 + tmp_81);
      real_t tmp_397 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_394*tmp_64 + tmp_395*tmp_77 + tmp_396*tmp_80 - 1.0/4.0) + tmp_59*(tmp_394*tmp_86 + tmp_395*tmp_87 + tmp_396*tmp_88 - 1.0/4.0) + tmp_62*(tmp_394*tmp_83 + tmp_395*tmp_84 + tmp_396*tmp_85 - 1.0/4.0));
      real_t tmp_398 = 0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_399 = tmp_19*(tmp_24 + tmp_398);
      real_t tmp_400 = tmp_399*tmp_8;
      real_t tmp_401 = tmp_27*tmp_399;
      real_t tmp_402 = 0.031405749086161582*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_403 = tmp_19*(tmp_34 + tmp_402);
      real_t tmp_404 = tmp_29*tmp_403;
      real_t tmp_405 = tmp_37*tmp_403;
      real_t tmp_406 = tmp_39*tmp_399;
      real_t tmp_407 = tmp_403*tmp_41;
      real_t tmp_408 = 0.031405749086161582*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_409 = tmp_19*(tmp_408 + tmp_48);
      real_t tmp_410 = tmp_409*tmp_43;
      real_t tmp_411 = tmp_409*tmp_51;
      real_t tmp_412 = tmp_409*tmp_53;
      real_t tmp_413 = tmp_74*(tmp_398 + tmp_75);
      real_t tmp_414 = tmp_74*(tmp_402 + tmp_78);
      real_t tmp_415 = tmp_74*(tmp_408 + tmp_81);
      real_t tmp_416 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_413*tmp_64 + tmp_414*tmp_77 + tmp_415*tmp_80 - 1.0/4.0) + tmp_59*(tmp_413*tmp_86 + tmp_414*tmp_87 + tmp_415*tmp_88 - 1.0/4.0) + tmp_62*(tmp_413*tmp_83 + tmp_414*tmp_84 + tmp_415*tmp_85 - 1.0/4.0));
      real_t tmp_417 = 0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_418 = tmp_19*(tmp_24 + tmp_417);
      real_t tmp_419 = tmp_418*tmp_8;
      real_t tmp_420 = tmp_27*tmp_418;
      real_t tmp_421 = 0.19601935860219369*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_422 = tmp_19*(tmp_34 + tmp_421);
      real_t tmp_423 = tmp_29*tmp_422;
      real_t tmp_424 = tmp_37*tmp_422;
      real_t tmp_425 = tmp_39*tmp_418;
      real_t tmp_426 = tmp_41*tmp_422;
      real_t tmp_427 = 0.19601935860219369*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_428 = tmp_19*(tmp_427 + tmp_48);
      real_t tmp_429 = tmp_428*tmp_43;
      real_t tmp_430 = tmp_428*tmp_51;
      real_t tmp_431 = tmp_428*tmp_53;
      real_t tmp_432 = tmp_74*(tmp_417 + tmp_75);
      real_t tmp_433 = tmp_74*(tmp_421 + tmp_78);
      real_t tmp_434 = tmp_74*(tmp_427 + tmp_81);
      real_t tmp_435 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_432*tmp_64 + tmp_433*tmp_77 + tmp_434*tmp_80 - 1.0/4.0) + tmp_59*(tmp_432*tmp_86 + tmp_433*tmp_87 + tmp_434*tmp_88 - 1.0/4.0) + tmp_62*(tmp_432*tmp_83 + tmp_433*tmp_84 + tmp_434*tmp_85 - 1.0/4.0));
      real_t tmp_436 = 0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_437 = tmp_19*(tmp_24 + tmp_436);
      real_t tmp_438 = tmp_437*tmp_8;
      real_t tmp_439 = tmp_27*tmp_437;
      real_t tmp_440 = 0.40446199974765351*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_441 = tmp_19*(tmp_34 + tmp_440);
      real_t tmp_442 = tmp_29*tmp_441;
      real_t tmp_443 = tmp_37*tmp_441;
      real_t tmp_444 = tmp_39*tmp_437;
      real_t tmp_445 = tmp_41*tmp_441;
      real_t tmp_446 = 0.40446199974765351*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_447 = tmp_19*(tmp_446 + tmp_48);
      real_t tmp_448 = tmp_43*tmp_447;
      real_t tmp_449 = tmp_447*tmp_51;
      real_t tmp_450 = tmp_447*tmp_53;
      real_t tmp_451 = tmp_74*(tmp_436 + tmp_75);
      real_t tmp_452 = tmp_74*(tmp_440 + tmp_78);
      real_t tmp_453 = tmp_74*(tmp_446 + tmp_81);
      real_t tmp_454 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_451*tmp_64 + tmp_452*tmp_77 + tmp_453*tmp_80 - 1.0/4.0) + tmp_59*(tmp_451*tmp_86 + tmp_452*tmp_87 + tmp_453*tmp_88 - 1.0/4.0) + tmp_62*(tmp_451*tmp_83 + tmp_452*tmp_84 + tmp_453*tmp_85 - 1.0/4.0));
      real_t tmp_455 = 0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_456 = tmp_19*(tmp_24 + tmp_455);
      real_t tmp_457 = tmp_456*tmp_8;
      real_t tmp_458 = tmp_27*tmp_456;
      real_t tmp_459 = 0.1711304259088916*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_460 = tmp_19*(tmp_34 + tmp_459);
      real_t tmp_461 = tmp_29*tmp_460;
      real_t tmp_462 = tmp_37*tmp_460;
      real_t tmp_463 = tmp_39*tmp_456;
      real_t tmp_464 = tmp_41*tmp_460;
      real_t tmp_465 = 0.1711304259088916*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_466 = tmp_19*(tmp_465 + tmp_48);
      real_t tmp_467 = tmp_43*tmp_466;
      real_t tmp_468 = tmp_466*tmp_51;
      real_t tmp_469 = tmp_466*tmp_53;
      real_t tmp_470 = tmp_74*(tmp_455 + tmp_75);
      real_t tmp_471 = tmp_74*(tmp_459 + tmp_78);
      real_t tmp_472 = tmp_74*(tmp_465 + tmp_81);
      real_t tmp_473 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_470*tmp_64 + tmp_471*tmp_77 + tmp_472*tmp_80 - 1.0/4.0) + tmp_59*(tmp_470*tmp_86 + tmp_471*tmp_87 + tmp_472*tmp_88 - 1.0/4.0) + tmp_62*(tmp_470*tmp_83 + tmp_471*tmp_84 + tmp_472*tmp_85 - 1.0/4.0));
      real_t a_0_0 = -tmp_112*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_106 - tmp_107 - tmp_108 - tmp_96 - tmp_97 + 1) - tmp_131*(-tmp_115 - tmp_116 - tmp_119 - tmp_120 - tmp_121 - tmp_122 - tmp_125 - tmp_126 - tmp_127 + 1) - tmp_150*(-tmp_134 - tmp_135 - tmp_138 - tmp_139 - tmp_140 - tmp_141 - tmp_144 - tmp_145 - tmp_146 + 1) - tmp_169*(-tmp_153 - tmp_154 - tmp_157 - tmp_158 - tmp_159 - tmp_160 - tmp_163 - tmp_164 - tmp_165 + 1) - tmp_188*(-tmp_172 - tmp_173 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_182 - tmp_183 - tmp_184 + 1) - tmp_207*(-tmp_191 - tmp_192 - tmp_195 - tmp_196 - tmp_197 - tmp_198 - tmp_201 - tmp_202 - tmp_203 + 1) - tmp_226*(-tmp_210 - tmp_211 - tmp_214 - tmp_215 - tmp_216 - tmp_217 - tmp_220 - tmp_221 - tmp_222 + 1) - tmp_245*(-tmp_229 - tmp_230 - tmp_233 - tmp_234 - tmp_235 - tmp_236 - tmp_239 - tmp_240 - tmp_241 + 1) - tmp_264*(-tmp_248 - tmp_249 - tmp_252 - tmp_253 - tmp_254 - tmp_255 - tmp_258 - tmp_259 - tmp_260 + 1) - tmp_283*(-tmp_267 - tmp_268 - tmp_271 - tmp_272 - tmp_273 - tmp_274 - tmp_277 - tmp_278 - tmp_279 + 1) - tmp_302*(-tmp_286 - tmp_287 - tmp_290 - tmp_291 - tmp_292 - tmp_293 - tmp_296 - tmp_297 - tmp_298 + 1) - tmp_321*(-tmp_305 - tmp_306 - tmp_309 - tmp_310 - tmp_311 - tmp_312 - tmp_315 - tmp_316 - tmp_317 + 1) - tmp_340*(-tmp_324 - tmp_325 - tmp_328 - tmp_329 - tmp_330 - tmp_331 - tmp_334 - tmp_335 - tmp_336 + 1) - tmp_359*(-tmp_343 - tmp_344 - tmp_347 - tmp_348 - tmp_349 - tmp_350 - tmp_353 - tmp_354 - tmp_355 + 1) - tmp_378*(-tmp_362 - tmp_363 - tmp_366 - tmp_367 - tmp_368 - tmp_369 - tmp_372 - tmp_373 - tmp_374 + 1) - tmp_397*(-tmp_381 - tmp_382 - tmp_385 - tmp_386 - tmp_387 - tmp_388 - tmp_391 - tmp_392 - tmp_393 + 1) - tmp_416*(-tmp_400 - tmp_401 - tmp_404 - tmp_405 - tmp_406 - tmp_407 - tmp_410 - tmp_411 - tmp_412 + 1) - tmp_435*(-tmp_419 - tmp_420 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_429 - tmp_430 - tmp_431 + 1) - tmp_454*(-tmp_438 - tmp_439 - tmp_442 - tmp_443 - tmp_444 - tmp_445 - tmp_448 - tmp_449 - tmp_450 + 1) - tmp_473*(-tmp_457 - tmp_458 - tmp_461 - tmp_462 - tmp_463 - tmp_464 - tmp_467 - tmp_468 - tmp_469 + 1) - tmp_93*(-tmp_26 - tmp_28 - tmp_36 - tmp_38 - tmp_40 - tmp_42 - tmp_50 - tmp_52 - tmp_54 + 1);
      real_t a_1_0 = -tmp_112*(tmp_102 + tmp_103 + tmp_108) - tmp_131*(tmp_121 + tmp_122 + tmp_127) - tmp_150*(tmp_140 + tmp_141 + tmp_146) - tmp_169*(tmp_159 + tmp_160 + tmp_165) - tmp_188*(tmp_178 + tmp_179 + tmp_184) - tmp_207*(tmp_197 + tmp_198 + tmp_203) - tmp_226*(tmp_216 + tmp_217 + tmp_222) - tmp_245*(tmp_235 + tmp_236 + tmp_241) - tmp_264*(tmp_254 + tmp_255 + tmp_260) - tmp_283*(tmp_273 + tmp_274 + tmp_279) - tmp_302*(tmp_292 + tmp_293 + tmp_298) - tmp_321*(tmp_311 + tmp_312 + tmp_317) - tmp_340*(tmp_330 + tmp_331 + tmp_336) - tmp_359*(tmp_349 + tmp_350 + tmp_355) - tmp_378*(tmp_368 + tmp_369 + tmp_374) - tmp_397*(tmp_387 + tmp_388 + tmp_393) - tmp_416*(tmp_406 + tmp_407 + tmp_412) - tmp_435*(tmp_425 + tmp_426 + tmp_431) - tmp_454*(tmp_444 + tmp_445 + tmp_450) - tmp_473*(tmp_463 + tmp_464 + tmp_469) - tmp_93*(tmp_40 + tmp_42 + tmp_54);
      real_t a_2_0 = -tmp_112*(tmp_101 + tmp_107 + tmp_97) - tmp_131*(tmp_116 + tmp_120 + tmp_126) - tmp_150*(tmp_135 + tmp_139 + tmp_145) - tmp_169*(tmp_154 + tmp_158 + tmp_164) - tmp_188*(tmp_173 + tmp_177 + tmp_183) - tmp_207*(tmp_192 + tmp_196 + tmp_202) - tmp_226*(tmp_211 + tmp_215 + tmp_221) - tmp_245*(tmp_230 + tmp_234 + tmp_240) - tmp_264*(tmp_249 + tmp_253 + tmp_259) - tmp_283*(tmp_268 + tmp_272 + tmp_278) - tmp_302*(tmp_287 + tmp_291 + tmp_297) - tmp_321*(tmp_306 + tmp_310 + tmp_316) - tmp_340*(tmp_325 + tmp_329 + tmp_335) - tmp_359*(tmp_344 + tmp_348 + tmp_354) - tmp_378*(tmp_363 + tmp_367 + tmp_373) - tmp_397*(tmp_382 + tmp_386 + tmp_392) - tmp_416*(tmp_401 + tmp_405 + tmp_411) - tmp_435*(tmp_420 + tmp_424 + tmp_430) - tmp_454*(tmp_439 + tmp_443 + tmp_449) - tmp_473*(tmp_458 + tmp_462 + tmp_468) - tmp_93*(tmp_28 + tmp_38 + tmp_52);
      real_t a_3_0 = -tmp_112*(tmp_100 + tmp_106 + tmp_96) - tmp_131*(tmp_115 + tmp_119 + tmp_125) - tmp_150*(tmp_134 + tmp_138 + tmp_144) - tmp_169*(tmp_153 + tmp_157 + tmp_163) - tmp_188*(tmp_172 + tmp_176 + tmp_182) - tmp_207*(tmp_191 + tmp_195 + tmp_201) - tmp_226*(tmp_210 + tmp_214 + tmp_220) - tmp_245*(tmp_229 + tmp_233 + tmp_239) - tmp_264*(tmp_248 + tmp_252 + tmp_258) - tmp_283*(tmp_267 + tmp_271 + tmp_277) - tmp_302*(tmp_286 + tmp_290 + tmp_296) - tmp_321*(tmp_305 + tmp_309 + tmp_315) - tmp_340*(tmp_324 + tmp_328 + tmp_334) - tmp_359*(tmp_343 + tmp_347 + tmp_353) - tmp_378*(tmp_362 + tmp_366 + tmp_372) - tmp_397*(tmp_381 + tmp_385 + tmp_391) - tmp_416*(tmp_400 + tmp_404 + tmp_410) - tmp_435*(tmp_419 + tmp_423 + tmp_429) - tmp_454*(tmp_438 + tmp_442 + tmp_448) - tmp_473*(tmp_457 + tmp_461 + tmp_467) - tmp_93*(tmp_26 + tmp_36 + tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 3.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_54 - 1.0/4.0) + tmp_3*(tmp_53 - 1.0/4.0) + tmp_5*(tmp_52 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_74 - 1.0/4.0) + tmp_3*(tmp_73 - 1.0/4.0) + tmp_5*(tmp_72 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_90 - 1.0/4.0) + tmp_3*(tmp_89 - 1.0/4.0) + tmp_5*(tmp_88 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_106 - 1.0/4.0) + tmp_3*(tmp_105 - 1.0/4.0) + tmp_5*(tmp_104 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_122 - 1.0/4.0) + tmp_3*(tmp_121 - 1.0/4.0) + tmp_5*(tmp_120 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_138 - 1.0/4.0) + tmp_3*(tmp_137 - 1.0/4.0) + tmp_5*(tmp_136 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_154 - 1.0/4.0) + tmp_3*(tmp_153 - 1.0/4.0) + tmp_5*(tmp_152 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_170 - 1.0/4.0) + tmp_3*(tmp_169 - 1.0/4.0) + tmp_5*(tmp_168 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_186 - 1.0/4.0) + tmp_3*(tmp_185 - 1.0/4.0) + tmp_5*(tmp_184 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_202 - 1.0/4.0) + tmp_3*(tmp_201 - 1.0/4.0) + tmp_5*(tmp_200 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_218 - 1.0/4.0) + tmp_3*(tmp_217 - 1.0/4.0) + tmp_5*(tmp_216 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_234 - 1.0/4.0) + tmp_3*(tmp_233 - 1.0/4.0) + tmp_5*(tmp_232 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_250 - 1.0/4.0) + tmp_3*(tmp_249 - 1.0/4.0) + tmp_5*(tmp_248 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_266 - 1.0/4.0) + tmp_3*(tmp_265 - 1.0/4.0) + tmp_5*(tmp_264 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_282 - 1.0/4.0) + tmp_3*(tmp_281 - 1.0/4.0) + tmp_5*(tmp_280 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_298 - 1.0/4.0) + tmp_3*(tmp_297 - 1.0/4.0) + tmp_5*(tmp_296 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_314 - 1.0/4.0) + tmp_3*(tmp_313 - 1.0/4.0) + tmp_5*(tmp_312 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_330 - 1.0/4.0) + tmp_3*(tmp_329 - 1.0/4.0) + tmp_5*(tmp_328 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_346 - 1.0/4.0) + tmp_3*(tmp_345 - 1.0/4.0) + tmp_5*(tmp_344 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_362 - 1.0/4.0) + tmp_3*(tmp_361 - 1.0/4.0) + tmp_5*(tmp_360 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_378 - 1.0/4.0) + tmp_3*(tmp_377 - 1.0/4.0) + tmp_5*(tmp_376 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_1_0 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_2_0 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_3_0 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }


};




class EGNIPGVectorLaplaceFormP1E_2 : public hyteg::dg::DGForm
{
 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
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

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_2;
      real_t tmp_9 = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_11;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_16 + tmp_4*tmp_9 - tmp_7*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_21 = tmp_18*(tmp_12 - tmp_17);
      real_t tmp_22 = tmp_10*tmp_21 + tmp_13*tmp_20 + tmp_19*tmp_9;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_10*tmp_25 + tmp_13*tmp_24 + tmp_23*tmp_9;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_10*tmp_29 + tmp_13*tmp_28 + tmp_27*tmp_9;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_1_0 = 0.1666666666666668*tmp_45;
      real_t a_2_0 = 0.1666666666666668*tmp_46;
      real_t a_3_0 = 0.1666666666666668*tmp_47;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 3.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_52 - 1.0/4.0) + tmp_14*(tmp_53 - 1.0/4.0) + tmp_8*(tmp_54 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_72 - 1.0/4.0) + tmp_14*(tmp_73 - 1.0/4.0) + tmp_8*(tmp_74 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_88 - 1.0/4.0) + tmp_14*(tmp_89 - 1.0/4.0) + tmp_8*(tmp_90 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_104 - 1.0/4.0) + tmp_14*(tmp_105 - 1.0/4.0) + tmp_8*(tmp_106 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_120 - 1.0/4.0) + tmp_14*(tmp_121 - 1.0/4.0) + tmp_8*(tmp_122 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_136 - 1.0/4.0) + tmp_14*(tmp_137 - 1.0/4.0) + tmp_8*(tmp_138 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_152 - 1.0/4.0) + tmp_14*(tmp_153 - 1.0/4.0) + tmp_8*(tmp_154 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_168 - 1.0/4.0) + tmp_14*(tmp_169 - 1.0/4.0) + tmp_8*(tmp_170 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_184 - 1.0/4.0) + tmp_14*(tmp_185 - 1.0/4.0) + tmp_8*(tmp_186 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_200 - 1.0/4.0) + tmp_14*(tmp_201 - 1.0/4.0) + tmp_8*(tmp_202 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_216 - 1.0/4.0) + tmp_14*(tmp_217 - 1.0/4.0) + tmp_8*(tmp_218 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_232 - 1.0/4.0) + tmp_14*(tmp_233 - 1.0/4.0) + tmp_8*(tmp_234 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_248 - 1.0/4.0) + tmp_14*(tmp_249 - 1.0/4.0) + tmp_8*(tmp_250 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_264 - 1.0/4.0) + tmp_14*(tmp_265 - 1.0/4.0) + tmp_8*(tmp_266 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_280 - 1.0/4.0) + tmp_14*(tmp_281 - 1.0/4.0) + tmp_8*(tmp_282 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_296 - 1.0/4.0) + tmp_14*(tmp_297 - 1.0/4.0) + tmp_8*(tmp_298 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_312 - 1.0/4.0) + tmp_14*(tmp_313 - 1.0/4.0) + tmp_8*(tmp_314 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_328 - 1.0/4.0) + tmp_14*(tmp_329 - 1.0/4.0) + tmp_8*(tmp_330 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_344 - 1.0/4.0) + tmp_14*(tmp_345 - 1.0/4.0) + tmp_8*(tmp_346 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_360 - 1.0/4.0) + tmp_14*(tmp_361 - 1.0/4.0) + tmp_8*(tmp_362 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_376 - 1.0/4.0) + tmp_14*(tmp_377 - 1.0/4.0) + tmp_8*(tmp_378 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_1_0 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_2_0 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_3_0 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = 0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22;
      real_t tmp_24 = p_affine_8_2 + tmp_9;
      real_t tmp_25 = tmp_19*(tmp_23 + tmp_24);
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_28 = tmp_25*tmp_27;
      real_t tmp_29 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_30 = -p_affine_8_1;
      real_t tmp_31 = p_affine_9_1 + tmp_30;
      real_t tmp_32 = p_affine_10_1 + tmp_30;
      real_t tmp_33 = 0.031405749086161582*tmp_31 + 0.93718850182767688*tmp_32;
      real_t tmp_34 = p_affine_8_1 + tmp_2;
      real_t tmp_35 = tmp_19*(tmp_33 + tmp_34);
      real_t tmp_36 = tmp_29*tmp_35;
      real_t tmp_37 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_40 = tmp_25*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_42 = tmp_35*tmp_41;
      real_t tmp_43 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_44 = -p_affine_8_0;
      real_t tmp_45 = p_affine_9_0 + tmp_44;
      real_t tmp_46 = p_affine_10_0 + tmp_44;
      real_t tmp_47 = 0.031405749086161582*tmp_45 + 0.93718850182767688*tmp_46;
      real_t tmp_48 = p_affine_8_0 + tmp_0;
      real_t tmp_49 = tmp_19*(tmp_47 + tmp_48);
      real_t tmp_50 = tmp_43*tmp_49;
      real_t tmp_51 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_52 = tmp_49*tmp_51;
      real_t tmp_53 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_54 = tmp_49*tmp_53;
      real_t tmp_55 = -p_affine_4_2;
      real_t tmp_56 = p_affine_5_2 + tmp_55;
      real_t tmp_57 = -p_affine_4_0;
      real_t tmp_58 = p_affine_6_0 + tmp_57;
      real_t tmp_59 = -p_affine_4_1;
      real_t tmp_60 = p_affine_7_1 + tmp_59;
      real_t tmp_61 = tmp_58*tmp_60;
      real_t tmp_62 = p_affine_7_0 + tmp_57;
      real_t tmp_63 = p_affine_6_1 + tmp_59;
      real_t tmp_64 = tmp_62*tmp_63;
      real_t tmp_65 = tmp_61 - tmp_64;
      real_t tmp_66 = p_affine_5_0 + tmp_57;
      real_t tmp_67 = p_affine_7_2 + tmp_55;
      real_t tmp_68 = tmp_63*tmp_67;
      real_t tmp_69 = p_affine_5_1 + tmp_59;
      real_t tmp_70 = p_affine_6_2 + tmp_55;
      real_t tmp_71 = tmp_62*tmp_70;
      real_t tmp_72 = tmp_60*tmp_70;
      real_t tmp_73 = tmp_58*tmp_67;
      real_t tmp_74 = 1.0 / (tmp_56*tmp_61 - tmp_56*tmp_64 + tmp_66*tmp_68 - tmp_66*tmp_72 + tmp_69*tmp_71 - tmp_69*tmp_73);
      real_t tmp_75 = p_affine_8_2 + tmp_55;
      real_t tmp_76 = tmp_74*(tmp_23 + tmp_75);
      real_t tmp_77 = tmp_71 - tmp_73;
      real_t tmp_78 = p_affine_8_1 + tmp_59;
      real_t tmp_79 = tmp_74*(tmp_33 + tmp_78);
      real_t tmp_80 = tmp_68 - tmp_72;
      real_t tmp_81 = p_affine_8_0 + tmp_57;
      real_t tmp_82 = tmp_74*(tmp_47 + tmp_81);
      real_t tmp_83 = -tmp_60*tmp_66 + tmp_62*tmp_69;
      real_t tmp_84 = -tmp_56*tmp_62 + tmp_66*tmp_67;
      real_t tmp_85 = tmp_56*tmp_60 - tmp_67*tmp_69;
      real_t tmp_86 = -tmp_58*tmp_69 + tmp_63*tmp_66;
      real_t tmp_87 = tmp_56*tmp_58 - tmp_66*tmp_70;
      real_t tmp_88 = -tmp_56*tmp_63 + tmp_69*tmp_70;
      real_t tmp_89 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_90 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_91 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_92 = 3.0*std::pow((std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)*std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)) + (std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)*std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)) + (std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)*std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)), 0.25);
      real_t tmp_93 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_65*tmp_76 + tmp_77*tmp_79 + tmp_80*tmp_82 - 1.0/4.0) + tmp_67*(tmp_76*tmp_86 + tmp_79*tmp_87 + tmp_82*tmp_88 - 1.0/4.0) + tmp_70*(tmp_76*tmp_83 + tmp_79*tmp_84 + tmp_82*tmp_85 - 1.0/4.0));
      real_t tmp_94 = 0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22;
      real_t tmp_95 = tmp_19*(tmp_24 + tmp_94);
      real_t tmp_96 = tmp_8*tmp_95;
      real_t tmp_97 = tmp_27*tmp_95;
      real_t tmp_98 = 0.19601935860219369*tmp_31 + 0.60796128279561268*tmp_32;
      real_t tmp_99 = tmp_19*(tmp_34 + tmp_98);
      real_t tmp_100 = tmp_29*tmp_99;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = tmp_39*tmp_95;
      real_t tmp_103 = tmp_41*tmp_99;
      real_t tmp_104 = 0.19601935860219369*tmp_45 + 0.60796128279561268*tmp_46;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_48);
      real_t tmp_106 = tmp_105*tmp_43;
      real_t tmp_107 = tmp_105*tmp_51;
      real_t tmp_108 = tmp_105*tmp_53;
      real_t tmp_109 = tmp_74*(tmp_75 + tmp_94);
      real_t tmp_110 = tmp_74*(tmp_78 + tmp_98);
      real_t tmp_111 = tmp_74*(tmp_104 + tmp_81);
      real_t tmp_112 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_109*tmp_65 + tmp_110*tmp_77 + tmp_111*tmp_80 - 1.0/4.0) + tmp_67*(tmp_109*tmp_86 + tmp_110*tmp_87 + tmp_111*tmp_88 - 1.0/4.0) + tmp_70*(tmp_109*tmp_83 + tmp_110*tmp_84 + tmp_111*tmp_85 - 1.0/4.0));
      real_t tmp_113 = 0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_114 = tmp_19*(tmp_113 + tmp_24);
      real_t tmp_115 = tmp_114*tmp_8;
      real_t tmp_116 = tmp_114*tmp_27;
      real_t tmp_117 = 0.37605877282253791*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_118 = tmp_19*(tmp_117 + tmp_34);
      real_t tmp_119 = tmp_118*tmp_29;
      real_t tmp_120 = tmp_118*tmp_37;
      real_t tmp_121 = tmp_114*tmp_39;
      real_t tmp_122 = tmp_118*tmp_41;
      real_t tmp_123 = 0.37605877282253791*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_124 = tmp_19*(tmp_123 + tmp_48);
      real_t tmp_125 = tmp_124*tmp_43;
      real_t tmp_126 = tmp_124*tmp_51;
      real_t tmp_127 = tmp_124*tmp_53;
      real_t tmp_128 = tmp_74*(tmp_113 + tmp_75);
      real_t tmp_129 = tmp_74*(tmp_117 + tmp_78);
      real_t tmp_130 = tmp_74*(tmp_123 + tmp_81);
      real_t tmp_131 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_128*tmp_65 + tmp_129*tmp_77 + tmp_130*tmp_80 - 1.0/4.0) + tmp_67*(tmp_128*tmp_86 + tmp_129*tmp_87 + tmp_130*tmp_88 - 1.0/4.0) + tmp_70*(tmp_128*tmp_83 + tmp_129*tmp_84 + tmp_130*tmp_85 - 1.0/4.0));
      real_t tmp_132 = 0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_133 = tmp_19*(tmp_132 + tmp_24);
      real_t tmp_134 = tmp_133*tmp_8;
      real_t tmp_135 = tmp_133*tmp_27;
      real_t tmp_136 = 0.78764240869137092*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_137 = tmp_19*(tmp_136 + tmp_34);
      real_t tmp_138 = tmp_137*tmp_29;
      real_t tmp_139 = tmp_137*tmp_37;
      real_t tmp_140 = tmp_133*tmp_39;
      real_t tmp_141 = tmp_137*tmp_41;
      real_t tmp_142 = 0.78764240869137092*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_143 = tmp_19*(tmp_142 + tmp_48);
      real_t tmp_144 = tmp_143*tmp_43;
      real_t tmp_145 = tmp_143*tmp_51;
      real_t tmp_146 = tmp_143*tmp_53;
      real_t tmp_147 = tmp_74*(tmp_132 + tmp_75);
      real_t tmp_148 = tmp_74*(tmp_136 + tmp_78);
      real_t tmp_149 = tmp_74*(tmp_142 + tmp_81);
      real_t tmp_150 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_147*tmp_65 + tmp_148*tmp_77 + tmp_149*tmp_80 - 1.0/4.0) + tmp_67*(tmp_147*tmp_86 + tmp_148*tmp_87 + tmp_149*tmp_88 - 1.0/4.0) + tmp_70*(tmp_147*tmp_83 + tmp_148*tmp_84 + tmp_149*tmp_85 - 1.0/4.0));
      real_t tmp_151 = 0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_152 = tmp_19*(tmp_151 + tmp_24);
      real_t tmp_153 = tmp_152*tmp_8;
      real_t tmp_154 = tmp_152*tmp_27;
      real_t tmp_155 = 0.58463275527740355*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_156 = tmp_19*(tmp_155 + tmp_34);
      real_t tmp_157 = tmp_156*tmp_29;
      real_t tmp_158 = tmp_156*tmp_37;
      real_t tmp_159 = tmp_152*tmp_39;
      real_t tmp_160 = tmp_156*tmp_41;
      real_t tmp_161 = 0.58463275527740355*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_162 = tmp_19*(tmp_161 + tmp_48);
      real_t tmp_163 = tmp_162*tmp_43;
      real_t tmp_164 = tmp_162*tmp_51;
      real_t tmp_165 = tmp_162*tmp_53;
      real_t tmp_166 = tmp_74*(tmp_151 + tmp_75);
      real_t tmp_167 = tmp_74*(tmp_155 + tmp_78);
      real_t tmp_168 = tmp_74*(tmp_161 + tmp_81);
      real_t tmp_169 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_166*tmp_65 + tmp_167*tmp_77 + tmp_168*tmp_80 - 1.0/4.0) + tmp_67*(tmp_166*tmp_86 + tmp_167*tmp_87 + tmp_168*tmp_88 - 1.0/4.0) + tmp_70*(tmp_166*tmp_83 + tmp_167*tmp_84 + tmp_168*tmp_85 - 1.0/4.0));
      real_t tmp_170 = 0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_24);
      real_t tmp_172 = tmp_171*tmp_8;
      real_t tmp_173 = tmp_171*tmp_27;
      real_t tmp_174 = 0.041227165399737475*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_34);
      real_t tmp_176 = tmp_175*tmp_29;
      real_t tmp_177 = tmp_175*tmp_37;
      real_t tmp_178 = tmp_171*tmp_39;
      real_t tmp_179 = tmp_175*tmp_41;
      real_t tmp_180 = 0.041227165399737475*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_181 = tmp_19*(tmp_180 + tmp_48);
      real_t tmp_182 = tmp_181*tmp_43;
      real_t tmp_183 = tmp_181*tmp_51;
      real_t tmp_184 = tmp_181*tmp_53;
      real_t tmp_185 = tmp_74*(tmp_170 + tmp_75);
      real_t tmp_186 = tmp_74*(tmp_174 + tmp_78);
      real_t tmp_187 = tmp_74*(tmp_180 + tmp_81);
      real_t tmp_188 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_185*tmp_65 + tmp_186*tmp_77 + tmp_187*tmp_80 - 1.0/4.0) + tmp_67*(tmp_185*tmp_86 + tmp_186*tmp_87 + tmp_187*tmp_88 - 1.0/4.0) + tmp_70*(tmp_185*tmp_83 + tmp_186*tmp_84 + tmp_187*tmp_85 - 1.0/4.0));
      real_t tmp_189 = 0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_190 = tmp_19*(tmp_189 + tmp_24);
      real_t tmp_191 = tmp_190*tmp_8;
      real_t tmp_192 = tmp_190*tmp_27;
      real_t tmp_193 = 0.039308471900058539*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_194 = tmp_19*(tmp_193 + tmp_34);
      real_t tmp_195 = tmp_194*tmp_29;
      real_t tmp_196 = tmp_194*tmp_37;
      real_t tmp_197 = tmp_190*tmp_39;
      real_t tmp_198 = tmp_194*tmp_41;
      real_t tmp_199 = 0.039308471900058539*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_200 = tmp_19*(tmp_199 + tmp_48);
      real_t tmp_201 = tmp_200*tmp_43;
      real_t tmp_202 = tmp_200*tmp_51;
      real_t tmp_203 = tmp_200*tmp_53;
      real_t tmp_204 = tmp_74*(tmp_189 + tmp_75);
      real_t tmp_205 = tmp_74*(tmp_193 + tmp_78);
      real_t tmp_206 = tmp_74*(tmp_199 + tmp_81);
      real_t tmp_207 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_204*tmp_65 + tmp_205*tmp_77 + tmp_206*tmp_80 - 1.0/4.0) + tmp_67*(tmp_204*tmp_86 + tmp_205*tmp_87 + tmp_206*tmp_88 - 1.0/4.0) + tmp_70*(tmp_204*tmp_83 + tmp_205*tmp_84 + tmp_206*tmp_85 - 1.0/4.0));
      real_t tmp_208 = 0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_209 = tmp_19*(tmp_208 + tmp_24);
      real_t tmp_210 = tmp_209*tmp_8;
      real_t tmp_211 = tmp_209*tmp_27;
      real_t tmp_212 = 0.78764240869137092*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_213 = tmp_19*(tmp_212 + tmp_34);
      real_t tmp_214 = tmp_213*tmp_29;
      real_t tmp_215 = tmp_213*tmp_37;
      real_t tmp_216 = tmp_209*tmp_39;
      real_t tmp_217 = tmp_213*tmp_41;
      real_t tmp_218 = 0.78764240869137092*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_48);
      real_t tmp_220 = tmp_219*tmp_43;
      real_t tmp_221 = tmp_219*tmp_51;
      real_t tmp_222 = tmp_219*tmp_53;
      real_t tmp_223 = tmp_74*(tmp_208 + tmp_75);
      real_t tmp_224 = tmp_74*(tmp_212 + tmp_78);
      real_t tmp_225 = tmp_74*(tmp_218 + tmp_81);
      real_t tmp_226 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_223*tmp_65 + tmp_224*tmp_77 + tmp_225*tmp_80 - 1.0/4.0) + tmp_67*(tmp_223*tmp_86 + tmp_224*tmp_87 + tmp_225*tmp_88 - 1.0/4.0) + tmp_70*(tmp_223*tmp_83 + tmp_224*tmp_84 + tmp_225*tmp_85 - 1.0/4.0));
      real_t tmp_227 = 0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_228 = tmp_19*(tmp_227 + tmp_24);
      real_t tmp_229 = tmp_228*tmp_8;
      real_t tmp_230 = tmp_228*tmp_27;
      real_t tmp_231 = 0.58463275527740355*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_232 = tmp_19*(tmp_231 + tmp_34);
      real_t tmp_233 = tmp_232*tmp_29;
      real_t tmp_234 = tmp_232*tmp_37;
      real_t tmp_235 = tmp_228*tmp_39;
      real_t tmp_236 = tmp_232*tmp_41;
      real_t tmp_237 = 0.58463275527740355*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_238 = tmp_19*(tmp_237 + tmp_48);
      real_t tmp_239 = tmp_238*tmp_43;
      real_t tmp_240 = tmp_238*tmp_51;
      real_t tmp_241 = tmp_238*tmp_53;
      real_t tmp_242 = tmp_74*(tmp_227 + tmp_75);
      real_t tmp_243 = tmp_74*(tmp_231 + tmp_78);
      real_t tmp_244 = tmp_74*(tmp_237 + tmp_81);
      real_t tmp_245 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_242*tmp_65 + tmp_243*tmp_77 + tmp_244*tmp_80 - 1.0/4.0) + tmp_67*(tmp_242*tmp_86 + tmp_243*tmp_87 + tmp_244*tmp_88 - 1.0/4.0) + tmp_70*(tmp_242*tmp_83 + tmp_243*tmp_84 + tmp_244*tmp_85 - 1.0/4.0));
      real_t tmp_246 = 0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_247 = tmp_19*(tmp_24 + tmp_246);
      real_t tmp_248 = tmp_247*tmp_8;
      real_t tmp_249 = tmp_247*tmp_27;
      real_t tmp_250 = 0.1711304259088916*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_251 = tmp_19*(tmp_250 + tmp_34);
      real_t tmp_252 = tmp_251*tmp_29;
      real_t tmp_253 = tmp_251*tmp_37;
      real_t tmp_254 = tmp_247*tmp_39;
      real_t tmp_255 = tmp_251*tmp_41;
      real_t tmp_256 = 0.1711304259088916*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_257 = tmp_19*(tmp_256 + tmp_48);
      real_t tmp_258 = tmp_257*tmp_43;
      real_t tmp_259 = tmp_257*tmp_51;
      real_t tmp_260 = tmp_257*tmp_53;
      real_t tmp_261 = tmp_74*(tmp_246 + tmp_75);
      real_t tmp_262 = tmp_74*(tmp_250 + tmp_78);
      real_t tmp_263 = tmp_74*(tmp_256 + tmp_81);
      real_t tmp_264 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_261*tmp_65 + tmp_262*tmp_77 + tmp_263*tmp_80 - 1.0/4.0) + tmp_67*(tmp_261*tmp_86 + tmp_262*tmp_87 + tmp_263*tmp_88 - 1.0/4.0) + tmp_70*(tmp_261*tmp_83 + tmp_262*tmp_84 + tmp_263*tmp_85 - 1.0/4.0));
      real_t tmp_265 = 0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_266 = tmp_19*(tmp_24 + tmp_265);
      real_t tmp_267 = tmp_266*tmp_8;
      real_t tmp_268 = tmp_266*tmp_27;
      real_t tmp_269 = 0.37605877282253791*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_270 = tmp_19*(tmp_269 + tmp_34);
      real_t tmp_271 = tmp_270*tmp_29;
      real_t tmp_272 = tmp_270*tmp_37;
      real_t tmp_273 = tmp_266*tmp_39;
      real_t tmp_274 = tmp_270*tmp_41;
      real_t tmp_275 = 0.37605877282253791*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_276 = tmp_19*(tmp_275 + tmp_48);
      real_t tmp_277 = tmp_276*tmp_43;
      real_t tmp_278 = tmp_276*tmp_51;
      real_t tmp_279 = tmp_276*tmp_53;
      real_t tmp_280 = tmp_74*(tmp_265 + tmp_75);
      real_t tmp_281 = tmp_74*(tmp_269 + tmp_78);
      real_t tmp_282 = tmp_74*(tmp_275 + tmp_81);
      real_t tmp_283 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_280*tmp_65 + tmp_281*tmp_77 + tmp_282*tmp_80 - 1.0/4.0) + tmp_67*(tmp_280*tmp_86 + tmp_281*tmp_87 + tmp_282*tmp_88 - 1.0/4.0) + tmp_70*(tmp_280*tmp_83 + tmp_281*tmp_84 + tmp_282*tmp_85 - 1.0/4.0));
      real_t tmp_284 = 0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_285 = tmp_19*(tmp_24 + tmp_284);
      real_t tmp_286 = tmp_285*tmp_8;
      real_t tmp_287 = tmp_27*tmp_285;
      real_t tmp_288 = 0.041227165399737475*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_289 = tmp_19*(tmp_288 + tmp_34);
      real_t tmp_290 = tmp_289*tmp_29;
      real_t tmp_291 = tmp_289*tmp_37;
      real_t tmp_292 = tmp_285*tmp_39;
      real_t tmp_293 = tmp_289*tmp_41;
      real_t tmp_294 = 0.041227165399737475*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_295 = tmp_19*(tmp_294 + tmp_48);
      real_t tmp_296 = tmp_295*tmp_43;
      real_t tmp_297 = tmp_295*tmp_51;
      real_t tmp_298 = tmp_295*tmp_53;
      real_t tmp_299 = tmp_74*(tmp_284 + tmp_75);
      real_t tmp_300 = tmp_74*(tmp_288 + tmp_78);
      real_t tmp_301 = tmp_74*(tmp_294 + tmp_81);
      real_t tmp_302 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_299*tmp_65 + tmp_300*tmp_77 + tmp_301*tmp_80 - 1.0/4.0) + tmp_67*(tmp_299*tmp_86 + tmp_300*tmp_87 + tmp_301*tmp_88 - 1.0/4.0) + tmp_70*(tmp_299*tmp_83 + tmp_300*tmp_84 + tmp_301*tmp_85 - 1.0/4.0));
      real_t tmp_303 = 0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22;
      real_t tmp_304 = tmp_19*(tmp_24 + tmp_303);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_27*tmp_304;
      real_t tmp_307 = 0.40446199974765351*tmp_31 + 0.19107600050469298*tmp_32;
      real_t tmp_308 = tmp_19*(tmp_307 + tmp_34);
      real_t tmp_309 = tmp_29*tmp_308;
      real_t tmp_310 = tmp_308*tmp_37;
      real_t tmp_311 = tmp_304*tmp_39;
      real_t tmp_312 = tmp_308*tmp_41;
      real_t tmp_313 = 0.40446199974765351*tmp_45 + 0.19107600050469298*tmp_46;
      real_t tmp_314 = tmp_19*(tmp_313 + tmp_48);
      real_t tmp_315 = tmp_314*tmp_43;
      real_t tmp_316 = tmp_314*tmp_51;
      real_t tmp_317 = tmp_314*tmp_53;
      real_t tmp_318 = tmp_74*(tmp_303 + tmp_75);
      real_t tmp_319 = tmp_74*(tmp_307 + tmp_78);
      real_t tmp_320 = tmp_74*(tmp_313 + tmp_81);
      real_t tmp_321 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_318*tmp_65 + tmp_319*tmp_77 + tmp_320*tmp_80 - 1.0/4.0) + tmp_67*(tmp_318*tmp_86 + tmp_319*tmp_87 + tmp_320*tmp_88 - 1.0/4.0) + tmp_70*(tmp_318*tmp_83 + tmp_319*tmp_84 + tmp_320*tmp_85 - 1.0/4.0));
      real_t tmp_322 = 0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_323 = tmp_19*(tmp_24 + tmp_322);
      real_t tmp_324 = tmp_323*tmp_8;
      real_t tmp_325 = tmp_27*tmp_323;
      real_t tmp_326 = 0.039308471900058539*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_327 = tmp_19*(tmp_326 + tmp_34);
      real_t tmp_328 = tmp_29*tmp_327;
      real_t tmp_329 = tmp_327*tmp_37;
      real_t tmp_330 = tmp_323*tmp_39;
      real_t tmp_331 = tmp_327*tmp_41;
      real_t tmp_332 = 0.039308471900058539*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_333 = tmp_19*(tmp_332 + tmp_48);
      real_t tmp_334 = tmp_333*tmp_43;
      real_t tmp_335 = tmp_333*tmp_51;
      real_t tmp_336 = tmp_333*tmp_53;
      real_t tmp_337 = tmp_74*(tmp_322 + tmp_75);
      real_t tmp_338 = tmp_74*(tmp_326 + tmp_78);
      real_t tmp_339 = tmp_74*(tmp_332 + tmp_81);
      real_t tmp_340 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_337*tmp_65 + tmp_338*tmp_77 + tmp_339*tmp_80 - 1.0/4.0) + tmp_67*(tmp_337*tmp_86 + tmp_338*tmp_87 + tmp_339*tmp_88 - 1.0/4.0) + tmp_70*(tmp_337*tmp_83 + tmp_338*tmp_84 + tmp_339*tmp_85 - 1.0/4.0));
      real_t tmp_341 = 0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_342 = tmp_19*(tmp_24 + tmp_341);
      real_t tmp_343 = tmp_342*tmp_8;
      real_t tmp_344 = tmp_27*tmp_342;
      real_t tmp_345 = 0.93718850182767688*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_346 = tmp_19*(tmp_34 + tmp_345);
      real_t tmp_347 = tmp_29*tmp_346;
      real_t tmp_348 = tmp_346*tmp_37;
      real_t tmp_349 = tmp_342*tmp_39;
      real_t tmp_350 = tmp_346*tmp_41;
      real_t tmp_351 = 0.93718850182767688*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_352 = tmp_19*(tmp_351 + tmp_48);
      real_t tmp_353 = tmp_352*tmp_43;
      real_t tmp_354 = tmp_352*tmp_51;
      real_t tmp_355 = tmp_352*tmp_53;
      real_t tmp_356 = tmp_74*(tmp_341 + tmp_75);
      real_t tmp_357 = tmp_74*(tmp_345 + tmp_78);
      real_t tmp_358 = tmp_74*(tmp_351 + tmp_81);
      real_t tmp_359 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_356*tmp_65 + tmp_357*tmp_77 + tmp_358*tmp_80 - 1.0/4.0) + tmp_67*(tmp_356*tmp_86 + tmp_357*tmp_87 + tmp_358*tmp_88 - 1.0/4.0) + tmp_70*(tmp_356*tmp_83 + tmp_357*tmp_84 + tmp_358*tmp_85 - 1.0/4.0));
      real_t tmp_360 = 0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_361 = tmp_19*(tmp_24 + tmp_360);
      real_t tmp_362 = tmp_361*tmp_8;
      real_t tmp_363 = tmp_27*tmp_361;
      real_t tmp_364 = 0.60796128279561268*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_365 = tmp_19*(tmp_34 + tmp_364);
      real_t tmp_366 = tmp_29*tmp_365;
      real_t tmp_367 = tmp_365*tmp_37;
      real_t tmp_368 = tmp_361*tmp_39;
      real_t tmp_369 = tmp_365*tmp_41;
      real_t tmp_370 = 0.60796128279561268*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_371 = tmp_19*(tmp_370 + tmp_48);
      real_t tmp_372 = tmp_371*tmp_43;
      real_t tmp_373 = tmp_371*tmp_51;
      real_t tmp_374 = tmp_371*tmp_53;
      real_t tmp_375 = tmp_74*(tmp_360 + tmp_75);
      real_t tmp_376 = tmp_74*(tmp_364 + tmp_78);
      real_t tmp_377 = tmp_74*(tmp_370 + tmp_81);
      real_t tmp_378 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_375*tmp_65 + tmp_376*tmp_77 + tmp_377*tmp_80 - 1.0/4.0) + tmp_67*(tmp_375*tmp_86 + tmp_376*tmp_87 + tmp_377*tmp_88 - 1.0/4.0) + tmp_70*(tmp_375*tmp_83 + tmp_376*tmp_84 + tmp_377*tmp_85 - 1.0/4.0));
      real_t tmp_379 = 0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_380 = tmp_19*(tmp_24 + tmp_379);
      real_t tmp_381 = tmp_380*tmp_8;
      real_t tmp_382 = tmp_27*tmp_380;
      real_t tmp_383 = 0.19107600050469298*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_384 = tmp_19*(tmp_34 + tmp_383);
      real_t tmp_385 = tmp_29*tmp_384;
      real_t tmp_386 = tmp_37*tmp_384;
      real_t tmp_387 = tmp_380*tmp_39;
      real_t tmp_388 = tmp_384*tmp_41;
      real_t tmp_389 = 0.19107600050469298*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_390 = tmp_19*(tmp_389 + tmp_48);
      real_t tmp_391 = tmp_390*tmp_43;
      real_t tmp_392 = tmp_390*tmp_51;
      real_t tmp_393 = tmp_390*tmp_53;
      real_t tmp_394 = tmp_74*(tmp_379 + tmp_75);
      real_t tmp_395 = tmp_74*(tmp_383 + tmp_78);
      real_t tmp_396 = tmp_74*(tmp_389 + tmp_81);
      real_t tmp_397 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_394*tmp_65 + tmp_395*tmp_77 + tmp_396*tmp_80 - 1.0/4.0) + tmp_67*(tmp_394*tmp_86 + tmp_395*tmp_87 + tmp_396*tmp_88 - 1.0/4.0) + tmp_70*(tmp_394*tmp_83 + tmp_395*tmp_84 + tmp_396*tmp_85 - 1.0/4.0));
      real_t tmp_398 = 0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_399 = tmp_19*(tmp_24 + tmp_398);
      real_t tmp_400 = tmp_399*tmp_8;
      real_t tmp_401 = tmp_27*tmp_399;
      real_t tmp_402 = 0.031405749086161582*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_403 = tmp_19*(tmp_34 + tmp_402);
      real_t tmp_404 = tmp_29*tmp_403;
      real_t tmp_405 = tmp_37*tmp_403;
      real_t tmp_406 = tmp_39*tmp_399;
      real_t tmp_407 = tmp_403*tmp_41;
      real_t tmp_408 = 0.031405749086161582*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_409 = tmp_19*(tmp_408 + tmp_48);
      real_t tmp_410 = tmp_409*tmp_43;
      real_t tmp_411 = tmp_409*tmp_51;
      real_t tmp_412 = tmp_409*tmp_53;
      real_t tmp_413 = tmp_74*(tmp_398 + tmp_75);
      real_t tmp_414 = tmp_74*(tmp_402 + tmp_78);
      real_t tmp_415 = tmp_74*(tmp_408 + tmp_81);
      real_t tmp_416 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_413*tmp_65 + tmp_414*tmp_77 + tmp_415*tmp_80 - 1.0/4.0) + tmp_67*(tmp_413*tmp_86 + tmp_414*tmp_87 + tmp_415*tmp_88 - 1.0/4.0) + tmp_70*(tmp_413*tmp_83 + tmp_414*tmp_84 + tmp_415*tmp_85 - 1.0/4.0));
      real_t tmp_417 = 0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_418 = tmp_19*(tmp_24 + tmp_417);
      real_t tmp_419 = tmp_418*tmp_8;
      real_t tmp_420 = tmp_27*tmp_418;
      real_t tmp_421 = 0.19601935860219369*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_422 = tmp_19*(tmp_34 + tmp_421);
      real_t tmp_423 = tmp_29*tmp_422;
      real_t tmp_424 = tmp_37*tmp_422;
      real_t tmp_425 = tmp_39*tmp_418;
      real_t tmp_426 = tmp_41*tmp_422;
      real_t tmp_427 = 0.19601935860219369*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_428 = tmp_19*(tmp_427 + tmp_48);
      real_t tmp_429 = tmp_428*tmp_43;
      real_t tmp_430 = tmp_428*tmp_51;
      real_t tmp_431 = tmp_428*tmp_53;
      real_t tmp_432 = tmp_74*(tmp_417 + tmp_75);
      real_t tmp_433 = tmp_74*(tmp_421 + tmp_78);
      real_t tmp_434 = tmp_74*(tmp_427 + tmp_81);
      real_t tmp_435 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_432*tmp_65 + tmp_433*tmp_77 + tmp_434*tmp_80 - 1.0/4.0) + tmp_67*(tmp_432*tmp_86 + tmp_433*tmp_87 + tmp_434*tmp_88 - 1.0/4.0) + tmp_70*(tmp_432*tmp_83 + tmp_433*tmp_84 + tmp_434*tmp_85 - 1.0/4.0));
      real_t tmp_436 = 0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_437 = tmp_19*(tmp_24 + tmp_436);
      real_t tmp_438 = tmp_437*tmp_8;
      real_t tmp_439 = tmp_27*tmp_437;
      real_t tmp_440 = 0.40446199974765351*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_441 = tmp_19*(tmp_34 + tmp_440);
      real_t tmp_442 = tmp_29*tmp_441;
      real_t tmp_443 = tmp_37*tmp_441;
      real_t tmp_444 = tmp_39*tmp_437;
      real_t tmp_445 = tmp_41*tmp_441;
      real_t tmp_446 = 0.40446199974765351*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_447 = tmp_19*(tmp_446 + tmp_48);
      real_t tmp_448 = tmp_43*tmp_447;
      real_t tmp_449 = tmp_447*tmp_51;
      real_t tmp_450 = tmp_447*tmp_53;
      real_t tmp_451 = tmp_74*(tmp_436 + tmp_75);
      real_t tmp_452 = tmp_74*(tmp_440 + tmp_78);
      real_t tmp_453 = tmp_74*(tmp_446 + tmp_81);
      real_t tmp_454 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_451*tmp_65 + tmp_452*tmp_77 + tmp_453*tmp_80 - 1.0/4.0) + tmp_67*(tmp_451*tmp_86 + tmp_452*tmp_87 + tmp_453*tmp_88 - 1.0/4.0) + tmp_70*(tmp_451*tmp_83 + tmp_452*tmp_84 + tmp_453*tmp_85 - 1.0/4.0));
      real_t tmp_455 = 0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_456 = tmp_19*(tmp_24 + tmp_455);
      real_t tmp_457 = tmp_456*tmp_8;
      real_t tmp_458 = tmp_27*tmp_456;
      real_t tmp_459 = 0.1711304259088916*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_460 = tmp_19*(tmp_34 + tmp_459);
      real_t tmp_461 = tmp_29*tmp_460;
      real_t tmp_462 = tmp_37*tmp_460;
      real_t tmp_463 = tmp_39*tmp_456;
      real_t tmp_464 = tmp_41*tmp_460;
      real_t tmp_465 = 0.1711304259088916*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_466 = tmp_19*(tmp_465 + tmp_48);
      real_t tmp_467 = tmp_43*tmp_466;
      real_t tmp_468 = tmp_466*tmp_51;
      real_t tmp_469 = tmp_466*tmp_53;
      real_t tmp_470 = tmp_74*(tmp_455 + tmp_75);
      real_t tmp_471 = tmp_74*(tmp_459 + tmp_78);
      real_t tmp_472 = tmp_74*(tmp_465 + tmp_81);
      real_t tmp_473 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_470*tmp_65 + tmp_471*tmp_77 + tmp_472*tmp_80 - 1.0/4.0) + tmp_67*(tmp_470*tmp_86 + tmp_471*tmp_87 + tmp_472*tmp_88 - 1.0/4.0) + tmp_70*(tmp_470*tmp_83 + tmp_471*tmp_84 + tmp_472*tmp_85 - 1.0/4.0));
      real_t a_0_0 = -tmp_112*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_106 - tmp_107 - tmp_108 - tmp_96 - tmp_97 + 1) - tmp_131*(-tmp_115 - tmp_116 - tmp_119 - tmp_120 - tmp_121 - tmp_122 - tmp_125 - tmp_126 - tmp_127 + 1) - tmp_150*(-tmp_134 - tmp_135 - tmp_138 - tmp_139 - tmp_140 - tmp_141 - tmp_144 - tmp_145 - tmp_146 + 1) - tmp_169*(-tmp_153 - tmp_154 - tmp_157 - tmp_158 - tmp_159 - tmp_160 - tmp_163 - tmp_164 - tmp_165 + 1) - tmp_188*(-tmp_172 - tmp_173 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_182 - tmp_183 - tmp_184 + 1) - tmp_207*(-tmp_191 - tmp_192 - tmp_195 - tmp_196 - tmp_197 - tmp_198 - tmp_201 - tmp_202 - tmp_203 + 1) - tmp_226*(-tmp_210 - tmp_211 - tmp_214 - tmp_215 - tmp_216 - tmp_217 - tmp_220 - tmp_221 - tmp_222 + 1) - tmp_245*(-tmp_229 - tmp_230 - tmp_233 - tmp_234 - tmp_235 - tmp_236 - tmp_239 - tmp_240 - tmp_241 + 1) - tmp_264*(-tmp_248 - tmp_249 - tmp_252 - tmp_253 - tmp_254 - tmp_255 - tmp_258 - tmp_259 - tmp_260 + 1) - tmp_283*(-tmp_267 - tmp_268 - tmp_271 - tmp_272 - tmp_273 - tmp_274 - tmp_277 - tmp_278 - tmp_279 + 1) - tmp_302*(-tmp_286 - tmp_287 - tmp_290 - tmp_291 - tmp_292 - tmp_293 - tmp_296 - tmp_297 - tmp_298 + 1) - tmp_321*(-tmp_305 - tmp_306 - tmp_309 - tmp_310 - tmp_311 - tmp_312 - tmp_315 - tmp_316 - tmp_317 + 1) - tmp_340*(-tmp_324 - tmp_325 - tmp_328 - tmp_329 - tmp_330 - tmp_331 - tmp_334 - tmp_335 - tmp_336 + 1) - tmp_359*(-tmp_343 - tmp_344 - tmp_347 - tmp_348 - tmp_349 - tmp_350 - tmp_353 - tmp_354 - tmp_355 + 1) - tmp_378*(-tmp_362 - tmp_363 - tmp_366 - tmp_367 - tmp_368 - tmp_369 - tmp_372 - tmp_373 - tmp_374 + 1) - tmp_397*(-tmp_381 - tmp_382 - tmp_385 - tmp_386 - tmp_387 - tmp_388 - tmp_391 - tmp_392 - tmp_393 + 1) - tmp_416*(-tmp_400 - tmp_401 - tmp_404 - tmp_405 - tmp_406 - tmp_407 - tmp_410 - tmp_411 - tmp_412 + 1) - tmp_435*(-tmp_419 - tmp_420 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_429 - tmp_430 - tmp_431 + 1) - tmp_454*(-tmp_438 - tmp_439 - tmp_442 - tmp_443 - tmp_444 - tmp_445 - tmp_448 - tmp_449 - tmp_450 + 1) - tmp_473*(-tmp_457 - tmp_458 - tmp_461 - tmp_462 - tmp_463 - tmp_464 - tmp_467 - tmp_468 - tmp_469 + 1) - tmp_93*(-tmp_26 - tmp_28 - tmp_36 - tmp_38 - tmp_40 - tmp_42 - tmp_50 - tmp_52 - tmp_54 + 1);
      real_t a_1_0 = -tmp_112*(tmp_102 + tmp_103 + tmp_108) - tmp_131*(tmp_121 + tmp_122 + tmp_127) - tmp_150*(tmp_140 + tmp_141 + tmp_146) - tmp_169*(tmp_159 + tmp_160 + tmp_165) - tmp_188*(tmp_178 + tmp_179 + tmp_184) - tmp_207*(tmp_197 + tmp_198 + tmp_203) - tmp_226*(tmp_216 + tmp_217 + tmp_222) - tmp_245*(tmp_235 + tmp_236 + tmp_241) - tmp_264*(tmp_254 + tmp_255 + tmp_260) - tmp_283*(tmp_273 + tmp_274 + tmp_279) - tmp_302*(tmp_292 + tmp_293 + tmp_298) - tmp_321*(tmp_311 + tmp_312 + tmp_317) - tmp_340*(tmp_330 + tmp_331 + tmp_336) - tmp_359*(tmp_349 + tmp_350 + tmp_355) - tmp_378*(tmp_368 + tmp_369 + tmp_374) - tmp_397*(tmp_387 + tmp_388 + tmp_393) - tmp_416*(tmp_406 + tmp_407 + tmp_412) - tmp_435*(tmp_425 + tmp_426 + tmp_431) - tmp_454*(tmp_444 + tmp_445 + tmp_450) - tmp_473*(tmp_463 + tmp_464 + tmp_469) - tmp_93*(tmp_40 + tmp_42 + tmp_54);
      real_t a_2_0 = -tmp_112*(tmp_101 + tmp_107 + tmp_97) - tmp_131*(tmp_116 + tmp_120 + tmp_126) - tmp_150*(tmp_135 + tmp_139 + tmp_145) - tmp_169*(tmp_154 + tmp_158 + tmp_164) - tmp_188*(tmp_173 + tmp_177 + tmp_183) - tmp_207*(tmp_192 + tmp_196 + tmp_202) - tmp_226*(tmp_211 + tmp_215 + tmp_221) - tmp_245*(tmp_230 + tmp_234 + tmp_240) - tmp_264*(tmp_249 + tmp_253 + tmp_259) - tmp_283*(tmp_268 + tmp_272 + tmp_278) - tmp_302*(tmp_287 + tmp_291 + tmp_297) - tmp_321*(tmp_306 + tmp_310 + tmp_316) - tmp_340*(tmp_325 + tmp_329 + tmp_335) - tmp_359*(tmp_344 + tmp_348 + tmp_354) - tmp_378*(tmp_363 + tmp_367 + tmp_373) - tmp_397*(tmp_382 + tmp_386 + tmp_392) - tmp_416*(tmp_401 + tmp_405 + tmp_411) - tmp_435*(tmp_420 + tmp_424 + tmp_430) - tmp_454*(tmp_439 + tmp_443 + tmp_449) - tmp_473*(tmp_458 + tmp_462 + tmp_468) - tmp_93*(tmp_28 + tmp_38 + tmp_52);
      real_t a_3_0 = -tmp_112*(tmp_100 + tmp_106 + tmp_96) - tmp_131*(tmp_115 + tmp_119 + tmp_125) - tmp_150*(tmp_134 + tmp_138 + tmp_144) - tmp_169*(tmp_153 + tmp_157 + tmp_163) - tmp_188*(tmp_172 + tmp_176 + tmp_182) - tmp_207*(tmp_191 + tmp_195 + tmp_201) - tmp_226*(tmp_210 + tmp_214 + tmp_220) - tmp_245*(tmp_229 + tmp_233 + tmp_239) - tmp_264*(tmp_248 + tmp_252 + tmp_258) - tmp_283*(tmp_267 + tmp_271 + tmp_277) - tmp_302*(tmp_286 + tmp_290 + tmp_296) - tmp_321*(tmp_305 + tmp_309 + tmp_315) - tmp_340*(tmp_324 + tmp_328 + tmp_334) - tmp_359*(tmp_343 + tmp_347 + tmp_353) - tmp_378*(tmp_362 + tmp_366 + tmp_372) - tmp_397*(tmp_381 + tmp_385 + tmp_391) - tmp_416*(tmp_400 + tmp_404 + tmp_410) - tmp_435*(tmp_419 + tmp_423 + tmp_429) - tmp_454*(tmp_438 + tmp_442 + tmp_448) - tmp_473*(tmp_457 + tmp_461 + tmp_467) - tmp_93*(tmp_26 + tmp_36 + tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 3.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_52 - 1.0/4.0) + tmp_14*(tmp_53 - 1.0/4.0) + tmp_8*(tmp_54 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_72 - 1.0/4.0) + tmp_14*(tmp_73 - 1.0/4.0) + tmp_8*(tmp_74 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_88 - 1.0/4.0) + tmp_14*(tmp_89 - 1.0/4.0) + tmp_8*(tmp_90 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_104 - 1.0/4.0) + tmp_14*(tmp_105 - 1.0/4.0) + tmp_8*(tmp_106 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_120 - 1.0/4.0) + tmp_14*(tmp_121 - 1.0/4.0) + tmp_8*(tmp_122 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_136 - 1.0/4.0) + tmp_14*(tmp_137 - 1.0/4.0) + tmp_8*(tmp_138 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_152 - 1.0/4.0) + tmp_14*(tmp_153 - 1.0/4.0) + tmp_8*(tmp_154 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_168 - 1.0/4.0) + tmp_14*(tmp_169 - 1.0/4.0) + tmp_8*(tmp_170 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_184 - 1.0/4.0) + tmp_14*(tmp_185 - 1.0/4.0) + tmp_8*(tmp_186 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_200 - 1.0/4.0) + tmp_14*(tmp_201 - 1.0/4.0) + tmp_8*(tmp_202 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_216 - 1.0/4.0) + tmp_14*(tmp_217 - 1.0/4.0) + tmp_8*(tmp_218 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_232 - 1.0/4.0) + tmp_14*(tmp_233 - 1.0/4.0) + tmp_8*(tmp_234 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_248 - 1.0/4.0) + tmp_14*(tmp_249 - 1.0/4.0) + tmp_8*(tmp_250 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_264 - 1.0/4.0) + tmp_14*(tmp_265 - 1.0/4.0) + tmp_8*(tmp_266 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_280 - 1.0/4.0) + tmp_14*(tmp_281 - 1.0/4.0) + tmp_8*(tmp_282 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_296 - 1.0/4.0) + tmp_14*(tmp_297 - 1.0/4.0) + tmp_8*(tmp_298 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_312 - 1.0/4.0) + tmp_14*(tmp_313 - 1.0/4.0) + tmp_8*(tmp_314 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_328 - 1.0/4.0) + tmp_14*(tmp_329 - 1.0/4.0) + tmp_8*(tmp_330 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_344 - 1.0/4.0) + tmp_14*(tmp_345 - 1.0/4.0) + tmp_8*(tmp_346 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_360 - 1.0/4.0) + tmp_14*(tmp_361 - 1.0/4.0) + tmp_8*(tmp_362 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_376 - 1.0/4.0) + tmp_14*(tmp_377 - 1.0/4.0) + tmp_8*(tmp_378 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_1_0 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_2_0 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_3_0 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }


};




class EGNIPGVectorLaplaceFormEP1_0 : public hyteg::dg::DGForm
{
 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_1_1 + tmp_2));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = tmp_1*tmp_8 + tmp_5*tmp_7;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_11*tmp_5 + tmp_4*tmp_6;
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

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_1_1 + tmp_2));
      real_t tmp_6 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7 = p_affine_6_1 + tmp_2;
      real_t tmp_8 = tmp_5*(0.046910077030668018*tmp_6 + tmp_7);
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5*(0.046910077030668018*tmp_12 + tmp_13);
      real_t tmp_15 = tmp_14*tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 0.35539032758428413*tmp_1*(tmp_18 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_19 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 0.71794300574904923*tmp_1*(tmp_27 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_28 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 0.8533333333333335*tmp_1*(tmp_36 - 1.0/3.0) + 0.8533333333333335*tmp_4*(tmp_37 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 0.71794300574904923*tmp_1*(tmp_45 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_46 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.35539032758428413*tmp_1*(tmp_54 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_55 - 1.0/3.0);
      real_t a_0_0 = tmp_20*(-tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1) + tmp_29*(-tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1) + tmp_38*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_47*(-tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1) + tmp_56*(-tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1);
      real_t a_0_1 = tmp_18*tmp_20 + tmp_27*tmp_29 + tmp_36*tmp_38 + tmp_45*tmp_47 + tmp_54*tmp_56;
      real_t a_0_2 = tmp_19*tmp_20 + tmp_28*tmp_29 + tmp_37*tmp_38 + tmp_46*tmp_47 + tmp_55*tmp_56;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_3_0;
      real_t tmp_1 = p_affine_4_0 + tmp_0;
      real_t tmp_2 = -p_affine_3_1;
      real_t tmp_3 = p_affine_5_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_4_1 + tmp_2)*(p_affine_5_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + 0.046910077030668018*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -p_affine_0_0;
      real_t tmp_18 = p_affine_1_0 + tmp_17;
      real_t tmp_19 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_20 = -p_affine_0_1;
      real_t tmp_21 = p_affine_2_1 + tmp_20;
      real_t tmp_22 = p_affine_2_0 + tmp_17;
      real_t tmp_23 = 1.0 / (tmp_18*tmp_21 - tmp_22*(p_affine_1_1 + tmp_20));
      real_t tmp_24 = tmp_23*(tmp_20 + tmp_6);
      real_t tmp_25 = tmp_23*(tmp_12 + tmp_17);
      real_t tmp_26 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_27 = 0.35539032758428413*tmp_18*(tmp_19*tmp_24 + tmp_21*tmp_25 - 1.0/3.0) + 0.35539032758428413*tmp_22*(tmp_18*tmp_24 + tmp_25*tmp_26 - 1.0/3.0);
      real_t tmp_28 = p_affine_6_1 + 0.23076534494715845*tmp_5;
      real_t tmp_29 = tmp_4*(tmp_2 + tmp_28);
      real_t tmp_30 = tmp_1*tmp_29;
      real_t tmp_31 = tmp_29*tmp_9;
      real_t tmp_32 = p_affine_6_0 + 0.23076534494715845*tmp_11;
      real_t tmp_33 = tmp_4*(tmp_0 + tmp_32);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_15*tmp_33;
      real_t tmp_36 = tmp_23*(tmp_20 + tmp_28);
      real_t tmp_37 = tmp_23*(tmp_17 + tmp_32);
      real_t tmp_38 = 0.71794300574904923*tmp_18*(tmp_19*tmp_36 + tmp_21*tmp_37 - 1.0/3.0) + 0.71794300574904923*tmp_22*(tmp_18*tmp_36 + tmp_26*tmp_37 - 1.0/3.0);
      real_t tmp_39 = p_affine_6_1 + 0.5*tmp_5;
      real_t tmp_40 = tmp_4*(tmp_2 + tmp_39);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = p_affine_6_0 + 0.5*tmp_11;
      real_t tmp_44 = tmp_4*(tmp_0 + tmp_43);
      real_t tmp_45 = tmp_3*tmp_44;
      real_t tmp_46 = tmp_15*tmp_44;
      real_t tmp_47 = tmp_23*(tmp_20 + tmp_39);
      real_t tmp_48 = tmp_23*(tmp_17 + tmp_43);
      real_t tmp_49 = 0.8533333333333335*tmp_18*(tmp_19*tmp_47 + tmp_21*tmp_48 - 1.0/3.0) + 0.8533333333333335*tmp_22*(tmp_18*tmp_47 + tmp_26*tmp_48 - 1.0/3.0);
      real_t tmp_50 = p_affine_6_1 + 0.7692346550528415*tmp_5;
      real_t tmp_51 = tmp_4*(tmp_2 + tmp_50);
      real_t tmp_52 = tmp_1*tmp_51;
      real_t tmp_53 = tmp_51*tmp_9;
      real_t tmp_54 = p_affine_6_0 + 0.7692346550528415*tmp_11;
      real_t tmp_55 = tmp_4*(tmp_0 + tmp_54);
      real_t tmp_56 = tmp_3*tmp_55;
      real_t tmp_57 = tmp_15*tmp_55;
      real_t tmp_58 = tmp_23*(tmp_20 + tmp_50);
      real_t tmp_59 = tmp_23*(tmp_17 + tmp_54);
      real_t tmp_60 = 0.71794300574904923*tmp_18*(tmp_19*tmp_58 + tmp_21*tmp_59 - 1.0/3.0) + 0.71794300574904923*tmp_22*(tmp_18*tmp_58 + tmp_26*tmp_59 - 1.0/3.0);
      real_t tmp_61 = p_affine_6_1 + 0.95308992296933193*tmp_5;
      real_t tmp_62 = tmp_4*(tmp_2 + tmp_61);
      real_t tmp_63 = tmp_1*tmp_62;
      real_t tmp_64 = tmp_62*tmp_9;
      real_t tmp_65 = p_affine_6_0 + 0.95308992296933193*tmp_11;
      real_t tmp_66 = tmp_4*(tmp_0 + tmp_65);
      real_t tmp_67 = tmp_3*tmp_66;
      real_t tmp_68 = tmp_15*tmp_66;
      real_t tmp_69 = tmp_23*(tmp_20 + tmp_61);
      real_t tmp_70 = tmp_23*(tmp_17 + tmp_65);
      real_t tmp_71 = 0.35539032758428413*tmp_18*(tmp_19*tmp_69 + tmp_21*tmp_70 - 1.0/3.0) + 0.35539032758428413*tmp_22*(tmp_18*tmp_69 + tmp_26*tmp_70 - 1.0/3.0);
      real_t a_0_0 = -tmp_27*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) - tmp_38*(-tmp_30 - tmp_31 - tmp_34 - tmp_35 + 1) - tmp_49*(-tmp_41 - tmp_42 - tmp_45 - tmp_46 + 1) - tmp_60*(-tmp_52 - tmp_53 - tmp_56 - tmp_57 + 1) - tmp_71*(-tmp_63 - tmp_64 - tmp_67 - tmp_68 + 1);
      real_t a_0_1 = -tmp_27*(tmp_10 + tmp_14) - tmp_38*(tmp_31 + tmp_34) - tmp_49*(tmp_42 + tmp_45) - tmp_60*(tmp_53 + tmp_56) - tmp_71*(tmp_64 + tmp_67);
      real_t a_0_2 = -tmp_27*(tmp_16 + tmp_8) - tmp_38*(tmp_30 + tmp_35) - tmp_49*(tmp_41 + tmp_46) - tmp_60*(tmp_52 + tmp_57) - tmp_71*(tmp_63 + tmp_68);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_2;
      real_t tmp_9 = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_11;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_16 + tmp_4*tmp_9 - tmp_7*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_21 = tmp_18*(tmp_12 - tmp_17);
      real_t tmp_22 = tmp_1*tmp_21 + tmp_14*tmp_19 + tmp_20*tmp_5;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_1*tmp_25 + tmp_14*tmp_23 + tmp_24*tmp_5;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_1*tmp_29 + tmp_14*tmp_27 + tmp_28*tmp_5;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_0_1 = 0.1666666666666668*tmp_45;
      real_t a_0_2 = 0.1666666666666668*tmp_46;
      real_t a_0_3 = 0.1666666666666668*tmp_47;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 3.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_52 - 1.0/4.0) + tmp_13*(tmp_54 - 1.0/4.0) + tmp_4*(tmp_53 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_72 - 1.0/4.0) + tmp_13*(tmp_74 - 1.0/4.0) + tmp_4*(tmp_73 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_88 - 1.0/4.0) + tmp_13*(tmp_90 - 1.0/4.0) + tmp_4*(tmp_89 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_104 - 1.0/4.0) + tmp_13*(tmp_106 - 1.0/4.0) + tmp_4*(tmp_105 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_120 - 1.0/4.0) + tmp_13*(tmp_122 - 1.0/4.0) + tmp_4*(tmp_121 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_136 - 1.0/4.0) + tmp_13*(tmp_138 - 1.0/4.0) + tmp_4*(tmp_137 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_152 - 1.0/4.0) + tmp_13*(tmp_154 - 1.0/4.0) + tmp_4*(tmp_153 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_168 - 1.0/4.0) + tmp_13*(tmp_170 - 1.0/4.0) + tmp_4*(tmp_169 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_184 - 1.0/4.0) + tmp_13*(tmp_186 - 1.0/4.0) + tmp_4*(tmp_185 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_200 - 1.0/4.0) + tmp_13*(tmp_202 - 1.0/4.0) + tmp_4*(tmp_201 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_216 - 1.0/4.0) + tmp_13*(tmp_218 - 1.0/4.0) + tmp_4*(tmp_217 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_232 - 1.0/4.0) + tmp_13*(tmp_234 - 1.0/4.0) + tmp_4*(tmp_233 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_248 - 1.0/4.0) + tmp_13*(tmp_250 - 1.0/4.0) + tmp_4*(tmp_249 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_264 - 1.0/4.0) + tmp_13*(tmp_266 - 1.0/4.0) + tmp_4*(tmp_265 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_280 - 1.0/4.0) + tmp_13*(tmp_282 - 1.0/4.0) + tmp_4*(tmp_281 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_296 - 1.0/4.0) + tmp_13*(tmp_298 - 1.0/4.0) + tmp_4*(tmp_297 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_312 - 1.0/4.0) + tmp_13*(tmp_314 - 1.0/4.0) + tmp_4*(tmp_313 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_328 - 1.0/4.0) + tmp_13*(tmp_330 - 1.0/4.0) + tmp_4*(tmp_329 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_344 - 1.0/4.0) + tmp_13*(tmp_346 - 1.0/4.0) + tmp_4*(tmp_345 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_360 - 1.0/4.0) + tmp_13*(tmp_362 - 1.0/4.0) + tmp_4*(tmp_361 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_376 - 1.0/4.0) + tmp_13*(tmp_378 - 1.0/4.0) + tmp_4*(tmp_377 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_0_1 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_0_2 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_0_3 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_4_0;
      real_t tmp_1 = p_affine_5_0 + tmp_0;
      real_t tmp_2 = -p_affine_4_1;
      real_t tmp_3 = p_affine_6_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_6_0 + tmp_0;
      real_t tmp_6 = p_affine_5_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_4_2;
      real_t tmp_10 = p_affine_7_2 + tmp_9;
      real_t tmp_11 = p_affine_7_1 + tmp_2;
      real_t tmp_12 = p_affine_5_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_7_0 + tmp_0;
      real_t tmp_15 = p_affine_6_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = 0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22;
      real_t tmp_24 = p_affine_8_2 + tmp_9;
      real_t tmp_25 = tmp_19*(tmp_23 + tmp_24);
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_28 = tmp_25*tmp_27;
      real_t tmp_29 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_30 = -p_affine_8_1;
      real_t tmp_31 = p_affine_9_1 + tmp_30;
      real_t tmp_32 = p_affine_10_1 + tmp_30;
      real_t tmp_33 = 0.031405749086161582*tmp_31 + 0.93718850182767688*tmp_32;
      real_t tmp_34 = p_affine_8_1 + tmp_2;
      real_t tmp_35 = tmp_19*(tmp_33 + tmp_34);
      real_t tmp_36 = tmp_29*tmp_35;
      real_t tmp_37 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_40 = tmp_25*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_42 = tmp_35*tmp_41;
      real_t tmp_43 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_44 = -p_affine_8_0;
      real_t tmp_45 = p_affine_9_0 + tmp_44;
      real_t tmp_46 = p_affine_10_0 + tmp_44;
      real_t tmp_47 = 0.031405749086161582*tmp_45 + 0.93718850182767688*tmp_46;
      real_t tmp_48 = p_affine_8_0 + tmp_0;
      real_t tmp_49 = tmp_19*(tmp_47 + tmp_48);
      real_t tmp_50 = tmp_43*tmp_49;
      real_t tmp_51 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_52 = tmp_49*tmp_51;
      real_t tmp_53 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_54 = tmp_49*tmp_53;
      real_t tmp_55 = -p_affine_0_0;
      real_t tmp_56 = p_affine_1_0 + tmp_55;
      real_t tmp_57 = p_affine_2_0 + tmp_55;
      real_t tmp_58 = -p_affine_0_1;
      real_t tmp_59 = p_affine_3_1 + tmp_58;
      real_t tmp_60 = tmp_57*tmp_59;
      real_t tmp_61 = p_affine_3_0 + tmp_55;
      real_t tmp_62 = p_affine_2_1 + tmp_58;
      real_t tmp_63 = tmp_61*tmp_62;
      real_t tmp_64 = tmp_60 - tmp_63;
      real_t tmp_65 = -p_affine_0_2;
      real_t tmp_66 = p_affine_3_2 + tmp_65;
      real_t tmp_67 = tmp_62*tmp_66;
      real_t tmp_68 = p_affine_1_2 + tmp_65;
      real_t tmp_69 = p_affine_1_1 + tmp_58;
      real_t tmp_70 = p_affine_2_2 + tmp_65;
      real_t tmp_71 = tmp_61*tmp_70;
      real_t tmp_72 = tmp_59*tmp_70;
      real_t tmp_73 = tmp_57*tmp_66;
      real_t tmp_74 = 1.0 / (tmp_56*tmp_67 - tmp_56*tmp_72 + tmp_60*tmp_68 - tmp_63*tmp_68 + tmp_69*tmp_71 - tmp_69*tmp_73);
      real_t tmp_75 = p_affine_8_2 + tmp_65;
      real_t tmp_76 = tmp_74*(tmp_23 + tmp_75);
      real_t tmp_77 = tmp_71 - tmp_73;
      real_t tmp_78 = p_affine_8_1 + tmp_58;
      real_t tmp_79 = tmp_74*(tmp_33 + tmp_78);
      real_t tmp_80 = tmp_67 - tmp_72;
      real_t tmp_81 = p_affine_8_0 + tmp_55;
      real_t tmp_82 = tmp_74*(tmp_47 + tmp_81);
      real_t tmp_83 = -tmp_56*tmp_59 + tmp_61*tmp_69;
      real_t tmp_84 = tmp_56*tmp_66 - tmp_61*tmp_68;
      real_t tmp_85 = tmp_59*tmp_68 - tmp_66*tmp_69;
      real_t tmp_86 = tmp_56*tmp_62 - tmp_57*tmp_69;
      real_t tmp_87 = -tmp_56*tmp_70 + tmp_57*tmp_68;
      real_t tmp_88 = -tmp_62*tmp_68 + tmp_69*tmp_70;
      real_t tmp_89 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_90 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_91 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_92 = 3.0*std::pow((std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)*std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)) + (std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)*std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)) + (std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)*std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)), 0.25);
      real_t tmp_93 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_64*tmp_76 + tmp_77*tmp_79 + tmp_80*tmp_82 - 1.0/4.0) + tmp_57*(tmp_76*tmp_83 + tmp_79*tmp_84 + tmp_82*tmp_85 - 1.0/4.0) + tmp_61*(tmp_76*tmp_86 + tmp_79*tmp_87 + tmp_82*tmp_88 - 1.0/4.0));
      real_t tmp_94 = 0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22;
      real_t tmp_95 = tmp_19*(tmp_24 + tmp_94);
      real_t tmp_96 = tmp_8*tmp_95;
      real_t tmp_97 = tmp_27*tmp_95;
      real_t tmp_98 = 0.19601935860219369*tmp_31 + 0.60796128279561268*tmp_32;
      real_t tmp_99 = tmp_19*(tmp_34 + tmp_98);
      real_t tmp_100 = tmp_29*tmp_99;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = tmp_39*tmp_95;
      real_t tmp_103 = tmp_41*tmp_99;
      real_t tmp_104 = 0.19601935860219369*tmp_45 + 0.60796128279561268*tmp_46;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_48);
      real_t tmp_106 = tmp_105*tmp_43;
      real_t tmp_107 = tmp_105*tmp_51;
      real_t tmp_108 = tmp_105*tmp_53;
      real_t tmp_109 = tmp_74*(tmp_75 + tmp_94);
      real_t tmp_110 = tmp_74*(tmp_78 + tmp_98);
      real_t tmp_111 = tmp_74*(tmp_104 + tmp_81);
      real_t tmp_112 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_109*tmp_64 + tmp_110*tmp_77 + tmp_111*tmp_80 - 1.0/4.0) + tmp_57*(tmp_109*tmp_83 + tmp_110*tmp_84 + tmp_111*tmp_85 - 1.0/4.0) + tmp_61*(tmp_109*tmp_86 + tmp_110*tmp_87 + tmp_111*tmp_88 - 1.0/4.0));
      real_t tmp_113 = 0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_114 = tmp_19*(tmp_113 + tmp_24);
      real_t tmp_115 = tmp_114*tmp_8;
      real_t tmp_116 = tmp_114*tmp_27;
      real_t tmp_117 = 0.37605877282253791*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_118 = tmp_19*(tmp_117 + tmp_34);
      real_t tmp_119 = tmp_118*tmp_29;
      real_t tmp_120 = tmp_118*tmp_37;
      real_t tmp_121 = tmp_114*tmp_39;
      real_t tmp_122 = tmp_118*tmp_41;
      real_t tmp_123 = 0.37605877282253791*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_124 = tmp_19*(tmp_123 + tmp_48);
      real_t tmp_125 = tmp_124*tmp_43;
      real_t tmp_126 = tmp_124*tmp_51;
      real_t tmp_127 = tmp_124*tmp_53;
      real_t tmp_128 = tmp_74*(tmp_113 + tmp_75);
      real_t tmp_129 = tmp_74*(tmp_117 + tmp_78);
      real_t tmp_130 = tmp_74*(tmp_123 + tmp_81);
      real_t tmp_131 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_128*tmp_64 + tmp_129*tmp_77 + tmp_130*tmp_80 - 1.0/4.0) + tmp_57*(tmp_128*tmp_83 + tmp_129*tmp_84 + tmp_130*tmp_85 - 1.0/4.0) + tmp_61*(tmp_128*tmp_86 + tmp_129*tmp_87 + tmp_130*tmp_88 - 1.0/4.0));
      real_t tmp_132 = 0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_133 = tmp_19*(tmp_132 + tmp_24);
      real_t tmp_134 = tmp_133*tmp_8;
      real_t tmp_135 = tmp_133*tmp_27;
      real_t tmp_136 = 0.78764240869137092*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_137 = tmp_19*(tmp_136 + tmp_34);
      real_t tmp_138 = tmp_137*tmp_29;
      real_t tmp_139 = tmp_137*tmp_37;
      real_t tmp_140 = tmp_133*tmp_39;
      real_t tmp_141 = tmp_137*tmp_41;
      real_t tmp_142 = 0.78764240869137092*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_143 = tmp_19*(tmp_142 + tmp_48);
      real_t tmp_144 = tmp_143*tmp_43;
      real_t tmp_145 = tmp_143*tmp_51;
      real_t tmp_146 = tmp_143*tmp_53;
      real_t tmp_147 = tmp_74*(tmp_132 + tmp_75);
      real_t tmp_148 = tmp_74*(tmp_136 + tmp_78);
      real_t tmp_149 = tmp_74*(tmp_142 + tmp_81);
      real_t tmp_150 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_147*tmp_64 + tmp_148*tmp_77 + tmp_149*tmp_80 - 1.0/4.0) + tmp_57*(tmp_147*tmp_83 + tmp_148*tmp_84 + tmp_149*tmp_85 - 1.0/4.0) + tmp_61*(tmp_147*tmp_86 + tmp_148*tmp_87 + tmp_149*tmp_88 - 1.0/4.0));
      real_t tmp_151 = 0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_152 = tmp_19*(tmp_151 + tmp_24);
      real_t tmp_153 = tmp_152*tmp_8;
      real_t tmp_154 = tmp_152*tmp_27;
      real_t tmp_155 = 0.58463275527740355*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_156 = tmp_19*(tmp_155 + tmp_34);
      real_t tmp_157 = tmp_156*tmp_29;
      real_t tmp_158 = tmp_156*tmp_37;
      real_t tmp_159 = tmp_152*tmp_39;
      real_t tmp_160 = tmp_156*tmp_41;
      real_t tmp_161 = 0.58463275527740355*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_162 = tmp_19*(tmp_161 + tmp_48);
      real_t tmp_163 = tmp_162*tmp_43;
      real_t tmp_164 = tmp_162*tmp_51;
      real_t tmp_165 = tmp_162*tmp_53;
      real_t tmp_166 = tmp_74*(tmp_151 + tmp_75);
      real_t tmp_167 = tmp_74*(tmp_155 + tmp_78);
      real_t tmp_168 = tmp_74*(tmp_161 + tmp_81);
      real_t tmp_169 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_166*tmp_64 + tmp_167*tmp_77 + tmp_168*tmp_80 - 1.0/4.0) + tmp_57*(tmp_166*tmp_83 + tmp_167*tmp_84 + tmp_168*tmp_85 - 1.0/4.0) + tmp_61*(tmp_166*tmp_86 + tmp_167*tmp_87 + tmp_168*tmp_88 - 1.0/4.0));
      real_t tmp_170 = 0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_24);
      real_t tmp_172 = tmp_171*tmp_8;
      real_t tmp_173 = tmp_171*tmp_27;
      real_t tmp_174 = 0.041227165399737475*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_34);
      real_t tmp_176 = tmp_175*tmp_29;
      real_t tmp_177 = tmp_175*tmp_37;
      real_t tmp_178 = tmp_171*tmp_39;
      real_t tmp_179 = tmp_175*tmp_41;
      real_t tmp_180 = 0.041227165399737475*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_181 = tmp_19*(tmp_180 + tmp_48);
      real_t tmp_182 = tmp_181*tmp_43;
      real_t tmp_183 = tmp_181*tmp_51;
      real_t tmp_184 = tmp_181*tmp_53;
      real_t tmp_185 = tmp_74*(tmp_170 + tmp_75);
      real_t tmp_186 = tmp_74*(tmp_174 + tmp_78);
      real_t tmp_187 = tmp_74*(tmp_180 + tmp_81);
      real_t tmp_188 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_185*tmp_64 + tmp_186*tmp_77 + tmp_187*tmp_80 - 1.0/4.0) + tmp_57*(tmp_185*tmp_83 + tmp_186*tmp_84 + tmp_187*tmp_85 - 1.0/4.0) + tmp_61*(tmp_185*tmp_86 + tmp_186*tmp_87 + tmp_187*tmp_88 - 1.0/4.0));
      real_t tmp_189 = 0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_190 = tmp_19*(tmp_189 + tmp_24);
      real_t tmp_191 = tmp_190*tmp_8;
      real_t tmp_192 = tmp_190*tmp_27;
      real_t tmp_193 = 0.039308471900058539*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_194 = tmp_19*(tmp_193 + tmp_34);
      real_t tmp_195 = tmp_194*tmp_29;
      real_t tmp_196 = tmp_194*tmp_37;
      real_t tmp_197 = tmp_190*tmp_39;
      real_t tmp_198 = tmp_194*tmp_41;
      real_t tmp_199 = 0.039308471900058539*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_200 = tmp_19*(tmp_199 + tmp_48);
      real_t tmp_201 = tmp_200*tmp_43;
      real_t tmp_202 = tmp_200*tmp_51;
      real_t tmp_203 = tmp_200*tmp_53;
      real_t tmp_204 = tmp_74*(tmp_189 + tmp_75);
      real_t tmp_205 = tmp_74*(tmp_193 + tmp_78);
      real_t tmp_206 = tmp_74*(tmp_199 + tmp_81);
      real_t tmp_207 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_204*tmp_64 + tmp_205*tmp_77 + tmp_206*tmp_80 - 1.0/4.0) + tmp_57*(tmp_204*tmp_83 + tmp_205*tmp_84 + tmp_206*tmp_85 - 1.0/4.0) + tmp_61*(tmp_204*tmp_86 + tmp_205*tmp_87 + tmp_206*tmp_88 - 1.0/4.0));
      real_t tmp_208 = 0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_209 = tmp_19*(tmp_208 + tmp_24);
      real_t tmp_210 = tmp_209*tmp_8;
      real_t tmp_211 = tmp_209*tmp_27;
      real_t tmp_212 = 0.78764240869137092*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_213 = tmp_19*(tmp_212 + tmp_34);
      real_t tmp_214 = tmp_213*tmp_29;
      real_t tmp_215 = tmp_213*tmp_37;
      real_t tmp_216 = tmp_209*tmp_39;
      real_t tmp_217 = tmp_213*tmp_41;
      real_t tmp_218 = 0.78764240869137092*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_48);
      real_t tmp_220 = tmp_219*tmp_43;
      real_t tmp_221 = tmp_219*tmp_51;
      real_t tmp_222 = tmp_219*tmp_53;
      real_t tmp_223 = tmp_74*(tmp_208 + tmp_75);
      real_t tmp_224 = tmp_74*(tmp_212 + tmp_78);
      real_t tmp_225 = tmp_74*(tmp_218 + tmp_81);
      real_t tmp_226 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_223*tmp_64 + tmp_224*tmp_77 + tmp_225*tmp_80 - 1.0/4.0) + tmp_57*(tmp_223*tmp_83 + tmp_224*tmp_84 + tmp_225*tmp_85 - 1.0/4.0) + tmp_61*(tmp_223*tmp_86 + tmp_224*tmp_87 + tmp_225*tmp_88 - 1.0/4.0));
      real_t tmp_227 = 0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_228 = tmp_19*(tmp_227 + tmp_24);
      real_t tmp_229 = tmp_228*tmp_8;
      real_t tmp_230 = tmp_228*tmp_27;
      real_t tmp_231 = 0.58463275527740355*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_232 = tmp_19*(tmp_231 + tmp_34);
      real_t tmp_233 = tmp_232*tmp_29;
      real_t tmp_234 = tmp_232*tmp_37;
      real_t tmp_235 = tmp_228*tmp_39;
      real_t tmp_236 = tmp_232*tmp_41;
      real_t tmp_237 = 0.58463275527740355*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_238 = tmp_19*(tmp_237 + tmp_48);
      real_t tmp_239 = tmp_238*tmp_43;
      real_t tmp_240 = tmp_238*tmp_51;
      real_t tmp_241 = tmp_238*tmp_53;
      real_t tmp_242 = tmp_74*(tmp_227 + tmp_75);
      real_t tmp_243 = tmp_74*(tmp_231 + tmp_78);
      real_t tmp_244 = tmp_74*(tmp_237 + tmp_81);
      real_t tmp_245 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_242*tmp_64 + tmp_243*tmp_77 + tmp_244*tmp_80 - 1.0/4.0) + tmp_57*(tmp_242*tmp_83 + tmp_243*tmp_84 + tmp_244*tmp_85 - 1.0/4.0) + tmp_61*(tmp_242*tmp_86 + tmp_243*tmp_87 + tmp_244*tmp_88 - 1.0/4.0));
      real_t tmp_246 = 0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_247 = tmp_19*(tmp_24 + tmp_246);
      real_t tmp_248 = tmp_247*tmp_8;
      real_t tmp_249 = tmp_247*tmp_27;
      real_t tmp_250 = 0.1711304259088916*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_251 = tmp_19*(tmp_250 + tmp_34);
      real_t tmp_252 = tmp_251*tmp_29;
      real_t tmp_253 = tmp_251*tmp_37;
      real_t tmp_254 = tmp_247*tmp_39;
      real_t tmp_255 = tmp_251*tmp_41;
      real_t tmp_256 = 0.1711304259088916*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_257 = tmp_19*(tmp_256 + tmp_48);
      real_t tmp_258 = tmp_257*tmp_43;
      real_t tmp_259 = tmp_257*tmp_51;
      real_t tmp_260 = tmp_257*tmp_53;
      real_t tmp_261 = tmp_74*(tmp_246 + tmp_75);
      real_t tmp_262 = tmp_74*(tmp_250 + tmp_78);
      real_t tmp_263 = tmp_74*(tmp_256 + tmp_81);
      real_t tmp_264 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_261*tmp_64 + tmp_262*tmp_77 + tmp_263*tmp_80 - 1.0/4.0) + tmp_57*(tmp_261*tmp_83 + tmp_262*tmp_84 + tmp_263*tmp_85 - 1.0/4.0) + tmp_61*(tmp_261*tmp_86 + tmp_262*tmp_87 + tmp_263*tmp_88 - 1.0/4.0));
      real_t tmp_265 = 0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_266 = tmp_19*(tmp_24 + tmp_265);
      real_t tmp_267 = tmp_266*tmp_8;
      real_t tmp_268 = tmp_266*tmp_27;
      real_t tmp_269 = 0.37605877282253791*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_270 = tmp_19*(tmp_269 + tmp_34);
      real_t tmp_271 = tmp_270*tmp_29;
      real_t tmp_272 = tmp_270*tmp_37;
      real_t tmp_273 = tmp_266*tmp_39;
      real_t tmp_274 = tmp_270*tmp_41;
      real_t tmp_275 = 0.37605877282253791*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_276 = tmp_19*(tmp_275 + tmp_48);
      real_t tmp_277 = tmp_276*tmp_43;
      real_t tmp_278 = tmp_276*tmp_51;
      real_t tmp_279 = tmp_276*tmp_53;
      real_t tmp_280 = tmp_74*(tmp_265 + tmp_75);
      real_t tmp_281 = tmp_74*(tmp_269 + tmp_78);
      real_t tmp_282 = tmp_74*(tmp_275 + tmp_81);
      real_t tmp_283 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_280*tmp_64 + tmp_281*tmp_77 + tmp_282*tmp_80 - 1.0/4.0) + tmp_57*(tmp_280*tmp_83 + tmp_281*tmp_84 + tmp_282*tmp_85 - 1.0/4.0) + tmp_61*(tmp_280*tmp_86 + tmp_281*tmp_87 + tmp_282*tmp_88 - 1.0/4.0));
      real_t tmp_284 = 0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_285 = tmp_19*(tmp_24 + tmp_284);
      real_t tmp_286 = tmp_285*tmp_8;
      real_t tmp_287 = tmp_27*tmp_285;
      real_t tmp_288 = 0.041227165399737475*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_289 = tmp_19*(tmp_288 + tmp_34);
      real_t tmp_290 = tmp_289*tmp_29;
      real_t tmp_291 = tmp_289*tmp_37;
      real_t tmp_292 = tmp_285*tmp_39;
      real_t tmp_293 = tmp_289*tmp_41;
      real_t tmp_294 = 0.041227165399737475*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_295 = tmp_19*(tmp_294 + tmp_48);
      real_t tmp_296 = tmp_295*tmp_43;
      real_t tmp_297 = tmp_295*tmp_51;
      real_t tmp_298 = tmp_295*tmp_53;
      real_t tmp_299 = tmp_74*(tmp_284 + tmp_75);
      real_t tmp_300 = tmp_74*(tmp_288 + tmp_78);
      real_t tmp_301 = tmp_74*(tmp_294 + tmp_81);
      real_t tmp_302 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_299*tmp_64 + tmp_300*tmp_77 + tmp_301*tmp_80 - 1.0/4.0) + tmp_57*(tmp_299*tmp_83 + tmp_300*tmp_84 + tmp_301*tmp_85 - 1.0/4.0) + tmp_61*(tmp_299*tmp_86 + tmp_300*tmp_87 + tmp_301*tmp_88 - 1.0/4.0));
      real_t tmp_303 = 0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22;
      real_t tmp_304 = tmp_19*(tmp_24 + tmp_303);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_27*tmp_304;
      real_t tmp_307 = 0.40446199974765351*tmp_31 + 0.19107600050469298*tmp_32;
      real_t tmp_308 = tmp_19*(tmp_307 + tmp_34);
      real_t tmp_309 = tmp_29*tmp_308;
      real_t tmp_310 = tmp_308*tmp_37;
      real_t tmp_311 = tmp_304*tmp_39;
      real_t tmp_312 = tmp_308*tmp_41;
      real_t tmp_313 = 0.40446199974765351*tmp_45 + 0.19107600050469298*tmp_46;
      real_t tmp_314 = tmp_19*(tmp_313 + tmp_48);
      real_t tmp_315 = tmp_314*tmp_43;
      real_t tmp_316 = tmp_314*tmp_51;
      real_t tmp_317 = tmp_314*tmp_53;
      real_t tmp_318 = tmp_74*(tmp_303 + tmp_75);
      real_t tmp_319 = tmp_74*(tmp_307 + tmp_78);
      real_t tmp_320 = tmp_74*(tmp_313 + tmp_81);
      real_t tmp_321 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_318*tmp_64 + tmp_319*tmp_77 + tmp_320*tmp_80 - 1.0/4.0) + tmp_57*(tmp_318*tmp_83 + tmp_319*tmp_84 + tmp_320*tmp_85 - 1.0/4.0) + tmp_61*(tmp_318*tmp_86 + tmp_319*tmp_87 + tmp_320*tmp_88 - 1.0/4.0));
      real_t tmp_322 = 0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_323 = tmp_19*(tmp_24 + tmp_322);
      real_t tmp_324 = tmp_323*tmp_8;
      real_t tmp_325 = tmp_27*tmp_323;
      real_t tmp_326 = 0.039308471900058539*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_327 = tmp_19*(tmp_326 + tmp_34);
      real_t tmp_328 = tmp_29*tmp_327;
      real_t tmp_329 = tmp_327*tmp_37;
      real_t tmp_330 = tmp_323*tmp_39;
      real_t tmp_331 = tmp_327*tmp_41;
      real_t tmp_332 = 0.039308471900058539*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_333 = tmp_19*(tmp_332 + tmp_48);
      real_t tmp_334 = tmp_333*tmp_43;
      real_t tmp_335 = tmp_333*tmp_51;
      real_t tmp_336 = tmp_333*tmp_53;
      real_t tmp_337 = tmp_74*(tmp_322 + tmp_75);
      real_t tmp_338 = tmp_74*(tmp_326 + tmp_78);
      real_t tmp_339 = tmp_74*(tmp_332 + tmp_81);
      real_t tmp_340 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_337*tmp_64 + tmp_338*tmp_77 + tmp_339*tmp_80 - 1.0/4.0) + tmp_57*(tmp_337*tmp_83 + tmp_338*tmp_84 + tmp_339*tmp_85 - 1.0/4.0) + tmp_61*(tmp_337*tmp_86 + tmp_338*tmp_87 + tmp_339*tmp_88 - 1.0/4.0));
      real_t tmp_341 = 0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_342 = tmp_19*(tmp_24 + tmp_341);
      real_t tmp_343 = tmp_342*tmp_8;
      real_t tmp_344 = tmp_27*tmp_342;
      real_t tmp_345 = 0.93718850182767688*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_346 = tmp_19*(tmp_34 + tmp_345);
      real_t tmp_347 = tmp_29*tmp_346;
      real_t tmp_348 = tmp_346*tmp_37;
      real_t tmp_349 = tmp_342*tmp_39;
      real_t tmp_350 = tmp_346*tmp_41;
      real_t tmp_351 = 0.93718850182767688*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_352 = tmp_19*(tmp_351 + tmp_48);
      real_t tmp_353 = tmp_352*tmp_43;
      real_t tmp_354 = tmp_352*tmp_51;
      real_t tmp_355 = tmp_352*tmp_53;
      real_t tmp_356 = tmp_74*(tmp_341 + tmp_75);
      real_t tmp_357 = tmp_74*(tmp_345 + tmp_78);
      real_t tmp_358 = tmp_74*(tmp_351 + tmp_81);
      real_t tmp_359 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_356*tmp_64 + tmp_357*tmp_77 + tmp_358*tmp_80 - 1.0/4.0) + tmp_57*(tmp_356*tmp_83 + tmp_357*tmp_84 + tmp_358*tmp_85 - 1.0/4.0) + tmp_61*(tmp_356*tmp_86 + tmp_357*tmp_87 + tmp_358*tmp_88 - 1.0/4.0));
      real_t tmp_360 = 0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_361 = tmp_19*(tmp_24 + tmp_360);
      real_t tmp_362 = tmp_361*tmp_8;
      real_t tmp_363 = tmp_27*tmp_361;
      real_t tmp_364 = 0.60796128279561268*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_365 = tmp_19*(tmp_34 + tmp_364);
      real_t tmp_366 = tmp_29*tmp_365;
      real_t tmp_367 = tmp_365*tmp_37;
      real_t tmp_368 = tmp_361*tmp_39;
      real_t tmp_369 = tmp_365*tmp_41;
      real_t tmp_370 = 0.60796128279561268*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_371 = tmp_19*(tmp_370 + tmp_48);
      real_t tmp_372 = tmp_371*tmp_43;
      real_t tmp_373 = tmp_371*tmp_51;
      real_t tmp_374 = tmp_371*tmp_53;
      real_t tmp_375 = tmp_74*(tmp_360 + tmp_75);
      real_t tmp_376 = tmp_74*(tmp_364 + tmp_78);
      real_t tmp_377 = tmp_74*(tmp_370 + tmp_81);
      real_t tmp_378 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_375*tmp_64 + tmp_376*tmp_77 + tmp_377*tmp_80 - 1.0/4.0) + tmp_57*(tmp_375*tmp_83 + tmp_376*tmp_84 + tmp_377*tmp_85 - 1.0/4.0) + tmp_61*(tmp_375*tmp_86 + tmp_376*tmp_87 + tmp_377*tmp_88 - 1.0/4.0));
      real_t tmp_379 = 0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_380 = tmp_19*(tmp_24 + tmp_379);
      real_t tmp_381 = tmp_380*tmp_8;
      real_t tmp_382 = tmp_27*tmp_380;
      real_t tmp_383 = 0.19107600050469298*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_384 = tmp_19*(tmp_34 + tmp_383);
      real_t tmp_385 = tmp_29*tmp_384;
      real_t tmp_386 = tmp_37*tmp_384;
      real_t tmp_387 = tmp_380*tmp_39;
      real_t tmp_388 = tmp_384*tmp_41;
      real_t tmp_389 = 0.19107600050469298*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_390 = tmp_19*(tmp_389 + tmp_48);
      real_t tmp_391 = tmp_390*tmp_43;
      real_t tmp_392 = tmp_390*tmp_51;
      real_t tmp_393 = tmp_390*tmp_53;
      real_t tmp_394 = tmp_74*(tmp_379 + tmp_75);
      real_t tmp_395 = tmp_74*(tmp_383 + tmp_78);
      real_t tmp_396 = tmp_74*(tmp_389 + tmp_81);
      real_t tmp_397 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_394*tmp_64 + tmp_395*tmp_77 + tmp_396*tmp_80 - 1.0/4.0) + tmp_57*(tmp_394*tmp_83 + tmp_395*tmp_84 + tmp_396*tmp_85 - 1.0/4.0) + tmp_61*(tmp_394*tmp_86 + tmp_395*tmp_87 + tmp_396*tmp_88 - 1.0/4.0));
      real_t tmp_398 = 0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_399 = tmp_19*(tmp_24 + tmp_398);
      real_t tmp_400 = tmp_399*tmp_8;
      real_t tmp_401 = tmp_27*tmp_399;
      real_t tmp_402 = 0.031405749086161582*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_403 = tmp_19*(tmp_34 + tmp_402);
      real_t tmp_404 = tmp_29*tmp_403;
      real_t tmp_405 = tmp_37*tmp_403;
      real_t tmp_406 = tmp_39*tmp_399;
      real_t tmp_407 = tmp_403*tmp_41;
      real_t tmp_408 = 0.031405749086161582*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_409 = tmp_19*(tmp_408 + tmp_48);
      real_t tmp_410 = tmp_409*tmp_43;
      real_t tmp_411 = tmp_409*tmp_51;
      real_t tmp_412 = tmp_409*tmp_53;
      real_t tmp_413 = tmp_74*(tmp_398 + tmp_75);
      real_t tmp_414 = tmp_74*(tmp_402 + tmp_78);
      real_t tmp_415 = tmp_74*(tmp_408 + tmp_81);
      real_t tmp_416 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_413*tmp_64 + tmp_414*tmp_77 + tmp_415*tmp_80 - 1.0/4.0) + tmp_57*(tmp_413*tmp_83 + tmp_414*tmp_84 + tmp_415*tmp_85 - 1.0/4.0) + tmp_61*(tmp_413*tmp_86 + tmp_414*tmp_87 + tmp_415*tmp_88 - 1.0/4.0));
      real_t tmp_417 = 0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_418 = tmp_19*(tmp_24 + tmp_417);
      real_t tmp_419 = tmp_418*tmp_8;
      real_t tmp_420 = tmp_27*tmp_418;
      real_t tmp_421 = 0.19601935860219369*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_422 = tmp_19*(tmp_34 + tmp_421);
      real_t tmp_423 = tmp_29*tmp_422;
      real_t tmp_424 = tmp_37*tmp_422;
      real_t tmp_425 = tmp_39*tmp_418;
      real_t tmp_426 = tmp_41*tmp_422;
      real_t tmp_427 = 0.19601935860219369*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_428 = tmp_19*(tmp_427 + tmp_48);
      real_t tmp_429 = tmp_428*tmp_43;
      real_t tmp_430 = tmp_428*tmp_51;
      real_t tmp_431 = tmp_428*tmp_53;
      real_t tmp_432 = tmp_74*(tmp_417 + tmp_75);
      real_t tmp_433 = tmp_74*(tmp_421 + tmp_78);
      real_t tmp_434 = tmp_74*(tmp_427 + tmp_81);
      real_t tmp_435 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_432*tmp_64 + tmp_433*tmp_77 + tmp_434*tmp_80 - 1.0/4.0) + tmp_57*(tmp_432*tmp_83 + tmp_433*tmp_84 + tmp_434*tmp_85 - 1.0/4.0) + tmp_61*(tmp_432*tmp_86 + tmp_433*tmp_87 + tmp_434*tmp_88 - 1.0/4.0));
      real_t tmp_436 = 0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_437 = tmp_19*(tmp_24 + tmp_436);
      real_t tmp_438 = tmp_437*tmp_8;
      real_t tmp_439 = tmp_27*tmp_437;
      real_t tmp_440 = 0.40446199974765351*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_441 = tmp_19*(tmp_34 + tmp_440);
      real_t tmp_442 = tmp_29*tmp_441;
      real_t tmp_443 = tmp_37*tmp_441;
      real_t tmp_444 = tmp_39*tmp_437;
      real_t tmp_445 = tmp_41*tmp_441;
      real_t tmp_446 = 0.40446199974765351*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_447 = tmp_19*(tmp_446 + tmp_48);
      real_t tmp_448 = tmp_43*tmp_447;
      real_t tmp_449 = tmp_447*tmp_51;
      real_t tmp_450 = tmp_447*tmp_53;
      real_t tmp_451 = tmp_74*(tmp_436 + tmp_75);
      real_t tmp_452 = tmp_74*(tmp_440 + tmp_78);
      real_t tmp_453 = tmp_74*(tmp_446 + tmp_81);
      real_t tmp_454 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_451*tmp_64 + tmp_452*tmp_77 + tmp_453*tmp_80 - 1.0/4.0) + tmp_57*(tmp_451*tmp_83 + tmp_452*tmp_84 + tmp_453*tmp_85 - 1.0/4.0) + tmp_61*(tmp_451*tmp_86 + tmp_452*tmp_87 + tmp_453*tmp_88 - 1.0/4.0));
      real_t tmp_455 = 0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_456 = tmp_19*(tmp_24 + tmp_455);
      real_t tmp_457 = tmp_456*tmp_8;
      real_t tmp_458 = tmp_27*tmp_456;
      real_t tmp_459 = 0.1711304259088916*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_460 = tmp_19*(tmp_34 + tmp_459);
      real_t tmp_461 = tmp_29*tmp_460;
      real_t tmp_462 = tmp_37*tmp_460;
      real_t tmp_463 = tmp_39*tmp_456;
      real_t tmp_464 = tmp_41*tmp_460;
      real_t tmp_465 = 0.1711304259088916*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_466 = tmp_19*(tmp_465 + tmp_48);
      real_t tmp_467 = tmp_43*tmp_466;
      real_t tmp_468 = tmp_466*tmp_51;
      real_t tmp_469 = tmp_466*tmp_53;
      real_t tmp_470 = tmp_74*(tmp_455 + tmp_75);
      real_t tmp_471 = tmp_74*(tmp_459 + tmp_78);
      real_t tmp_472 = tmp_74*(tmp_465 + tmp_81);
      real_t tmp_473 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_470*tmp_64 + tmp_471*tmp_77 + tmp_472*tmp_80 - 1.0/4.0) + tmp_57*(tmp_470*tmp_83 + tmp_471*tmp_84 + tmp_472*tmp_85 - 1.0/4.0) + tmp_61*(tmp_470*tmp_86 + tmp_471*tmp_87 + tmp_472*tmp_88 - 1.0/4.0));
      real_t a_0_0 = -tmp_112*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_106 - tmp_107 - tmp_108 - tmp_96 - tmp_97 + 1) - tmp_131*(-tmp_115 - tmp_116 - tmp_119 - tmp_120 - tmp_121 - tmp_122 - tmp_125 - tmp_126 - tmp_127 + 1) - tmp_150*(-tmp_134 - tmp_135 - tmp_138 - tmp_139 - tmp_140 - tmp_141 - tmp_144 - tmp_145 - tmp_146 + 1) - tmp_169*(-tmp_153 - tmp_154 - tmp_157 - tmp_158 - tmp_159 - tmp_160 - tmp_163 - tmp_164 - tmp_165 + 1) - tmp_188*(-tmp_172 - tmp_173 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_182 - tmp_183 - tmp_184 + 1) - tmp_207*(-tmp_191 - tmp_192 - tmp_195 - tmp_196 - tmp_197 - tmp_198 - tmp_201 - tmp_202 - tmp_203 + 1) - tmp_226*(-tmp_210 - tmp_211 - tmp_214 - tmp_215 - tmp_216 - tmp_217 - tmp_220 - tmp_221 - tmp_222 + 1) - tmp_245*(-tmp_229 - tmp_230 - tmp_233 - tmp_234 - tmp_235 - tmp_236 - tmp_239 - tmp_240 - tmp_241 + 1) - tmp_264*(-tmp_248 - tmp_249 - tmp_252 - tmp_253 - tmp_254 - tmp_255 - tmp_258 - tmp_259 - tmp_260 + 1) - tmp_283*(-tmp_267 - tmp_268 - tmp_271 - tmp_272 - tmp_273 - tmp_274 - tmp_277 - tmp_278 - tmp_279 + 1) - tmp_302*(-tmp_286 - tmp_287 - tmp_290 - tmp_291 - tmp_292 - tmp_293 - tmp_296 - tmp_297 - tmp_298 + 1) - tmp_321*(-tmp_305 - tmp_306 - tmp_309 - tmp_310 - tmp_311 - tmp_312 - tmp_315 - tmp_316 - tmp_317 + 1) - tmp_340*(-tmp_324 - tmp_325 - tmp_328 - tmp_329 - tmp_330 - tmp_331 - tmp_334 - tmp_335 - tmp_336 + 1) - tmp_359*(-tmp_343 - tmp_344 - tmp_347 - tmp_348 - tmp_349 - tmp_350 - tmp_353 - tmp_354 - tmp_355 + 1) - tmp_378*(-tmp_362 - tmp_363 - tmp_366 - tmp_367 - tmp_368 - tmp_369 - tmp_372 - tmp_373 - tmp_374 + 1) - tmp_397*(-tmp_381 - tmp_382 - tmp_385 - tmp_386 - tmp_387 - tmp_388 - tmp_391 - tmp_392 - tmp_393 + 1) - tmp_416*(-tmp_400 - tmp_401 - tmp_404 - tmp_405 - tmp_406 - tmp_407 - tmp_410 - tmp_411 - tmp_412 + 1) - tmp_435*(-tmp_419 - tmp_420 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_429 - tmp_430 - tmp_431 + 1) - tmp_454*(-tmp_438 - tmp_439 - tmp_442 - tmp_443 - tmp_444 - tmp_445 - tmp_448 - tmp_449 - tmp_450 + 1) - tmp_473*(-tmp_457 - tmp_458 - tmp_461 - tmp_462 - tmp_463 - tmp_464 - tmp_467 - tmp_468 - tmp_469 + 1) - tmp_93*(-tmp_26 - tmp_28 - tmp_36 - tmp_38 - tmp_40 - tmp_42 - tmp_50 - tmp_52 - tmp_54 + 1);
      real_t a_0_1 = -tmp_112*(tmp_102 + tmp_103 + tmp_108) - tmp_131*(tmp_121 + tmp_122 + tmp_127) - tmp_150*(tmp_140 + tmp_141 + tmp_146) - tmp_169*(tmp_159 + tmp_160 + tmp_165) - tmp_188*(tmp_178 + tmp_179 + tmp_184) - tmp_207*(tmp_197 + tmp_198 + tmp_203) - tmp_226*(tmp_216 + tmp_217 + tmp_222) - tmp_245*(tmp_235 + tmp_236 + tmp_241) - tmp_264*(tmp_254 + tmp_255 + tmp_260) - tmp_283*(tmp_273 + tmp_274 + tmp_279) - tmp_302*(tmp_292 + tmp_293 + tmp_298) - tmp_321*(tmp_311 + tmp_312 + tmp_317) - tmp_340*(tmp_330 + tmp_331 + tmp_336) - tmp_359*(tmp_349 + tmp_350 + tmp_355) - tmp_378*(tmp_368 + tmp_369 + tmp_374) - tmp_397*(tmp_387 + tmp_388 + tmp_393) - tmp_416*(tmp_406 + tmp_407 + tmp_412) - tmp_435*(tmp_425 + tmp_426 + tmp_431) - tmp_454*(tmp_444 + tmp_445 + tmp_450) - tmp_473*(tmp_463 + tmp_464 + tmp_469) - tmp_93*(tmp_40 + tmp_42 + tmp_54);
      real_t a_0_2 = -tmp_112*(tmp_101 + tmp_107 + tmp_97) - tmp_131*(tmp_116 + tmp_120 + tmp_126) - tmp_150*(tmp_135 + tmp_139 + tmp_145) - tmp_169*(tmp_154 + tmp_158 + tmp_164) - tmp_188*(tmp_173 + tmp_177 + tmp_183) - tmp_207*(tmp_192 + tmp_196 + tmp_202) - tmp_226*(tmp_211 + tmp_215 + tmp_221) - tmp_245*(tmp_230 + tmp_234 + tmp_240) - tmp_264*(tmp_249 + tmp_253 + tmp_259) - tmp_283*(tmp_268 + tmp_272 + tmp_278) - tmp_302*(tmp_287 + tmp_291 + tmp_297) - tmp_321*(tmp_306 + tmp_310 + tmp_316) - tmp_340*(tmp_325 + tmp_329 + tmp_335) - tmp_359*(tmp_344 + tmp_348 + tmp_354) - tmp_378*(tmp_363 + tmp_367 + tmp_373) - tmp_397*(tmp_382 + tmp_386 + tmp_392) - tmp_416*(tmp_401 + tmp_405 + tmp_411) - tmp_435*(tmp_420 + tmp_424 + tmp_430) - tmp_454*(tmp_439 + tmp_443 + tmp_449) - tmp_473*(tmp_458 + tmp_462 + tmp_468) - tmp_93*(tmp_28 + tmp_38 + tmp_52);
      real_t a_0_3 = -tmp_112*(tmp_100 + tmp_106 + tmp_96) - tmp_131*(tmp_115 + tmp_119 + tmp_125) - tmp_150*(tmp_134 + tmp_138 + tmp_144) - tmp_169*(tmp_153 + tmp_157 + tmp_163) - tmp_188*(tmp_172 + tmp_176 + tmp_182) - tmp_207*(tmp_191 + tmp_195 + tmp_201) - tmp_226*(tmp_210 + tmp_214 + tmp_220) - tmp_245*(tmp_229 + tmp_233 + tmp_239) - tmp_264*(tmp_248 + tmp_252 + tmp_258) - tmp_283*(tmp_267 + tmp_271 + tmp_277) - tmp_302*(tmp_286 + tmp_290 + tmp_296) - tmp_321*(tmp_305 + tmp_309 + tmp_315) - tmp_340*(tmp_324 + tmp_328 + tmp_334) - tmp_359*(tmp_343 + tmp_347 + tmp_353) - tmp_378*(tmp_362 + tmp_366 + tmp_372) - tmp_397*(tmp_381 + tmp_385 + tmp_391) - tmp_416*(tmp_400 + tmp_404 + tmp_410) - tmp_435*(tmp_419 + tmp_423 + tmp_429) - tmp_454*(tmp_438 + tmp_442 + tmp_448) - tmp_473*(tmp_457 + tmp_461 + tmp_467) - tmp_93*(tmp_26 + tmp_36 + tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_0_3 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }


};




class EGNIPGVectorLaplaceFormEP1_1 : public hyteg::dg::DGForm
{
 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_2_0 + tmp_0));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = tmp_4*tmp_6 + tmp_5*tmp_8;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_10*tmp_5 + tmp_11*tmp_3;
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

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_1_1 + tmp_2;
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_2_0 + tmp_0));
      real_t tmp_6 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7 = p_affine_6_1 + tmp_2;
      real_t tmp_8 = tmp_5*(0.046910077030668018*tmp_6 + tmp_7);
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5*(0.046910077030668018*tmp_12 + tmp_13);
      real_t tmp_15 = tmp_14*tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 0.35539032758428413*tmp_3*(tmp_19 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_18 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 0.71794300574904923*tmp_3*(tmp_28 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_27 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 0.8533333333333335*tmp_3*(tmp_37 - 1.0/3.0) + 0.8533333333333335*tmp_4*(tmp_36 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 0.71794300574904923*tmp_3*(tmp_46 - 1.0/3.0) + 0.71794300574904923*tmp_4*(tmp_45 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.35539032758428413*tmp_3*(tmp_55 - 1.0/3.0) + 0.35539032758428413*tmp_4*(tmp_54 - 1.0/3.0);
      real_t a_0_0 = tmp_20*(-tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1) + tmp_29*(-tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1) + tmp_38*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_47*(-tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1) + tmp_56*(-tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1);
      real_t a_0_1 = tmp_18*tmp_20 + tmp_27*tmp_29 + tmp_36*tmp_38 + tmp_45*tmp_47 + tmp_54*tmp_56;
      real_t a_0_2 = tmp_19*tmp_20 + tmp_28*tmp_29 + tmp_37*tmp_38 + tmp_46*tmp_47 + tmp_55*tmp_56;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_3_0;
      real_t tmp_1 = p_affine_4_0 + tmp_0;
      real_t tmp_2 = -p_affine_3_1;
      real_t tmp_3 = p_affine_5_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_4_1 + tmp_2)*(p_affine_5_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + 0.046910077030668018*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -p_affine_0_1;
      real_t tmp_18 = p_affine_1_1 + tmp_17;
      real_t tmp_19 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_20 = -p_affine_0_0;
      real_t tmp_21 = p_affine_1_0 + tmp_20;
      real_t tmp_22 = p_affine_2_1 + tmp_17;
      real_t tmp_23 = 1.0 / (-tmp_18*(p_affine_2_0 + tmp_20) + tmp_21*tmp_22);
      real_t tmp_24 = tmp_23*(tmp_17 + tmp_6);
      real_t tmp_25 = tmp_23*(tmp_12 + tmp_20);
      real_t tmp_26 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_27 = 0.35539032758428413*tmp_18*(tmp_19*tmp_24 + tmp_22*tmp_25 - 1.0/3.0) + 0.35539032758428413*tmp_22*(tmp_21*tmp_24 + tmp_25*tmp_26 - 1.0/3.0);
      real_t tmp_28 = p_affine_6_1 + 0.23076534494715845*tmp_5;
      real_t tmp_29 = tmp_4*(tmp_2 + tmp_28);
      real_t tmp_30 = tmp_1*tmp_29;
      real_t tmp_31 = tmp_29*tmp_9;
      real_t tmp_32 = p_affine_6_0 + 0.23076534494715845*tmp_11;
      real_t tmp_33 = tmp_4*(tmp_0 + tmp_32);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_15*tmp_33;
      real_t tmp_36 = tmp_23*(tmp_17 + tmp_28);
      real_t tmp_37 = tmp_23*(tmp_20 + tmp_32);
      real_t tmp_38 = 0.71794300574904923*tmp_18*(tmp_19*tmp_36 + tmp_22*tmp_37 - 1.0/3.0) + 0.71794300574904923*tmp_22*(tmp_21*tmp_36 + tmp_26*tmp_37 - 1.0/3.0);
      real_t tmp_39 = p_affine_6_1 + 0.5*tmp_5;
      real_t tmp_40 = tmp_4*(tmp_2 + tmp_39);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = p_affine_6_0 + 0.5*tmp_11;
      real_t tmp_44 = tmp_4*(tmp_0 + tmp_43);
      real_t tmp_45 = tmp_3*tmp_44;
      real_t tmp_46 = tmp_15*tmp_44;
      real_t tmp_47 = tmp_23*(tmp_17 + tmp_39);
      real_t tmp_48 = tmp_23*(tmp_20 + tmp_43);
      real_t tmp_49 = 0.8533333333333335*tmp_18*(tmp_19*tmp_47 + tmp_22*tmp_48 - 1.0/3.0) + 0.8533333333333335*tmp_22*(tmp_21*tmp_47 + tmp_26*tmp_48 - 1.0/3.0);
      real_t tmp_50 = p_affine_6_1 + 0.7692346550528415*tmp_5;
      real_t tmp_51 = tmp_4*(tmp_2 + tmp_50);
      real_t tmp_52 = tmp_1*tmp_51;
      real_t tmp_53 = tmp_51*tmp_9;
      real_t tmp_54 = p_affine_6_0 + 0.7692346550528415*tmp_11;
      real_t tmp_55 = tmp_4*(tmp_0 + tmp_54);
      real_t tmp_56 = tmp_3*tmp_55;
      real_t tmp_57 = tmp_15*tmp_55;
      real_t tmp_58 = tmp_23*(tmp_17 + tmp_50);
      real_t tmp_59 = tmp_23*(tmp_20 + tmp_54);
      real_t tmp_60 = 0.71794300574904923*tmp_18*(tmp_19*tmp_58 + tmp_22*tmp_59 - 1.0/3.0) + 0.71794300574904923*tmp_22*(tmp_21*tmp_58 + tmp_26*tmp_59 - 1.0/3.0);
      real_t tmp_61 = p_affine_6_1 + 0.95308992296933193*tmp_5;
      real_t tmp_62 = tmp_4*(tmp_2 + tmp_61);
      real_t tmp_63 = tmp_1*tmp_62;
      real_t tmp_64 = tmp_62*tmp_9;
      real_t tmp_65 = p_affine_6_0 + 0.95308992296933193*tmp_11;
      real_t tmp_66 = tmp_4*(tmp_0 + tmp_65);
      real_t tmp_67 = tmp_3*tmp_66;
      real_t tmp_68 = tmp_15*tmp_66;
      real_t tmp_69 = tmp_23*(tmp_17 + tmp_61);
      real_t tmp_70 = tmp_23*(tmp_20 + tmp_65);
      real_t tmp_71 = 0.35539032758428413*tmp_18*(tmp_19*tmp_69 + tmp_22*tmp_70 - 1.0/3.0) + 0.35539032758428413*tmp_22*(tmp_21*tmp_69 + tmp_26*tmp_70 - 1.0/3.0);
      real_t a_0_0 = -tmp_27*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) - tmp_38*(-tmp_30 - tmp_31 - tmp_34 - tmp_35 + 1) - tmp_49*(-tmp_41 - tmp_42 - tmp_45 - tmp_46 + 1) - tmp_60*(-tmp_52 - tmp_53 - tmp_56 - tmp_57 + 1) - tmp_71*(-tmp_63 - tmp_64 - tmp_67 - tmp_68 + 1);
      real_t a_0_1 = -tmp_27*(tmp_10 + tmp_14) - tmp_38*(tmp_31 + tmp_34) - tmp_49*(tmp_42 + tmp_45) - tmp_60*(tmp_53 + tmp_56) - tmp_71*(tmp_64 + tmp_67);
      real_t a_0_2 = -tmp_27*(tmp_16 + tmp_8) - tmp_38*(tmp_30 + tmp_35) - tmp_49*(tmp_41 + tmp_46) - tmp_60*(tmp_52 + tmp_57) - tmp_71*(tmp_63 + tmp_68);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_2;
      real_t tmp_9 = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_11;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_16 + tmp_4*tmp_9 - tmp_7*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_21 = tmp_18*(tmp_12 - tmp_17);
      real_t tmp_22 = tmp_11*tmp_19 + tmp_20*tmp_3 + tmp_21*tmp_6;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_11*tmp_23 + tmp_24*tmp_3 + tmp_25*tmp_6;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_11*tmp_27 + tmp_28*tmp_3 + tmp_29*tmp_6;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_0_1 = 0.1666666666666668*tmp_45;
      real_t a_0_2 = 0.1666666666666668*tmp_46;
      real_t a_0_3 = 0.1666666666666668*tmp_47;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 3.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_54 - 1.0/4.0) + tmp_3*(tmp_53 - 1.0/4.0) + tmp_5*(tmp_52 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_74 - 1.0/4.0) + tmp_3*(tmp_73 - 1.0/4.0) + tmp_5*(tmp_72 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_90 - 1.0/4.0) + tmp_3*(tmp_89 - 1.0/4.0) + tmp_5*(tmp_88 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_106 - 1.0/4.0) + tmp_3*(tmp_105 - 1.0/4.0) + tmp_5*(tmp_104 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_122 - 1.0/4.0) + tmp_3*(tmp_121 - 1.0/4.0) + tmp_5*(tmp_120 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_138 - 1.0/4.0) + tmp_3*(tmp_137 - 1.0/4.0) + tmp_5*(tmp_136 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_154 - 1.0/4.0) + tmp_3*(tmp_153 - 1.0/4.0) + tmp_5*(tmp_152 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_170 - 1.0/4.0) + tmp_3*(tmp_169 - 1.0/4.0) + tmp_5*(tmp_168 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_186 - 1.0/4.0) + tmp_3*(tmp_185 - 1.0/4.0) + tmp_5*(tmp_184 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_202 - 1.0/4.0) + tmp_3*(tmp_201 - 1.0/4.0) + tmp_5*(tmp_200 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_218 - 1.0/4.0) + tmp_3*(tmp_217 - 1.0/4.0) + tmp_5*(tmp_216 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_234 - 1.0/4.0) + tmp_3*(tmp_233 - 1.0/4.0) + tmp_5*(tmp_232 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_250 - 1.0/4.0) + tmp_3*(tmp_249 - 1.0/4.0) + tmp_5*(tmp_248 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_266 - 1.0/4.0) + tmp_3*(tmp_265 - 1.0/4.0) + tmp_5*(tmp_264 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_282 - 1.0/4.0) + tmp_3*(tmp_281 - 1.0/4.0) + tmp_5*(tmp_280 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_298 - 1.0/4.0) + tmp_3*(tmp_297 - 1.0/4.0) + tmp_5*(tmp_296 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_314 - 1.0/4.0) + tmp_3*(tmp_313 - 1.0/4.0) + tmp_5*(tmp_312 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_330 - 1.0/4.0) + tmp_3*(tmp_329 - 1.0/4.0) + tmp_5*(tmp_328 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_346 - 1.0/4.0) + tmp_3*(tmp_345 - 1.0/4.0) + tmp_5*(tmp_344 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_362 - 1.0/4.0) + tmp_3*(tmp_361 - 1.0/4.0) + tmp_5*(tmp_360 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_378 - 1.0/4.0) + tmp_3*(tmp_377 - 1.0/4.0) + tmp_5*(tmp_376 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_0_1 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_0_2 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_0_3 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_4_0;
      real_t tmp_1 = p_affine_5_0 + tmp_0;
      real_t tmp_2 = -p_affine_4_1;
      real_t tmp_3 = p_affine_6_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_6_0 + tmp_0;
      real_t tmp_6 = p_affine_5_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_4_2;
      real_t tmp_10 = p_affine_7_2 + tmp_9;
      real_t tmp_11 = p_affine_7_1 + tmp_2;
      real_t tmp_12 = p_affine_5_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_7_0 + tmp_0;
      real_t tmp_15 = p_affine_6_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = 0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22;
      real_t tmp_24 = p_affine_8_2 + tmp_9;
      real_t tmp_25 = tmp_19*(tmp_23 + tmp_24);
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_28 = tmp_25*tmp_27;
      real_t tmp_29 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_30 = -p_affine_8_1;
      real_t tmp_31 = p_affine_9_1 + tmp_30;
      real_t tmp_32 = p_affine_10_1 + tmp_30;
      real_t tmp_33 = 0.031405749086161582*tmp_31 + 0.93718850182767688*tmp_32;
      real_t tmp_34 = p_affine_8_1 + tmp_2;
      real_t tmp_35 = tmp_19*(tmp_33 + tmp_34);
      real_t tmp_36 = tmp_29*tmp_35;
      real_t tmp_37 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_40 = tmp_25*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_42 = tmp_35*tmp_41;
      real_t tmp_43 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_44 = -p_affine_8_0;
      real_t tmp_45 = p_affine_9_0 + tmp_44;
      real_t tmp_46 = p_affine_10_0 + tmp_44;
      real_t tmp_47 = 0.031405749086161582*tmp_45 + 0.93718850182767688*tmp_46;
      real_t tmp_48 = p_affine_8_0 + tmp_0;
      real_t tmp_49 = tmp_19*(tmp_47 + tmp_48);
      real_t tmp_50 = tmp_43*tmp_49;
      real_t tmp_51 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_52 = tmp_49*tmp_51;
      real_t tmp_53 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_54 = tmp_49*tmp_53;
      real_t tmp_55 = -p_affine_0_1;
      real_t tmp_56 = p_affine_1_1 + tmp_55;
      real_t tmp_57 = -p_affine_0_0;
      real_t tmp_58 = p_affine_2_0 + tmp_57;
      real_t tmp_59 = p_affine_3_1 + tmp_55;
      real_t tmp_60 = tmp_58*tmp_59;
      real_t tmp_61 = p_affine_3_0 + tmp_57;
      real_t tmp_62 = p_affine_2_1 + tmp_55;
      real_t tmp_63 = tmp_61*tmp_62;
      real_t tmp_64 = tmp_60 - tmp_63;
      real_t tmp_65 = p_affine_1_0 + tmp_57;
      real_t tmp_66 = -p_affine_0_2;
      real_t tmp_67 = p_affine_3_2 + tmp_66;
      real_t tmp_68 = tmp_62*tmp_67;
      real_t tmp_69 = p_affine_1_2 + tmp_66;
      real_t tmp_70 = p_affine_2_2 + tmp_66;
      real_t tmp_71 = tmp_61*tmp_70;
      real_t tmp_72 = tmp_59*tmp_70;
      real_t tmp_73 = tmp_58*tmp_67;
      real_t tmp_74 = 1.0 / (tmp_56*tmp_71 - tmp_56*tmp_73 + tmp_60*tmp_69 - tmp_63*tmp_69 + tmp_65*tmp_68 - tmp_65*tmp_72);
      real_t tmp_75 = p_affine_8_2 + tmp_66;
      real_t tmp_76 = tmp_74*(tmp_23 + tmp_75);
      real_t tmp_77 = tmp_71 - tmp_73;
      real_t tmp_78 = p_affine_8_1 + tmp_55;
      real_t tmp_79 = tmp_74*(tmp_33 + tmp_78);
      real_t tmp_80 = tmp_68 - tmp_72;
      real_t tmp_81 = p_affine_8_0 + tmp_57;
      real_t tmp_82 = tmp_74*(tmp_47 + tmp_81);
      real_t tmp_83 = tmp_56*tmp_61 - tmp_59*tmp_65;
      real_t tmp_84 = -tmp_61*tmp_69 + tmp_65*tmp_67;
      real_t tmp_85 = -tmp_56*tmp_67 + tmp_59*tmp_69;
      real_t tmp_86 = -tmp_56*tmp_58 + tmp_62*tmp_65;
      real_t tmp_87 = tmp_58*tmp_69 - tmp_65*tmp_70;
      real_t tmp_88 = tmp_56*tmp_70 - tmp_62*tmp_69;
      real_t tmp_89 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_90 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_91 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_92 = 3.0*std::pow((std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)*std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)) + (std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)*std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)) + (std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)*std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)), 0.25);
      real_t tmp_93 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_64*tmp_76 + tmp_77*tmp_79 + tmp_80*tmp_82 - 1.0/4.0) + tmp_59*(tmp_76*tmp_86 + tmp_79*tmp_87 + tmp_82*tmp_88 - 1.0/4.0) + tmp_62*(tmp_76*tmp_83 + tmp_79*tmp_84 + tmp_82*tmp_85 - 1.0/4.0));
      real_t tmp_94 = 0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22;
      real_t tmp_95 = tmp_19*(tmp_24 + tmp_94);
      real_t tmp_96 = tmp_8*tmp_95;
      real_t tmp_97 = tmp_27*tmp_95;
      real_t tmp_98 = 0.19601935860219369*tmp_31 + 0.60796128279561268*tmp_32;
      real_t tmp_99 = tmp_19*(tmp_34 + tmp_98);
      real_t tmp_100 = tmp_29*tmp_99;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = tmp_39*tmp_95;
      real_t tmp_103 = tmp_41*tmp_99;
      real_t tmp_104 = 0.19601935860219369*tmp_45 + 0.60796128279561268*tmp_46;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_48);
      real_t tmp_106 = tmp_105*tmp_43;
      real_t tmp_107 = tmp_105*tmp_51;
      real_t tmp_108 = tmp_105*tmp_53;
      real_t tmp_109 = tmp_74*(tmp_75 + tmp_94);
      real_t tmp_110 = tmp_74*(tmp_78 + tmp_98);
      real_t tmp_111 = tmp_74*(tmp_104 + tmp_81);
      real_t tmp_112 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_109*tmp_64 + tmp_110*tmp_77 + tmp_111*tmp_80 - 1.0/4.0) + tmp_59*(tmp_109*tmp_86 + tmp_110*tmp_87 + tmp_111*tmp_88 - 1.0/4.0) + tmp_62*(tmp_109*tmp_83 + tmp_110*tmp_84 + tmp_111*tmp_85 - 1.0/4.0));
      real_t tmp_113 = 0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_114 = tmp_19*(tmp_113 + tmp_24);
      real_t tmp_115 = tmp_114*tmp_8;
      real_t tmp_116 = tmp_114*tmp_27;
      real_t tmp_117 = 0.37605877282253791*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_118 = tmp_19*(tmp_117 + tmp_34);
      real_t tmp_119 = tmp_118*tmp_29;
      real_t tmp_120 = tmp_118*tmp_37;
      real_t tmp_121 = tmp_114*tmp_39;
      real_t tmp_122 = tmp_118*tmp_41;
      real_t tmp_123 = 0.37605877282253791*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_124 = tmp_19*(tmp_123 + tmp_48);
      real_t tmp_125 = tmp_124*tmp_43;
      real_t tmp_126 = tmp_124*tmp_51;
      real_t tmp_127 = tmp_124*tmp_53;
      real_t tmp_128 = tmp_74*(tmp_113 + tmp_75);
      real_t tmp_129 = tmp_74*(tmp_117 + tmp_78);
      real_t tmp_130 = tmp_74*(tmp_123 + tmp_81);
      real_t tmp_131 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_128*tmp_64 + tmp_129*tmp_77 + tmp_130*tmp_80 - 1.0/4.0) + tmp_59*(tmp_128*tmp_86 + tmp_129*tmp_87 + tmp_130*tmp_88 - 1.0/4.0) + tmp_62*(tmp_128*tmp_83 + tmp_129*tmp_84 + tmp_130*tmp_85 - 1.0/4.0));
      real_t tmp_132 = 0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_133 = tmp_19*(tmp_132 + tmp_24);
      real_t tmp_134 = tmp_133*tmp_8;
      real_t tmp_135 = tmp_133*tmp_27;
      real_t tmp_136 = 0.78764240869137092*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_137 = tmp_19*(tmp_136 + tmp_34);
      real_t tmp_138 = tmp_137*tmp_29;
      real_t tmp_139 = tmp_137*tmp_37;
      real_t tmp_140 = tmp_133*tmp_39;
      real_t tmp_141 = tmp_137*tmp_41;
      real_t tmp_142 = 0.78764240869137092*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_143 = tmp_19*(tmp_142 + tmp_48);
      real_t tmp_144 = tmp_143*tmp_43;
      real_t tmp_145 = tmp_143*tmp_51;
      real_t tmp_146 = tmp_143*tmp_53;
      real_t tmp_147 = tmp_74*(tmp_132 + tmp_75);
      real_t tmp_148 = tmp_74*(tmp_136 + tmp_78);
      real_t tmp_149 = tmp_74*(tmp_142 + tmp_81);
      real_t tmp_150 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_147*tmp_64 + tmp_148*tmp_77 + tmp_149*tmp_80 - 1.0/4.0) + tmp_59*(tmp_147*tmp_86 + tmp_148*tmp_87 + tmp_149*tmp_88 - 1.0/4.0) + tmp_62*(tmp_147*tmp_83 + tmp_148*tmp_84 + tmp_149*tmp_85 - 1.0/4.0));
      real_t tmp_151 = 0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_152 = tmp_19*(tmp_151 + tmp_24);
      real_t tmp_153 = tmp_152*tmp_8;
      real_t tmp_154 = tmp_152*tmp_27;
      real_t tmp_155 = 0.58463275527740355*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_156 = tmp_19*(tmp_155 + tmp_34);
      real_t tmp_157 = tmp_156*tmp_29;
      real_t tmp_158 = tmp_156*tmp_37;
      real_t tmp_159 = tmp_152*tmp_39;
      real_t tmp_160 = tmp_156*tmp_41;
      real_t tmp_161 = 0.58463275527740355*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_162 = tmp_19*(tmp_161 + tmp_48);
      real_t tmp_163 = tmp_162*tmp_43;
      real_t tmp_164 = tmp_162*tmp_51;
      real_t tmp_165 = tmp_162*tmp_53;
      real_t tmp_166 = tmp_74*(tmp_151 + tmp_75);
      real_t tmp_167 = tmp_74*(tmp_155 + tmp_78);
      real_t tmp_168 = tmp_74*(tmp_161 + tmp_81);
      real_t tmp_169 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_166*tmp_64 + tmp_167*tmp_77 + tmp_168*tmp_80 - 1.0/4.0) + tmp_59*(tmp_166*tmp_86 + tmp_167*tmp_87 + tmp_168*tmp_88 - 1.0/4.0) + tmp_62*(tmp_166*tmp_83 + tmp_167*tmp_84 + tmp_168*tmp_85 - 1.0/4.0));
      real_t tmp_170 = 0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_24);
      real_t tmp_172 = tmp_171*tmp_8;
      real_t tmp_173 = tmp_171*tmp_27;
      real_t tmp_174 = 0.041227165399737475*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_34);
      real_t tmp_176 = tmp_175*tmp_29;
      real_t tmp_177 = tmp_175*tmp_37;
      real_t tmp_178 = tmp_171*tmp_39;
      real_t tmp_179 = tmp_175*tmp_41;
      real_t tmp_180 = 0.041227165399737475*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_181 = tmp_19*(tmp_180 + tmp_48);
      real_t tmp_182 = tmp_181*tmp_43;
      real_t tmp_183 = tmp_181*tmp_51;
      real_t tmp_184 = tmp_181*tmp_53;
      real_t tmp_185 = tmp_74*(tmp_170 + tmp_75);
      real_t tmp_186 = tmp_74*(tmp_174 + tmp_78);
      real_t tmp_187 = tmp_74*(tmp_180 + tmp_81);
      real_t tmp_188 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_185*tmp_64 + tmp_186*tmp_77 + tmp_187*tmp_80 - 1.0/4.0) + tmp_59*(tmp_185*tmp_86 + tmp_186*tmp_87 + tmp_187*tmp_88 - 1.0/4.0) + tmp_62*(tmp_185*tmp_83 + tmp_186*tmp_84 + tmp_187*tmp_85 - 1.0/4.0));
      real_t tmp_189 = 0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_190 = tmp_19*(tmp_189 + tmp_24);
      real_t tmp_191 = tmp_190*tmp_8;
      real_t tmp_192 = tmp_190*tmp_27;
      real_t tmp_193 = 0.039308471900058539*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_194 = tmp_19*(tmp_193 + tmp_34);
      real_t tmp_195 = tmp_194*tmp_29;
      real_t tmp_196 = tmp_194*tmp_37;
      real_t tmp_197 = tmp_190*tmp_39;
      real_t tmp_198 = tmp_194*tmp_41;
      real_t tmp_199 = 0.039308471900058539*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_200 = tmp_19*(tmp_199 + tmp_48);
      real_t tmp_201 = tmp_200*tmp_43;
      real_t tmp_202 = tmp_200*tmp_51;
      real_t tmp_203 = tmp_200*tmp_53;
      real_t tmp_204 = tmp_74*(tmp_189 + tmp_75);
      real_t tmp_205 = tmp_74*(tmp_193 + tmp_78);
      real_t tmp_206 = tmp_74*(tmp_199 + tmp_81);
      real_t tmp_207 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_204*tmp_64 + tmp_205*tmp_77 + tmp_206*tmp_80 - 1.0/4.0) + tmp_59*(tmp_204*tmp_86 + tmp_205*tmp_87 + tmp_206*tmp_88 - 1.0/4.0) + tmp_62*(tmp_204*tmp_83 + tmp_205*tmp_84 + tmp_206*tmp_85 - 1.0/4.0));
      real_t tmp_208 = 0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_209 = tmp_19*(tmp_208 + tmp_24);
      real_t tmp_210 = tmp_209*tmp_8;
      real_t tmp_211 = tmp_209*tmp_27;
      real_t tmp_212 = 0.78764240869137092*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_213 = tmp_19*(tmp_212 + tmp_34);
      real_t tmp_214 = tmp_213*tmp_29;
      real_t tmp_215 = tmp_213*tmp_37;
      real_t tmp_216 = tmp_209*tmp_39;
      real_t tmp_217 = tmp_213*tmp_41;
      real_t tmp_218 = 0.78764240869137092*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_48);
      real_t tmp_220 = tmp_219*tmp_43;
      real_t tmp_221 = tmp_219*tmp_51;
      real_t tmp_222 = tmp_219*tmp_53;
      real_t tmp_223 = tmp_74*(tmp_208 + tmp_75);
      real_t tmp_224 = tmp_74*(tmp_212 + tmp_78);
      real_t tmp_225 = tmp_74*(tmp_218 + tmp_81);
      real_t tmp_226 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_223*tmp_64 + tmp_224*tmp_77 + tmp_225*tmp_80 - 1.0/4.0) + tmp_59*(tmp_223*tmp_86 + tmp_224*tmp_87 + tmp_225*tmp_88 - 1.0/4.0) + tmp_62*(tmp_223*tmp_83 + tmp_224*tmp_84 + tmp_225*tmp_85 - 1.0/4.0));
      real_t tmp_227 = 0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_228 = tmp_19*(tmp_227 + tmp_24);
      real_t tmp_229 = tmp_228*tmp_8;
      real_t tmp_230 = tmp_228*tmp_27;
      real_t tmp_231 = 0.58463275527740355*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_232 = tmp_19*(tmp_231 + tmp_34);
      real_t tmp_233 = tmp_232*tmp_29;
      real_t tmp_234 = tmp_232*tmp_37;
      real_t tmp_235 = tmp_228*tmp_39;
      real_t tmp_236 = tmp_232*tmp_41;
      real_t tmp_237 = 0.58463275527740355*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_238 = tmp_19*(tmp_237 + tmp_48);
      real_t tmp_239 = tmp_238*tmp_43;
      real_t tmp_240 = tmp_238*tmp_51;
      real_t tmp_241 = tmp_238*tmp_53;
      real_t tmp_242 = tmp_74*(tmp_227 + tmp_75);
      real_t tmp_243 = tmp_74*(tmp_231 + tmp_78);
      real_t tmp_244 = tmp_74*(tmp_237 + tmp_81);
      real_t tmp_245 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_242*tmp_64 + tmp_243*tmp_77 + tmp_244*tmp_80 - 1.0/4.0) + tmp_59*(tmp_242*tmp_86 + tmp_243*tmp_87 + tmp_244*tmp_88 - 1.0/4.0) + tmp_62*(tmp_242*tmp_83 + tmp_243*tmp_84 + tmp_244*tmp_85 - 1.0/4.0));
      real_t tmp_246 = 0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_247 = tmp_19*(tmp_24 + tmp_246);
      real_t tmp_248 = tmp_247*tmp_8;
      real_t tmp_249 = tmp_247*tmp_27;
      real_t tmp_250 = 0.1711304259088916*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_251 = tmp_19*(tmp_250 + tmp_34);
      real_t tmp_252 = tmp_251*tmp_29;
      real_t tmp_253 = tmp_251*tmp_37;
      real_t tmp_254 = tmp_247*tmp_39;
      real_t tmp_255 = tmp_251*tmp_41;
      real_t tmp_256 = 0.1711304259088916*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_257 = tmp_19*(tmp_256 + tmp_48);
      real_t tmp_258 = tmp_257*tmp_43;
      real_t tmp_259 = tmp_257*tmp_51;
      real_t tmp_260 = tmp_257*tmp_53;
      real_t tmp_261 = tmp_74*(tmp_246 + tmp_75);
      real_t tmp_262 = tmp_74*(tmp_250 + tmp_78);
      real_t tmp_263 = tmp_74*(tmp_256 + tmp_81);
      real_t tmp_264 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_261*tmp_64 + tmp_262*tmp_77 + tmp_263*tmp_80 - 1.0/4.0) + tmp_59*(tmp_261*tmp_86 + tmp_262*tmp_87 + tmp_263*tmp_88 - 1.0/4.0) + tmp_62*(tmp_261*tmp_83 + tmp_262*tmp_84 + tmp_263*tmp_85 - 1.0/4.0));
      real_t tmp_265 = 0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_266 = tmp_19*(tmp_24 + tmp_265);
      real_t tmp_267 = tmp_266*tmp_8;
      real_t tmp_268 = tmp_266*tmp_27;
      real_t tmp_269 = 0.37605877282253791*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_270 = tmp_19*(tmp_269 + tmp_34);
      real_t tmp_271 = tmp_270*tmp_29;
      real_t tmp_272 = tmp_270*tmp_37;
      real_t tmp_273 = tmp_266*tmp_39;
      real_t tmp_274 = tmp_270*tmp_41;
      real_t tmp_275 = 0.37605877282253791*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_276 = tmp_19*(tmp_275 + tmp_48);
      real_t tmp_277 = tmp_276*tmp_43;
      real_t tmp_278 = tmp_276*tmp_51;
      real_t tmp_279 = tmp_276*tmp_53;
      real_t tmp_280 = tmp_74*(tmp_265 + tmp_75);
      real_t tmp_281 = tmp_74*(tmp_269 + tmp_78);
      real_t tmp_282 = tmp_74*(tmp_275 + tmp_81);
      real_t tmp_283 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_280*tmp_64 + tmp_281*tmp_77 + tmp_282*tmp_80 - 1.0/4.0) + tmp_59*(tmp_280*tmp_86 + tmp_281*tmp_87 + tmp_282*tmp_88 - 1.0/4.0) + tmp_62*(tmp_280*tmp_83 + tmp_281*tmp_84 + tmp_282*tmp_85 - 1.0/4.0));
      real_t tmp_284 = 0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_285 = tmp_19*(tmp_24 + tmp_284);
      real_t tmp_286 = tmp_285*tmp_8;
      real_t tmp_287 = tmp_27*tmp_285;
      real_t tmp_288 = 0.041227165399737475*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_289 = tmp_19*(tmp_288 + tmp_34);
      real_t tmp_290 = tmp_289*tmp_29;
      real_t tmp_291 = tmp_289*tmp_37;
      real_t tmp_292 = tmp_285*tmp_39;
      real_t tmp_293 = tmp_289*tmp_41;
      real_t tmp_294 = 0.041227165399737475*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_295 = tmp_19*(tmp_294 + tmp_48);
      real_t tmp_296 = tmp_295*tmp_43;
      real_t tmp_297 = tmp_295*tmp_51;
      real_t tmp_298 = tmp_295*tmp_53;
      real_t tmp_299 = tmp_74*(tmp_284 + tmp_75);
      real_t tmp_300 = tmp_74*(tmp_288 + tmp_78);
      real_t tmp_301 = tmp_74*(tmp_294 + tmp_81);
      real_t tmp_302 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_299*tmp_64 + tmp_300*tmp_77 + tmp_301*tmp_80 - 1.0/4.0) + tmp_59*(tmp_299*tmp_86 + tmp_300*tmp_87 + tmp_301*tmp_88 - 1.0/4.0) + tmp_62*(tmp_299*tmp_83 + tmp_300*tmp_84 + tmp_301*tmp_85 - 1.0/4.0));
      real_t tmp_303 = 0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22;
      real_t tmp_304 = tmp_19*(tmp_24 + tmp_303);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_27*tmp_304;
      real_t tmp_307 = 0.40446199974765351*tmp_31 + 0.19107600050469298*tmp_32;
      real_t tmp_308 = tmp_19*(tmp_307 + tmp_34);
      real_t tmp_309 = tmp_29*tmp_308;
      real_t tmp_310 = tmp_308*tmp_37;
      real_t tmp_311 = tmp_304*tmp_39;
      real_t tmp_312 = tmp_308*tmp_41;
      real_t tmp_313 = 0.40446199974765351*tmp_45 + 0.19107600050469298*tmp_46;
      real_t tmp_314 = tmp_19*(tmp_313 + tmp_48);
      real_t tmp_315 = tmp_314*tmp_43;
      real_t tmp_316 = tmp_314*tmp_51;
      real_t tmp_317 = tmp_314*tmp_53;
      real_t tmp_318 = tmp_74*(tmp_303 + tmp_75);
      real_t tmp_319 = tmp_74*(tmp_307 + tmp_78);
      real_t tmp_320 = tmp_74*(tmp_313 + tmp_81);
      real_t tmp_321 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_318*tmp_64 + tmp_319*tmp_77 + tmp_320*tmp_80 - 1.0/4.0) + tmp_59*(tmp_318*tmp_86 + tmp_319*tmp_87 + tmp_320*tmp_88 - 1.0/4.0) + tmp_62*(tmp_318*tmp_83 + tmp_319*tmp_84 + tmp_320*tmp_85 - 1.0/4.0));
      real_t tmp_322 = 0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_323 = tmp_19*(tmp_24 + tmp_322);
      real_t tmp_324 = tmp_323*tmp_8;
      real_t tmp_325 = tmp_27*tmp_323;
      real_t tmp_326 = 0.039308471900058539*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_327 = tmp_19*(tmp_326 + tmp_34);
      real_t tmp_328 = tmp_29*tmp_327;
      real_t tmp_329 = tmp_327*tmp_37;
      real_t tmp_330 = tmp_323*tmp_39;
      real_t tmp_331 = tmp_327*tmp_41;
      real_t tmp_332 = 0.039308471900058539*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_333 = tmp_19*(tmp_332 + tmp_48);
      real_t tmp_334 = tmp_333*tmp_43;
      real_t tmp_335 = tmp_333*tmp_51;
      real_t tmp_336 = tmp_333*tmp_53;
      real_t tmp_337 = tmp_74*(tmp_322 + tmp_75);
      real_t tmp_338 = tmp_74*(tmp_326 + tmp_78);
      real_t tmp_339 = tmp_74*(tmp_332 + tmp_81);
      real_t tmp_340 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_337*tmp_64 + tmp_338*tmp_77 + tmp_339*tmp_80 - 1.0/4.0) + tmp_59*(tmp_337*tmp_86 + tmp_338*tmp_87 + tmp_339*tmp_88 - 1.0/4.0) + tmp_62*(tmp_337*tmp_83 + tmp_338*tmp_84 + tmp_339*tmp_85 - 1.0/4.0));
      real_t tmp_341 = 0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_342 = tmp_19*(tmp_24 + tmp_341);
      real_t tmp_343 = tmp_342*tmp_8;
      real_t tmp_344 = tmp_27*tmp_342;
      real_t tmp_345 = 0.93718850182767688*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_346 = tmp_19*(tmp_34 + tmp_345);
      real_t tmp_347 = tmp_29*tmp_346;
      real_t tmp_348 = tmp_346*tmp_37;
      real_t tmp_349 = tmp_342*tmp_39;
      real_t tmp_350 = tmp_346*tmp_41;
      real_t tmp_351 = 0.93718850182767688*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_352 = tmp_19*(tmp_351 + tmp_48);
      real_t tmp_353 = tmp_352*tmp_43;
      real_t tmp_354 = tmp_352*tmp_51;
      real_t tmp_355 = tmp_352*tmp_53;
      real_t tmp_356 = tmp_74*(tmp_341 + tmp_75);
      real_t tmp_357 = tmp_74*(tmp_345 + tmp_78);
      real_t tmp_358 = tmp_74*(tmp_351 + tmp_81);
      real_t tmp_359 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_356*tmp_64 + tmp_357*tmp_77 + tmp_358*tmp_80 - 1.0/4.0) + tmp_59*(tmp_356*tmp_86 + tmp_357*tmp_87 + tmp_358*tmp_88 - 1.0/4.0) + tmp_62*(tmp_356*tmp_83 + tmp_357*tmp_84 + tmp_358*tmp_85 - 1.0/4.0));
      real_t tmp_360 = 0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_361 = tmp_19*(tmp_24 + tmp_360);
      real_t tmp_362 = tmp_361*tmp_8;
      real_t tmp_363 = tmp_27*tmp_361;
      real_t tmp_364 = 0.60796128279561268*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_365 = tmp_19*(tmp_34 + tmp_364);
      real_t tmp_366 = tmp_29*tmp_365;
      real_t tmp_367 = tmp_365*tmp_37;
      real_t tmp_368 = tmp_361*tmp_39;
      real_t tmp_369 = tmp_365*tmp_41;
      real_t tmp_370 = 0.60796128279561268*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_371 = tmp_19*(tmp_370 + tmp_48);
      real_t tmp_372 = tmp_371*tmp_43;
      real_t tmp_373 = tmp_371*tmp_51;
      real_t tmp_374 = tmp_371*tmp_53;
      real_t tmp_375 = tmp_74*(tmp_360 + tmp_75);
      real_t tmp_376 = tmp_74*(tmp_364 + tmp_78);
      real_t tmp_377 = tmp_74*(tmp_370 + tmp_81);
      real_t tmp_378 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_375*tmp_64 + tmp_376*tmp_77 + tmp_377*tmp_80 - 1.0/4.0) + tmp_59*(tmp_375*tmp_86 + tmp_376*tmp_87 + tmp_377*tmp_88 - 1.0/4.0) + tmp_62*(tmp_375*tmp_83 + tmp_376*tmp_84 + tmp_377*tmp_85 - 1.0/4.0));
      real_t tmp_379 = 0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_380 = tmp_19*(tmp_24 + tmp_379);
      real_t tmp_381 = tmp_380*tmp_8;
      real_t tmp_382 = tmp_27*tmp_380;
      real_t tmp_383 = 0.19107600050469298*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_384 = tmp_19*(tmp_34 + tmp_383);
      real_t tmp_385 = tmp_29*tmp_384;
      real_t tmp_386 = tmp_37*tmp_384;
      real_t tmp_387 = tmp_380*tmp_39;
      real_t tmp_388 = tmp_384*tmp_41;
      real_t tmp_389 = 0.19107600050469298*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_390 = tmp_19*(tmp_389 + tmp_48);
      real_t tmp_391 = tmp_390*tmp_43;
      real_t tmp_392 = tmp_390*tmp_51;
      real_t tmp_393 = tmp_390*tmp_53;
      real_t tmp_394 = tmp_74*(tmp_379 + tmp_75);
      real_t tmp_395 = tmp_74*(tmp_383 + tmp_78);
      real_t tmp_396 = tmp_74*(tmp_389 + tmp_81);
      real_t tmp_397 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_394*tmp_64 + tmp_395*tmp_77 + tmp_396*tmp_80 - 1.0/4.0) + tmp_59*(tmp_394*tmp_86 + tmp_395*tmp_87 + tmp_396*tmp_88 - 1.0/4.0) + tmp_62*(tmp_394*tmp_83 + tmp_395*tmp_84 + tmp_396*tmp_85 - 1.0/4.0));
      real_t tmp_398 = 0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_399 = tmp_19*(tmp_24 + tmp_398);
      real_t tmp_400 = tmp_399*tmp_8;
      real_t tmp_401 = tmp_27*tmp_399;
      real_t tmp_402 = 0.031405749086161582*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_403 = tmp_19*(tmp_34 + tmp_402);
      real_t tmp_404 = tmp_29*tmp_403;
      real_t tmp_405 = tmp_37*tmp_403;
      real_t tmp_406 = tmp_39*tmp_399;
      real_t tmp_407 = tmp_403*tmp_41;
      real_t tmp_408 = 0.031405749086161582*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_409 = tmp_19*(tmp_408 + tmp_48);
      real_t tmp_410 = tmp_409*tmp_43;
      real_t tmp_411 = tmp_409*tmp_51;
      real_t tmp_412 = tmp_409*tmp_53;
      real_t tmp_413 = tmp_74*(tmp_398 + tmp_75);
      real_t tmp_414 = tmp_74*(tmp_402 + tmp_78);
      real_t tmp_415 = tmp_74*(tmp_408 + tmp_81);
      real_t tmp_416 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_413*tmp_64 + tmp_414*tmp_77 + tmp_415*tmp_80 - 1.0/4.0) + tmp_59*(tmp_413*tmp_86 + tmp_414*tmp_87 + tmp_415*tmp_88 - 1.0/4.0) + tmp_62*(tmp_413*tmp_83 + tmp_414*tmp_84 + tmp_415*tmp_85 - 1.0/4.0));
      real_t tmp_417 = 0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_418 = tmp_19*(tmp_24 + tmp_417);
      real_t tmp_419 = tmp_418*tmp_8;
      real_t tmp_420 = tmp_27*tmp_418;
      real_t tmp_421 = 0.19601935860219369*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_422 = tmp_19*(tmp_34 + tmp_421);
      real_t tmp_423 = tmp_29*tmp_422;
      real_t tmp_424 = tmp_37*tmp_422;
      real_t tmp_425 = tmp_39*tmp_418;
      real_t tmp_426 = tmp_41*tmp_422;
      real_t tmp_427 = 0.19601935860219369*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_428 = tmp_19*(tmp_427 + tmp_48);
      real_t tmp_429 = tmp_428*tmp_43;
      real_t tmp_430 = tmp_428*tmp_51;
      real_t tmp_431 = tmp_428*tmp_53;
      real_t tmp_432 = tmp_74*(tmp_417 + tmp_75);
      real_t tmp_433 = tmp_74*(tmp_421 + tmp_78);
      real_t tmp_434 = tmp_74*(tmp_427 + tmp_81);
      real_t tmp_435 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_432*tmp_64 + tmp_433*tmp_77 + tmp_434*tmp_80 - 1.0/4.0) + tmp_59*(tmp_432*tmp_86 + tmp_433*tmp_87 + tmp_434*tmp_88 - 1.0/4.0) + tmp_62*(tmp_432*tmp_83 + tmp_433*tmp_84 + tmp_434*tmp_85 - 1.0/4.0));
      real_t tmp_436 = 0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_437 = tmp_19*(tmp_24 + tmp_436);
      real_t tmp_438 = tmp_437*tmp_8;
      real_t tmp_439 = tmp_27*tmp_437;
      real_t tmp_440 = 0.40446199974765351*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_441 = tmp_19*(tmp_34 + tmp_440);
      real_t tmp_442 = tmp_29*tmp_441;
      real_t tmp_443 = tmp_37*tmp_441;
      real_t tmp_444 = tmp_39*tmp_437;
      real_t tmp_445 = tmp_41*tmp_441;
      real_t tmp_446 = 0.40446199974765351*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_447 = tmp_19*(tmp_446 + tmp_48);
      real_t tmp_448 = tmp_43*tmp_447;
      real_t tmp_449 = tmp_447*tmp_51;
      real_t tmp_450 = tmp_447*tmp_53;
      real_t tmp_451 = tmp_74*(tmp_436 + tmp_75);
      real_t tmp_452 = tmp_74*(tmp_440 + tmp_78);
      real_t tmp_453 = tmp_74*(tmp_446 + tmp_81);
      real_t tmp_454 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_451*tmp_64 + tmp_452*tmp_77 + tmp_453*tmp_80 - 1.0/4.0) + tmp_59*(tmp_451*tmp_86 + tmp_452*tmp_87 + tmp_453*tmp_88 - 1.0/4.0) + tmp_62*(tmp_451*tmp_83 + tmp_452*tmp_84 + tmp_453*tmp_85 - 1.0/4.0));
      real_t tmp_455 = 0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_456 = tmp_19*(tmp_24 + tmp_455);
      real_t tmp_457 = tmp_456*tmp_8;
      real_t tmp_458 = tmp_27*tmp_456;
      real_t tmp_459 = 0.1711304259088916*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_460 = tmp_19*(tmp_34 + tmp_459);
      real_t tmp_461 = tmp_29*tmp_460;
      real_t tmp_462 = tmp_37*tmp_460;
      real_t tmp_463 = tmp_39*tmp_456;
      real_t tmp_464 = tmp_41*tmp_460;
      real_t tmp_465 = 0.1711304259088916*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_466 = tmp_19*(tmp_465 + tmp_48);
      real_t tmp_467 = tmp_43*tmp_466;
      real_t tmp_468 = tmp_466*tmp_51;
      real_t tmp_469 = tmp_466*tmp_53;
      real_t tmp_470 = tmp_74*(tmp_455 + tmp_75);
      real_t tmp_471 = tmp_74*(tmp_459 + tmp_78);
      real_t tmp_472 = tmp_74*(tmp_465 + tmp_81);
      real_t tmp_473 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_470*tmp_64 + tmp_471*tmp_77 + tmp_472*tmp_80 - 1.0/4.0) + tmp_59*(tmp_470*tmp_86 + tmp_471*tmp_87 + tmp_472*tmp_88 - 1.0/4.0) + tmp_62*(tmp_470*tmp_83 + tmp_471*tmp_84 + tmp_472*tmp_85 - 1.0/4.0));
      real_t a_0_0 = -tmp_112*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_106 - tmp_107 - tmp_108 - tmp_96 - tmp_97 + 1) - tmp_131*(-tmp_115 - tmp_116 - tmp_119 - tmp_120 - tmp_121 - tmp_122 - tmp_125 - tmp_126 - tmp_127 + 1) - tmp_150*(-tmp_134 - tmp_135 - tmp_138 - tmp_139 - tmp_140 - tmp_141 - tmp_144 - tmp_145 - tmp_146 + 1) - tmp_169*(-tmp_153 - tmp_154 - tmp_157 - tmp_158 - tmp_159 - tmp_160 - tmp_163 - tmp_164 - tmp_165 + 1) - tmp_188*(-tmp_172 - tmp_173 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_182 - tmp_183 - tmp_184 + 1) - tmp_207*(-tmp_191 - tmp_192 - tmp_195 - tmp_196 - tmp_197 - tmp_198 - tmp_201 - tmp_202 - tmp_203 + 1) - tmp_226*(-tmp_210 - tmp_211 - tmp_214 - tmp_215 - tmp_216 - tmp_217 - tmp_220 - tmp_221 - tmp_222 + 1) - tmp_245*(-tmp_229 - tmp_230 - tmp_233 - tmp_234 - tmp_235 - tmp_236 - tmp_239 - tmp_240 - tmp_241 + 1) - tmp_264*(-tmp_248 - tmp_249 - tmp_252 - tmp_253 - tmp_254 - tmp_255 - tmp_258 - tmp_259 - tmp_260 + 1) - tmp_283*(-tmp_267 - tmp_268 - tmp_271 - tmp_272 - tmp_273 - tmp_274 - tmp_277 - tmp_278 - tmp_279 + 1) - tmp_302*(-tmp_286 - tmp_287 - tmp_290 - tmp_291 - tmp_292 - tmp_293 - tmp_296 - tmp_297 - tmp_298 + 1) - tmp_321*(-tmp_305 - tmp_306 - tmp_309 - tmp_310 - tmp_311 - tmp_312 - tmp_315 - tmp_316 - tmp_317 + 1) - tmp_340*(-tmp_324 - tmp_325 - tmp_328 - tmp_329 - tmp_330 - tmp_331 - tmp_334 - tmp_335 - tmp_336 + 1) - tmp_359*(-tmp_343 - tmp_344 - tmp_347 - tmp_348 - tmp_349 - tmp_350 - tmp_353 - tmp_354 - tmp_355 + 1) - tmp_378*(-tmp_362 - tmp_363 - tmp_366 - tmp_367 - tmp_368 - tmp_369 - tmp_372 - tmp_373 - tmp_374 + 1) - tmp_397*(-tmp_381 - tmp_382 - tmp_385 - tmp_386 - tmp_387 - tmp_388 - tmp_391 - tmp_392 - tmp_393 + 1) - tmp_416*(-tmp_400 - tmp_401 - tmp_404 - tmp_405 - tmp_406 - tmp_407 - tmp_410 - tmp_411 - tmp_412 + 1) - tmp_435*(-tmp_419 - tmp_420 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_429 - tmp_430 - tmp_431 + 1) - tmp_454*(-tmp_438 - tmp_439 - tmp_442 - tmp_443 - tmp_444 - tmp_445 - tmp_448 - tmp_449 - tmp_450 + 1) - tmp_473*(-tmp_457 - tmp_458 - tmp_461 - tmp_462 - tmp_463 - tmp_464 - tmp_467 - tmp_468 - tmp_469 + 1) - tmp_93*(-tmp_26 - tmp_28 - tmp_36 - tmp_38 - tmp_40 - tmp_42 - tmp_50 - tmp_52 - tmp_54 + 1);
      real_t a_0_1 = -tmp_112*(tmp_102 + tmp_103 + tmp_108) - tmp_131*(tmp_121 + tmp_122 + tmp_127) - tmp_150*(tmp_140 + tmp_141 + tmp_146) - tmp_169*(tmp_159 + tmp_160 + tmp_165) - tmp_188*(tmp_178 + tmp_179 + tmp_184) - tmp_207*(tmp_197 + tmp_198 + tmp_203) - tmp_226*(tmp_216 + tmp_217 + tmp_222) - tmp_245*(tmp_235 + tmp_236 + tmp_241) - tmp_264*(tmp_254 + tmp_255 + tmp_260) - tmp_283*(tmp_273 + tmp_274 + tmp_279) - tmp_302*(tmp_292 + tmp_293 + tmp_298) - tmp_321*(tmp_311 + tmp_312 + tmp_317) - tmp_340*(tmp_330 + tmp_331 + tmp_336) - tmp_359*(tmp_349 + tmp_350 + tmp_355) - tmp_378*(tmp_368 + tmp_369 + tmp_374) - tmp_397*(tmp_387 + tmp_388 + tmp_393) - tmp_416*(tmp_406 + tmp_407 + tmp_412) - tmp_435*(tmp_425 + tmp_426 + tmp_431) - tmp_454*(tmp_444 + tmp_445 + tmp_450) - tmp_473*(tmp_463 + tmp_464 + tmp_469) - tmp_93*(tmp_40 + tmp_42 + tmp_54);
      real_t a_0_2 = -tmp_112*(tmp_101 + tmp_107 + tmp_97) - tmp_131*(tmp_116 + tmp_120 + tmp_126) - tmp_150*(tmp_135 + tmp_139 + tmp_145) - tmp_169*(tmp_154 + tmp_158 + tmp_164) - tmp_188*(tmp_173 + tmp_177 + tmp_183) - tmp_207*(tmp_192 + tmp_196 + tmp_202) - tmp_226*(tmp_211 + tmp_215 + tmp_221) - tmp_245*(tmp_230 + tmp_234 + tmp_240) - tmp_264*(tmp_249 + tmp_253 + tmp_259) - tmp_283*(tmp_268 + tmp_272 + tmp_278) - tmp_302*(tmp_287 + tmp_291 + tmp_297) - tmp_321*(tmp_306 + tmp_310 + tmp_316) - tmp_340*(tmp_325 + tmp_329 + tmp_335) - tmp_359*(tmp_344 + tmp_348 + tmp_354) - tmp_378*(tmp_363 + tmp_367 + tmp_373) - tmp_397*(tmp_382 + tmp_386 + tmp_392) - tmp_416*(tmp_401 + tmp_405 + tmp_411) - tmp_435*(tmp_420 + tmp_424 + tmp_430) - tmp_454*(tmp_439 + tmp_443 + tmp_449) - tmp_473*(tmp_458 + tmp_462 + tmp_468) - tmp_93*(tmp_28 + tmp_38 + tmp_52);
      real_t a_0_3 = -tmp_112*(tmp_100 + tmp_106 + tmp_96) - tmp_131*(tmp_115 + tmp_119 + tmp_125) - tmp_150*(tmp_134 + tmp_138 + tmp_144) - tmp_169*(tmp_153 + tmp_157 + tmp_163) - tmp_188*(tmp_172 + tmp_176 + tmp_182) - tmp_207*(tmp_191 + tmp_195 + tmp_201) - tmp_226*(tmp_210 + tmp_214 + tmp_220) - tmp_245*(tmp_229 + tmp_233 + tmp_239) - tmp_264*(tmp_248 + tmp_252 + tmp_258) - tmp_283*(tmp_267 + tmp_271 + tmp_277) - tmp_302*(tmp_286 + tmp_290 + tmp_296) - tmp_321*(tmp_305 + tmp_309 + tmp_315) - tmp_340*(tmp_324 + tmp_328 + tmp_334) - tmp_359*(tmp_343 + tmp_347 + tmp_353) - tmp_378*(tmp_362 + tmp_366 + tmp_372) - tmp_397*(tmp_381 + tmp_385 + tmp_391) - tmp_416*(tmp_400 + tmp_404 + tmp_410) - tmp_435*(tmp_419 + tmp_423 + tmp_429) - tmp_454*(tmp_438 + tmp_442 + tmp_448) - tmp_473*(tmp_457 + tmp_461 + tmp_467) - tmp_93*(tmp_26 + tmp_36 + tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_0_3 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }


};




class EGNIPGVectorLaplaceFormEP1_2 : public hyteg::dg::DGForm
{
 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
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

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_2;
      real_t tmp_9 = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_11;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_16 + tmp_4*tmp_9 - tmp_7*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_21 = tmp_18*(tmp_12 - tmp_17);
      real_t tmp_22 = tmp_10*tmp_21 + tmp_13*tmp_20 + tmp_19*tmp_9;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_10*tmp_25 + tmp_13*tmp_24 + tmp_23*tmp_9;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_10*tmp_29 + tmp_13*tmp_28 + tmp_27*tmp_9;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_0_1 = 0.1666666666666668*tmp_45;
      real_t a_0_2 = 0.1666666666666668*tmp_46;
      real_t a_0_3 = 0.1666666666666668*tmp_47;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 3.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_52 - 1.0/4.0) + tmp_14*(tmp_53 - 1.0/4.0) + tmp_8*(tmp_54 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_72 - 1.0/4.0) + tmp_14*(tmp_73 - 1.0/4.0) + tmp_8*(tmp_74 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_88 - 1.0/4.0) + tmp_14*(tmp_89 - 1.0/4.0) + tmp_8*(tmp_90 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_104 - 1.0/4.0) + tmp_14*(tmp_105 - 1.0/4.0) + tmp_8*(tmp_106 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_120 - 1.0/4.0) + tmp_14*(tmp_121 - 1.0/4.0) + tmp_8*(tmp_122 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_136 - 1.0/4.0) + tmp_14*(tmp_137 - 1.0/4.0) + tmp_8*(tmp_138 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_152 - 1.0/4.0) + tmp_14*(tmp_153 - 1.0/4.0) + tmp_8*(tmp_154 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_168 - 1.0/4.0) + tmp_14*(tmp_169 - 1.0/4.0) + tmp_8*(tmp_170 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_184 - 1.0/4.0) + tmp_14*(tmp_185 - 1.0/4.0) + tmp_8*(tmp_186 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_200 - 1.0/4.0) + tmp_14*(tmp_201 - 1.0/4.0) + tmp_8*(tmp_202 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_216 - 1.0/4.0) + tmp_14*(tmp_217 - 1.0/4.0) + tmp_8*(tmp_218 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_232 - 1.0/4.0) + tmp_14*(tmp_233 - 1.0/4.0) + tmp_8*(tmp_234 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_248 - 1.0/4.0) + tmp_14*(tmp_249 - 1.0/4.0) + tmp_8*(tmp_250 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_264 - 1.0/4.0) + tmp_14*(tmp_265 - 1.0/4.0) + tmp_8*(tmp_266 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_280 - 1.0/4.0) + tmp_14*(tmp_281 - 1.0/4.0) + tmp_8*(tmp_282 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_296 - 1.0/4.0) + tmp_14*(tmp_297 - 1.0/4.0) + tmp_8*(tmp_298 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_312 - 1.0/4.0) + tmp_14*(tmp_313 - 1.0/4.0) + tmp_8*(tmp_314 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_328 - 1.0/4.0) + tmp_14*(tmp_329 - 1.0/4.0) + tmp_8*(tmp_330 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_344 - 1.0/4.0) + tmp_14*(tmp_345 - 1.0/4.0) + tmp_8*(tmp_346 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_360 - 1.0/4.0) + tmp_14*(tmp_361 - 1.0/4.0) + tmp_8*(tmp_362 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_376 - 1.0/4.0) + tmp_14*(tmp_377 - 1.0/4.0) + tmp_8*(tmp_378 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_0_1 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_0_2 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_0_3 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_4_0;
      real_t tmp_1 = p_affine_5_0 + tmp_0;
      real_t tmp_2 = -p_affine_4_1;
      real_t tmp_3 = p_affine_6_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_6_0 + tmp_0;
      real_t tmp_6 = p_affine_5_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_4_2;
      real_t tmp_10 = p_affine_7_2 + tmp_9;
      real_t tmp_11 = p_affine_7_1 + tmp_2;
      real_t tmp_12 = p_affine_5_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_7_0 + tmp_0;
      real_t tmp_15 = p_affine_6_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = 0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22;
      real_t tmp_24 = p_affine_8_2 + tmp_9;
      real_t tmp_25 = tmp_19*(tmp_23 + tmp_24);
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_28 = tmp_25*tmp_27;
      real_t tmp_29 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_30 = -p_affine_8_1;
      real_t tmp_31 = p_affine_9_1 + tmp_30;
      real_t tmp_32 = p_affine_10_1 + tmp_30;
      real_t tmp_33 = 0.031405749086161582*tmp_31 + 0.93718850182767688*tmp_32;
      real_t tmp_34 = p_affine_8_1 + tmp_2;
      real_t tmp_35 = tmp_19*(tmp_33 + tmp_34);
      real_t tmp_36 = tmp_29*tmp_35;
      real_t tmp_37 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_40 = tmp_25*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_42 = tmp_35*tmp_41;
      real_t tmp_43 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_44 = -p_affine_8_0;
      real_t tmp_45 = p_affine_9_0 + tmp_44;
      real_t tmp_46 = p_affine_10_0 + tmp_44;
      real_t tmp_47 = 0.031405749086161582*tmp_45 + 0.93718850182767688*tmp_46;
      real_t tmp_48 = p_affine_8_0 + tmp_0;
      real_t tmp_49 = tmp_19*(tmp_47 + tmp_48);
      real_t tmp_50 = tmp_43*tmp_49;
      real_t tmp_51 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_52 = tmp_49*tmp_51;
      real_t tmp_53 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_54 = tmp_49*tmp_53;
      real_t tmp_55 = -p_affine_0_2;
      real_t tmp_56 = p_affine_1_2 + tmp_55;
      real_t tmp_57 = -p_affine_0_0;
      real_t tmp_58 = p_affine_2_0 + tmp_57;
      real_t tmp_59 = -p_affine_0_1;
      real_t tmp_60 = p_affine_3_1 + tmp_59;
      real_t tmp_61 = tmp_58*tmp_60;
      real_t tmp_62 = p_affine_3_0 + tmp_57;
      real_t tmp_63 = p_affine_2_1 + tmp_59;
      real_t tmp_64 = tmp_62*tmp_63;
      real_t tmp_65 = tmp_61 - tmp_64;
      real_t tmp_66 = p_affine_1_0 + tmp_57;
      real_t tmp_67 = p_affine_3_2 + tmp_55;
      real_t tmp_68 = tmp_63*tmp_67;
      real_t tmp_69 = p_affine_1_1 + tmp_59;
      real_t tmp_70 = p_affine_2_2 + tmp_55;
      real_t tmp_71 = tmp_62*tmp_70;
      real_t tmp_72 = tmp_60*tmp_70;
      real_t tmp_73 = tmp_58*tmp_67;
      real_t tmp_74 = 1.0 / (tmp_56*tmp_61 - tmp_56*tmp_64 + tmp_66*tmp_68 - tmp_66*tmp_72 + tmp_69*tmp_71 - tmp_69*tmp_73);
      real_t tmp_75 = p_affine_8_2 + tmp_55;
      real_t tmp_76 = tmp_74*(tmp_23 + tmp_75);
      real_t tmp_77 = tmp_71 - tmp_73;
      real_t tmp_78 = p_affine_8_1 + tmp_59;
      real_t tmp_79 = tmp_74*(tmp_33 + tmp_78);
      real_t tmp_80 = tmp_68 - tmp_72;
      real_t tmp_81 = p_affine_8_0 + tmp_57;
      real_t tmp_82 = tmp_74*(tmp_47 + tmp_81);
      real_t tmp_83 = -tmp_60*tmp_66 + tmp_62*tmp_69;
      real_t tmp_84 = -tmp_56*tmp_62 + tmp_66*tmp_67;
      real_t tmp_85 = tmp_56*tmp_60 - tmp_67*tmp_69;
      real_t tmp_86 = -tmp_58*tmp_69 + tmp_63*tmp_66;
      real_t tmp_87 = tmp_56*tmp_58 - tmp_66*tmp_70;
      real_t tmp_88 = -tmp_56*tmp_63 + tmp_69*tmp_70;
      real_t tmp_89 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_90 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_91 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_92 = 3.0*std::pow((std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)*std::abs(tmp_22*tmp_89 - tmp_32*tmp_91)) + (std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)*std::abs(tmp_22*tmp_90 - tmp_46*tmp_91)) + (std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)*std::abs(tmp_32*tmp_90 - tmp_46*tmp_89)), 0.25);
      real_t tmp_93 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_65*tmp_76 + tmp_77*tmp_79 + tmp_80*tmp_82 - 1.0/4.0) + tmp_67*(tmp_76*tmp_86 + tmp_79*tmp_87 + tmp_82*tmp_88 - 1.0/4.0) + tmp_70*(tmp_76*tmp_83 + tmp_79*tmp_84 + tmp_82*tmp_85 - 1.0/4.0));
      real_t tmp_94 = 0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22;
      real_t tmp_95 = tmp_19*(tmp_24 + tmp_94);
      real_t tmp_96 = tmp_8*tmp_95;
      real_t tmp_97 = tmp_27*tmp_95;
      real_t tmp_98 = 0.19601935860219369*tmp_31 + 0.60796128279561268*tmp_32;
      real_t tmp_99 = tmp_19*(tmp_34 + tmp_98);
      real_t tmp_100 = tmp_29*tmp_99;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = tmp_39*tmp_95;
      real_t tmp_103 = tmp_41*tmp_99;
      real_t tmp_104 = 0.19601935860219369*tmp_45 + 0.60796128279561268*tmp_46;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_48);
      real_t tmp_106 = tmp_105*tmp_43;
      real_t tmp_107 = tmp_105*tmp_51;
      real_t tmp_108 = tmp_105*tmp_53;
      real_t tmp_109 = tmp_74*(tmp_75 + tmp_94);
      real_t tmp_110 = tmp_74*(tmp_78 + tmp_98);
      real_t tmp_111 = tmp_74*(tmp_104 + tmp_81);
      real_t tmp_112 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_109*tmp_65 + tmp_110*tmp_77 + tmp_111*tmp_80 - 1.0/4.0) + tmp_67*(tmp_109*tmp_86 + tmp_110*tmp_87 + tmp_111*tmp_88 - 1.0/4.0) + tmp_70*(tmp_109*tmp_83 + tmp_110*tmp_84 + tmp_111*tmp_85 - 1.0/4.0));
      real_t tmp_113 = 0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_114 = tmp_19*(tmp_113 + tmp_24);
      real_t tmp_115 = tmp_114*tmp_8;
      real_t tmp_116 = tmp_114*tmp_27;
      real_t tmp_117 = 0.37605877282253791*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_118 = tmp_19*(tmp_117 + tmp_34);
      real_t tmp_119 = tmp_118*tmp_29;
      real_t tmp_120 = tmp_118*tmp_37;
      real_t tmp_121 = tmp_114*tmp_39;
      real_t tmp_122 = tmp_118*tmp_41;
      real_t tmp_123 = 0.37605877282253791*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_124 = tmp_19*(tmp_123 + tmp_48);
      real_t tmp_125 = tmp_124*tmp_43;
      real_t tmp_126 = tmp_124*tmp_51;
      real_t tmp_127 = tmp_124*tmp_53;
      real_t tmp_128 = tmp_74*(tmp_113 + tmp_75);
      real_t tmp_129 = tmp_74*(tmp_117 + tmp_78);
      real_t tmp_130 = tmp_74*(tmp_123 + tmp_81);
      real_t tmp_131 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_128*tmp_65 + tmp_129*tmp_77 + tmp_130*tmp_80 - 1.0/4.0) + tmp_67*(tmp_128*tmp_86 + tmp_129*tmp_87 + tmp_130*tmp_88 - 1.0/4.0) + tmp_70*(tmp_128*tmp_83 + tmp_129*tmp_84 + tmp_130*tmp_85 - 1.0/4.0));
      real_t tmp_132 = 0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_133 = tmp_19*(tmp_132 + tmp_24);
      real_t tmp_134 = tmp_133*tmp_8;
      real_t tmp_135 = tmp_133*tmp_27;
      real_t tmp_136 = 0.78764240869137092*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_137 = tmp_19*(tmp_136 + tmp_34);
      real_t tmp_138 = tmp_137*tmp_29;
      real_t tmp_139 = tmp_137*tmp_37;
      real_t tmp_140 = tmp_133*tmp_39;
      real_t tmp_141 = tmp_137*tmp_41;
      real_t tmp_142 = 0.78764240869137092*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_143 = tmp_19*(tmp_142 + tmp_48);
      real_t tmp_144 = tmp_143*tmp_43;
      real_t tmp_145 = tmp_143*tmp_51;
      real_t tmp_146 = tmp_143*tmp_53;
      real_t tmp_147 = tmp_74*(tmp_132 + tmp_75);
      real_t tmp_148 = tmp_74*(tmp_136 + tmp_78);
      real_t tmp_149 = tmp_74*(tmp_142 + tmp_81);
      real_t tmp_150 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_147*tmp_65 + tmp_148*tmp_77 + tmp_149*tmp_80 - 1.0/4.0) + tmp_67*(tmp_147*tmp_86 + tmp_148*tmp_87 + tmp_149*tmp_88 - 1.0/4.0) + tmp_70*(tmp_147*tmp_83 + tmp_148*tmp_84 + tmp_149*tmp_85 - 1.0/4.0));
      real_t tmp_151 = 0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_152 = tmp_19*(tmp_151 + tmp_24);
      real_t tmp_153 = tmp_152*tmp_8;
      real_t tmp_154 = tmp_152*tmp_27;
      real_t tmp_155 = 0.58463275527740355*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_156 = tmp_19*(tmp_155 + tmp_34);
      real_t tmp_157 = tmp_156*tmp_29;
      real_t tmp_158 = tmp_156*tmp_37;
      real_t tmp_159 = tmp_152*tmp_39;
      real_t tmp_160 = tmp_156*tmp_41;
      real_t tmp_161 = 0.58463275527740355*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_162 = tmp_19*(tmp_161 + tmp_48);
      real_t tmp_163 = tmp_162*tmp_43;
      real_t tmp_164 = tmp_162*tmp_51;
      real_t tmp_165 = tmp_162*tmp_53;
      real_t tmp_166 = tmp_74*(tmp_151 + tmp_75);
      real_t tmp_167 = tmp_74*(tmp_155 + tmp_78);
      real_t tmp_168 = tmp_74*(tmp_161 + tmp_81);
      real_t tmp_169 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_166*tmp_65 + tmp_167*tmp_77 + tmp_168*tmp_80 - 1.0/4.0) + tmp_67*(tmp_166*tmp_86 + tmp_167*tmp_87 + tmp_168*tmp_88 - 1.0/4.0) + tmp_70*(tmp_166*tmp_83 + tmp_167*tmp_84 + tmp_168*tmp_85 - 1.0/4.0));
      real_t tmp_170 = 0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_24);
      real_t tmp_172 = tmp_171*tmp_8;
      real_t tmp_173 = tmp_171*tmp_27;
      real_t tmp_174 = 0.041227165399737475*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_34);
      real_t tmp_176 = tmp_175*tmp_29;
      real_t tmp_177 = tmp_175*tmp_37;
      real_t tmp_178 = tmp_171*tmp_39;
      real_t tmp_179 = tmp_175*tmp_41;
      real_t tmp_180 = 0.041227165399737475*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_181 = tmp_19*(tmp_180 + tmp_48);
      real_t tmp_182 = tmp_181*tmp_43;
      real_t tmp_183 = tmp_181*tmp_51;
      real_t tmp_184 = tmp_181*tmp_53;
      real_t tmp_185 = tmp_74*(tmp_170 + tmp_75);
      real_t tmp_186 = tmp_74*(tmp_174 + tmp_78);
      real_t tmp_187 = tmp_74*(tmp_180 + tmp_81);
      real_t tmp_188 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_185*tmp_65 + tmp_186*tmp_77 + tmp_187*tmp_80 - 1.0/4.0) + tmp_67*(tmp_185*tmp_86 + tmp_186*tmp_87 + tmp_187*tmp_88 - 1.0/4.0) + tmp_70*(tmp_185*tmp_83 + tmp_186*tmp_84 + tmp_187*tmp_85 - 1.0/4.0));
      real_t tmp_189 = 0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_190 = tmp_19*(tmp_189 + tmp_24);
      real_t tmp_191 = tmp_190*tmp_8;
      real_t tmp_192 = tmp_190*tmp_27;
      real_t tmp_193 = 0.039308471900058539*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_194 = tmp_19*(tmp_193 + tmp_34);
      real_t tmp_195 = tmp_194*tmp_29;
      real_t tmp_196 = tmp_194*tmp_37;
      real_t tmp_197 = tmp_190*tmp_39;
      real_t tmp_198 = tmp_194*tmp_41;
      real_t tmp_199 = 0.039308471900058539*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_200 = tmp_19*(tmp_199 + tmp_48);
      real_t tmp_201 = tmp_200*tmp_43;
      real_t tmp_202 = tmp_200*tmp_51;
      real_t tmp_203 = tmp_200*tmp_53;
      real_t tmp_204 = tmp_74*(tmp_189 + tmp_75);
      real_t tmp_205 = tmp_74*(tmp_193 + tmp_78);
      real_t tmp_206 = tmp_74*(tmp_199 + tmp_81);
      real_t tmp_207 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_204*tmp_65 + tmp_205*tmp_77 + tmp_206*tmp_80 - 1.0/4.0) + tmp_67*(tmp_204*tmp_86 + tmp_205*tmp_87 + tmp_206*tmp_88 - 1.0/4.0) + tmp_70*(tmp_204*tmp_83 + tmp_205*tmp_84 + tmp_206*tmp_85 - 1.0/4.0));
      real_t tmp_208 = 0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_209 = tmp_19*(tmp_208 + tmp_24);
      real_t tmp_210 = tmp_209*tmp_8;
      real_t tmp_211 = tmp_209*tmp_27;
      real_t tmp_212 = 0.78764240869137092*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_213 = tmp_19*(tmp_212 + tmp_34);
      real_t tmp_214 = tmp_213*tmp_29;
      real_t tmp_215 = tmp_213*tmp_37;
      real_t tmp_216 = tmp_209*tmp_39;
      real_t tmp_217 = tmp_213*tmp_41;
      real_t tmp_218 = 0.78764240869137092*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_48);
      real_t tmp_220 = tmp_219*tmp_43;
      real_t tmp_221 = tmp_219*tmp_51;
      real_t tmp_222 = tmp_219*tmp_53;
      real_t tmp_223 = tmp_74*(tmp_208 + tmp_75);
      real_t tmp_224 = tmp_74*(tmp_212 + tmp_78);
      real_t tmp_225 = tmp_74*(tmp_218 + tmp_81);
      real_t tmp_226 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_223*tmp_65 + tmp_224*tmp_77 + tmp_225*tmp_80 - 1.0/4.0) + tmp_67*(tmp_223*tmp_86 + tmp_224*tmp_87 + tmp_225*tmp_88 - 1.0/4.0) + tmp_70*(tmp_223*tmp_83 + tmp_224*tmp_84 + tmp_225*tmp_85 - 1.0/4.0));
      real_t tmp_227 = 0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_228 = tmp_19*(tmp_227 + tmp_24);
      real_t tmp_229 = tmp_228*tmp_8;
      real_t tmp_230 = tmp_228*tmp_27;
      real_t tmp_231 = 0.58463275527740355*tmp_31 + 0.039308471900058539*tmp_32;
      real_t tmp_232 = tmp_19*(tmp_231 + tmp_34);
      real_t tmp_233 = tmp_232*tmp_29;
      real_t tmp_234 = tmp_232*tmp_37;
      real_t tmp_235 = tmp_228*tmp_39;
      real_t tmp_236 = tmp_232*tmp_41;
      real_t tmp_237 = 0.58463275527740355*tmp_45 + 0.039308471900058539*tmp_46;
      real_t tmp_238 = tmp_19*(tmp_237 + tmp_48);
      real_t tmp_239 = tmp_238*tmp_43;
      real_t tmp_240 = tmp_238*tmp_51;
      real_t tmp_241 = tmp_238*tmp_53;
      real_t tmp_242 = tmp_74*(tmp_227 + tmp_75);
      real_t tmp_243 = tmp_74*(tmp_231 + tmp_78);
      real_t tmp_244 = tmp_74*(tmp_237 + tmp_81);
      real_t tmp_245 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_242*tmp_65 + tmp_243*tmp_77 + tmp_244*tmp_80 - 1.0/4.0) + tmp_67*(tmp_242*tmp_86 + tmp_243*tmp_87 + tmp_244*tmp_88 - 1.0/4.0) + tmp_70*(tmp_242*tmp_83 + tmp_243*tmp_84 + tmp_244*tmp_85 - 1.0/4.0));
      real_t tmp_246 = 0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_247 = tmp_19*(tmp_24 + tmp_246);
      real_t tmp_248 = tmp_247*tmp_8;
      real_t tmp_249 = tmp_247*tmp_27;
      real_t tmp_250 = 0.1711304259088916*tmp_31 + 0.78764240869137092*tmp_32;
      real_t tmp_251 = tmp_19*(tmp_250 + tmp_34);
      real_t tmp_252 = tmp_251*tmp_29;
      real_t tmp_253 = tmp_251*tmp_37;
      real_t tmp_254 = tmp_247*tmp_39;
      real_t tmp_255 = tmp_251*tmp_41;
      real_t tmp_256 = 0.1711304259088916*tmp_45 + 0.78764240869137092*tmp_46;
      real_t tmp_257 = tmp_19*(tmp_256 + tmp_48);
      real_t tmp_258 = tmp_257*tmp_43;
      real_t tmp_259 = tmp_257*tmp_51;
      real_t tmp_260 = tmp_257*tmp_53;
      real_t tmp_261 = tmp_74*(tmp_246 + tmp_75);
      real_t tmp_262 = tmp_74*(tmp_250 + tmp_78);
      real_t tmp_263 = tmp_74*(tmp_256 + tmp_81);
      real_t tmp_264 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_261*tmp_65 + tmp_262*tmp_77 + tmp_263*tmp_80 - 1.0/4.0) + tmp_67*(tmp_261*tmp_86 + tmp_262*tmp_87 + tmp_263*tmp_88 - 1.0/4.0) + tmp_70*(tmp_261*tmp_83 + tmp_262*tmp_84 + tmp_263*tmp_85 - 1.0/4.0));
      real_t tmp_265 = 0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_266 = tmp_19*(tmp_24 + tmp_265);
      real_t tmp_267 = tmp_266*tmp_8;
      real_t tmp_268 = tmp_266*tmp_27;
      real_t tmp_269 = 0.37605877282253791*tmp_31 + 0.58463275527740355*tmp_32;
      real_t tmp_270 = tmp_19*(tmp_269 + tmp_34);
      real_t tmp_271 = tmp_270*tmp_29;
      real_t tmp_272 = tmp_270*tmp_37;
      real_t tmp_273 = tmp_266*tmp_39;
      real_t tmp_274 = tmp_270*tmp_41;
      real_t tmp_275 = 0.37605877282253791*tmp_45 + 0.58463275527740355*tmp_46;
      real_t tmp_276 = tmp_19*(tmp_275 + tmp_48);
      real_t tmp_277 = tmp_276*tmp_43;
      real_t tmp_278 = tmp_276*tmp_51;
      real_t tmp_279 = tmp_276*tmp_53;
      real_t tmp_280 = tmp_74*(tmp_265 + tmp_75);
      real_t tmp_281 = tmp_74*(tmp_269 + tmp_78);
      real_t tmp_282 = tmp_74*(tmp_275 + tmp_81);
      real_t tmp_283 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_280*tmp_65 + tmp_281*tmp_77 + tmp_282*tmp_80 - 1.0/4.0) + tmp_67*(tmp_280*tmp_86 + tmp_281*tmp_87 + tmp_282*tmp_88 - 1.0/4.0) + tmp_70*(tmp_280*tmp_83 + tmp_281*tmp_84 + tmp_282*tmp_85 - 1.0/4.0));
      real_t tmp_284 = 0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_285 = tmp_19*(tmp_24 + tmp_284);
      real_t tmp_286 = tmp_285*tmp_8;
      real_t tmp_287 = tmp_27*tmp_285;
      real_t tmp_288 = 0.041227165399737475*tmp_31 + 0.1711304259088916*tmp_32;
      real_t tmp_289 = tmp_19*(tmp_288 + tmp_34);
      real_t tmp_290 = tmp_289*tmp_29;
      real_t tmp_291 = tmp_289*tmp_37;
      real_t tmp_292 = tmp_285*tmp_39;
      real_t tmp_293 = tmp_289*tmp_41;
      real_t tmp_294 = 0.041227165399737475*tmp_45 + 0.1711304259088916*tmp_46;
      real_t tmp_295 = tmp_19*(tmp_294 + tmp_48);
      real_t tmp_296 = tmp_295*tmp_43;
      real_t tmp_297 = tmp_295*tmp_51;
      real_t tmp_298 = tmp_295*tmp_53;
      real_t tmp_299 = tmp_74*(tmp_284 + tmp_75);
      real_t tmp_300 = tmp_74*(tmp_288 + tmp_78);
      real_t tmp_301 = tmp_74*(tmp_294 + tmp_81);
      real_t tmp_302 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_299*tmp_65 + tmp_300*tmp_77 + tmp_301*tmp_80 - 1.0/4.0) + tmp_67*(tmp_299*tmp_86 + tmp_300*tmp_87 + tmp_301*tmp_88 - 1.0/4.0) + tmp_70*(tmp_299*tmp_83 + tmp_300*tmp_84 + tmp_301*tmp_85 - 1.0/4.0));
      real_t tmp_303 = 0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22;
      real_t tmp_304 = tmp_19*(tmp_24 + tmp_303);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_27*tmp_304;
      real_t tmp_307 = 0.40446199974765351*tmp_31 + 0.19107600050469298*tmp_32;
      real_t tmp_308 = tmp_19*(tmp_307 + tmp_34);
      real_t tmp_309 = tmp_29*tmp_308;
      real_t tmp_310 = tmp_308*tmp_37;
      real_t tmp_311 = tmp_304*tmp_39;
      real_t tmp_312 = tmp_308*tmp_41;
      real_t tmp_313 = 0.40446199974765351*tmp_45 + 0.19107600050469298*tmp_46;
      real_t tmp_314 = tmp_19*(tmp_313 + tmp_48);
      real_t tmp_315 = tmp_314*tmp_43;
      real_t tmp_316 = tmp_314*tmp_51;
      real_t tmp_317 = tmp_314*tmp_53;
      real_t tmp_318 = tmp_74*(tmp_303 + tmp_75);
      real_t tmp_319 = tmp_74*(tmp_307 + tmp_78);
      real_t tmp_320 = tmp_74*(tmp_313 + tmp_81);
      real_t tmp_321 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_318*tmp_65 + tmp_319*tmp_77 + tmp_320*tmp_80 - 1.0/4.0) + tmp_67*(tmp_318*tmp_86 + tmp_319*tmp_87 + tmp_320*tmp_88 - 1.0/4.0) + tmp_70*(tmp_318*tmp_83 + tmp_319*tmp_84 + tmp_320*tmp_85 - 1.0/4.0));
      real_t tmp_322 = 0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_323 = tmp_19*(tmp_24 + tmp_322);
      real_t tmp_324 = tmp_323*tmp_8;
      real_t tmp_325 = tmp_27*tmp_323;
      real_t tmp_326 = 0.039308471900058539*tmp_31 + 0.37605877282253791*tmp_32;
      real_t tmp_327 = tmp_19*(tmp_326 + tmp_34);
      real_t tmp_328 = tmp_29*tmp_327;
      real_t tmp_329 = tmp_327*tmp_37;
      real_t tmp_330 = tmp_323*tmp_39;
      real_t tmp_331 = tmp_327*tmp_41;
      real_t tmp_332 = 0.039308471900058539*tmp_45 + 0.37605877282253791*tmp_46;
      real_t tmp_333 = tmp_19*(tmp_332 + tmp_48);
      real_t tmp_334 = tmp_333*tmp_43;
      real_t tmp_335 = tmp_333*tmp_51;
      real_t tmp_336 = tmp_333*tmp_53;
      real_t tmp_337 = tmp_74*(tmp_322 + tmp_75);
      real_t tmp_338 = tmp_74*(tmp_326 + tmp_78);
      real_t tmp_339 = tmp_74*(tmp_332 + tmp_81);
      real_t tmp_340 = 0.020848748529055869*tmp_92*(tmp_56*(tmp_337*tmp_65 + tmp_338*tmp_77 + tmp_339*tmp_80 - 1.0/4.0) + tmp_67*(tmp_337*tmp_86 + tmp_338*tmp_87 + tmp_339*tmp_88 - 1.0/4.0) + tmp_70*(tmp_337*tmp_83 + tmp_338*tmp_84 + tmp_339*tmp_85 - 1.0/4.0));
      real_t tmp_341 = 0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_342 = tmp_19*(tmp_24 + tmp_341);
      real_t tmp_343 = tmp_342*tmp_8;
      real_t tmp_344 = tmp_27*tmp_342;
      real_t tmp_345 = 0.93718850182767688*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_346 = tmp_19*(tmp_34 + tmp_345);
      real_t tmp_347 = tmp_29*tmp_346;
      real_t tmp_348 = tmp_346*tmp_37;
      real_t tmp_349 = tmp_342*tmp_39;
      real_t tmp_350 = tmp_346*tmp_41;
      real_t tmp_351 = 0.93718850182767688*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_352 = tmp_19*(tmp_351 + tmp_48);
      real_t tmp_353 = tmp_352*tmp_43;
      real_t tmp_354 = tmp_352*tmp_51;
      real_t tmp_355 = tmp_352*tmp_53;
      real_t tmp_356 = tmp_74*(tmp_341 + tmp_75);
      real_t tmp_357 = tmp_74*(tmp_345 + tmp_78);
      real_t tmp_358 = tmp_74*(tmp_351 + tmp_81);
      real_t tmp_359 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_356*tmp_65 + tmp_357*tmp_77 + tmp_358*tmp_80 - 1.0/4.0) + tmp_67*(tmp_356*tmp_86 + tmp_357*tmp_87 + tmp_358*tmp_88 - 1.0/4.0) + tmp_70*(tmp_356*tmp_83 + tmp_357*tmp_84 + tmp_358*tmp_85 - 1.0/4.0));
      real_t tmp_360 = 0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_361 = tmp_19*(tmp_24 + tmp_360);
      real_t tmp_362 = tmp_361*tmp_8;
      real_t tmp_363 = tmp_27*tmp_361;
      real_t tmp_364 = 0.60796128279561268*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_365 = tmp_19*(tmp_34 + tmp_364);
      real_t tmp_366 = tmp_29*tmp_365;
      real_t tmp_367 = tmp_365*tmp_37;
      real_t tmp_368 = tmp_361*tmp_39;
      real_t tmp_369 = tmp_365*tmp_41;
      real_t tmp_370 = 0.60796128279561268*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_371 = tmp_19*(tmp_370 + tmp_48);
      real_t tmp_372 = tmp_371*tmp_43;
      real_t tmp_373 = tmp_371*tmp_51;
      real_t tmp_374 = tmp_371*tmp_53;
      real_t tmp_375 = tmp_74*(tmp_360 + tmp_75);
      real_t tmp_376 = tmp_74*(tmp_364 + tmp_78);
      real_t tmp_377 = tmp_74*(tmp_370 + tmp_81);
      real_t tmp_378 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_375*tmp_65 + tmp_376*tmp_77 + tmp_377*tmp_80 - 1.0/4.0) + tmp_67*(tmp_375*tmp_86 + tmp_376*tmp_87 + tmp_377*tmp_88 - 1.0/4.0) + tmp_70*(tmp_375*tmp_83 + tmp_376*tmp_84 + tmp_377*tmp_85 - 1.0/4.0));
      real_t tmp_379 = 0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_380 = tmp_19*(tmp_24 + tmp_379);
      real_t tmp_381 = tmp_380*tmp_8;
      real_t tmp_382 = tmp_27*tmp_380;
      real_t tmp_383 = 0.19107600050469298*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_384 = tmp_19*(tmp_34 + tmp_383);
      real_t tmp_385 = tmp_29*tmp_384;
      real_t tmp_386 = tmp_37*tmp_384;
      real_t tmp_387 = tmp_380*tmp_39;
      real_t tmp_388 = tmp_384*tmp_41;
      real_t tmp_389 = 0.19107600050469298*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_390 = tmp_19*(tmp_389 + tmp_48);
      real_t tmp_391 = tmp_390*tmp_43;
      real_t tmp_392 = tmp_390*tmp_51;
      real_t tmp_393 = tmp_390*tmp_53;
      real_t tmp_394 = tmp_74*(tmp_379 + tmp_75);
      real_t tmp_395 = tmp_74*(tmp_383 + tmp_78);
      real_t tmp_396 = tmp_74*(tmp_389 + tmp_81);
      real_t tmp_397 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_394*tmp_65 + tmp_395*tmp_77 + tmp_396*tmp_80 - 1.0/4.0) + tmp_67*(tmp_394*tmp_86 + tmp_395*tmp_87 + tmp_396*tmp_88 - 1.0/4.0) + tmp_70*(tmp_394*tmp_83 + tmp_395*tmp_84 + tmp_396*tmp_85 - 1.0/4.0));
      real_t tmp_398 = 0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_399 = tmp_19*(tmp_24 + tmp_398);
      real_t tmp_400 = tmp_399*tmp_8;
      real_t tmp_401 = tmp_27*tmp_399;
      real_t tmp_402 = 0.031405749086161582*tmp_31 + 0.031405749086161582*tmp_32;
      real_t tmp_403 = tmp_19*(tmp_34 + tmp_402);
      real_t tmp_404 = tmp_29*tmp_403;
      real_t tmp_405 = tmp_37*tmp_403;
      real_t tmp_406 = tmp_39*tmp_399;
      real_t tmp_407 = tmp_403*tmp_41;
      real_t tmp_408 = 0.031405749086161582*tmp_45 + 0.031405749086161582*tmp_46;
      real_t tmp_409 = tmp_19*(tmp_408 + tmp_48);
      real_t tmp_410 = tmp_409*tmp_43;
      real_t tmp_411 = tmp_409*tmp_51;
      real_t tmp_412 = tmp_409*tmp_53;
      real_t tmp_413 = tmp_74*(tmp_398 + tmp_75);
      real_t tmp_414 = tmp_74*(tmp_402 + tmp_78);
      real_t tmp_415 = tmp_74*(tmp_408 + tmp_81);
      real_t tmp_416 = 0.0068572537431980923*tmp_92*(tmp_56*(tmp_413*tmp_65 + tmp_414*tmp_77 + tmp_415*tmp_80 - 1.0/4.0) + tmp_67*(tmp_413*tmp_86 + tmp_414*tmp_87 + tmp_415*tmp_88 - 1.0/4.0) + tmp_70*(tmp_413*tmp_83 + tmp_414*tmp_84 + tmp_415*tmp_85 - 1.0/4.0));
      real_t tmp_417 = 0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_418 = tmp_19*(tmp_24 + tmp_417);
      real_t tmp_419 = tmp_418*tmp_8;
      real_t tmp_420 = tmp_27*tmp_418;
      real_t tmp_421 = 0.19601935860219369*tmp_31 + 0.19601935860219369*tmp_32;
      real_t tmp_422 = tmp_19*(tmp_34 + tmp_421);
      real_t tmp_423 = tmp_29*tmp_422;
      real_t tmp_424 = tmp_37*tmp_422;
      real_t tmp_425 = tmp_39*tmp_418;
      real_t tmp_426 = tmp_41*tmp_422;
      real_t tmp_427 = 0.19601935860219369*tmp_45 + 0.19601935860219369*tmp_46;
      real_t tmp_428 = tmp_19*(tmp_427 + tmp_48);
      real_t tmp_429 = tmp_428*tmp_43;
      real_t tmp_430 = tmp_428*tmp_51;
      real_t tmp_431 = tmp_428*tmp_53;
      real_t tmp_432 = tmp_74*(tmp_417 + tmp_75);
      real_t tmp_433 = tmp_74*(tmp_421 + tmp_78);
      real_t tmp_434 = tmp_74*(tmp_427 + tmp_81);
      real_t tmp_435 = 0.037198804536718075*tmp_92*(tmp_56*(tmp_432*tmp_65 + tmp_433*tmp_77 + tmp_434*tmp_80 - 1.0/4.0) + tmp_67*(tmp_432*tmp_86 + tmp_433*tmp_87 + tmp_434*tmp_88 - 1.0/4.0) + tmp_70*(tmp_432*tmp_83 + tmp_433*tmp_84 + tmp_434*tmp_85 - 1.0/4.0));
      real_t tmp_436 = 0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_437 = tmp_19*(tmp_24 + tmp_436);
      real_t tmp_438 = tmp_437*tmp_8;
      real_t tmp_439 = tmp_27*tmp_437;
      real_t tmp_440 = 0.40446199974765351*tmp_31 + 0.40446199974765351*tmp_32;
      real_t tmp_441 = tmp_19*(tmp_34 + tmp_440);
      real_t tmp_442 = tmp_29*tmp_441;
      real_t tmp_443 = tmp_37*tmp_441;
      real_t tmp_444 = tmp_39*tmp_437;
      real_t tmp_445 = tmp_41*tmp_441;
      real_t tmp_446 = 0.40446199974765351*tmp_45 + 0.40446199974765351*tmp_46;
      real_t tmp_447 = tmp_19*(tmp_446 + tmp_48);
      real_t tmp_448 = tmp_43*tmp_447;
      real_t tmp_449 = tmp_447*tmp_51;
      real_t tmp_450 = tmp_447*tmp_53;
      real_t tmp_451 = tmp_74*(tmp_436 + tmp_75);
      real_t tmp_452 = tmp_74*(tmp_440 + tmp_78);
      real_t tmp_453 = tmp_74*(tmp_446 + tmp_81);
      real_t tmp_454 = 0.042507265838595799*tmp_92*(tmp_56*(tmp_451*tmp_65 + tmp_452*tmp_77 + tmp_453*tmp_80 - 1.0/4.0) + tmp_67*(tmp_451*tmp_86 + tmp_452*tmp_87 + tmp_453*tmp_88 - 1.0/4.0) + tmp_70*(tmp_451*tmp_83 + tmp_452*tmp_84 + tmp_453*tmp_85 - 1.0/4.0));
      real_t tmp_455 = 0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_456 = tmp_19*(tmp_24 + tmp_455);
      real_t tmp_457 = tmp_456*tmp_8;
      real_t tmp_458 = tmp_27*tmp_456;
      real_t tmp_459 = 0.1711304259088916*tmp_31 + 0.041227165399737475*tmp_32;
      real_t tmp_460 = tmp_19*(tmp_34 + tmp_459);
      real_t tmp_461 = tmp_29*tmp_460;
      real_t tmp_462 = tmp_37*tmp_460;
      real_t tmp_463 = tmp_39*tmp_456;
      real_t tmp_464 = tmp_41*tmp_460;
      real_t tmp_465 = 0.1711304259088916*tmp_45 + 0.041227165399737475*tmp_46;
      real_t tmp_466 = tmp_19*(tmp_465 + tmp_48);
      real_t tmp_467 = tmp_43*tmp_466;
      real_t tmp_468 = tmp_466*tmp_51;
      real_t tmp_469 = tmp_466*tmp_53;
      real_t tmp_470 = tmp_74*(tmp_455 + tmp_75);
      real_t tmp_471 = tmp_74*(tmp_459 + tmp_78);
      real_t tmp_472 = tmp_74*(tmp_465 + tmp_81);
      real_t tmp_473 = 0.019202922745021479*tmp_92*(tmp_56*(tmp_470*tmp_65 + tmp_471*tmp_77 + tmp_472*tmp_80 - 1.0/4.0) + tmp_67*(tmp_470*tmp_86 + tmp_471*tmp_87 + tmp_472*tmp_88 - 1.0/4.0) + tmp_70*(tmp_470*tmp_83 + tmp_471*tmp_84 + tmp_472*tmp_85 - 1.0/4.0));
      real_t a_0_0 = -tmp_112*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_106 - tmp_107 - tmp_108 - tmp_96 - tmp_97 + 1) - tmp_131*(-tmp_115 - tmp_116 - tmp_119 - tmp_120 - tmp_121 - tmp_122 - tmp_125 - tmp_126 - tmp_127 + 1) - tmp_150*(-tmp_134 - tmp_135 - tmp_138 - tmp_139 - tmp_140 - tmp_141 - tmp_144 - tmp_145 - tmp_146 + 1) - tmp_169*(-tmp_153 - tmp_154 - tmp_157 - tmp_158 - tmp_159 - tmp_160 - tmp_163 - tmp_164 - tmp_165 + 1) - tmp_188*(-tmp_172 - tmp_173 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_182 - tmp_183 - tmp_184 + 1) - tmp_207*(-tmp_191 - tmp_192 - tmp_195 - tmp_196 - tmp_197 - tmp_198 - tmp_201 - tmp_202 - tmp_203 + 1) - tmp_226*(-tmp_210 - tmp_211 - tmp_214 - tmp_215 - tmp_216 - tmp_217 - tmp_220 - tmp_221 - tmp_222 + 1) - tmp_245*(-tmp_229 - tmp_230 - tmp_233 - tmp_234 - tmp_235 - tmp_236 - tmp_239 - tmp_240 - tmp_241 + 1) - tmp_264*(-tmp_248 - tmp_249 - tmp_252 - tmp_253 - tmp_254 - tmp_255 - tmp_258 - tmp_259 - tmp_260 + 1) - tmp_283*(-tmp_267 - tmp_268 - tmp_271 - tmp_272 - tmp_273 - tmp_274 - tmp_277 - tmp_278 - tmp_279 + 1) - tmp_302*(-tmp_286 - tmp_287 - tmp_290 - tmp_291 - tmp_292 - tmp_293 - tmp_296 - tmp_297 - tmp_298 + 1) - tmp_321*(-tmp_305 - tmp_306 - tmp_309 - tmp_310 - tmp_311 - tmp_312 - tmp_315 - tmp_316 - tmp_317 + 1) - tmp_340*(-tmp_324 - tmp_325 - tmp_328 - tmp_329 - tmp_330 - tmp_331 - tmp_334 - tmp_335 - tmp_336 + 1) - tmp_359*(-tmp_343 - tmp_344 - tmp_347 - tmp_348 - tmp_349 - tmp_350 - tmp_353 - tmp_354 - tmp_355 + 1) - tmp_378*(-tmp_362 - tmp_363 - tmp_366 - tmp_367 - tmp_368 - tmp_369 - tmp_372 - tmp_373 - tmp_374 + 1) - tmp_397*(-tmp_381 - tmp_382 - tmp_385 - tmp_386 - tmp_387 - tmp_388 - tmp_391 - tmp_392 - tmp_393 + 1) - tmp_416*(-tmp_400 - tmp_401 - tmp_404 - tmp_405 - tmp_406 - tmp_407 - tmp_410 - tmp_411 - tmp_412 + 1) - tmp_435*(-tmp_419 - tmp_420 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_429 - tmp_430 - tmp_431 + 1) - tmp_454*(-tmp_438 - tmp_439 - tmp_442 - tmp_443 - tmp_444 - tmp_445 - tmp_448 - tmp_449 - tmp_450 + 1) - tmp_473*(-tmp_457 - tmp_458 - tmp_461 - tmp_462 - tmp_463 - tmp_464 - tmp_467 - tmp_468 - tmp_469 + 1) - tmp_93*(-tmp_26 - tmp_28 - tmp_36 - tmp_38 - tmp_40 - tmp_42 - tmp_50 - tmp_52 - tmp_54 + 1);
      real_t a_0_1 = -tmp_112*(tmp_102 + tmp_103 + tmp_108) - tmp_131*(tmp_121 + tmp_122 + tmp_127) - tmp_150*(tmp_140 + tmp_141 + tmp_146) - tmp_169*(tmp_159 + tmp_160 + tmp_165) - tmp_188*(tmp_178 + tmp_179 + tmp_184) - tmp_207*(tmp_197 + tmp_198 + tmp_203) - tmp_226*(tmp_216 + tmp_217 + tmp_222) - tmp_245*(tmp_235 + tmp_236 + tmp_241) - tmp_264*(tmp_254 + tmp_255 + tmp_260) - tmp_283*(tmp_273 + tmp_274 + tmp_279) - tmp_302*(tmp_292 + tmp_293 + tmp_298) - tmp_321*(tmp_311 + tmp_312 + tmp_317) - tmp_340*(tmp_330 + tmp_331 + tmp_336) - tmp_359*(tmp_349 + tmp_350 + tmp_355) - tmp_378*(tmp_368 + tmp_369 + tmp_374) - tmp_397*(tmp_387 + tmp_388 + tmp_393) - tmp_416*(tmp_406 + tmp_407 + tmp_412) - tmp_435*(tmp_425 + tmp_426 + tmp_431) - tmp_454*(tmp_444 + tmp_445 + tmp_450) - tmp_473*(tmp_463 + tmp_464 + tmp_469) - tmp_93*(tmp_40 + tmp_42 + tmp_54);
      real_t a_0_2 = -tmp_112*(tmp_101 + tmp_107 + tmp_97) - tmp_131*(tmp_116 + tmp_120 + tmp_126) - tmp_150*(tmp_135 + tmp_139 + tmp_145) - tmp_169*(tmp_154 + tmp_158 + tmp_164) - tmp_188*(tmp_173 + tmp_177 + tmp_183) - tmp_207*(tmp_192 + tmp_196 + tmp_202) - tmp_226*(tmp_211 + tmp_215 + tmp_221) - tmp_245*(tmp_230 + tmp_234 + tmp_240) - tmp_264*(tmp_249 + tmp_253 + tmp_259) - tmp_283*(tmp_268 + tmp_272 + tmp_278) - tmp_302*(tmp_287 + tmp_291 + tmp_297) - tmp_321*(tmp_306 + tmp_310 + tmp_316) - tmp_340*(tmp_325 + tmp_329 + tmp_335) - tmp_359*(tmp_344 + tmp_348 + tmp_354) - tmp_378*(tmp_363 + tmp_367 + tmp_373) - tmp_397*(tmp_382 + tmp_386 + tmp_392) - tmp_416*(tmp_401 + tmp_405 + tmp_411) - tmp_435*(tmp_420 + tmp_424 + tmp_430) - tmp_454*(tmp_439 + tmp_443 + tmp_449) - tmp_473*(tmp_458 + tmp_462 + tmp_468) - tmp_93*(tmp_28 + tmp_38 + tmp_52);
      real_t a_0_3 = -tmp_112*(tmp_100 + tmp_106 + tmp_96) - tmp_131*(tmp_115 + tmp_119 + tmp_125) - tmp_150*(tmp_134 + tmp_138 + tmp_144) - tmp_169*(tmp_153 + tmp_157 + tmp_163) - tmp_188*(tmp_172 + tmp_176 + tmp_182) - tmp_207*(tmp_191 + tmp_195 + tmp_201) - tmp_226*(tmp_210 + tmp_214 + tmp_220) - tmp_245*(tmp_229 + tmp_233 + tmp_239) - tmp_264*(tmp_248 + tmp_252 + tmp_258) - tmp_283*(tmp_267 + tmp_271 + tmp_277) - tmp_302*(tmp_286 + tmp_290 + tmp_296) - tmp_321*(tmp_305 + tmp_309 + tmp_315) - tmp_340*(tmp_324 + tmp_328 + tmp_334) - tmp_359*(tmp_343 + tmp_347 + tmp_353) - tmp_378*(tmp_362 + tmp_366 + tmp_372) - tmp_397*(tmp_381 + tmp_385 + tmp_391) - tmp_416*(tmp_400 + tmp_404 + tmp_410) - tmp_435*(tmp_419 + tmp_423 + tmp_429) - tmp_454*(tmp_438 + tmp_442 + tmp_448) - tmp_473*(tmp_457 + tmp_461 + tmp_467) - tmp_93*(tmp_26 + tmp_36 + tmp_50);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_0_3 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }


};




class EGNIPGVectorLaplaceFormEE : public hyteg::dg::DGForm
{
 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_2_0 + tmp_0;
      real_t tmp_2 = p_affine_1_0 + tmp_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_1_1 + tmp_3;
      real_t tmp_7 = 1.0 / (-tmp_1*tmp_6 + tmp_5);
      real_t tmp_8 = tmp_2*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_5*tmp_7;
      real_t tmp_11 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_6*tmp_7;
      real_t tmp_13 = (((tmp_10 + tmp_12*tmp_9)*(tmp_10 + tmp_12*tmp_9)) + ((tmp_1*tmp_11 + tmp_10)*(tmp_1*tmp_11 + tmp_10)) + ((tmp_1*tmp_8 + tmp_8*tmp_9)*(tmp_1*tmp_8 + tmp_8*tmp_9)) + ((tmp_11*tmp_4 + tmp_12*tmp_4)*(tmp_11*tmp_4 + tmp_12*tmp_4)))*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.5*tmp_13;
      elMat( 0, 0) = a_0_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_3;
      real_t tmp_7 = 1.0 / (tmp_1*tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = p_affine_6_1 + tmp_3;
      real_t tmp_10 = tmp_7*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_7*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_10*tmp_2 + tmp_13*tmp_4 - 1.0/3.0;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_1*tmp_10 + tmp_13*tmp_15 - 1.0/3.0;
      real_t tmp_17 = tmp_7*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_18 = tmp_7*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_19 = tmp_17*tmp_2 + tmp_18*tmp_4 - 1.0/3.0;
      real_t tmp_20 = tmp_1*tmp_17 + tmp_15*tmp_18 - 1.0/3.0;
      real_t tmp_21 = tmp_7*(0.5*tmp_8 + tmp_9);
      real_t tmp_22 = tmp_7*(0.5*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_2*tmp_21 + tmp_22*tmp_4 - 1.0/3.0;
      real_t tmp_24 = tmp_1*tmp_21 + tmp_15*tmp_22 - 1.0/3.0;
      real_t tmp_25 = tmp_7*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_26 = tmp_7*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_27 = tmp_2*tmp_25 + tmp_26*tmp_4 - 1.0/3.0;
      real_t tmp_28 = tmp_1*tmp_25 + tmp_15*tmp_26 - 1.0/3.0;
      real_t tmp_29 = tmp_7*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_30 = tmp_7*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_31 = tmp_2*tmp_29 + tmp_30*tmp_4 - 1.0/3.0;
      real_t tmp_32 = tmp_1*tmp_29 + tmp_15*tmp_30 - 1.0/3.0;
      real_t a_0_0 = 0.35539032758428413*((tmp_1*tmp_14 + tmp_16*tmp_5)*(tmp_1*tmp_14 + tmp_16*tmp_5)) + 0.71794300574904923*((tmp_1*tmp_19 + tmp_20*tmp_5)*(tmp_1*tmp_19 + tmp_20*tmp_5)) + 0.8533333333333335*((tmp_1*tmp_23 + tmp_24*tmp_5)*(tmp_1*tmp_23 + tmp_24*tmp_5)) + 0.71794300574904923*((tmp_1*tmp_27 + tmp_28*tmp_5)*(tmp_1*tmp_27 + tmp_28*tmp_5)) + 0.35539032758428413*((tmp_1*tmp_31 + tmp_32*tmp_5)*(tmp_1*tmp_31 + tmp_32*tmp_5)) + 0.35539032758428413*((tmp_14*tmp_6 + tmp_16*tmp_4)*(tmp_14*tmp_6 + tmp_16*tmp_4)) + 0.71794300574904923*((tmp_19*tmp_6 + tmp_20*tmp_4)*(tmp_19*tmp_6 + tmp_20*tmp_4)) + 0.8533333333333335*((tmp_23*tmp_6 + tmp_24*tmp_4)*(tmp_23*tmp_6 + tmp_24*tmp_4)) + 0.71794300574904923*((tmp_27*tmp_6 + tmp_28*tmp_4)*(tmp_27*tmp_6 + tmp_28*tmp_4)) + 0.35539032758428413*((tmp_31*tmp_6 + tmp_32*tmp_4)*(tmp_31*tmp_6 + tmp_32*tmp_4));
      elMat( 0, 0) = a_0_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_3;
      real_t tmp_7 = 1.0 / (tmp_1*tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = p_affine_6_1 + 0.046910077030668018*tmp_8;
      real_t tmp_10 = tmp_7*(tmp_3 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018*tmp_11;
      real_t tmp_13 = tmp_7*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_10*tmp_2 + tmp_13*tmp_4 - 1.0/3.0;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_1*tmp_10 + tmp_13*tmp_15 - 1.0/3.0;
      real_t tmp_17 = -p_affine_3_0;
      real_t tmp_18 = p_affine_4_0 + tmp_17;
      real_t tmp_19 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_20 = -p_affine_3_1;
      real_t tmp_21 = p_affine_5_1 + tmp_20;
      real_t tmp_22 = p_affine_5_0 + tmp_17;
      real_t tmp_23 = p_affine_4_1 + tmp_20;
      real_t tmp_24 = 1.0 / (tmp_18*tmp_21 - tmp_22*tmp_23);
      real_t tmp_25 = tmp_24*(tmp_20 + tmp_9);
      real_t tmp_26 = tmp_24*(tmp_12 + tmp_17);
      real_t tmp_27 = tmp_19*tmp_25 + tmp_21*tmp_26 - 1.0/3.0;
      real_t tmp_28 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_29 = tmp_18*tmp_25 + tmp_26*tmp_28 - 1.0/3.0;
      real_t tmp_30 = p_affine_6_1 + 0.23076534494715845*tmp_8;
      real_t tmp_31 = tmp_7*(tmp_3 + tmp_30);
      real_t tmp_32 = p_affine_6_0 + 0.23076534494715845*tmp_11;
      real_t tmp_33 = tmp_7*(tmp_0 + tmp_32);
      real_t tmp_34 = tmp_2*tmp_31 + tmp_33*tmp_4 - 1.0/3.0;
      real_t tmp_35 = tmp_1*tmp_31 + tmp_15*tmp_33 - 1.0/3.0;
      real_t tmp_36 = tmp_24*(tmp_20 + tmp_30);
      real_t tmp_37 = tmp_24*(tmp_17 + tmp_32);
      real_t tmp_38 = tmp_19*tmp_36 + tmp_21*tmp_37 - 1.0/3.0;
      real_t tmp_39 = tmp_18*tmp_36 + tmp_28*tmp_37 - 1.0/3.0;
      real_t tmp_40 = p_affine_6_1 + 0.5*tmp_8;
      real_t tmp_41 = tmp_7*(tmp_3 + tmp_40);
      real_t tmp_42 = p_affine_6_0 + 0.5*tmp_11;
      real_t tmp_43 = tmp_7*(tmp_0 + tmp_42);
      real_t tmp_44 = tmp_2*tmp_41 + tmp_4*tmp_43 - 1.0/3.0;
      real_t tmp_45 = tmp_1*tmp_41 + tmp_15*tmp_43 - 1.0/3.0;
      real_t tmp_46 = tmp_24*(tmp_20 + tmp_40);
      real_t tmp_47 = tmp_24*(tmp_17 + tmp_42);
      real_t tmp_48 = tmp_19*tmp_46 + tmp_21*tmp_47 - 1.0/3.0;
      real_t tmp_49 = tmp_18*tmp_46 + tmp_28*tmp_47 - 1.0/3.0;
      real_t tmp_50 = p_affine_6_1 + 0.7692346550528415*tmp_8;
      real_t tmp_51 = tmp_7*(tmp_3 + tmp_50);
      real_t tmp_52 = p_affine_6_0 + 0.7692346550528415*tmp_11;
      real_t tmp_53 = tmp_7*(tmp_0 + tmp_52);
      real_t tmp_54 = tmp_2*tmp_51 + tmp_4*tmp_53 - 1.0/3.0;
      real_t tmp_55 = tmp_1*tmp_51 + tmp_15*tmp_53 - 1.0/3.0;
      real_t tmp_56 = tmp_24*(tmp_20 + tmp_50);
      real_t tmp_57 = tmp_24*(tmp_17 + tmp_52);
      real_t tmp_58 = tmp_19*tmp_56 + tmp_21*tmp_57 - 1.0/3.0;
      real_t tmp_59 = tmp_18*tmp_56 + tmp_28*tmp_57 - 1.0/3.0;
      real_t tmp_60 = p_affine_6_1 + 0.95308992296933193*tmp_8;
      real_t tmp_61 = tmp_7*(tmp_3 + tmp_60);
      real_t tmp_62 = p_affine_6_0 + 0.95308992296933193*tmp_11;
      real_t tmp_63 = tmp_7*(tmp_0 + tmp_62);
      real_t tmp_64 = tmp_2*tmp_61 + tmp_4*tmp_63 - 1.0/3.0;
      real_t tmp_65 = tmp_1*tmp_61 + tmp_15*tmp_63 - 1.0/3.0;
      real_t tmp_66 = tmp_24*(tmp_20 + tmp_60);
      real_t tmp_67 = tmp_24*(tmp_17 + tmp_62);
      real_t tmp_68 = tmp_19*tmp_66 + tmp_21*tmp_67 - 1.0/3.0;
      real_t tmp_69 = tmp_18*tmp_66 + tmp_28*tmp_67 - 1.0/3.0;
      real_t a_0_0 = -0.35539032758428413*(tmp_1*tmp_14 + tmp_16*tmp_5)*(tmp_18*tmp_27 + tmp_22*tmp_29) - 0.71794300574904923*(tmp_1*tmp_34 + tmp_35*tmp_5)*(tmp_18*tmp_38 + tmp_22*tmp_39) - 0.8533333333333335*(tmp_1*tmp_44 + tmp_45*tmp_5)*(tmp_18*tmp_48 + tmp_22*tmp_49) - 0.71794300574904923*(tmp_1*tmp_54 + tmp_5*tmp_55)*(tmp_18*tmp_58 + tmp_22*tmp_59) - 0.35539032758428413*(tmp_1*tmp_64 + tmp_5*tmp_65)*(tmp_18*tmp_68 + tmp_22*tmp_69) - 0.35539032758428413*(tmp_14*tmp_6 + tmp_16*tmp_4)*(tmp_21*tmp_29 + tmp_23*tmp_27) - 0.71794300574904923*(tmp_21*tmp_39 + tmp_23*tmp_38)*(tmp_34*tmp_6 + tmp_35*tmp_4) - 0.8533333333333335*(tmp_21*tmp_49 + tmp_23*tmp_48)*(tmp_4*tmp_45 + tmp_44*tmp_6) - 0.71794300574904923*(tmp_21*tmp_59 + tmp_23*tmp_58)*(tmp_4*tmp_55 + tmp_54*tmp_6) - 0.35539032758428413*(tmp_21*tmp_69 + tmp_23*tmp_68)*(tmp_4*tmp_65 + tmp_6*tmp_64);
      elMat( 0, 0) = a_0_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_3;
      real_t tmp_7 = 1.0 / (tmp_1*tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = p_affine_6_1 + tmp_3;
      real_t tmp_10 = tmp_7*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_7*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_10*tmp_2 + tmp_13*tmp_4 - 1.0/3.0;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_1*tmp_10 + tmp_13*tmp_15 - 1.0/3.0;
      real_t tmp_17 = tmp_7*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_18 = tmp_7*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_19 = tmp_17*tmp_2 + tmp_18*tmp_4 - 1.0/3.0;
      real_t tmp_20 = tmp_1*tmp_17 + tmp_15*tmp_18 - 1.0/3.0;
      real_t tmp_21 = tmp_7*(0.5*tmp_8 + tmp_9);
      real_t tmp_22 = tmp_7*(0.5*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_2*tmp_21 + tmp_22*tmp_4 - 1.0/3.0;
      real_t tmp_24 = tmp_1*tmp_21 + tmp_15*tmp_22 - 1.0/3.0;
      real_t tmp_25 = tmp_7*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_26 = tmp_7*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_27 = tmp_2*tmp_25 + tmp_26*tmp_4 - 1.0/3.0;
      real_t tmp_28 = tmp_1*tmp_25 + tmp_15*tmp_26 - 1.0/3.0;
      real_t tmp_29 = tmp_7*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_30 = tmp_7*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_31 = tmp_2*tmp_29 + tmp_30*tmp_4 - 1.0/3.0;
      real_t tmp_32 = tmp_1*tmp_29 + tmp_15*tmp_30 - 1.0/3.0;
      real_t a_0_0 = 0.35539032758428413*((tmp_1*tmp_14 + tmp_16*tmp_5)*(tmp_1*tmp_14 + tmp_16*tmp_5)) + 0.71794300574904923*((tmp_1*tmp_19 + tmp_20*tmp_5)*(tmp_1*tmp_19 + tmp_20*tmp_5)) + 0.8533333333333335*((tmp_1*tmp_23 + tmp_24*tmp_5)*(tmp_1*tmp_23 + tmp_24*tmp_5)) + 0.71794300574904923*((tmp_1*tmp_27 + tmp_28*tmp_5)*(tmp_1*tmp_27 + tmp_28*tmp_5)) + 0.35539032758428413*((tmp_1*tmp_31 + tmp_32*tmp_5)*(tmp_1*tmp_31 + tmp_32*tmp_5)) + 0.35539032758428413*((tmp_14*tmp_6 + tmp_16*tmp_4)*(tmp_14*tmp_6 + tmp_16*tmp_4)) + 0.71794300574904923*((tmp_19*tmp_6 + tmp_20*tmp_4)*(tmp_19*tmp_6 + tmp_20*tmp_4)) + 0.8533333333333335*((tmp_23*tmp_6 + tmp_24*tmp_4)*(tmp_23*tmp_6 + tmp_24*tmp_4)) + 0.71794300574904923*((tmp_27*tmp_6 + tmp_28*tmp_4)*(tmp_27*tmp_6 + tmp_28*tmp_4)) + 0.35539032758428413*((tmp_31*tmp_6 + tmp_32*tmp_4)*(tmp_31*tmp_6 + tmp_32*tmp_4));
      elMat( 0, 0) = a_0_0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_2_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_3_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_3_0 + tmp_0;
      real_t tmp_6 = p_affine_2_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = p_affine_1_0 + tmp_0;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_6*tmp_9;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_2_2 + tmp_10;
      real_t tmp_15 = p_affine_1_1 + tmp_2;
      real_t tmp_16 = tmp_15*tmp_5;
      real_t tmp_17 = tmp_3*tmp_9;
      real_t tmp_18 = tmp_1*tmp_15;
      real_t tmp_19 = 1.0 / (tmp_11*tmp_12 - tmp_11*tmp_18 + tmp_13*tmp_4 - tmp_13*tmp_7 + tmp_14*tmp_16 - tmp_14*tmp_17);
      real_t tmp_20 = tmp_19*tmp_9;
      real_t tmp_21 = tmp_16 - tmp_17;
      real_t tmp_22 = tmp_1*tmp_19;
      real_t tmp_23 = tmp_12 - tmp_18;
      real_t tmp_24 = tmp_19*tmp_5;
      real_t tmp_25 = -tmp_1*tmp_11 + tmp_14*tmp_5;
      real_t tmp_26 = tmp_11*tmp_9 - tmp_13*tmp_5;
      real_t tmp_27 = tmp_1*tmp_13 - tmp_14*tmp_9;
      real_t tmp_28 = tmp_11*tmp_6 - tmp_14*tmp_3;
      real_t tmp_29 = -tmp_11*tmp_15 + tmp_13*tmp_3;
      real_t tmp_30 = -tmp_13*tmp_6 + tmp_14*tmp_15;
      real_t tmp_31 = tmp_15*tmp_19;
      real_t tmp_32 = tmp_19*tmp_6;
      real_t tmp_33 = tmp_19*tmp_3;
      real_t tmp_34 = tmp_13*tmp_19;
      real_t tmp_35 = tmp_14*tmp_19;
      real_t tmp_36 = tmp_11*tmp_19;
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
      real_t tmp_49 = (((tmp_20*tmp_25 + tmp_22*tmp_26 + tmp_24*tmp_27)*(tmp_20*tmp_25 + tmp_22*tmp_26 + tmp_24*tmp_27)) + ((tmp_20*tmp_28 + tmp_22*tmp_29 + tmp_24*tmp_30)*(tmp_20*tmp_28 + tmp_22*tmp_29 + tmp_24*tmp_30)) + ((tmp_20*tmp_8 + tmp_21*tmp_22 + tmp_23*tmp_24)*(tmp_20*tmp_8 + tmp_21*tmp_22 + tmp_23*tmp_24)) + ((tmp_21*tmp_32 + tmp_23*tmp_33 + tmp_31*tmp_8)*(tmp_21*tmp_32 + tmp_23*tmp_33 + tmp_31*tmp_8)) + ((tmp_21*tmp_35 + tmp_23*tmp_36 + tmp_34*tmp_8)*(tmp_21*tmp_35 + tmp_23*tmp_36 + tmp_34*tmp_8)) + ((tmp_25*tmp_31 + tmp_26*tmp_32 + tmp_27*tmp_33)*(tmp_25*tmp_31 + tmp_26*tmp_32 + tmp_27*tmp_33)) + ((tmp_25*tmp_34 + tmp_26*tmp_35 + tmp_27*tmp_36)*(tmp_25*tmp_34 + tmp_26*tmp_35 + tmp_27*tmp_36)) + ((tmp_28*tmp_31 + tmp_29*tmp_32 + tmp_30*tmp_33)*(tmp_28*tmp_31 + tmp_29*tmp_32 + tmp_30*tmp_33)) + ((tmp_28*tmp_34 + tmp_29*tmp_35 + tmp_30*tmp_36)*(tmp_28*tmp_34 + tmp_29*tmp_35 + tmp_30*tmp_36)))*std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t a_0_0 = 0.1666666666666668*tmp_49;
      elMat( 0, 0) = a_0_0;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_2_0 + tmp_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_3_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_3_0 + tmp_0;
      real_t tmp_7 = p_affine_2_1 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11*tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_15*tmp_4;
      real_t tmp_18 = tmp_11*tmp_2;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_12 - tmp_1*tmp_17 + tmp_13*tmp_5 - tmp_13*tmp_8 + tmp_14*tmp_16 - tmp_14*tmp_18);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_10;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_16 - tmp_18;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_3;
      real_t tmp_30 = tmp_19*(0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29);
      real_t tmp_31 = tmp_12 - tmp_17;
      real_t tmp_32 = -p_affine_8_0;
      real_t tmp_33 = p_affine_9_0 + tmp_32;
      real_t tmp_34 = p_affine_10_0 + tmp_32;
      real_t tmp_35 = p_affine_8_0 + tmp_0;
      real_t tmp_36 = tmp_19*(0.031405749086161582*tmp_33 + 0.93718850182767688*tmp_34 + tmp_35);
      real_t tmp_37 = tmp_24*tmp_9 + tmp_25*tmp_30 + tmp_31*tmp_36 - 1.0/4.0;
      real_t tmp_38 = -tmp_1*tmp_4 + tmp_14*tmp_6;
      real_t tmp_39 = tmp_1*tmp_11 - tmp_13*tmp_6;
      real_t tmp_40 = -tmp_11*tmp_14 + tmp_13*tmp_4;
      real_t tmp_41 = tmp_24*tmp_38 + tmp_30*tmp_39 + tmp_36*tmp_40 - 1.0/4.0;
      real_t tmp_42 = tmp_1*tmp_7 - tmp_14*tmp_2;
      real_t tmp_43 = -tmp_1*tmp_15 + tmp_13*tmp_2;
      real_t tmp_44 = -tmp_13*tmp_7 + tmp_14*tmp_15;
      real_t tmp_45 = tmp_24*tmp_42 + tmp_30*tmp_43 + tmp_36*tmp_44 - 1.0/4.0;
      real_t tmp_46 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_47 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_48 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_49 = 3.0*std::pow((std::abs(tmp_22*tmp_46 - tmp_28*tmp_48)*std::abs(tmp_22*tmp_46 - tmp_28*tmp_48)) + (std::abs(tmp_22*tmp_47 - tmp_34*tmp_48)*std::abs(tmp_22*tmp_47 - tmp_34*tmp_48)) + (std::abs(tmp_28*tmp_47 - tmp_34*tmp_46)*std::abs(tmp_28*tmp_47 - tmp_34*tmp_46)), 0.25);
      real_t tmp_50 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_51 = tmp_19*(0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29);
      real_t tmp_52 = tmp_19*(0.19601935860219369*tmp_33 + 0.60796128279561268*tmp_34 + tmp_35);
      real_t tmp_53 = tmp_25*tmp_51 + tmp_31*tmp_52 + tmp_50*tmp_9 - 1.0/4.0;
      real_t tmp_54 = tmp_38*tmp_50 + tmp_39*tmp_51 + tmp_40*tmp_52 - 1.0/4.0;
      real_t tmp_55 = tmp_42*tmp_50 + tmp_43*tmp_51 + tmp_44*tmp_52 - 1.0/4.0;
      real_t tmp_56 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_57 = tmp_19*(0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_58 = tmp_19*(0.37605877282253791*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_59 = tmp_25*tmp_57 + tmp_31*tmp_58 + tmp_56*tmp_9 - 1.0/4.0;
      real_t tmp_60 = tmp_38*tmp_56 + tmp_39*tmp_57 + tmp_40*tmp_58 - 1.0/4.0;
      real_t tmp_61 = tmp_42*tmp_56 + tmp_43*tmp_57 + tmp_44*tmp_58 - 1.0/4.0;
      real_t tmp_62 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_63 = tmp_19*(0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_64 = tmp_19*(0.78764240869137092*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_65 = tmp_25*tmp_63 + tmp_31*tmp_64 + tmp_62*tmp_9 - 1.0/4.0;
      real_t tmp_66 = tmp_38*tmp_62 + tmp_39*tmp_63 + tmp_40*tmp_64 - 1.0/4.0;
      real_t tmp_67 = tmp_42*tmp_62 + tmp_43*tmp_63 + tmp_44*tmp_64 - 1.0/4.0;
      real_t tmp_68 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_69 = tmp_19*(0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_70 = tmp_19*(0.58463275527740355*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_71 = tmp_25*tmp_69 + tmp_31*tmp_70 + tmp_68*tmp_9 - 1.0/4.0;
      real_t tmp_72 = tmp_38*tmp_68 + tmp_39*tmp_69 + tmp_40*tmp_70 - 1.0/4.0;
      real_t tmp_73 = tmp_42*tmp_68 + tmp_43*tmp_69 + tmp_44*tmp_70 - 1.0/4.0;
      real_t tmp_74 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_75 = tmp_19*(0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_76 = tmp_19*(0.041227165399737475*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_77 = tmp_25*tmp_75 + tmp_31*tmp_76 + tmp_74*tmp_9 - 1.0/4.0;
      real_t tmp_78 = tmp_38*tmp_74 + tmp_39*tmp_75 + tmp_40*tmp_76 - 1.0/4.0;
      real_t tmp_79 = tmp_42*tmp_74 + tmp_43*tmp_75 + tmp_44*tmp_76 - 1.0/4.0;
      real_t tmp_80 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_81 = tmp_19*(0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_82 = tmp_19*(0.039308471900058539*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_83 = tmp_25*tmp_81 + tmp_31*tmp_82 + tmp_80*tmp_9 - 1.0/4.0;
      real_t tmp_84 = tmp_38*tmp_80 + tmp_39*tmp_81 + tmp_40*tmp_82 - 1.0/4.0;
      real_t tmp_85 = tmp_42*tmp_80 + tmp_43*tmp_81 + tmp_44*tmp_82 - 1.0/4.0;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_87 = tmp_19*(0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_88 = tmp_19*(0.78764240869137092*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_89 = tmp_25*tmp_87 + tmp_31*tmp_88 + tmp_86*tmp_9 - 1.0/4.0;
      real_t tmp_90 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88 - 1.0/4.0;
      real_t tmp_91 = tmp_42*tmp_86 + tmp_43*tmp_87 + tmp_44*tmp_88 - 1.0/4.0;
      real_t tmp_92 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_19*(0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_94 = tmp_19*(0.58463275527740355*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_95 = tmp_25*tmp_93 + tmp_31*tmp_94 + tmp_9*tmp_92 - 1.0/4.0;
      real_t tmp_96 = tmp_38*tmp_92 + tmp_39*tmp_93 + tmp_40*tmp_94 - 1.0/4.0;
      real_t tmp_97 = tmp_42*tmp_92 + tmp_43*tmp_93 + tmp_44*tmp_94 - 1.0/4.0;
      real_t tmp_98 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_99 = tmp_19*(0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_100 = tmp_19*(0.1711304259088916*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_101 = tmp_100*tmp_31 + tmp_25*tmp_99 + tmp_9*tmp_98 - 1.0/4.0;
      real_t tmp_102 = tmp_100*tmp_40 + tmp_38*tmp_98 + tmp_39*tmp_99 - 1.0/4.0;
      real_t tmp_103 = tmp_100*tmp_44 + tmp_42*tmp_98 + tmp_43*tmp_99 - 1.0/4.0;
      real_t tmp_104 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_105 = tmp_19*(0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_106 = tmp_19*(0.37605877282253791*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_107 = tmp_104*tmp_9 + tmp_105*tmp_25 + tmp_106*tmp_31 - 1.0/4.0;
      real_t tmp_108 = tmp_104*tmp_38 + tmp_105*tmp_39 + tmp_106*tmp_40 - 1.0/4.0;
      real_t tmp_109 = tmp_104*tmp_42 + tmp_105*tmp_43 + tmp_106*tmp_44 - 1.0/4.0;
      real_t tmp_110 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_111 = tmp_19*(0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_113 = tmp_110*tmp_9 + tmp_111*tmp_25 + tmp_112*tmp_31 - 1.0/4.0;
      real_t tmp_114 = tmp_110*tmp_38 + tmp_111*tmp_39 + tmp_112*tmp_40 - 1.0/4.0;
      real_t tmp_115 = tmp_110*tmp_42 + tmp_111*tmp_43 + tmp_112*tmp_44 - 1.0/4.0;
      real_t tmp_116 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_117 = tmp_19*(0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29);
      real_t tmp_118 = tmp_19*(0.40446199974765351*tmp_33 + 0.19107600050469298*tmp_34 + tmp_35);
      real_t tmp_119 = tmp_116*tmp_9 + tmp_117*tmp_25 + tmp_118*tmp_31 - 1.0/4.0;
      real_t tmp_120 = tmp_116*tmp_38 + tmp_117*tmp_39 + tmp_118*tmp_40 - 1.0/4.0;
      real_t tmp_121 = tmp_116*tmp_42 + tmp_117*tmp_43 + tmp_118*tmp_44 - 1.0/4.0;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_19*(0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_124 = tmp_19*(0.039308471900058539*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_125 = tmp_122*tmp_9 + tmp_123*tmp_25 + tmp_124*tmp_31 - 1.0/4.0;
      real_t tmp_126 = tmp_122*tmp_38 + tmp_123*tmp_39 + tmp_124*tmp_40 - 1.0/4.0;
      real_t tmp_127 = tmp_122*tmp_42 + tmp_123*tmp_43 + tmp_124*tmp_44 - 1.0/4.0;
      real_t tmp_128 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_129 = tmp_19*(0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_130 = tmp_19*(0.93718850182767688*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_131 = tmp_128*tmp_9 + tmp_129*tmp_25 + tmp_130*tmp_31 - 1.0/4.0;
      real_t tmp_132 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40 - 1.0/4.0;
      real_t tmp_133 = tmp_128*tmp_42 + tmp_129*tmp_43 + tmp_130*tmp_44 - 1.0/4.0;
      real_t tmp_134 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_135 = tmp_19*(0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_136 = tmp_19*(0.60796128279561268*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_137 = tmp_134*tmp_9 + tmp_135*tmp_25 + tmp_136*tmp_31 - 1.0/4.0;
      real_t tmp_138 = tmp_134*tmp_38 + tmp_135*tmp_39 + tmp_136*tmp_40 - 1.0/4.0;
      real_t tmp_139 = tmp_134*tmp_42 + tmp_135*tmp_43 + tmp_136*tmp_44 - 1.0/4.0;
      real_t tmp_140 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_19*(0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_142 = tmp_19*(0.19107600050469298*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_143 = tmp_140*tmp_9 + tmp_141*tmp_25 + tmp_142*tmp_31 - 1.0/4.0;
      real_t tmp_144 = tmp_140*tmp_38 + tmp_141*tmp_39 + tmp_142*tmp_40 - 1.0/4.0;
      real_t tmp_145 = tmp_140*tmp_42 + tmp_141*tmp_43 + tmp_142*tmp_44 - 1.0/4.0;
      real_t tmp_146 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_147 = tmp_19*(0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_148 = tmp_19*(0.031405749086161582*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_149 = tmp_146*tmp_9 + tmp_147*tmp_25 + tmp_148*tmp_31 - 1.0/4.0;
      real_t tmp_150 = tmp_146*tmp_38 + tmp_147*tmp_39 + tmp_148*tmp_40 - 1.0/4.0;
      real_t tmp_151 = tmp_146*tmp_42 + tmp_147*tmp_43 + tmp_148*tmp_44 - 1.0/4.0;
      real_t tmp_152 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_153 = tmp_19*(0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_154 = tmp_19*(0.19601935860219369*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_155 = tmp_152*tmp_9 + tmp_153*tmp_25 + tmp_154*tmp_31 - 1.0/4.0;
      real_t tmp_156 = tmp_152*tmp_38 + tmp_153*tmp_39 + tmp_154*tmp_40 - 1.0/4.0;
      real_t tmp_157 = tmp_152*tmp_42 + tmp_153*tmp_43 + tmp_154*tmp_44 - 1.0/4.0;
      real_t tmp_158 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_159 = tmp_19*(0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_160 = tmp_19*(0.40446199974765351*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_161 = tmp_158*tmp_9 + tmp_159*tmp_25 + tmp_160*tmp_31 - 1.0/4.0;
      real_t tmp_162 = tmp_158*tmp_38 + tmp_159*tmp_39 + tmp_160*tmp_40 - 1.0/4.0;
      real_t tmp_163 = tmp_158*tmp_42 + tmp_159*tmp_43 + tmp_160*tmp_44 - 1.0/4.0;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_165 = tmp_19*(0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_166 = tmp_19*(0.1711304259088916*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_167 = tmp_164*tmp_9 + tmp_165*tmp_25 + tmp_166*tmp_31 - 1.0/4.0;
      real_t tmp_168 = tmp_164*tmp_38 + tmp_165*tmp_39 + tmp_166*tmp_40 - 1.0/4.0;
      real_t tmp_169 = tmp_164*tmp_42 + tmp_165*tmp_43 + tmp_166*tmp_44 - 1.0/4.0;
      real_t a_0_0 = 0.019202922745021479*tmp_49*(((tmp_1*tmp_101 + tmp_102*tmp_2 + tmp_103*tmp_6)*(tmp_1*tmp_101 + tmp_102*tmp_2 + tmp_103*tmp_6)) + ((tmp_101*tmp_13 + tmp_102*tmp_15 + tmp_103*tmp_11)*(tmp_101*tmp_13 + tmp_102*tmp_15 + tmp_103*tmp_11)) + ((tmp_101*tmp_14 + tmp_102*tmp_7 + tmp_103*tmp_4)*(tmp_101*tmp_14 + tmp_102*tmp_7 + tmp_103*tmp_4))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_107 + tmp_108*tmp_2 + tmp_109*tmp_6)*(tmp_1*tmp_107 + tmp_108*tmp_2 + tmp_109*tmp_6)) + ((tmp_107*tmp_13 + tmp_108*tmp_15 + tmp_109*tmp_11)*(tmp_107*tmp_13 + tmp_108*tmp_15 + tmp_109*tmp_11)) + ((tmp_107*tmp_14 + tmp_108*tmp_7 + tmp_109*tmp_4)*(tmp_107*tmp_14 + tmp_108*tmp_7 + tmp_109*tmp_4))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_113 + tmp_114*tmp_2 + tmp_115*tmp_6)*(tmp_1*tmp_113 + tmp_114*tmp_2 + tmp_115*tmp_6)) + ((tmp_11*tmp_115 + tmp_113*tmp_13 + tmp_114*tmp_15)*(tmp_11*tmp_115 + tmp_113*tmp_13 + tmp_114*tmp_15)) + ((tmp_113*tmp_14 + tmp_114*tmp_7 + tmp_115*tmp_4)*(tmp_113*tmp_14 + tmp_114*tmp_7 + tmp_115*tmp_4))) + 0.042507265838595799*tmp_49*(((tmp_1*tmp_119 + tmp_120*tmp_2 + tmp_121*tmp_6)*(tmp_1*tmp_119 + tmp_120*tmp_2 + tmp_121*tmp_6)) + ((tmp_11*tmp_121 + tmp_119*tmp_13 + tmp_120*tmp_15)*(tmp_11*tmp_121 + tmp_119*tmp_13 + tmp_120*tmp_15)) + ((tmp_119*tmp_14 + tmp_120*tmp_7 + tmp_121*tmp_4)*(tmp_119*tmp_14 + tmp_120*tmp_7 + tmp_121*tmp_4))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_125 + tmp_126*tmp_2 + tmp_127*tmp_6)*(tmp_1*tmp_125 + tmp_126*tmp_2 + tmp_127*tmp_6)) + ((tmp_11*tmp_127 + tmp_125*tmp_13 + tmp_126*tmp_15)*(tmp_11*tmp_127 + tmp_125*tmp_13 + tmp_126*tmp_15)) + ((tmp_125*tmp_14 + tmp_126*tmp_7 + tmp_127*tmp_4)*(tmp_125*tmp_14 + tmp_126*tmp_7 + tmp_127*tmp_4))) + 0.0068572537431980923*tmp_49*(((tmp_1*tmp_131 + tmp_132*tmp_2 + tmp_133*tmp_6)*(tmp_1*tmp_131 + tmp_132*tmp_2 + tmp_133*tmp_6)) + ((tmp_11*tmp_133 + tmp_13*tmp_131 + tmp_132*tmp_15)*(tmp_11*tmp_133 + tmp_13*tmp_131 + tmp_132*tmp_15)) + ((tmp_131*tmp_14 + tmp_132*tmp_7 + tmp_133*tmp_4)*(tmp_131*tmp_14 + tmp_132*tmp_7 + tmp_133*tmp_4))) + 0.037198804536718075*tmp_49*(((tmp_1*tmp_137 + tmp_138*tmp_2 + tmp_139*tmp_6)*(tmp_1*tmp_137 + tmp_138*tmp_2 + tmp_139*tmp_6)) + ((tmp_11*tmp_139 + tmp_13*tmp_137 + tmp_138*tmp_15)*(tmp_11*tmp_139 + tmp_13*tmp_137 + tmp_138*tmp_15)) + ((tmp_137*tmp_14 + tmp_138*tmp_7 + tmp_139*tmp_4)*(tmp_137*tmp_14 + tmp_138*tmp_7 + tmp_139*tmp_4))) + 0.042507265838595799*tmp_49*(((tmp_1*tmp_143 + tmp_144*tmp_2 + tmp_145*tmp_6)*(tmp_1*tmp_143 + tmp_144*tmp_2 + tmp_145*tmp_6)) + ((tmp_11*tmp_145 + tmp_13*tmp_143 + tmp_144*tmp_15)*(tmp_11*tmp_145 + tmp_13*tmp_143 + tmp_144*tmp_15)) + ((tmp_14*tmp_143 + tmp_144*tmp_7 + tmp_145*tmp_4)*(tmp_14*tmp_143 + tmp_144*tmp_7 + tmp_145*tmp_4))) + 0.0068572537431980923*tmp_49*(((tmp_1*tmp_149 + tmp_150*tmp_2 + tmp_151*tmp_6)*(tmp_1*tmp_149 + tmp_150*tmp_2 + tmp_151*tmp_6)) + ((tmp_11*tmp_151 + tmp_13*tmp_149 + tmp_15*tmp_150)*(tmp_11*tmp_151 + tmp_13*tmp_149 + tmp_15*tmp_150)) + ((tmp_14*tmp_149 + tmp_150*tmp_7 + tmp_151*tmp_4)*(tmp_14*tmp_149 + tmp_150*tmp_7 + tmp_151*tmp_4))) + 0.037198804536718075*tmp_49*(((tmp_1*tmp_155 + tmp_156*tmp_2 + tmp_157*tmp_6)*(tmp_1*tmp_155 + tmp_156*tmp_2 + tmp_157*tmp_6)) + ((tmp_11*tmp_157 + tmp_13*tmp_155 + tmp_15*tmp_156)*(tmp_11*tmp_157 + tmp_13*tmp_155 + tmp_15*tmp_156)) + ((tmp_14*tmp_155 + tmp_156*tmp_7 + tmp_157*tmp_4)*(tmp_14*tmp_155 + tmp_156*tmp_7 + tmp_157*tmp_4))) + 0.042507265838595799*tmp_49*(((tmp_1*tmp_161 + tmp_162*tmp_2 + tmp_163*tmp_6)*(tmp_1*tmp_161 + tmp_162*tmp_2 + tmp_163*tmp_6)) + ((tmp_11*tmp_163 + tmp_13*tmp_161 + tmp_15*tmp_162)*(tmp_11*tmp_163 + tmp_13*tmp_161 + tmp_15*tmp_162)) + ((tmp_14*tmp_161 + tmp_162*tmp_7 + tmp_163*tmp_4)*(tmp_14*tmp_161 + tmp_162*tmp_7 + tmp_163*tmp_4))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_167 + tmp_168*tmp_2 + tmp_169*tmp_6)*(tmp_1*tmp_167 + tmp_168*tmp_2 + tmp_169*tmp_6)) + ((tmp_11*tmp_169 + tmp_13*tmp_167 + tmp_15*tmp_168)*(tmp_11*tmp_169 + tmp_13*tmp_167 + tmp_15*tmp_168)) + ((tmp_14*tmp_167 + tmp_168*tmp_7 + tmp_169*tmp_4)*(tmp_14*tmp_167 + tmp_168*tmp_7 + tmp_169*tmp_4))) + 0.0068572537431980923*tmp_49*(((tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6)*(tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6)) + ((tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41)*(tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41)) + ((tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7)*(tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7))) + 0.037198804536718075*tmp_49*(((tmp_1*tmp_53 + tmp_2*tmp_54 + tmp_55*tmp_6)*(tmp_1*tmp_53 + tmp_2*tmp_54 + tmp_55*tmp_6)) + ((tmp_11*tmp_55 + tmp_13*tmp_53 + tmp_15*tmp_54)*(tmp_11*tmp_55 + tmp_13*tmp_53 + tmp_15*tmp_54)) + ((tmp_14*tmp_53 + tmp_4*tmp_55 + tmp_54*tmp_7)*(tmp_14*tmp_53 + tmp_4*tmp_55 + tmp_54*tmp_7))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_59 + tmp_2*tmp_60 + tmp_6*tmp_61)*(tmp_1*tmp_59 + tmp_2*tmp_60 + tmp_6*tmp_61)) + ((tmp_11*tmp_61 + tmp_13*tmp_59 + tmp_15*tmp_60)*(tmp_11*tmp_61 + tmp_13*tmp_59 + tmp_15*tmp_60)) + ((tmp_14*tmp_59 + tmp_4*tmp_61 + tmp_60*tmp_7)*(tmp_14*tmp_59 + tmp_4*tmp_61 + tmp_60*tmp_7))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_65 + tmp_2*tmp_66 + tmp_6*tmp_67)*(tmp_1*tmp_65 + tmp_2*tmp_66 + tmp_6*tmp_67)) + ((tmp_11*tmp_67 + tmp_13*tmp_65 + tmp_15*tmp_66)*(tmp_11*tmp_67 + tmp_13*tmp_65 + tmp_15*tmp_66)) + ((tmp_14*tmp_65 + tmp_4*tmp_67 + tmp_66*tmp_7)*(tmp_14*tmp_65 + tmp_4*tmp_67 + tmp_66*tmp_7))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_71 + tmp_2*tmp_72 + tmp_6*tmp_73)*(tmp_1*tmp_71 + tmp_2*tmp_72 + tmp_6*tmp_73)) + ((tmp_11*tmp_73 + tmp_13*tmp_71 + tmp_15*tmp_72)*(tmp_11*tmp_73 + tmp_13*tmp_71 + tmp_15*tmp_72)) + ((tmp_14*tmp_71 + tmp_4*tmp_73 + tmp_7*tmp_72)*(tmp_14*tmp_71 + tmp_4*tmp_73 + tmp_7*tmp_72))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_77 + tmp_2*tmp_78 + tmp_6*tmp_79)*(tmp_1*tmp_77 + tmp_2*tmp_78 + tmp_6*tmp_79)) + ((tmp_11*tmp_79 + tmp_13*tmp_77 + tmp_15*tmp_78)*(tmp_11*tmp_79 + tmp_13*tmp_77 + tmp_15*tmp_78)) + ((tmp_14*tmp_77 + tmp_4*tmp_79 + tmp_7*tmp_78)*(tmp_14*tmp_77 + tmp_4*tmp_79 + tmp_7*tmp_78))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_83 + tmp_2*tmp_84 + tmp_6*tmp_85)*(tmp_1*tmp_83 + tmp_2*tmp_84 + tmp_6*tmp_85)) + ((tmp_11*tmp_85 + tmp_13*tmp_83 + tmp_15*tmp_84)*(tmp_11*tmp_85 + tmp_13*tmp_83 + tmp_15*tmp_84)) + ((tmp_14*tmp_83 + tmp_4*tmp_85 + tmp_7*tmp_84)*(tmp_14*tmp_83 + tmp_4*tmp_85 + tmp_7*tmp_84))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_89 + tmp_2*tmp_90 + tmp_6*tmp_91)*(tmp_1*tmp_89 + tmp_2*tmp_90 + tmp_6*tmp_91)) + ((tmp_11*tmp_91 + tmp_13*tmp_89 + tmp_15*tmp_90)*(tmp_11*tmp_91 + tmp_13*tmp_89 + tmp_15*tmp_90)) + ((tmp_14*tmp_89 + tmp_4*tmp_91 + tmp_7*tmp_90)*(tmp_14*tmp_89 + tmp_4*tmp_91 + tmp_7*tmp_90))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_95 + tmp_2*tmp_96 + tmp_6*tmp_97)*(tmp_1*tmp_95 + tmp_2*tmp_96 + tmp_6*tmp_97)) + ((tmp_11*tmp_97 + tmp_13*tmp_95 + tmp_15*tmp_96)*(tmp_11*tmp_97 + tmp_13*tmp_95 + tmp_15*tmp_96)) + ((tmp_14*tmp_95 + tmp_4*tmp_97 + tmp_7*tmp_96)*(tmp_14*tmp_95 + tmp_4*tmp_97 + tmp_7*tmp_96)));
      elMat( 0, 0) = a_0_0;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_2_0 + tmp_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_3_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_3_0 + tmp_0;
      real_t tmp_7 = p_affine_2_1 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11*tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_15*tmp_4;
      real_t tmp_18 = tmp_11*tmp_2;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_12 - tmp_1*tmp_17 + tmp_13*tmp_5 - tmp_13*tmp_8 + tmp_14*tmp_16 - tmp_14*tmp_18);
      real_t tmp_20 = p_affine_8_2 + tmp_10;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582*tmp_22 + 0.93718850182767688*tmp_23;
      real_t tmp_25 = tmp_19*(tmp_20 + tmp_24);
      real_t tmp_26 = tmp_16 - tmp_18;
      real_t tmp_27 = p_affine_8_1 + tmp_3;
      real_t tmp_28 = -p_affine_8_1;
      real_t tmp_29 = p_affine_9_1 + tmp_28;
      real_t tmp_30 = p_affine_10_1 + tmp_28;
      real_t tmp_31 = 0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30;
      real_t tmp_32 = tmp_19*(tmp_27 + tmp_31);
      real_t tmp_33 = tmp_12 - tmp_17;
      real_t tmp_34 = p_affine_8_0 + tmp_0;
      real_t tmp_35 = -p_affine_8_0;
      real_t tmp_36 = p_affine_9_0 + tmp_35;
      real_t tmp_37 = p_affine_10_0 + tmp_35;
      real_t tmp_38 = 0.031405749086161582*tmp_36 + 0.93718850182767688*tmp_37;
      real_t tmp_39 = tmp_19*(tmp_34 + tmp_38);
      real_t tmp_40 = tmp_25*tmp_9 + tmp_26*tmp_32 + tmp_33*tmp_39 - 1.0/4.0;
      real_t tmp_41 = -tmp_1*tmp_4 + tmp_14*tmp_6;
      real_t tmp_42 = tmp_1*tmp_11 - tmp_13*tmp_6;
      real_t tmp_43 = -tmp_11*tmp_14 + tmp_13*tmp_4;
      real_t tmp_44 = tmp_25*tmp_41 + tmp_32*tmp_42 + tmp_39*tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_1*tmp_7 - tmp_14*tmp_2;
      real_t tmp_46 = -tmp_1*tmp_15 + tmp_13*tmp_2;
      real_t tmp_47 = -tmp_13*tmp_7 + tmp_14*tmp_15;
      real_t tmp_48 = tmp_25*tmp_45 + tmp_32*tmp_46 + tmp_39*tmp_47 - 1.0/4.0;
      real_t tmp_49 = -p_affine_4_0;
      real_t tmp_50 = p_affine_5_0 + tmp_49;
      real_t tmp_51 = p_affine_6_0 + tmp_49;
      real_t tmp_52 = -p_affine_4_1;
      real_t tmp_53 = p_affine_7_1 + tmp_52;
      real_t tmp_54 = tmp_51*tmp_53;
      real_t tmp_55 = p_affine_7_0 + tmp_49;
      real_t tmp_56 = p_affine_6_1 + tmp_52;
      real_t tmp_57 = tmp_55*tmp_56;
      real_t tmp_58 = tmp_54 - tmp_57;
      real_t tmp_59 = -p_affine_4_2;
      real_t tmp_60 = p_affine_7_2 + tmp_59;
      real_t tmp_61 = tmp_56*tmp_60;
      real_t tmp_62 = p_affine_5_2 + tmp_59;
      real_t tmp_63 = p_affine_5_1 + tmp_52;
      real_t tmp_64 = p_affine_6_2 + tmp_59;
      real_t tmp_65 = tmp_55*tmp_64;
      real_t tmp_66 = tmp_53*tmp_64;
      real_t tmp_67 = tmp_51*tmp_60;
      real_t tmp_68 = 1.0 / (tmp_50*tmp_61 - tmp_50*tmp_66 + tmp_54*tmp_62 - tmp_57*tmp_62 + tmp_63*tmp_65 - tmp_63*tmp_67);
      real_t tmp_69 = p_affine_8_2 + tmp_59;
      real_t tmp_70 = tmp_68*(tmp_24 + tmp_69);
      real_t tmp_71 = tmp_65 - tmp_67;
      real_t tmp_72 = p_affine_8_1 + tmp_52;
      real_t tmp_73 = tmp_68*(tmp_31 + tmp_72);
      real_t tmp_74 = tmp_61 - tmp_66;
      real_t tmp_75 = p_affine_8_0 + tmp_49;
      real_t tmp_76 = tmp_68*(tmp_38 + tmp_75);
      real_t tmp_77 = tmp_58*tmp_70 + tmp_71*tmp_73 + tmp_74*tmp_76 - 1.0/4.0;
      real_t tmp_78 = -tmp_50*tmp_53 + tmp_55*tmp_63;
      real_t tmp_79 = tmp_50*tmp_60 - tmp_55*tmp_62;
      real_t tmp_80 = tmp_53*tmp_62 - tmp_60*tmp_63;
      real_t tmp_81 = tmp_70*tmp_78 + tmp_73*tmp_79 + tmp_76*tmp_80 - 1.0/4.0;
      real_t tmp_82 = tmp_50*tmp_56 - tmp_51*tmp_63;
      real_t tmp_83 = -tmp_50*tmp_64 + tmp_51*tmp_62;
      real_t tmp_84 = -tmp_56*tmp_62 + tmp_63*tmp_64;
      real_t tmp_85 = tmp_70*tmp_82 + tmp_73*tmp_83 + tmp_76*tmp_84 - 1.0/4.0;
      real_t tmp_86 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_87 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_88 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_89 = 3.0*std::pow((std::abs(tmp_23*tmp_86 - tmp_30*tmp_88)*std::abs(tmp_23*tmp_86 - tmp_30*tmp_88)) + (std::abs(tmp_23*tmp_87 - tmp_37*tmp_88)*std::abs(tmp_23*tmp_87 - tmp_37*tmp_88)) + (std::abs(tmp_30*tmp_87 - tmp_37*tmp_86)*std::abs(tmp_30*tmp_87 - tmp_37*tmp_86)), 0.25);
      real_t tmp_90 = 0.19601935860219369*tmp_22 + 0.60796128279561268*tmp_23;
      real_t tmp_91 = tmp_19*(tmp_20 + tmp_90);
      real_t tmp_92 = 0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30;
      real_t tmp_93 = tmp_19*(tmp_27 + tmp_92);
      real_t tmp_94 = 0.19601935860219369*tmp_36 + 0.60796128279561268*tmp_37;
      real_t tmp_95 = tmp_19*(tmp_34 + tmp_94);
      real_t tmp_96 = tmp_26*tmp_93 + tmp_33*tmp_95 + tmp_9*tmp_91 - 1.0/4.0;
      real_t tmp_97 = tmp_41*tmp_91 + tmp_42*tmp_93 + tmp_43*tmp_95 - 1.0/4.0;
      real_t tmp_98 = tmp_45*tmp_91 + tmp_46*tmp_93 + tmp_47*tmp_95 - 1.0/4.0;
      real_t tmp_99 = tmp_68*(tmp_69 + tmp_90);
      real_t tmp_100 = tmp_68*(tmp_72 + tmp_92);
      real_t tmp_101 = tmp_68*(tmp_75 + tmp_94);
      real_t tmp_102 = tmp_100*tmp_71 + tmp_101*tmp_74 + tmp_58*tmp_99 - 1.0/4.0;
      real_t tmp_103 = tmp_100*tmp_79 + tmp_101*tmp_80 + tmp_78*tmp_99 - 1.0/4.0;
      real_t tmp_104 = tmp_100*tmp_83 + tmp_101*tmp_84 + tmp_82*tmp_99 - 1.0/4.0;
      real_t tmp_105 = 0.37605877282253791*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_106 = tmp_19*(tmp_105 + tmp_20);
      real_t tmp_107 = 0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_108 = tmp_19*(tmp_107 + tmp_27);
      real_t tmp_109 = 0.37605877282253791*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_110 = tmp_19*(tmp_109 + tmp_34);
      real_t tmp_111 = tmp_106*tmp_9 + tmp_108*tmp_26 + tmp_110*tmp_33 - 1.0/4.0;
      real_t tmp_112 = tmp_106*tmp_41 + tmp_108*tmp_42 + tmp_110*tmp_43 - 1.0/4.0;
      real_t tmp_113 = tmp_106*tmp_45 + tmp_108*tmp_46 + tmp_110*tmp_47 - 1.0/4.0;
      real_t tmp_114 = tmp_68*(tmp_105 + tmp_69);
      real_t tmp_115 = tmp_68*(tmp_107 + tmp_72);
      real_t tmp_116 = tmp_68*(tmp_109 + tmp_75);
      real_t tmp_117 = tmp_114*tmp_58 + tmp_115*tmp_71 + tmp_116*tmp_74 - 1.0/4.0;
      real_t tmp_118 = tmp_114*tmp_78 + tmp_115*tmp_79 + tmp_116*tmp_80 - 1.0/4.0;
      real_t tmp_119 = tmp_114*tmp_82 + tmp_115*tmp_83 + tmp_116*tmp_84 - 1.0/4.0;
      real_t tmp_120 = 0.78764240869137092*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_121 = tmp_19*(tmp_120 + tmp_20);
      real_t tmp_122 = 0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_123 = tmp_19*(tmp_122 + tmp_27);
      real_t tmp_124 = 0.78764240869137092*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_125 = tmp_19*(tmp_124 + tmp_34);
      real_t tmp_126 = tmp_121*tmp_9 + tmp_123*tmp_26 + tmp_125*tmp_33 - 1.0/4.0;
      real_t tmp_127 = tmp_121*tmp_41 + tmp_123*tmp_42 + tmp_125*tmp_43 - 1.0/4.0;
      real_t tmp_128 = tmp_121*tmp_45 + tmp_123*tmp_46 + tmp_125*tmp_47 - 1.0/4.0;
      real_t tmp_129 = tmp_68*(tmp_120 + tmp_69);
      real_t tmp_130 = tmp_68*(tmp_122 + tmp_72);
      real_t tmp_131 = tmp_68*(tmp_124 + tmp_75);
      real_t tmp_132 = tmp_129*tmp_58 + tmp_130*tmp_71 + tmp_131*tmp_74 - 1.0/4.0;
      real_t tmp_133 = tmp_129*tmp_78 + tmp_130*tmp_79 + tmp_131*tmp_80 - 1.0/4.0;
      real_t tmp_134 = tmp_129*tmp_82 + tmp_130*tmp_83 + tmp_131*tmp_84 - 1.0/4.0;
      real_t tmp_135 = 0.58463275527740355*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_136 = tmp_19*(tmp_135 + tmp_20);
      real_t tmp_137 = 0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_138 = tmp_19*(tmp_137 + tmp_27);
      real_t tmp_139 = 0.58463275527740355*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_140 = tmp_19*(tmp_139 + tmp_34);
      real_t tmp_141 = tmp_136*tmp_9 + tmp_138*tmp_26 + tmp_140*tmp_33 - 1.0/4.0;
      real_t tmp_142 = tmp_136*tmp_41 + tmp_138*tmp_42 + tmp_140*tmp_43 - 1.0/4.0;
      real_t tmp_143 = tmp_136*tmp_45 + tmp_138*tmp_46 + tmp_140*tmp_47 - 1.0/4.0;
      real_t tmp_144 = tmp_68*(tmp_135 + tmp_69);
      real_t tmp_145 = tmp_68*(tmp_137 + tmp_72);
      real_t tmp_146 = tmp_68*(tmp_139 + tmp_75);
      real_t tmp_147 = tmp_144*tmp_58 + tmp_145*tmp_71 + tmp_146*tmp_74 - 1.0/4.0;
      real_t tmp_148 = tmp_144*tmp_78 + tmp_145*tmp_79 + tmp_146*tmp_80 - 1.0/4.0;
      real_t tmp_149 = tmp_144*tmp_82 + tmp_145*tmp_83 + tmp_146*tmp_84 - 1.0/4.0;
      real_t tmp_150 = 0.041227165399737475*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_151 = tmp_19*(tmp_150 + tmp_20);
      real_t tmp_152 = 0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_153 = tmp_19*(tmp_152 + tmp_27);
      real_t tmp_154 = 0.041227165399737475*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_155 = tmp_19*(tmp_154 + tmp_34);
      real_t tmp_156 = tmp_151*tmp_9 + tmp_153*tmp_26 + tmp_155*tmp_33 - 1.0/4.0;
      real_t tmp_157 = tmp_151*tmp_41 + tmp_153*tmp_42 + tmp_155*tmp_43 - 1.0/4.0;
      real_t tmp_158 = tmp_151*tmp_45 + tmp_153*tmp_46 + tmp_155*tmp_47 - 1.0/4.0;
      real_t tmp_159 = tmp_68*(tmp_150 + tmp_69);
      real_t tmp_160 = tmp_68*(tmp_152 + tmp_72);
      real_t tmp_161 = tmp_68*(tmp_154 + tmp_75);
      real_t tmp_162 = tmp_159*tmp_58 + tmp_160*tmp_71 + tmp_161*tmp_74 - 1.0/4.0;
      real_t tmp_163 = tmp_159*tmp_78 + tmp_160*tmp_79 + tmp_161*tmp_80 - 1.0/4.0;
      real_t tmp_164 = tmp_159*tmp_82 + tmp_160*tmp_83 + tmp_161*tmp_84 - 1.0/4.0;
      real_t tmp_165 = 0.039308471900058539*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_166 = tmp_19*(tmp_165 + tmp_20);
      real_t tmp_167 = 0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_168 = tmp_19*(tmp_167 + tmp_27);
      real_t tmp_169 = 0.039308471900058539*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_170 = tmp_19*(tmp_169 + tmp_34);
      real_t tmp_171 = tmp_166*tmp_9 + tmp_168*tmp_26 + tmp_170*tmp_33 - 1.0/4.0;
      real_t tmp_172 = tmp_166*tmp_41 + tmp_168*tmp_42 + tmp_170*tmp_43 - 1.0/4.0;
      real_t tmp_173 = tmp_166*tmp_45 + tmp_168*tmp_46 + tmp_170*tmp_47 - 1.0/4.0;
      real_t tmp_174 = tmp_68*(tmp_165 + tmp_69);
      real_t tmp_175 = tmp_68*(tmp_167 + tmp_72);
      real_t tmp_176 = tmp_68*(tmp_169 + tmp_75);
      real_t tmp_177 = tmp_174*tmp_58 + tmp_175*tmp_71 + tmp_176*tmp_74 - 1.0/4.0;
      real_t tmp_178 = tmp_174*tmp_78 + tmp_175*tmp_79 + tmp_176*tmp_80 - 1.0/4.0;
      real_t tmp_179 = tmp_174*tmp_82 + tmp_175*tmp_83 + tmp_176*tmp_84 - 1.0/4.0;
      real_t tmp_180 = 0.78764240869137092*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_181 = tmp_19*(tmp_180 + tmp_20);
      real_t tmp_182 = 0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_183 = tmp_19*(tmp_182 + tmp_27);
      real_t tmp_184 = 0.78764240869137092*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_185 = tmp_19*(tmp_184 + tmp_34);
      real_t tmp_186 = tmp_181*tmp_9 + tmp_183*tmp_26 + tmp_185*tmp_33 - 1.0/4.0;
      real_t tmp_187 = tmp_181*tmp_41 + tmp_183*tmp_42 + tmp_185*tmp_43 - 1.0/4.0;
      real_t tmp_188 = tmp_181*tmp_45 + tmp_183*tmp_46 + tmp_185*tmp_47 - 1.0/4.0;
      real_t tmp_189 = tmp_68*(tmp_180 + tmp_69);
      real_t tmp_190 = tmp_68*(tmp_182 + tmp_72);
      real_t tmp_191 = tmp_68*(tmp_184 + tmp_75);
      real_t tmp_192 = tmp_189*tmp_58 + tmp_190*tmp_71 + tmp_191*tmp_74 - 1.0/4.0;
      real_t tmp_193 = tmp_189*tmp_78 + tmp_190*tmp_79 + tmp_191*tmp_80 - 1.0/4.0;
      real_t tmp_194 = tmp_189*tmp_82 + tmp_190*tmp_83 + tmp_191*tmp_84 - 1.0/4.0;
      real_t tmp_195 = 0.58463275527740355*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_196 = tmp_19*(tmp_195 + tmp_20);
      real_t tmp_197 = 0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_198 = tmp_19*(tmp_197 + tmp_27);
      real_t tmp_199 = 0.58463275527740355*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_200 = tmp_19*(tmp_199 + tmp_34);
      real_t tmp_201 = tmp_196*tmp_9 + tmp_198*tmp_26 + tmp_200*tmp_33 - 1.0/4.0;
      real_t tmp_202 = tmp_196*tmp_41 + tmp_198*tmp_42 + tmp_200*tmp_43 - 1.0/4.0;
      real_t tmp_203 = tmp_196*tmp_45 + tmp_198*tmp_46 + tmp_200*tmp_47 - 1.0/4.0;
      real_t tmp_204 = tmp_68*(tmp_195 + tmp_69);
      real_t tmp_205 = tmp_68*(tmp_197 + tmp_72);
      real_t tmp_206 = tmp_68*(tmp_199 + tmp_75);
      real_t tmp_207 = tmp_204*tmp_58 + tmp_205*tmp_71 + tmp_206*tmp_74 - 1.0/4.0;
      real_t tmp_208 = tmp_204*tmp_78 + tmp_205*tmp_79 + tmp_206*tmp_80 - 1.0/4.0;
      real_t tmp_209 = tmp_204*tmp_82 + tmp_205*tmp_83 + tmp_206*tmp_84 - 1.0/4.0;
      real_t tmp_210 = 0.1711304259088916*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_211 = tmp_19*(tmp_20 + tmp_210);
      real_t tmp_212 = 0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_213 = tmp_19*(tmp_212 + tmp_27);
      real_t tmp_214 = 0.1711304259088916*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_215 = tmp_19*(tmp_214 + tmp_34);
      real_t tmp_216 = tmp_211*tmp_9 + tmp_213*tmp_26 + tmp_215*tmp_33 - 1.0/4.0;
      real_t tmp_217 = tmp_211*tmp_41 + tmp_213*tmp_42 + tmp_215*tmp_43 - 1.0/4.0;
      real_t tmp_218 = tmp_211*tmp_45 + tmp_213*tmp_46 + tmp_215*tmp_47 - 1.0/4.0;
      real_t tmp_219 = tmp_68*(tmp_210 + tmp_69);
      real_t tmp_220 = tmp_68*(tmp_212 + tmp_72);
      real_t tmp_221 = tmp_68*(tmp_214 + tmp_75);
      real_t tmp_222 = tmp_219*tmp_58 + tmp_220*tmp_71 + tmp_221*tmp_74 - 1.0/4.0;
      real_t tmp_223 = tmp_219*tmp_78 + tmp_220*tmp_79 + tmp_221*tmp_80 - 1.0/4.0;
      real_t tmp_224 = tmp_219*tmp_82 + tmp_220*tmp_83 + tmp_221*tmp_84 - 1.0/4.0;
      real_t tmp_225 = 0.37605877282253791*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_226 = tmp_19*(tmp_20 + tmp_225);
      real_t tmp_227 = 0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_228 = tmp_19*(tmp_227 + tmp_27);
      real_t tmp_229 = 0.37605877282253791*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_230 = tmp_19*(tmp_229 + tmp_34);
      real_t tmp_231 = tmp_226*tmp_9 + tmp_228*tmp_26 + tmp_230*tmp_33 - 1.0/4.0;
      real_t tmp_232 = tmp_226*tmp_41 + tmp_228*tmp_42 + tmp_230*tmp_43 - 1.0/4.0;
      real_t tmp_233 = tmp_226*tmp_45 + tmp_228*tmp_46 + tmp_230*tmp_47 - 1.0/4.0;
      real_t tmp_234 = tmp_68*(tmp_225 + tmp_69);
      real_t tmp_235 = tmp_68*(tmp_227 + tmp_72);
      real_t tmp_236 = tmp_68*(tmp_229 + tmp_75);
      real_t tmp_237 = tmp_234*tmp_58 + tmp_235*tmp_71 + tmp_236*tmp_74 - 1.0/4.0;
      real_t tmp_238 = tmp_234*tmp_78 + tmp_235*tmp_79 + tmp_236*tmp_80 - 1.0/4.0;
      real_t tmp_239 = tmp_234*tmp_82 + tmp_235*tmp_83 + tmp_236*tmp_84 - 1.0/4.0;
      real_t tmp_240 = 0.041227165399737475*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_241 = tmp_19*(tmp_20 + tmp_240);
      real_t tmp_242 = 0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_243 = tmp_19*(tmp_242 + tmp_27);
      real_t tmp_244 = 0.041227165399737475*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_245 = tmp_19*(tmp_244 + tmp_34);
      real_t tmp_246 = tmp_241*tmp_9 + tmp_243*tmp_26 + tmp_245*tmp_33 - 1.0/4.0;
      real_t tmp_247 = tmp_241*tmp_41 + tmp_243*tmp_42 + tmp_245*tmp_43 - 1.0/4.0;
      real_t tmp_248 = tmp_241*tmp_45 + tmp_243*tmp_46 + tmp_245*tmp_47 - 1.0/4.0;
      real_t tmp_249 = tmp_68*(tmp_240 + tmp_69);
      real_t tmp_250 = tmp_68*(tmp_242 + tmp_72);
      real_t tmp_251 = tmp_68*(tmp_244 + tmp_75);
      real_t tmp_252 = tmp_249*tmp_58 + tmp_250*tmp_71 + tmp_251*tmp_74 - 1.0/4.0;
      real_t tmp_253 = tmp_249*tmp_78 + tmp_250*tmp_79 + tmp_251*tmp_80 - 1.0/4.0;
      real_t tmp_254 = tmp_249*tmp_82 + tmp_250*tmp_83 + tmp_251*tmp_84 - 1.0/4.0;
      real_t tmp_255 = 0.40446199974765351*tmp_22 + 0.19107600050469298*tmp_23;
      real_t tmp_256 = tmp_19*(tmp_20 + tmp_255);
      real_t tmp_257 = 0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30;
      real_t tmp_258 = tmp_19*(tmp_257 + tmp_27);
      real_t tmp_259 = 0.40446199974765351*tmp_36 + 0.19107600050469298*tmp_37;
      real_t tmp_260 = tmp_19*(tmp_259 + tmp_34);
      real_t tmp_261 = tmp_256*tmp_9 + tmp_258*tmp_26 + tmp_260*tmp_33 - 1.0/4.0;
      real_t tmp_262 = tmp_256*tmp_41 + tmp_258*tmp_42 + tmp_260*tmp_43 - 1.0/4.0;
      real_t tmp_263 = tmp_256*tmp_45 + tmp_258*tmp_46 + tmp_260*tmp_47 - 1.0/4.0;
      real_t tmp_264 = tmp_68*(tmp_255 + tmp_69);
      real_t tmp_265 = tmp_68*(tmp_257 + tmp_72);
      real_t tmp_266 = tmp_68*(tmp_259 + tmp_75);
      real_t tmp_267 = tmp_264*tmp_58 + tmp_265*tmp_71 + tmp_266*tmp_74 - 1.0/4.0;
      real_t tmp_268 = tmp_264*tmp_78 + tmp_265*tmp_79 + tmp_266*tmp_80 - 1.0/4.0;
      real_t tmp_269 = tmp_264*tmp_82 + tmp_265*tmp_83 + tmp_266*tmp_84 - 1.0/4.0;
      real_t tmp_270 = 0.039308471900058539*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_271 = tmp_19*(tmp_20 + tmp_270);
      real_t tmp_272 = 0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_273 = tmp_19*(tmp_27 + tmp_272);
      real_t tmp_274 = 0.039308471900058539*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_275 = tmp_19*(tmp_274 + tmp_34);
      real_t tmp_276 = tmp_26*tmp_273 + tmp_271*tmp_9 + tmp_275*tmp_33 - 1.0/4.0;
      real_t tmp_277 = tmp_271*tmp_41 + tmp_273*tmp_42 + tmp_275*tmp_43 - 1.0/4.0;
      real_t tmp_278 = tmp_271*tmp_45 + tmp_273*tmp_46 + tmp_275*tmp_47 - 1.0/4.0;
      real_t tmp_279 = tmp_68*(tmp_270 + tmp_69);
      real_t tmp_280 = tmp_68*(tmp_272 + tmp_72);
      real_t tmp_281 = tmp_68*(tmp_274 + tmp_75);
      real_t tmp_282 = tmp_279*tmp_58 + tmp_280*tmp_71 + tmp_281*tmp_74 - 1.0/4.0;
      real_t tmp_283 = tmp_279*tmp_78 + tmp_280*tmp_79 + tmp_281*tmp_80 - 1.0/4.0;
      real_t tmp_284 = tmp_279*tmp_82 + tmp_280*tmp_83 + tmp_281*tmp_84 - 1.0/4.0;
      real_t tmp_285 = 0.93718850182767688*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_286 = tmp_19*(tmp_20 + tmp_285);
      real_t tmp_287 = 0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_288 = tmp_19*(tmp_27 + tmp_287);
      real_t tmp_289 = 0.93718850182767688*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_290 = tmp_19*(tmp_289 + tmp_34);
      real_t tmp_291 = tmp_26*tmp_288 + tmp_286*tmp_9 + tmp_290*tmp_33 - 1.0/4.0;
      real_t tmp_292 = tmp_286*tmp_41 + tmp_288*tmp_42 + tmp_290*tmp_43 - 1.0/4.0;
      real_t tmp_293 = tmp_286*tmp_45 + tmp_288*tmp_46 + tmp_290*tmp_47 - 1.0/4.0;
      real_t tmp_294 = tmp_68*(tmp_285 + tmp_69);
      real_t tmp_295 = tmp_68*(tmp_287 + tmp_72);
      real_t tmp_296 = tmp_68*(tmp_289 + tmp_75);
      real_t tmp_297 = tmp_294*tmp_58 + tmp_295*tmp_71 + tmp_296*tmp_74 - 1.0/4.0;
      real_t tmp_298 = tmp_294*tmp_78 + tmp_295*tmp_79 + tmp_296*tmp_80 - 1.0/4.0;
      real_t tmp_299 = tmp_294*tmp_82 + tmp_295*tmp_83 + tmp_296*tmp_84 - 1.0/4.0;
      real_t tmp_300 = 0.60796128279561268*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_301 = tmp_19*(tmp_20 + tmp_300);
      real_t tmp_302 = 0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_303 = tmp_19*(tmp_27 + tmp_302);
      real_t tmp_304 = 0.60796128279561268*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_305 = tmp_19*(tmp_304 + tmp_34);
      real_t tmp_306 = tmp_26*tmp_303 + tmp_301*tmp_9 + tmp_305*tmp_33 - 1.0/4.0;
      real_t tmp_307 = tmp_301*tmp_41 + tmp_303*tmp_42 + tmp_305*tmp_43 - 1.0/4.0;
      real_t tmp_308 = tmp_301*tmp_45 + tmp_303*tmp_46 + tmp_305*tmp_47 - 1.0/4.0;
      real_t tmp_309 = tmp_68*(tmp_300 + tmp_69);
      real_t tmp_310 = tmp_68*(tmp_302 + tmp_72);
      real_t tmp_311 = tmp_68*(tmp_304 + tmp_75);
      real_t tmp_312 = tmp_309*tmp_58 + tmp_310*tmp_71 + tmp_311*tmp_74 - 1.0/4.0;
      real_t tmp_313 = tmp_309*tmp_78 + tmp_310*tmp_79 + tmp_311*tmp_80 - 1.0/4.0;
      real_t tmp_314 = tmp_309*tmp_82 + tmp_310*tmp_83 + tmp_311*tmp_84 - 1.0/4.0;
      real_t tmp_315 = 0.19107600050469298*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_316 = tmp_19*(tmp_20 + tmp_315);
      real_t tmp_317 = 0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_318 = tmp_19*(tmp_27 + tmp_317);
      real_t tmp_319 = 0.19107600050469298*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_320 = tmp_19*(tmp_319 + tmp_34);
      real_t tmp_321 = tmp_26*tmp_318 + tmp_316*tmp_9 + tmp_320*tmp_33 - 1.0/4.0;
      real_t tmp_322 = tmp_316*tmp_41 + tmp_318*tmp_42 + tmp_320*tmp_43 - 1.0/4.0;
      real_t tmp_323 = tmp_316*tmp_45 + tmp_318*tmp_46 + tmp_320*tmp_47 - 1.0/4.0;
      real_t tmp_324 = tmp_68*(tmp_315 + tmp_69);
      real_t tmp_325 = tmp_68*(tmp_317 + tmp_72);
      real_t tmp_326 = tmp_68*(tmp_319 + tmp_75);
      real_t tmp_327 = tmp_324*tmp_58 + tmp_325*tmp_71 + tmp_326*tmp_74 - 1.0/4.0;
      real_t tmp_328 = tmp_324*tmp_78 + tmp_325*tmp_79 + tmp_326*tmp_80 - 1.0/4.0;
      real_t tmp_329 = tmp_324*tmp_82 + tmp_325*tmp_83 + tmp_326*tmp_84 - 1.0/4.0;
      real_t tmp_330 = 0.031405749086161582*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_331 = tmp_19*(tmp_20 + tmp_330);
      real_t tmp_332 = 0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_333 = tmp_19*(tmp_27 + tmp_332);
      real_t tmp_334 = 0.031405749086161582*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_335 = tmp_19*(tmp_334 + tmp_34);
      real_t tmp_336 = tmp_26*tmp_333 + tmp_33*tmp_335 + tmp_331*tmp_9 - 1.0/4.0;
      real_t tmp_337 = tmp_331*tmp_41 + tmp_333*tmp_42 + tmp_335*tmp_43 - 1.0/4.0;
      real_t tmp_338 = tmp_331*tmp_45 + tmp_333*tmp_46 + tmp_335*tmp_47 - 1.0/4.0;
      real_t tmp_339 = tmp_68*(tmp_330 + tmp_69);
      real_t tmp_340 = tmp_68*(tmp_332 + tmp_72);
      real_t tmp_341 = tmp_68*(tmp_334 + tmp_75);
      real_t tmp_342 = tmp_339*tmp_58 + tmp_340*tmp_71 + tmp_341*tmp_74 - 1.0/4.0;
      real_t tmp_343 = tmp_339*tmp_78 + tmp_340*tmp_79 + tmp_341*tmp_80 - 1.0/4.0;
      real_t tmp_344 = tmp_339*tmp_82 + tmp_340*tmp_83 + tmp_341*tmp_84 - 1.0/4.0;
      real_t tmp_345 = 0.19601935860219369*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_346 = tmp_19*(tmp_20 + tmp_345);
      real_t tmp_347 = 0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_348 = tmp_19*(tmp_27 + tmp_347);
      real_t tmp_349 = 0.19601935860219369*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_350 = tmp_19*(tmp_34 + tmp_349);
      real_t tmp_351 = tmp_26*tmp_348 + tmp_33*tmp_350 + tmp_346*tmp_9 - 1.0/4.0;
      real_t tmp_352 = tmp_346*tmp_41 + tmp_348*tmp_42 + tmp_350*tmp_43 - 1.0/4.0;
      real_t tmp_353 = tmp_346*tmp_45 + tmp_348*tmp_46 + tmp_350*tmp_47 - 1.0/4.0;
      real_t tmp_354 = tmp_68*(tmp_345 + tmp_69);
      real_t tmp_355 = tmp_68*(tmp_347 + tmp_72);
      real_t tmp_356 = tmp_68*(tmp_349 + tmp_75);
      real_t tmp_357 = tmp_354*tmp_58 + tmp_355*tmp_71 + tmp_356*tmp_74 - 1.0/4.0;
      real_t tmp_358 = tmp_354*tmp_78 + tmp_355*tmp_79 + tmp_356*tmp_80 - 1.0/4.0;
      real_t tmp_359 = tmp_354*tmp_82 + tmp_355*tmp_83 + tmp_356*tmp_84 - 1.0/4.0;
      real_t tmp_360 = 0.40446199974765351*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_361 = tmp_19*(tmp_20 + tmp_360);
      real_t tmp_362 = 0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_363 = tmp_19*(tmp_27 + tmp_362);
      real_t tmp_364 = 0.40446199974765351*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_365 = tmp_19*(tmp_34 + tmp_364);
      real_t tmp_366 = tmp_26*tmp_363 + tmp_33*tmp_365 + tmp_361*tmp_9 - 1.0/4.0;
      real_t tmp_367 = tmp_361*tmp_41 + tmp_363*tmp_42 + tmp_365*tmp_43 - 1.0/4.0;
      real_t tmp_368 = tmp_361*tmp_45 + tmp_363*tmp_46 + tmp_365*tmp_47 - 1.0/4.0;
      real_t tmp_369 = tmp_68*(tmp_360 + tmp_69);
      real_t tmp_370 = tmp_68*(tmp_362 + tmp_72);
      real_t tmp_371 = tmp_68*(tmp_364 + tmp_75);
      real_t tmp_372 = tmp_369*tmp_58 + tmp_370*tmp_71 + tmp_371*tmp_74 - 1.0/4.0;
      real_t tmp_373 = tmp_369*tmp_78 + tmp_370*tmp_79 + tmp_371*tmp_80 - 1.0/4.0;
      real_t tmp_374 = tmp_369*tmp_82 + tmp_370*tmp_83 + tmp_371*tmp_84 - 1.0/4.0;
      real_t tmp_375 = 0.1711304259088916*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_376 = tmp_19*(tmp_20 + tmp_375);
      real_t tmp_377 = 0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_378 = tmp_19*(tmp_27 + tmp_377);
      real_t tmp_379 = 0.1711304259088916*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_380 = tmp_19*(tmp_34 + tmp_379);
      real_t tmp_381 = tmp_26*tmp_378 + tmp_33*tmp_380 + tmp_376*tmp_9 - 1.0/4.0;
      real_t tmp_382 = tmp_376*tmp_41 + tmp_378*tmp_42 + tmp_380*tmp_43 - 1.0/4.0;
      real_t tmp_383 = tmp_376*tmp_45 + tmp_378*tmp_46 + tmp_380*tmp_47 - 1.0/4.0;
      real_t tmp_384 = tmp_68*(tmp_375 + tmp_69);
      real_t tmp_385 = tmp_68*(tmp_377 + tmp_72);
      real_t tmp_386 = tmp_68*(tmp_379 + tmp_75);
      real_t tmp_387 = tmp_384*tmp_58 + tmp_385*tmp_71 + tmp_386*tmp_74 - 1.0/4.0;
      real_t tmp_388 = tmp_384*tmp_78 + tmp_385*tmp_79 + tmp_386*tmp_80 - 1.0/4.0;
      real_t tmp_389 = tmp_384*tmp_82 + tmp_385*tmp_83 + tmp_386*tmp_84 - 1.0/4.0;
      real_t a_0_0 = -0.020848748529055869*tmp_89*((tmp_1*tmp_111 + tmp_112*tmp_2 + tmp_113*tmp_6)*(tmp_117*tmp_50 + tmp_118*tmp_51 + tmp_119*tmp_55) + (tmp_11*tmp_113 + tmp_111*tmp_13 + tmp_112*tmp_15)*(tmp_117*tmp_62 + tmp_118*tmp_64 + tmp_119*tmp_60) + (tmp_111*tmp_14 + tmp_112*tmp_7 + tmp_113*tmp_4)*(tmp_117*tmp_63 + tmp_118*tmp_56 + tmp_119*tmp_53)) - 0.019202922745021479*tmp_89*((tmp_1*tmp_126 + tmp_127*tmp_2 + tmp_128*tmp_6)*(tmp_132*tmp_50 + tmp_133*tmp_51 + tmp_134*tmp_55) + (tmp_11*tmp_128 + tmp_126*tmp_13 + tmp_127*tmp_15)*(tmp_132*tmp_62 + tmp_133*tmp_64 + tmp_134*tmp_60) + (tmp_126*tmp_14 + tmp_127*tmp_7 + tmp_128*tmp_4)*(tmp_132*tmp_63 + tmp_133*tmp_56 + tmp_134*tmp_53)) - 0.020848748529055869*tmp_89*((tmp_1*tmp_141 + tmp_142*tmp_2 + tmp_143*tmp_6)*(tmp_147*tmp_50 + tmp_148*tmp_51 + tmp_149*tmp_55) + (tmp_11*tmp_143 + tmp_13*tmp_141 + tmp_142*tmp_15)*(tmp_147*tmp_62 + tmp_148*tmp_64 + tmp_149*tmp_60) + (tmp_14*tmp_141 + tmp_142*tmp_7 + tmp_143*tmp_4)*(tmp_147*tmp_63 + tmp_148*tmp_56 + tmp_149*tmp_53)) - 0.019202922745021479*tmp_89*((tmp_1*tmp_156 + tmp_157*tmp_2 + tmp_158*tmp_6)*(tmp_162*tmp_50 + tmp_163*tmp_51 + tmp_164*tmp_55) + (tmp_11*tmp_158 + tmp_13*tmp_156 + tmp_15*tmp_157)*(tmp_162*tmp_62 + tmp_163*tmp_64 + tmp_164*tmp_60) + (tmp_14*tmp_156 + tmp_157*tmp_7 + tmp_158*tmp_4)*(tmp_162*tmp_63 + tmp_163*tmp_56 + tmp_164*tmp_53)) - 0.020848748529055869*tmp_89*((tmp_1*tmp_171 + tmp_172*tmp_2 + tmp_173*tmp_6)*(tmp_177*tmp_50 + tmp_178*tmp_51 + tmp_179*tmp_55) + (tmp_11*tmp_173 + tmp_13*tmp_171 + tmp_15*tmp_172)*(tmp_177*tmp_62 + tmp_178*tmp_64 + tmp_179*tmp_60) + (tmp_14*tmp_171 + tmp_172*tmp_7 + tmp_173*tmp_4)*(tmp_177*tmp_63 + tmp_178*tmp_56 + tmp_179*tmp_53)) - 0.019202922745021479*tmp_89*((tmp_1*tmp_186 + tmp_187*tmp_2 + tmp_188*tmp_6)*(tmp_192*tmp_50 + tmp_193*tmp_51 + tmp_194*tmp_55) + (tmp_11*tmp_188 + tmp_13*tmp_186 + tmp_15*tmp_187)*(tmp_192*tmp_62 + tmp_193*tmp_64 + tmp_194*tmp_60) + (tmp_14*tmp_186 + tmp_187*tmp_7 + tmp_188*tmp_4)*(tmp_192*tmp_63 + tmp_193*tmp_56 + tmp_194*tmp_53)) - 0.020848748529055869*tmp_89*((tmp_1*tmp_201 + tmp_2*tmp_202 + tmp_203*tmp_6)*(tmp_207*tmp_50 + tmp_208*tmp_51 + tmp_209*tmp_55) + (tmp_11*tmp_203 + tmp_13*tmp_201 + tmp_15*tmp_202)*(tmp_207*tmp_62 + tmp_208*tmp_64 + tmp_209*tmp_60) + (tmp_14*tmp_201 + tmp_202*tmp_7 + tmp_203*tmp_4)*(tmp_207*tmp_63 + tmp_208*tmp_56 + tmp_209*tmp_53)) - 0.019202922745021479*tmp_89*((tmp_1*tmp_216 + tmp_2*tmp_217 + tmp_218*tmp_6)*(tmp_222*tmp_50 + tmp_223*tmp_51 + tmp_224*tmp_55) + (tmp_11*tmp_218 + tmp_13*tmp_216 + tmp_15*tmp_217)*(tmp_222*tmp_62 + tmp_223*tmp_64 + tmp_224*tmp_60) + (tmp_14*tmp_216 + tmp_217*tmp_7 + tmp_218*tmp_4)*(tmp_222*tmp_63 + tmp_223*tmp_56 + tmp_224*tmp_53)) - 0.020848748529055869*tmp_89*((tmp_1*tmp_231 + tmp_2*tmp_232 + tmp_233*tmp_6)*(tmp_237*tmp_50 + tmp_238*tmp_51 + tmp_239*tmp_55) + (tmp_11*tmp_233 + tmp_13*tmp_231 + tmp_15*tmp_232)*(tmp_237*tmp_62 + tmp_238*tmp_64 + tmp_239*tmp_60) + (tmp_14*tmp_231 + tmp_232*tmp_7 + tmp_233*tmp_4)*(tmp_237*tmp_63 + tmp_238*tmp_56 + tmp_239*tmp_53)) - 0.019202922745021479*tmp_89*((tmp_1*tmp_246 + tmp_2*tmp_247 + tmp_248*tmp_6)*(tmp_252*tmp_50 + tmp_253*tmp_51 + tmp_254*tmp_55) + (tmp_11*tmp_248 + tmp_13*tmp_246 + tmp_15*tmp_247)*(tmp_252*tmp_62 + tmp_253*tmp_64 + tmp_254*tmp_60) + (tmp_14*tmp_246 + tmp_247*tmp_7 + tmp_248*tmp_4)*(tmp_252*tmp_63 + tmp_253*tmp_56 + tmp_254*tmp_53)) - 0.042507265838595799*tmp_89*((tmp_1*tmp_261 + tmp_2*tmp_262 + tmp_263*tmp_6)*(tmp_267*tmp_50 + tmp_268*tmp_51 + tmp_269*tmp_55) + (tmp_11*tmp_263 + tmp_13*tmp_261 + tmp_15*tmp_262)*(tmp_267*tmp_62 + tmp_268*tmp_64 + tmp_269*tmp_60) + (tmp_14*tmp_261 + tmp_262*tmp_7 + tmp_263*tmp_4)*(tmp_267*tmp_63 + tmp_268*tmp_56 + tmp_269*tmp_53)) - 0.020848748529055869*tmp_89*((tmp_1*tmp_276 + tmp_2*tmp_277 + tmp_278*tmp_6)*(tmp_282*tmp_50 + tmp_283*tmp_51 + tmp_284*tmp_55) + (tmp_11*tmp_278 + tmp_13*tmp_276 + tmp_15*tmp_277)*(tmp_282*tmp_62 + tmp_283*tmp_64 + tmp_284*tmp_60) + (tmp_14*tmp_276 + tmp_277*tmp_7 + tmp_278*tmp_4)*(tmp_282*tmp_63 + tmp_283*tmp_56 + tmp_284*tmp_53)) - 0.0068572537431980923*tmp_89*((tmp_1*tmp_291 + tmp_2*tmp_292 + tmp_293*tmp_6)*(tmp_297*tmp_50 + tmp_298*tmp_51 + tmp_299*tmp_55) + (tmp_11*tmp_293 + tmp_13*tmp_291 + tmp_15*tmp_292)*(tmp_297*tmp_62 + tmp_298*tmp_64 + tmp_299*tmp_60) + (tmp_14*tmp_291 + tmp_292*tmp_7 + tmp_293*tmp_4)*(tmp_297*tmp_63 + tmp_298*tmp_56 + tmp_299*tmp_53)) - 0.037198804536718075*tmp_89*((tmp_1*tmp_306 + tmp_2*tmp_307 + tmp_308*tmp_6)*(tmp_312*tmp_50 + tmp_313*tmp_51 + tmp_314*tmp_55) + (tmp_11*tmp_308 + tmp_13*tmp_306 + tmp_15*tmp_307)*(tmp_312*tmp_62 + tmp_313*tmp_64 + tmp_314*tmp_60) + (tmp_14*tmp_306 + tmp_307*tmp_7 + tmp_308*tmp_4)*(tmp_312*tmp_63 + tmp_313*tmp_56 + tmp_314*tmp_53)) - 0.042507265838595799*tmp_89*((tmp_1*tmp_321 + tmp_2*tmp_322 + tmp_323*tmp_6)*(tmp_327*tmp_50 + tmp_328*tmp_51 + tmp_329*tmp_55) + (tmp_11*tmp_323 + tmp_13*tmp_321 + tmp_15*tmp_322)*(tmp_327*tmp_62 + tmp_328*tmp_64 + tmp_329*tmp_60) + (tmp_14*tmp_321 + tmp_322*tmp_7 + tmp_323*tmp_4)*(tmp_327*tmp_63 + tmp_328*tmp_56 + tmp_329*tmp_53)) - 0.0068572537431980923*tmp_89*((tmp_1*tmp_336 + tmp_2*tmp_337 + tmp_338*tmp_6)*(tmp_342*tmp_50 + tmp_343*tmp_51 + tmp_344*tmp_55) + (tmp_11*tmp_338 + tmp_13*tmp_336 + tmp_15*tmp_337)*(tmp_342*tmp_62 + tmp_343*tmp_64 + tmp_344*tmp_60) + (tmp_14*tmp_336 + tmp_337*tmp_7 + tmp_338*tmp_4)*(tmp_342*tmp_63 + tmp_343*tmp_56 + tmp_344*tmp_53)) - 0.037198804536718075*tmp_89*((tmp_1*tmp_351 + tmp_2*tmp_352 + tmp_353*tmp_6)*(tmp_357*tmp_50 + tmp_358*tmp_51 + tmp_359*tmp_55) + (tmp_11*tmp_353 + tmp_13*tmp_351 + tmp_15*tmp_352)*(tmp_357*tmp_62 + tmp_358*tmp_64 + tmp_359*tmp_60) + (tmp_14*tmp_351 + tmp_352*tmp_7 + tmp_353*tmp_4)*(tmp_357*tmp_63 + tmp_358*tmp_56 + tmp_359*tmp_53)) - 0.042507265838595799*tmp_89*((tmp_1*tmp_366 + tmp_2*tmp_367 + tmp_368*tmp_6)*(tmp_372*tmp_50 + tmp_373*tmp_51 + tmp_374*tmp_55) + (tmp_11*tmp_368 + tmp_13*tmp_366 + tmp_15*tmp_367)*(tmp_372*tmp_62 + tmp_373*tmp_64 + tmp_374*tmp_60) + (tmp_14*tmp_366 + tmp_367*tmp_7 + tmp_368*tmp_4)*(tmp_372*tmp_63 + tmp_373*tmp_56 + tmp_374*tmp_53)) - 0.019202922745021479*tmp_89*((tmp_1*tmp_381 + tmp_2*tmp_382 + tmp_383*tmp_6)*(tmp_387*tmp_50 + tmp_388*tmp_51 + tmp_389*tmp_55) + (tmp_11*tmp_383 + tmp_13*tmp_381 + tmp_15*tmp_382)*(tmp_387*tmp_62 + tmp_388*tmp_64 + tmp_389*tmp_60) + (tmp_14*tmp_381 + tmp_382*tmp_7 + tmp_383*tmp_4)*(tmp_387*tmp_63 + tmp_388*tmp_56 + tmp_389*tmp_53)) - 0.0068572537431980923*tmp_89*((tmp_1*tmp_40 + tmp_2*tmp_44 + tmp_48*tmp_6)*(tmp_50*tmp_77 + tmp_51*tmp_81 + tmp_55*tmp_85) + (tmp_11*tmp_48 + tmp_13*tmp_40 + tmp_15*tmp_44)*(tmp_60*tmp_85 + tmp_62*tmp_77 + tmp_64*tmp_81) + (tmp_14*tmp_40 + tmp_4*tmp_48 + tmp_44*tmp_7)*(tmp_53*tmp_85 + tmp_56*tmp_81 + tmp_63*tmp_77)) - 0.037198804536718075*tmp_89*((tmp_1*tmp_96 + tmp_2*tmp_97 + tmp_6*tmp_98)*(tmp_102*tmp_50 + tmp_103*tmp_51 + tmp_104*tmp_55) + (tmp_102*tmp_62 + tmp_103*tmp_64 + tmp_104*tmp_60)*(tmp_11*tmp_98 + tmp_13*tmp_96 + tmp_15*tmp_97) + (tmp_102*tmp_63 + tmp_103*tmp_56 + tmp_104*tmp_53)*(tmp_14*tmp_96 + tmp_4*tmp_98 + tmp_7*tmp_97));
      elMat( 0, 0) = a_0_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_2_0 + tmp_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_3_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_3_0 + tmp_0;
      real_t tmp_7 = p_affine_2_1 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11*tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_15*tmp_4;
      real_t tmp_18 = tmp_11*tmp_2;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_12 - tmp_1*tmp_17 + tmp_13*tmp_5 - tmp_13*tmp_8 + tmp_14*tmp_16 - tmp_14*tmp_18);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_10;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_16 - tmp_18;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_3;
      real_t tmp_30 = tmp_19*(0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29);
      real_t tmp_31 = tmp_12 - tmp_17;
      real_t tmp_32 = -p_affine_8_0;
      real_t tmp_33 = p_affine_9_0 + tmp_32;
      real_t tmp_34 = p_affine_10_0 + tmp_32;
      real_t tmp_35 = p_affine_8_0 + tmp_0;
      real_t tmp_36 = tmp_19*(0.031405749086161582*tmp_33 + 0.93718850182767688*tmp_34 + tmp_35);
      real_t tmp_37 = tmp_24*tmp_9 + tmp_25*tmp_30 + tmp_31*tmp_36 - 1.0/4.0;
      real_t tmp_38 = -tmp_1*tmp_4 + tmp_14*tmp_6;
      real_t tmp_39 = tmp_1*tmp_11 - tmp_13*tmp_6;
      real_t tmp_40 = -tmp_11*tmp_14 + tmp_13*tmp_4;
      real_t tmp_41 = tmp_24*tmp_38 + tmp_30*tmp_39 + tmp_36*tmp_40 - 1.0/4.0;
      real_t tmp_42 = tmp_1*tmp_7 - tmp_14*tmp_2;
      real_t tmp_43 = -tmp_1*tmp_15 + tmp_13*tmp_2;
      real_t tmp_44 = -tmp_13*tmp_7 + tmp_14*tmp_15;
      real_t tmp_45 = tmp_24*tmp_42 + tmp_30*tmp_43 + tmp_36*tmp_44 - 1.0/4.0;
      real_t tmp_46 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_47 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_48 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_49 = 3.0*std::pow((std::abs(tmp_22*tmp_46 - tmp_28*tmp_48)*std::abs(tmp_22*tmp_46 - tmp_28*tmp_48)) + (std::abs(tmp_22*tmp_47 - tmp_34*tmp_48)*std::abs(tmp_22*tmp_47 - tmp_34*tmp_48)) + (std::abs(tmp_28*tmp_47 - tmp_34*tmp_46)*std::abs(tmp_28*tmp_47 - tmp_34*tmp_46)), 0.25);
      real_t tmp_50 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_51 = tmp_19*(0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29);
      real_t tmp_52 = tmp_19*(0.19601935860219369*tmp_33 + 0.60796128279561268*tmp_34 + tmp_35);
      real_t tmp_53 = tmp_25*tmp_51 + tmp_31*tmp_52 + tmp_50*tmp_9 - 1.0/4.0;
      real_t tmp_54 = tmp_38*tmp_50 + tmp_39*tmp_51 + tmp_40*tmp_52 - 1.0/4.0;
      real_t tmp_55 = tmp_42*tmp_50 + tmp_43*tmp_51 + tmp_44*tmp_52 - 1.0/4.0;
      real_t tmp_56 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_57 = tmp_19*(0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_58 = tmp_19*(0.37605877282253791*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_59 = tmp_25*tmp_57 + tmp_31*tmp_58 + tmp_56*tmp_9 - 1.0/4.0;
      real_t tmp_60 = tmp_38*tmp_56 + tmp_39*tmp_57 + tmp_40*tmp_58 - 1.0/4.0;
      real_t tmp_61 = tmp_42*tmp_56 + tmp_43*tmp_57 + tmp_44*tmp_58 - 1.0/4.0;
      real_t tmp_62 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_63 = tmp_19*(0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_64 = tmp_19*(0.78764240869137092*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_65 = tmp_25*tmp_63 + tmp_31*tmp_64 + tmp_62*tmp_9 - 1.0/4.0;
      real_t tmp_66 = tmp_38*tmp_62 + tmp_39*tmp_63 + tmp_40*tmp_64 - 1.0/4.0;
      real_t tmp_67 = tmp_42*tmp_62 + tmp_43*tmp_63 + tmp_44*tmp_64 - 1.0/4.0;
      real_t tmp_68 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_69 = tmp_19*(0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_70 = tmp_19*(0.58463275527740355*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_71 = tmp_25*tmp_69 + tmp_31*tmp_70 + tmp_68*tmp_9 - 1.0/4.0;
      real_t tmp_72 = tmp_38*tmp_68 + tmp_39*tmp_69 + tmp_40*tmp_70 - 1.0/4.0;
      real_t tmp_73 = tmp_42*tmp_68 + tmp_43*tmp_69 + tmp_44*tmp_70 - 1.0/4.0;
      real_t tmp_74 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_75 = tmp_19*(0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_76 = tmp_19*(0.041227165399737475*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_77 = tmp_25*tmp_75 + tmp_31*tmp_76 + tmp_74*tmp_9 - 1.0/4.0;
      real_t tmp_78 = tmp_38*tmp_74 + tmp_39*tmp_75 + tmp_40*tmp_76 - 1.0/4.0;
      real_t tmp_79 = tmp_42*tmp_74 + tmp_43*tmp_75 + tmp_44*tmp_76 - 1.0/4.0;
      real_t tmp_80 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_81 = tmp_19*(0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_82 = tmp_19*(0.039308471900058539*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_83 = tmp_25*tmp_81 + tmp_31*tmp_82 + tmp_80*tmp_9 - 1.0/4.0;
      real_t tmp_84 = tmp_38*tmp_80 + tmp_39*tmp_81 + tmp_40*tmp_82 - 1.0/4.0;
      real_t tmp_85 = tmp_42*tmp_80 + tmp_43*tmp_81 + tmp_44*tmp_82 - 1.0/4.0;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_87 = tmp_19*(0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_88 = tmp_19*(0.78764240869137092*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_89 = tmp_25*tmp_87 + tmp_31*tmp_88 + tmp_86*tmp_9 - 1.0/4.0;
      real_t tmp_90 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88 - 1.0/4.0;
      real_t tmp_91 = tmp_42*tmp_86 + tmp_43*tmp_87 + tmp_44*tmp_88 - 1.0/4.0;
      real_t tmp_92 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_19*(0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_94 = tmp_19*(0.58463275527740355*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_95 = tmp_25*tmp_93 + tmp_31*tmp_94 + tmp_9*tmp_92 - 1.0/4.0;
      real_t tmp_96 = tmp_38*tmp_92 + tmp_39*tmp_93 + tmp_40*tmp_94 - 1.0/4.0;
      real_t tmp_97 = tmp_42*tmp_92 + tmp_43*tmp_93 + tmp_44*tmp_94 - 1.0/4.0;
      real_t tmp_98 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_99 = tmp_19*(0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_100 = tmp_19*(0.1711304259088916*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_101 = tmp_100*tmp_31 + tmp_25*tmp_99 + tmp_9*tmp_98 - 1.0/4.0;
      real_t tmp_102 = tmp_100*tmp_40 + tmp_38*tmp_98 + tmp_39*tmp_99 - 1.0/4.0;
      real_t tmp_103 = tmp_100*tmp_44 + tmp_42*tmp_98 + tmp_43*tmp_99 - 1.0/4.0;
      real_t tmp_104 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_105 = tmp_19*(0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_106 = tmp_19*(0.37605877282253791*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_107 = tmp_104*tmp_9 + tmp_105*tmp_25 + tmp_106*tmp_31 - 1.0/4.0;
      real_t tmp_108 = tmp_104*tmp_38 + tmp_105*tmp_39 + tmp_106*tmp_40 - 1.0/4.0;
      real_t tmp_109 = tmp_104*tmp_42 + tmp_105*tmp_43 + tmp_106*tmp_44 - 1.0/4.0;
      real_t tmp_110 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_111 = tmp_19*(0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_113 = tmp_110*tmp_9 + tmp_111*tmp_25 + tmp_112*tmp_31 - 1.0/4.0;
      real_t tmp_114 = tmp_110*tmp_38 + tmp_111*tmp_39 + tmp_112*tmp_40 - 1.0/4.0;
      real_t tmp_115 = tmp_110*tmp_42 + tmp_111*tmp_43 + tmp_112*tmp_44 - 1.0/4.0;
      real_t tmp_116 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_117 = tmp_19*(0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29);
      real_t tmp_118 = tmp_19*(0.40446199974765351*tmp_33 + 0.19107600050469298*tmp_34 + tmp_35);
      real_t tmp_119 = tmp_116*tmp_9 + tmp_117*tmp_25 + tmp_118*tmp_31 - 1.0/4.0;
      real_t tmp_120 = tmp_116*tmp_38 + tmp_117*tmp_39 + tmp_118*tmp_40 - 1.0/4.0;
      real_t tmp_121 = tmp_116*tmp_42 + tmp_117*tmp_43 + tmp_118*tmp_44 - 1.0/4.0;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_19*(0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_124 = tmp_19*(0.039308471900058539*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_125 = tmp_122*tmp_9 + tmp_123*tmp_25 + tmp_124*tmp_31 - 1.0/4.0;
      real_t tmp_126 = tmp_122*tmp_38 + tmp_123*tmp_39 + tmp_124*tmp_40 - 1.0/4.0;
      real_t tmp_127 = tmp_122*tmp_42 + tmp_123*tmp_43 + tmp_124*tmp_44 - 1.0/4.0;
      real_t tmp_128 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_129 = tmp_19*(0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_130 = tmp_19*(0.93718850182767688*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_131 = tmp_128*tmp_9 + tmp_129*tmp_25 + tmp_130*tmp_31 - 1.0/4.0;
      real_t tmp_132 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40 - 1.0/4.0;
      real_t tmp_133 = tmp_128*tmp_42 + tmp_129*tmp_43 + tmp_130*tmp_44 - 1.0/4.0;
      real_t tmp_134 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_135 = tmp_19*(0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_136 = tmp_19*(0.60796128279561268*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_137 = tmp_134*tmp_9 + tmp_135*tmp_25 + tmp_136*tmp_31 - 1.0/4.0;
      real_t tmp_138 = tmp_134*tmp_38 + tmp_135*tmp_39 + tmp_136*tmp_40 - 1.0/4.0;
      real_t tmp_139 = tmp_134*tmp_42 + tmp_135*tmp_43 + tmp_136*tmp_44 - 1.0/4.0;
      real_t tmp_140 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_19*(0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_142 = tmp_19*(0.19107600050469298*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_143 = tmp_140*tmp_9 + tmp_141*tmp_25 + tmp_142*tmp_31 - 1.0/4.0;
      real_t tmp_144 = tmp_140*tmp_38 + tmp_141*tmp_39 + tmp_142*tmp_40 - 1.0/4.0;
      real_t tmp_145 = tmp_140*tmp_42 + tmp_141*tmp_43 + tmp_142*tmp_44 - 1.0/4.0;
      real_t tmp_146 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_147 = tmp_19*(0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_148 = tmp_19*(0.031405749086161582*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_149 = tmp_146*tmp_9 + tmp_147*tmp_25 + tmp_148*tmp_31 - 1.0/4.0;
      real_t tmp_150 = tmp_146*tmp_38 + tmp_147*tmp_39 + tmp_148*tmp_40 - 1.0/4.0;
      real_t tmp_151 = tmp_146*tmp_42 + tmp_147*tmp_43 + tmp_148*tmp_44 - 1.0/4.0;
      real_t tmp_152 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_153 = tmp_19*(0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_154 = tmp_19*(0.19601935860219369*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_155 = tmp_152*tmp_9 + tmp_153*tmp_25 + tmp_154*tmp_31 - 1.0/4.0;
      real_t tmp_156 = tmp_152*tmp_38 + tmp_153*tmp_39 + tmp_154*tmp_40 - 1.0/4.0;
      real_t tmp_157 = tmp_152*tmp_42 + tmp_153*tmp_43 + tmp_154*tmp_44 - 1.0/4.0;
      real_t tmp_158 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_159 = tmp_19*(0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_160 = tmp_19*(0.40446199974765351*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_161 = tmp_158*tmp_9 + tmp_159*tmp_25 + tmp_160*tmp_31 - 1.0/4.0;
      real_t tmp_162 = tmp_158*tmp_38 + tmp_159*tmp_39 + tmp_160*tmp_40 - 1.0/4.0;
      real_t tmp_163 = tmp_158*tmp_42 + tmp_159*tmp_43 + tmp_160*tmp_44 - 1.0/4.0;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_165 = tmp_19*(0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_166 = tmp_19*(0.1711304259088916*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_167 = tmp_164*tmp_9 + tmp_165*tmp_25 + tmp_166*tmp_31 - 1.0/4.0;
      real_t tmp_168 = tmp_164*tmp_38 + tmp_165*tmp_39 + tmp_166*tmp_40 - 1.0/4.0;
      real_t tmp_169 = tmp_164*tmp_42 + tmp_165*tmp_43 + tmp_166*tmp_44 - 1.0/4.0;
      real_t a_0_0 = 0.019202922745021479*tmp_49*(((tmp_1*tmp_101 + tmp_102*tmp_2 + tmp_103*tmp_6)*(tmp_1*tmp_101 + tmp_102*tmp_2 + tmp_103*tmp_6)) + ((tmp_101*tmp_13 + tmp_102*tmp_15 + tmp_103*tmp_11)*(tmp_101*tmp_13 + tmp_102*tmp_15 + tmp_103*tmp_11)) + ((tmp_101*tmp_14 + tmp_102*tmp_7 + tmp_103*tmp_4)*(tmp_101*tmp_14 + tmp_102*tmp_7 + tmp_103*tmp_4))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_107 + tmp_108*tmp_2 + tmp_109*tmp_6)*(tmp_1*tmp_107 + tmp_108*tmp_2 + tmp_109*tmp_6)) + ((tmp_107*tmp_13 + tmp_108*tmp_15 + tmp_109*tmp_11)*(tmp_107*tmp_13 + tmp_108*tmp_15 + tmp_109*tmp_11)) + ((tmp_107*tmp_14 + tmp_108*tmp_7 + tmp_109*tmp_4)*(tmp_107*tmp_14 + tmp_108*tmp_7 + tmp_109*tmp_4))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_113 + tmp_114*tmp_2 + tmp_115*tmp_6)*(tmp_1*tmp_113 + tmp_114*tmp_2 + tmp_115*tmp_6)) + ((tmp_11*tmp_115 + tmp_113*tmp_13 + tmp_114*tmp_15)*(tmp_11*tmp_115 + tmp_113*tmp_13 + tmp_114*tmp_15)) + ((tmp_113*tmp_14 + tmp_114*tmp_7 + tmp_115*tmp_4)*(tmp_113*tmp_14 + tmp_114*tmp_7 + tmp_115*tmp_4))) + 0.042507265838595799*tmp_49*(((tmp_1*tmp_119 + tmp_120*tmp_2 + tmp_121*tmp_6)*(tmp_1*tmp_119 + tmp_120*tmp_2 + tmp_121*tmp_6)) + ((tmp_11*tmp_121 + tmp_119*tmp_13 + tmp_120*tmp_15)*(tmp_11*tmp_121 + tmp_119*tmp_13 + tmp_120*tmp_15)) + ((tmp_119*tmp_14 + tmp_120*tmp_7 + tmp_121*tmp_4)*(tmp_119*tmp_14 + tmp_120*tmp_7 + tmp_121*tmp_4))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_125 + tmp_126*tmp_2 + tmp_127*tmp_6)*(tmp_1*tmp_125 + tmp_126*tmp_2 + tmp_127*tmp_6)) + ((tmp_11*tmp_127 + tmp_125*tmp_13 + tmp_126*tmp_15)*(tmp_11*tmp_127 + tmp_125*tmp_13 + tmp_126*tmp_15)) + ((tmp_125*tmp_14 + tmp_126*tmp_7 + tmp_127*tmp_4)*(tmp_125*tmp_14 + tmp_126*tmp_7 + tmp_127*tmp_4))) + 0.0068572537431980923*tmp_49*(((tmp_1*tmp_131 + tmp_132*tmp_2 + tmp_133*tmp_6)*(tmp_1*tmp_131 + tmp_132*tmp_2 + tmp_133*tmp_6)) + ((tmp_11*tmp_133 + tmp_13*tmp_131 + tmp_132*tmp_15)*(tmp_11*tmp_133 + tmp_13*tmp_131 + tmp_132*tmp_15)) + ((tmp_131*tmp_14 + tmp_132*tmp_7 + tmp_133*tmp_4)*(tmp_131*tmp_14 + tmp_132*tmp_7 + tmp_133*tmp_4))) + 0.037198804536718075*tmp_49*(((tmp_1*tmp_137 + tmp_138*tmp_2 + tmp_139*tmp_6)*(tmp_1*tmp_137 + tmp_138*tmp_2 + tmp_139*tmp_6)) + ((tmp_11*tmp_139 + tmp_13*tmp_137 + tmp_138*tmp_15)*(tmp_11*tmp_139 + tmp_13*tmp_137 + tmp_138*tmp_15)) + ((tmp_137*tmp_14 + tmp_138*tmp_7 + tmp_139*tmp_4)*(tmp_137*tmp_14 + tmp_138*tmp_7 + tmp_139*tmp_4))) + 0.042507265838595799*tmp_49*(((tmp_1*tmp_143 + tmp_144*tmp_2 + tmp_145*tmp_6)*(tmp_1*tmp_143 + tmp_144*tmp_2 + tmp_145*tmp_6)) + ((tmp_11*tmp_145 + tmp_13*tmp_143 + tmp_144*tmp_15)*(tmp_11*tmp_145 + tmp_13*tmp_143 + tmp_144*tmp_15)) + ((tmp_14*tmp_143 + tmp_144*tmp_7 + tmp_145*tmp_4)*(tmp_14*tmp_143 + tmp_144*tmp_7 + tmp_145*tmp_4))) + 0.0068572537431980923*tmp_49*(((tmp_1*tmp_149 + tmp_150*tmp_2 + tmp_151*tmp_6)*(tmp_1*tmp_149 + tmp_150*tmp_2 + tmp_151*tmp_6)) + ((tmp_11*tmp_151 + tmp_13*tmp_149 + tmp_15*tmp_150)*(tmp_11*tmp_151 + tmp_13*tmp_149 + tmp_15*tmp_150)) + ((tmp_14*tmp_149 + tmp_150*tmp_7 + tmp_151*tmp_4)*(tmp_14*tmp_149 + tmp_150*tmp_7 + tmp_151*tmp_4))) + 0.037198804536718075*tmp_49*(((tmp_1*tmp_155 + tmp_156*tmp_2 + tmp_157*tmp_6)*(tmp_1*tmp_155 + tmp_156*tmp_2 + tmp_157*tmp_6)) + ((tmp_11*tmp_157 + tmp_13*tmp_155 + tmp_15*tmp_156)*(tmp_11*tmp_157 + tmp_13*tmp_155 + tmp_15*tmp_156)) + ((tmp_14*tmp_155 + tmp_156*tmp_7 + tmp_157*tmp_4)*(tmp_14*tmp_155 + tmp_156*tmp_7 + tmp_157*tmp_4))) + 0.042507265838595799*tmp_49*(((tmp_1*tmp_161 + tmp_162*tmp_2 + tmp_163*tmp_6)*(tmp_1*tmp_161 + tmp_162*tmp_2 + tmp_163*tmp_6)) + ((tmp_11*tmp_163 + tmp_13*tmp_161 + tmp_15*tmp_162)*(tmp_11*tmp_163 + tmp_13*tmp_161 + tmp_15*tmp_162)) + ((tmp_14*tmp_161 + tmp_162*tmp_7 + tmp_163*tmp_4)*(tmp_14*tmp_161 + tmp_162*tmp_7 + tmp_163*tmp_4))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_167 + tmp_168*tmp_2 + tmp_169*tmp_6)*(tmp_1*tmp_167 + tmp_168*tmp_2 + tmp_169*tmp_6)) + ((tmp_11*tmp_169 + tmp_13*tmp_167 + tmp_15*tmp_168)*(tmp_11*tmp_169 + tmp_13*tmp_167 + tmp_15*tmp_168)) + ((tmp_14*tmp_167 + tmp_168*tmp_7 + tmp_169*tmp_4)*(tmp_14*tmp_167 + tmp_168*tmp_7 + tmp_169*tmp_4))) + 0.0068572537431980923*tmp_49*(((tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6)*(tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6)) + ((tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41)*(tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41)) + ((tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7)*(tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7))) + 0.037198804536718075*tmp_49*(((tmp_1*tmp_53 + tmp_2*tmp_54 + tmp_55*tmp_6)*(tmp_1*tmp_53 + tmp_2*tmp_54 + tmp_55*tmp_6)) + ((tmp_11*tmp_55 + tmp_13*tmp_53 + tmp_15*tmp_54)*(tmp_11*tmp_55 + tmp_13*tmp_53 + tmp_15*tmp_54)) + ((tmp_14*tmp_53 + tmp_4*tmp_55 + tmp_54*tmp_7)*(tmp_14*tmp_53 + tmp_4*tmp_55 + tmp_54*tmp_7))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_59 + tmp_2*tmp_60 + tmp_6*tmp_61)*(tmp_1*tmp_59 + tmp_2*tmp_60 + tmp_6*tmp_61)) + ((tmp_11*tmp_61 + tmp_13*tmp_59 + tmp_15*tmp_60)*(tmp_11*tmp_61 + tmp_13*tmp_59 + tmp_15*tmp_60)) + ((tmp_14*tmp_59 + tmp_4*tmp_61 + tmp_60*tmp_7)*(tmp_14*tmp_59 + tmp_4*tmp_61 + tmp_60*tmp_7))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_65 + tmp_2*tmp_66 + tmp_6*tmp_67)*(tmp_1*tmp_65 + tmp_2*tmp_66 + tmp_6*tmp_67)) + ((tmp_11*tmp_67 + tmp_13*tmp_65 + tmp_15*tmp_66)*(tmp_11*tmp_67 + tmp_13*tmp_65 + tmp_15*tmp_66)) + ((tmp_14*tmp_65 + tmp_4*tmp_67 + tmp_66*tmp_7)*(tmp_14*tmp_65 + tmp_4*tmp_67 + tmp_66*tmp_7))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_71 + tmp_2*tmp_72 + tmp_6*tmp_73)*(tmp_1*tmp_71 + tmp_2*tmp_72 + tmp_6*tmp_73)) + ((tmp_11*tmp_73 + tmp_13*tmp_71 + tmp_15*tmp_72)*(tmp_11*tmp_73 + tmp_13*tmp_71 + tmp_15*tmp_72)) + ((tmp_14*tmp_71 + tmp_4*tmp_73 + tmp_7*tmp_72)*(tmp_14*tmp_71 + tmp_4*tmp_73 + tmp_7*tmp_72))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_77 + tmp_2*tmp_78 + tmp_6*tmp_79)*(tmp_1*tmp_77 + tmp_2*tmp_78 + tmp_6*tmp_79)) + ((tmp_11*tmp_79 + tmp_13*tmp_77 + tmp_15*tmp_78)*(tmp_11*tmp_79 + tmp_13*tmp_77 + tmp_15*tmp_78)) + ((tmp_14*tmp_77 + tmp_4*tmp_79 + tmp_7*tmp_78)*(tmp_14*tmp_77 + tmp_4*tmp_79 + tmp_7*tmp_78))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_83 + tmp_2*tmp_84 + tmp_6*tmp_85)*(tmp_1*tmp_83 + tmp_2*tmp_84 + tmp_6*tmp_85)) + ((tmp_11*tmp_85 + tmp_13*tmp_83 + tmp_15*tmp_84)*(tmp_11*tmp_85 + tmp_13*tmp_83 + tmp_15*tmp_84)) + ((tmp_14*tmp_83 + tmp_4*tmp_85 + tmp_7*tmp_84)*(tmp_14*tmp_83 + tmp_4*tmp_85 + tmp_7*tmp_84))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_89 + tmp_2*tmp_90 + tmp_6*tmp_91)*(tmp_1*tmp_89 + tmp_2*tmp_90 + tmp_6*tmp_91)) + ((tmp_11*tmp_91 + tmp_13*tmp_89 + tmp_15*tmp_90)*(tmp_11*tmp_91 + tmp_13*tmp_89 + tmp_15*tmp_90)) + ((tmp_14*tmp_89 + tmp_4*tmp_91 + tmp_7*tmp_90)*(tmp_14*tmp_89 + tmp_4*tmp_91 + tmp_7*tmp_90))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_95 + tmp_2*tmp_96 + tmp_6*tmp_97)*(tmp_1*tmp_95 + tmp_2*tmp_96 + tmp_6*tmp_97)) + ((tmp_11*tmp_97 + tmp_13*tmp_95 + tmp_15*tmp_96)*(tmp_11*tmp_97 + tmp_13*tmp_95 + tmp_15*tmp_96)) + ((tmp_14*tmp_95 + tmp_4*tmp_97 + tmp_7*tmp_96)*(tmp_14*tmp_95 + tmp_4*tmp_97 + tmp_7*tmp_96)));
      elMat( 0, 0) = a_0_0;
   }


};


} //eg
} // dg
} // hyteg
