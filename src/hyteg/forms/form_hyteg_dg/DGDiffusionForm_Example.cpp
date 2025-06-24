/*
* Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr.
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

#include "hyteg/forms/form_hyteg_dg/DGDiffusionForm_Example.hpp"

namespace hyteg {
namespace dg {

using walberla::real_c;
using walberla::uint_c;

void DGDiffusionForm_Example::integrateVolume2D( const std::vector< Point3D >& coords,
                                                 const DGBasisInfo&            trialBasis,
                                                 const DGBasisInfo&            testBasis,
                                                 int                           trialDegree,
                                                 int                           testDegree,
                                                 MatrixXr&                     elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 2, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 2, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coords[0]( 0 );
   const auto p_affine_0_1 = coords[0]( 1 );

   const auto p_affine_1_0 = coords[1]( 0 );
   const auto p_affine_1_1 = coords[1]( 1 );

   const auto p_affine_2_0 = coords[2]( 0 );
   const auto p_affine_2_1 = coords[2]( 1 );

   real_t tmp_0  = p_affine_0_1 * p_affine_1_0;
   real_t tmp_1  = p_affine_0_0 * p_affine_1_1;
   real_t tmp_2  = 2 * tmp_1;
   real_t tmp_3  = 2 * p_affine_0_0;
   real_t tmp_4  = p_affine_1_0 * p_affine_2_1;
   real_t tmp_5  = p_affine_0_1 * tmp_4;
   real_t tmp_6  = p_affine_0_1 * p_affine_2_0;
   real_t tmp_7  = p_affine_0_0 * p_affine_2_1;
   real_t tmp_8  = 2 * tmp_6;
   real_t tmp_9  = p_affine_2_0 * p_affine_2_1;
   real_t tmp_10 = p_affine_1_0 * p_affine_1_1;
   real_t tmp_11 = p_affine_1_1 * p_affine_2_0;
   real_t tmp_12 = tmp_11 * tmp_4;
   real_t tmp_13 = ( p_affine_0_0 * p_affine_0_0 );
   real_t tmp_14 = ( p_affine_1_1 * p_affine_1_1 );
   real_t tmp_15 = tmp_13 * tmp_14;
   real_t tmp_16 = ( p_affine_2_1 * p_affine_2_1 );
   real_t tmp_17 = tmp_13 * tmp_16;
   real_t tmp_18 = ( p_affine_0_1 * p_affine_0_1 );
   real_t tmp_19 = ( p_affine_1_0 * p_affine_1_0 );
   real_t tmp_20 = tmp_18 * tmp_19;
   real_t tmp_21 = ( p_affine_2_0 * p_affine_2_0 );
   real_t tmp_22 = tmp_18 * tmp_21;
   real_t tmp_23 = tmp_16 * tmp_19;
   real_t tmp_24 = tmp_14 * tmp_21;
   real_t tmp_25 = p_affine_1_0 * tmp_16;
   real_t tmp_26 = p_affine_2_0 * tmp_14;
   real_t tmp_27 = p_affine_1_1 * p_affine_2_1;
   real_t tmp_28 = tmp_13 * tmp_27;
   real_t tmp_29 = 2 * p_affine_0_1;
   real_t tmp_30 = p_affine_2_1 * tmp_19;
   real_t tmp_31 = p_affine_1_1 * tmp_21;
   real_t tmp_32 = p_affine_1_0 * p_affine_2_0;
   real_t tmp_33 = tmp_18 * tmp_32;
   real_t tmp_34 = std::abs( tmp_0 - tmp_1 + tmp_11 - tmp_4 - tmp_6 + tmp_7 );
   real_t tmp_35 = tmp_34 / ( -tmp_0 * tmp_2 + tmp_10 * tmp_8 - 2 * tmp_12 + tmp_15 + tmp_17 + tmp_2 * tmp_4 + tmp_2 * tmp_6 +
                              tmp_2 * tmp_9 + tmp_20 + tmp_22 + tmp_23 + tmp_24 - tmp_25 * tmp_3 - tmp_26 * tmp_3 - 2 * tmp_28 -
                              tmp_29 * tmp_30 - tmp_29 * tmp_31 + tmp_3 * tmp_5 - 2 * tmp_33 + tmp_4 * tmp_8 - tmp_7 * tmp_8 );
   real_t tmp_36 = tmp_32 * tmp_35;
   real_t tmp_37 = tmp_27 * tmp_35;
   real_t tmp_38 = 4 * tmp_1;
   real_t tmp_39 = 4 * p_affine_0_0;
   real_t tmp_40 = 4 * tmp_6;
   real_t tmp_41 = 4 * p_affine_0_1;
   real_t tmp_42 = tmp_34 / ( -tmp_0 * tmp_38 + tmp_10 * tmp_40 - 4 * tmp_12 + 2 * tmp_15 + 2 * tmp_17 + 2 * tmp_20 + 2 * tmp_22 +
                              2 * tmp_23 + 2 * tmp_24 - tmp_25 * tmp_39 - tmp_26 * tmp_39 - 4 * tmp_28 - tmp_30 * tmp_41 -
                              tmp_31 * tmp_41 - 4 * tmp_33 + tmp_38 * tmp_4 + tmp_38 * tmp_6 + tmp_38 * tmp_9 + tmp_39 * tmp_5 +
                              tmp_4 * tmp_40 - tmp_40 * tmp_7 );
   real_t tmp_43 = tmp_32 * tmp_42;
   real_t tmp_44 = tmp_27 * tmp_42;
   real_t tmp_45 = tmp_21 * tmp_35;
   real_t tmp_46 = tmp_16 * tmp_35;
   real_t tmp_47 = tmp_21 * tmp_42;
   real_t tmp_48 = tmp_16 * tmp_42;
   real_t tmp_49 = tmp_45 + tmp_46 - tmp_47 - tmp_48;
   real_t tmp_50 = tmp_19 * tmp_35;
   real_t tmp_51 = tmp_14 * tmp_35;
   real_t tmp_52 = tmp_19 * tmp_42;
   real_t tmp_53 = tmp_14 * tmp_42;
   real_t tmp_54 = tmp_50 + tmp_51 - tmp_52 - tmp_53;
   real_t tmp_55 = p_affine_0_0 * p_affine_1_0;
   real_t tmp_56 = tmp_42 * tmp_55;
   real_t tmp_57 = p_affine_0_1 * p_affine_1_1;
   real_t tmp_58 = tmp_42 * tmp_57;
   real_t tmp_59 = tmp_35 * tmp_55;
   real_t tmp_60 = tmp_35 * tmp_57;
   real_t tmp_61 = tmp_36 + tmp_37 - tmp_43 - tmp_44;
   real_t tmp_62 = p_affine_0_0 * p_affine_2_0;
   real_t tmp_63 = tmp_35 * tmp_62;
   real_t tmp_64 = p_affine_0_1 * p_affine_2_1;
   real_t tmp_65 = tmp_35 * tmp_64;
   real_t tmp_66 = tmp_42 * tmp_62;
   real_t tmp_67 = tmp_42 * tmp_64;
   real_t tmp_68 = tmp_63 + tmp_65 - tmp_66 - tmp_67;
   real_t tmp_69 = -tmp_45 - tmp_46 + tmp_47 + tmp_48 + tmp_56 + tmp_58 - tmp_59 - tmp_60 + tmp_61 + tmp_68;
   real_t tmp_70 = -tmp_56 - tmp_58 + tmp_59 + tmp_60;
   real_t tmp_71 = -tmp_50 - tmp_51 + tmp_52 + tmp_53 + tmp_61 - tmp_63 - tmp_65 + tmp_66 + tmp_67 + tmp_70;
   real_t tmp_72 = p_affine_2_0 * tmp_3;
   real_t tmp_73 = p_affine_2_1 * tmp_29;
   real_t tmp_74 = tmp_13 * tmp_35;
   real_t tmp_75 = tmp_18 * tmp_35;
   real_t tmp_76 = tmp_13 * tmp_42;
   real_t tmp_77 = tmp_18 * tmp_42;
   real_t tmp_78 = tmp_74 + tmp_75 - tmp_76 - tmp_77;
   real_t tmp_79 = -tmp_36 - tmp_37 + tmp_43 + tmp_44 + tmp_68 + tmp_70 - tmp_74 - tmp_75 + tmp_76 + tmp_77;
   real_t tmp_80 = p_affine_1_0 * tmp_3;
   real_t tmp_81 = p_affine_1_1 * tmp_29;
   real_t a_0_0  = -2 * tmp_36 - 2 * tmp_37 + 2 * tmp_43 + 2 * tmp_44 + tmp_49 + tmp_54;
   real_t a_0_1  = tmp_69;
   real_t a_0_2  = tmp_71;
   real_t a_1_0  = tmp_69;
   real_t a_1_1  = -tmp_35 * tmp_72 - tmp_35 * tmp_73 + tmp_42 * tmp_72 + tmp_42 * tmp_73 + tmp_49 + tmp_78;
   real_t a_1_2  = tmp_79;
   real_t a_2_0  = tmp_71;
   real_t a_2_1  = tmp_79;
   real_t a_2_2  = -tmp_35 * tmp_80 - tmp_35 * tmp_81 + tmp_42 * tmp_80 + tmp_42 * tmp_81 + tmp_54 + tmp_78;

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

void DGDiffusionForm_Example::integrateVolume3D( const std::vector< Point3D >& coords,
                                                 const DGBasisInfo&            trialBasis,
                                                 const DGBasisInfo&            testBasis,
                                                 int                           trialDegree,
                                                 int                           testDegree,
                                                 MatrixXr&                     elMat ) const
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

   real_t tmp_0  = -p_affine_0_0;
   real_t tmp_1  = p_affine_1_0 + tmp_0;
   real_t tmp_2  = -p_affine_0_1;
   real_t tmp_3  = p_affine_2_1 + tmp_2;
   real_t tmp_4  = tmp_1 * tmp_3;
   real_t tmp_5  = p_affine_2_0 + tmp_0;
   real_t tmp_6  = p_affine_1_1 + tmp_2;
   real_t tmp_7  = tmp_5 * tmp_6;
   real_t tmp_8  = tmp_4 - tmp_7;
   real_t tmp_9  = -p_affine_0_2;
   real_t tmp_10 = p_affine_3_2 + tmp_9;
   real_t tmp_11 = p_affine_1_2 + tmp_9;
   real_t tmp_12 = p_affine_3_1 + tmp_2;
   real_t tmp_13 = tmp_12 * tmp_5;
   real_t tmp_14 = p_affine_2_2 + tmp_9;
   real_t tmp_15 = p_affine_3_0 + tmp_0;
   real_t tmp_16 = tmp_15 * tmp_6;
   real_t tmp_17 = tmp_1 * tmp_12;
   real_t tmp_18 = tmp_15 * tmp_3;
   real_t tmp_19 = tmp_10 * tmp_4 - tmp_10 * tmp_7 + tmp_11 * tmp_13 - tmp_11 * tmp_18 + tmp_14 * tmp_16 - tmp_14 * tmp_17;
   real_t tmp_20 = 1.0 / ( tmp_19 );
   real_t tmp_21 = tmp_20 * tmp_8;
   real_t tmp_22 = tmp_16 - tmp_17;
   real_t tmp_23 = tmp_20 * tmp_22;
   real_t tmp_24 = tmp_13 - tmp_18;
   real_t tmp_25 = tmp_20 * tmp_24;
   real_t tmp_26 = -tmp_21 - tmp_23 - tmp_25;
   real_t tmp_27 = -tmp_1 * tmp_14 + tmp_11 * tmp_5;
   real_t tmp_28 = tmp_20 * tmp_27;
   real_t tmp_29 = tmp_1 * tmp_10 - tmp_11 * tmp_15;
   real_t tmp_30 = tmp_20 * tmp_29;
   real_t tmp_31 = -tmp_10 * tmp_5 + tmp_14 * tmp_15;
   real_t tmp_32 = tmp_20 * tmp_31;
   real_t tmp_33 = -tmp_28 - tmp_30 - tmp_32;
   real_t tmp_34 = -tmp_11 * tmp_3 + tmp_14 * tmp_6;
   real_t tmp_35 = tmp_20 * tmp_34;
   real_t tmp_36 = -tmp_10 * tmp_6 + tmp_11 * tmp_12;
   real_t tmp_37 = tmp_20 * tmp_36;
   real_t tmp_38 = tmp_10 * tmp_3 - tmp_12 * tmp_14;
   real_t tmp_39 = tmp_20 * tmp_38;
   real_t tmp_40 = -tmp_35 - tmp_37 - tmp_39;
   real_t tmp_41 = p_affine_0_0 * p_affine_1_1;
   real_t tmp_42 = p_affine_0_0 * p_affine_1_2;
   real_t tmp_43 = p_affine_2_1 * p_affine_3_2;
   real_t tmp_44 = p_affine_0_1 * p_affine_1_0;
   real_t tmp_45 = p_affine_0_1 * p_affine_1_2;
   real_t tmp_46 = p_affine_2_2 * p_affine_3_0;
   real_t tmp_47 = p_affine_0_2 * p_affine_1_0;
   real_t tmp_48 = p_affine_0_2 * p_affine_1_1;
   real_t tmp_49 = p_affine_2_0 * p_affine_3_1;
   real_t tmp_50 = p_affine_2_2 * p_affine_3_1;
   real_t tmp_51 = p_affine_2_0 * p_affine_3_2;
   real_t tmp_52 = p_affine_2_1 * p_affine_3_0;
   real_t tmp_53 = std::abs( p_affine_0_0 * tmp_43 - p_affine_0_0 * tmp_50 + p_affine_0_1 * tmp_46 - p_affine_0_1 * tmp_51 +
                             p_affine_0_2 * tmp_49 - p_affine_0_2 * tmp_52 - p_affine_1_0 * tmp_43 + p_affine_1_0 * tmp_50 -
                             p_affine_1_1 * tmp_46 + p_affine_1_1 * tmp_51 - p_affine_1_2 * tmp_49 + p_affine_1_2 * tmp_52 +
                             p_affine_2_0 * tmp_45 - p_affine_2_0 * tmp_48 - p_affine_2_1 * tmp_42 + p_affine_2_1 * tmp_47 +
                             p_affine_2_2 * tmp_41 - p_affine_2_2 * tmp_44 - p_affine_3_0 * tmp_45 + p_affine_3_0 * tmp_48 +
                             p_affine_3_1 * tmp_42 - p_affine_3_1 * tmp_47 - p_affine_3_2 * tmp_41 + p_affine_3_2 * tmp_44 );
   real_t tmp_54 = tmp_53 * ( ( tmp_26 * tmp_26 ) + ( tmp_33 * tmp_33 ) + ( tmp_40 * tmp_40 ) );
   real_t tmp_55 = tmp_53 * ( tmp_25 * tmp_26 + tmp_32 * tmp_33 + tmp_39 * tmp_40 );
   real_t tmp_56 = 0.16666666666666666 * tmp_55;
   real_t tmp_57 = tmp_53 * ( tmp_23 * tmp_26 + tmp_30 * tmp_33 + tmp_37 * tmp_40 );
   real_t tmp_58 = 0.16666666666666666 * tmp_57;
   real_t tmp_59 = tmp_53 * ( tmp_21 * tmp_26 + tmp_28 * tmp_33 + tmp_35 * tmp_40 );
   real_t tmp_60 = 0.16666666666666666 * tmp_59;
   real_t tmp_61 = 1.0 / ( tmp_19 * tmp_19 );
   real_t tmp_62 = tmp_53 * ( ( tmp_24 * tmp_24 ) * tmp_61 + ( tmp_31 * tmp_31 ) * tmp_61 + ( tmp_38 * tmp_38 ) * tmp_61 );
   real_t tmp_63 = tmp_24 * tmp_61;
   real_t tmp_64 = tmp_31 * tmp_61;
   real_t tmp_65 = tmp_38 * tmp_61;
   real_t tmp_66 = tmp_53 * ( tmp_22 * tmp_63 + tmp_29 * tmp_64 + tmp_36 * tmp_65 );
   real_t tmp_67 = 0.16666666666666666 * tmp_66;
   real_t tmp_68 = tmp_53 * ( tmp_27 * tmp_64 + tmp_34 * tmp_65 + tmp_63 * tmp_8 );
   real_t tmp_69 = 0.16666666666666666 * tmp_68;
   real_t tmp_70 = tmp_53 * ( ( tmp_22 * tmp_22 ) * tmp_61 + ( tmp_29 * tmp_29 ) * tmp_61 + ( tmp_36 * tmp_36 ) * tmp_61 );
   real_t tmp_71 = tmp_53 * ( tmp_22 * tmp_61 * tmp_8 + tmp_27 * tmp_29 * tmp_61 + tmp_34 * tmp_36 * tmp_61 );
   real_t tmp_72 = 0.16666666666666666 * tmp_71;
   real_t tmp_73 = tmp_53 * ( ( tmp_27 * tmp_27 ) * tmp_61 + ( tmp_34 * tmp_34 ) * tmp_61 + tmp_61 * ( tmp_8 * tmp_8 ) );
   real_t a_0_0  = 0.16666666666666666 * tmp_54;
   real_t a_0_1  = tmp_56;
   real_t a_0_2  = tmp_58;
   real_t a_0_3  = tmp_60;
   real_t a_1_0  = tmp_56;
   real_t a_1_1  = 0.16666666666666666 * tmp_62;
   real_t a_1_2  = tmp_67;
   real_t a_1_3  = tmp_69;
   real_t a_2_0  = tmp_58;
   real_t a_2_1  = tmp_67;
   real_t a_2_2  = 0.16666666666666666 * tmp_70;
   real_t a_2_3  = tmp_72;
   real_t a_3_0  = tmp_60;
   real_t a_3_1  = tmp_69;
   real_t a_3_2  = tmp_72;
   real_t a_3_3  = 0.16666666666666666 * tmp_73;

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

void DGDiffusionForm_Example::integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&     outwardNormal,
                                                     const DGBasisInfo& trialBasis,
                                                     const DGBasisInfo& testBasis,
                                                     int                trialDegree,
                                                     int                testDegree,
                                                     MatrixXr&          elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 2, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 2, uint_c( trialDegree ) ) ) );

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

   real_t tmp_0  = -p_affine_6_1 + p_affine_7_1;
   real_t tmp_1  = -p_affine_0_1;
   real_t tmp_2  = p_affine_6_1 + tmp_1;
   real_t tmp_3  = 0.21132486540518713 * tmp_0 + tmp_2;
   real_t tmp_4  = -p_affine_0_0;
   real_t tmp_5  = p_affine_1_0 + tmp_4;
   real_t tmp_6  = p_affine_2_1 + tmp_1;
   real_t tmp_7  = 1.0 / ( tmp_5 * tmp_6 - ( p_affine_1_1 + tmp_1 ) * ( p_affine_2_0 + tmp_4 ) );
   real_t tmp_8  = tmp_5 * tmp_7;
   real_t tmp_9  = tmp_3 * tmp_8;
   real_t tmp_10 = tmp_7 * ( p_affine_0_0 - p_affine_2_0 );
   real_t tmp_11 = tmp_10 * tmp_3;
   real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
   real_t tmp_13 = p_affine_6_0 + tmp_4;
   real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
   real_t tmp_15 = tmp_6 * tmp_7;
   real_t tmp_16 = tmp_14 * tmp_15;
   real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
   real_t tmp_18 = tmp_14 * tmp_17;
   real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
   real_t tmp_20 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
   real_t tmp_21 = sigma_0 * std::pow( tmp_20, -beta_0 );
   real_t tmp_22 = 0.5 * p_affine_10_0;
   real_t tmp_23 = tmp_22 * ( -tmp_15 - tmp_17 );
   real_t tmp_24 = 0.5 * p_affine_10_1;
   real_t tmp_25 = tmp_24 * ( -tmp_10 - tmp_8 );
   real_t tmp_26 = -tmp_23 - tmp_25;
   real_t tmp_27 = tmp_23 + tmp_25;
   real_t tmp_28 = 0.5 * tmp_20;
   real_t tmp_29 = 0.78867513459481287 * tmp_0 + tmp_2;
   real_t tmp_30 = tmp_29 * tmp_8;
   real_t tmp_31 = tmp_10 * tmp_29;
   real_t tmp_32 = 0.78867513459481287 * tmp_12 + tmp_13;
   real_t tmp_33 = tmp_15 * tmp_32;
   real_t tmp_34 = tmp_17 * tmp_32;
   real_t tmp_35 = -tmp_30 - tmp_31 - tmp_33 - tmp_34 + 1;
   real_t tmp_36 = 0.5 * tmp_20;
   real_t tmp_37 = tmp_11 + tmp_16;
   real_t tmp_38 = tmp_15 * tmp_22;
   real_t tmp_39 = tmp_10 * tmp_24;
   real_t tmp_40 = tmp_38 + tmp_39;
   real_t tmp_41 = tmp_19 * tmp_21;
   real_t tmp_42 = tmp_37 * tmp_41;
   real_t tmp_43 = tmp_31 + tmp_33;
   real_t tmp_44 = tmp_21 * tmp_35;
   real_t tmp_45 = tmp_43 * tmp_44;
   real_t tmp_46 = tmp_18 + tmp_9;
   real_t tmp_47 = tmp_17 * tmp_22;
   real_t tmp_48 = tmp_24 * tmp_8;
   real_t tmp_49 = tmp_47 + tmp_48;
   real_t tmp_50 = tmp_41 * tmp_46;
   real_t tmp_51 = tmp_30 + tmp_34;
   real_t tmp_52 = tmp_44 * tmp_51;
   real_t tmp_53 = -tmp_38 - tmp_39;
   real_t tmp_54 = tmp_21 * tmp_37 * tmp_46;
   real_t tmp_55 = tmp_21 * tmp_43 * tmp_51;
   real_t tmp_56 = -tmp_47 - tmp_48;
   real_t a_0_0  = tmp_28 * ( ( tmp_19 * tmp_19 ) * tmp_21 + tmp_19 * tmp_26 - tmp_19 * tmp_27 ) +
                  tmp_36 * ( tmp_21 * ( tmp_35 * tmp_35 ) + tmp_26 * tmp_35 - tmp_27 * tmp_35 );
   real_t a_0_1 =
       tmp_28 * ( -tmp_19 * tmp_40 + tmp_26 * tmp_37 + tmp_42 ) + tmp_36 * ( tmp_26 * tmp_43 - tmp_35 * tmp_40 + tmp_45 );
   real_t a_0_2 =
       tmp_28 * ( -tmp_19 * tmp_49 + tmp_26 * tmp_46 + tmp_50 ) + tmp_36 * ( tmp_26 * tmp_51 - tmp_35 * tmp_49 + tmp_52 );
   real_t a_1_0 =
       tmp_28 * ( tmp_19 * tmp_53 - tmp_27 * tmp_37 + tmp_42 ) + tmp_36 * ( -tmp_27 * tmp_43 + tmp_35 * tmp_53 + tmp_45 );
   real_t a_1_1 = tmp_28 * ( tmp_21 * ( tmp_37 * tmp_37 ) - tmp_37 * tmp_40 + tmp_37 * tmp_53 ) +
                  tmp_36 * ( tmp_21 * ( tmp_43 * tmp_43 ) - tmp_40 * tmp_43 + tmp_43 * tmp_53 );
   real_t a_1_2 =
       tmp_28 * ( -tmp_37 * tmp_49 + tmp_46 * tmp_53 + tmp_54 ) + tmp_36 * ( -tmp_43 * tmp_49 + tmp_51 * tmp_53 + tmp_55 );
   real_t a_2_0 =
       tmp_28 * ( tmp_19 * tmp_56 - tmp_27 * tmp_46 + tmp_50 ) + tmp_36 * ( -tmp_27 * tmp_51 + tmp_35 * tmp_56 + tmp_52 );
   real_t a_2_1 =
       tmp_28 * ( tmp_37 * tmp_56 - tmp_40 * tmp_46 + tmp_54 ) + tmp_36 * ( -tmp_40 * tmp_51 + tmp_43 * tmp_56 + tmp_55 );
   real_t a_2_2 = tmp_28 * ( tmp_21 * ( tmp_46 * tmp_46 ) - tmp_46 * tmp_49 + tmp_46 * tmp_56 ) +
                  tmp_36 * ( tmp_21 * ( tmp_51 * tmp_51 ) - tmp_49 * tmp_51 + tmp_51 * tmp_56 );

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

void DGDiffusionForm_Example::integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&     outwardNormal,
                                                     const DGBasisInfo& trialBasis,
                                                     const DGBasisInfo& testBasis,
                                                     int                trialDegree,
                                                     int                testDegree,
                                                     MatrixXr&          elMat ) const
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

   real_t tmp_0  = -p_affine_8_2;
   real_t tmp_1  = p_affine_9_2 + tmp_0;
   real_t tmp_2  = p_affine_10_2 + tmp_0;
   real_t tmp_3  = -p_affine_0_2;
   real_t tmp_4  = p_affine_8_2 + tmp_3;
   real_t tmp_5  = 0.091576213509770743 * tmp_1 + 0.81684757298045851 * tmp_2 + tmp_4;
   real_t tmp_6  = -p_affine_0_0;
   real_t tmp_7  = p_affine_1_0 + tmp_6;
   real_t tmp_8  = -p_affine_0_1;
   real_t tmp_9  = p_affine_2_1 + tmp_8;
   real_t tmp_10 = p_affine_2_0 + tmp_6;
   real_t tmp_11 = p_affine_1_1 + tmp_8;
   real_t tmp_12 = p_affine_3_2 + tmp_3;
   real_t tmp_13 = tmp_12 * tmp_9;
   real_t tmp_14 = p_affine_3_1 + tmp_8;
   real_t tmp_15 = p_affine_1_2 + tmp_3;
   real_t tmp_16 = tmp_14 * tmp_15;
   real_t tmp_17 = p_affine_3_0 + tmp_6;
   real_t tmp_18 = p_affine_2_2 + tmp_3;
   real_t tmp_19 = tmp_11 * tmp_18;
   real_t tmp_20 = tmp_14 * tmp_18;
   real_t tmp_21 = tmp_11 * tmp_12;
   real_t tmp_22 = tmp_15 * tmp_9;
   real_t tmp_23 =
       1.0 / ( tmp_10 * tmp_16 - tmp_10 * tmp_21 + tmp_13 * tmp_7 + tmp_17 * tmp_19 - tmp_17 * tmp_22 - tmp_20 * tmp_7 );
   real_t tmp_24 = tmp_23 * ( -tmp_10 * tmp_11 + tmp_7 * tmp_9 );
   real_t tmp_25 = tmp_24 * tmp_5;
   real_t tmp_26 = tmp_23 * ( tmp_11 * tmp_17 - tmp_14 * tmp_7 );
   real_t tmp_27 = tmp_26 * tmp_5;
   real_t tmp_28 = -p_affine_8_1;
   real_t tmp_29 = p_affine_9_1 + tmp_28;
   real_t tmp_30 = p_affine_10_1 + tmp_28;
   real_t tmp_31 = p_affine_8_1 + tmp_8;
   real_t tmp_32 = 0.091576213509770743 * tmp_29 + 0.81684757298045851 * tmp_30 + tmp_31;
   real_t tmp_33 = tmp_23 * ( tmp_10 * tmp_15 - tmp_18 * tmp_7 );
   real_t tmp_34 = tmp_32 * tmp_33;
   real_t tmp_35 = tmp_23 * ( tmp_12 * tmp_7 - tmp_15 * tmp_17 );
   real_t tmp_36 = tmp_32 * tmp_35;
   real_t tmp_37 = tmp_23 * ( tmp_10 * tmp_14 - tmp_17 * tmp_9 );
   real_t tmp_38 = tmp_37 * tmp_5;
   real_t tmp_39 = tmp_23 * ( -tmp_10 * tmp_12 + tmp_17 * tmp_18 );
   real_t tmp_40 = tmp_32 * tmp_39;
   real_t tmp_41 = -p_affine_8_0;
   real_t tmp_42 = p_affine_9_0 + tmp_41;
   real_t tmp_43 = p_affine_10_0 + tmp_41;
   real_t tmp_44 = p_affine_8_0 + tmp_6;
   real_t tmp_45 = 0.091576213509770743 * tmp_42 + 0.81684757298045851 * tmp_43 + tmp_44;
   real_t tmp_46 = tmp_23 * ( tmp_19 - tmp_22 );
   real_t tmp_47 = tmp_45 * tmp_46;
   real_t tmp_48 = tmp_23 * ( tmp_16 - tmp_21 );
   real_t tmp_49 = tmp_45 * tmp_48;
   real_t tmp_50 = tmp_23 * ( tmp_13 - tmp_20 );
   real_t tmp_51 = tmp_45 * tmp_50;
   real_t tmp_52 = -tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1;
   real_t tmp_53 = p_affine_8_1 - p_affine_9_1;
   real_t tmp_54 = p_affine_8_0 - p_affine_9_0;
   real_t tmp_55 = p_affine_8_2 - p_affine_9_2;
   real_t tmp_56 =
       std::pow( ( std::abs( tmp_2 * tmp_53 - tmp_30 * tmp_55 ) * std::abs( tmp_2 * tmp_53 - tmp_30 * tmp_55 ) ) +
                     ( std::abs( tmp_2 * tmp_54 - tmp_43 * tmp_55 ) * std::abs( tmp_2 * tmp_54 - tmp_43 * tmp_55 ) ) +
                     ( std::abs( tmp_30 * tmp_54 - tmp_43 * tmp_53 ) * std::abs( tmp_30 * tmp_54 - tmp_43 * tmp_53 ) ),
                 1.0 / 2.0 );
   real_t tmp_57  = sigma_0 * std::pow( 0.5 * tmp_56, -beta_0 );
   real_t tmp_58  = 0.5 * p_affine_13_0;
   real_t tmp_59  = tmp_58 * ( -tmp_46 - tmp_48 - tmp_50 );
   real_t tmp_60  = 0.5 * p_affine_13_1;
   real_t tmp_61  = tmp_60 * ( -tmp_33 - tmp_35 - tmp_39 );
   real_t tmp_62  = 0.5 * p_affine_13_2;
   real_t tmp_63  = tmp_62 * ( -tmp_24 - tmp_26 - tmp_37 );
   real_t tmp_64  = -tmp_59 - tmp_61 - tmp_63;
   real_t tmp_65  = tmp_59 + tmp_61 + tmp_63;
   real_t tmp_66  = 1.0 * tmp_56;
   real_t tmp_67  = 0.054975871827660928 * tmp_66;
   real_t tmp_68  = 0.44594849091596489 * tmp_1 + 0.10810301816807022 * tmp_2 + tmp_4;
   real_t tmp_69  = tmp_24 * tmp_68;
   real_t tmp_70  = tmp_26 * tmp_68;
   real_t tmp_71  = 0.44594849091596489 * tmp_29 + 0.10810301816807022 * tmp_30 + tmp_31;
   real_t tmp_72  = tmp_33 * tmp_71;
   real_t tmp_73  = tmp_35 * tmp_71;
   real_t tmp_74  = tmp_37 * tmp_68;
   real_t tmp_75  = tmp_39 * tmp_71;
   real_t tmp_76  = 0.44594849091596489 * tmp_42 + 0.10810301816807022 * tmp_43 + tmp_44;
   real_t tmp_77  = tmp_46 * tmp_76;
   real_t tmp_78  = tmp_48 * tmp_76;
   real_t tmp_79  = tmp_50 * tmp_76;
   real_t tmp_80  = -tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74 - tmp_75 - tmp_77 - tmp_78 - tmp_79 + 1;
   real_t tmp_81  = 0.11169079483900572 * tmp_66;
   real_t tmp_82  = 0.81684757298045851 * tmp_1 + 0.091576213509770743 * tmp_2 + tmp_4;
   real_t tmp_83  = tmp_24 * tmp_82;
   real_t tmp_84  = tmp_26 * tmp_82;
   real_t tmp_85  = 0.81684757298045851 * tmp_29 + 0.091576213509770743 * tmp_30 + tmp_31;
   real_t tmp_86  = tmp_33 * tmp_85;
   real_t tmp_87  = tmp_35 * tmp_85;
   real_t tmp_88  = tmp_37 * tmp_82;
   real_t tmp_89  = tmp_39 * tmp_85;
   real_t tmp_90  = 0.81684757298045851 * tmp_42 + 0.091576213509770743 * tmp_43 + tmp_44;
   real_t tmp_91  = tmp_46 * tmp_90;
   real_t tmp_92  = tmp_48 * tmp_90;
   real_t tmp_93  = tmp_50 * tmp_90;
   real_t tmp_94  = -tmp_83 - tmp_84 - tmp_86 - tmp_87 - tmp_88 - tmp_89 - tmp_91 - tmp_92 - tmp_93 + 1;
   real_t tmp_95  = 0.054975871827660928 * tmp_66;
   real_t tmp_96  = 0.10810301816807022 * tmp_1 + 0.44594849091596489 * tmp_2 + tmp_4;
   real_t tmp_97  = tmp_24 * tmp_96;
   real_t tmp_98  = tmp_26 * tmp_96;
   real_t tmp_99  = 0.10810301816807022 * tmp_29 + 0.44594849091596489 * tmp_30 + tmp_31;
   real_t tmp_100 = tmp_33 * tmp_99;
   real_t tmp_101 = tmp_35 * tmp_99;
   real_t tmp_102 = tmp_37 * tmp_96;
   real_t tmp_103 = tmp_39 * tmp_99;
   real_t tmp_104 = 0.10810301816807022 * tmp_42 + 0.44594849091596489 * tmp_43 + tmp_44;
   real_t tmp_105 = tmp_104 * tmp_46;
   real_t tmp_106 = tmp_104 * tmp_48;
   real_t tmp_107 = tmp_104 * tmp_50;
   real_t tmp_108 = -tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1;
   real_t tmp_109 = 0.11169079483900572 * tmp_66;
   real_t tmp_110 = 0.091576213509770743 * tmp_1 + 0.091576213509770743 * tmp_2 + tmp_4;
   real_t tmp_111 = tmp_110 * tmp_24;
   real_t tmp_112 = tmp_110 * tmp_26;
   real_t tmp_113 = 0.091576213509770743 * tmp_29 + 0.091576213509770743 * tmp_30 + tmp_31;
   real_t tmp_114 = tmp_113 * tmp_33;
   real_t tmp_115 = tmp_113 * tmp_35;
   real_t tmp_116 = tmp_110 * tmp_37;
   real_t tmp_117 = tmp_113 * tmp_39;
   real_t tmp_118 = 0.091576213509770743 * tmp_42 + 0.091576213509770743 * tmp_43 + tmp_44;
   real_t tmp_119 = tmp_118 * tmp_46;
   real_t tmp_120 = tmp_118 * tmp_48;
   real_t tmp_121 = tmp_118 * tmp_50;
   real_t tmp_122 = -tmp_111 - tmp_112 - tmp_114 - tmp_115 - tmp_116 - tmp_117 - tmp_119 - tmp_120 - tmp_121 + 1;
   real_t tmp_123 = 0.054975871827660928 * tmp_66;
   real_t tmp_124 = 0.44594849091596489 * tmp_1 + 0.44594849091596489 * tmp_2 + tmp_4;
   real_t tmp_125 = tmp_124 * tmp_24;
   real_t tmp_126 = tmp_124 * tmp_26;
   real_t tmp_127 = 0.44594849091596489 * tmp_29 + 0.44594849091596489 * tmp_30 + tmp_31;
   real_t tmp_128 = tmp_127 * tmp_33;
   real_t tmp_129 = tmp_127 * tmp_35;
   real_t tmp_130 = tmp_124 * tmp_37;
   real_t tmp_131 = tmp_127 * tmp_39;
   real_t tmp_132 = 0.44594849091596489 * tmp_42 + 0.44594849091596489 * tmp_43 + tmp_44;
   real_t tmp_133 = tmp_132 * tmp_46;
   real_t tmp_134 = tmp_132 * tmp_48;
   real_t tmp_135 = tmp_132 * tmp_50;
   real_t tmp_136 = -tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1;
   real_t tmp_137 = 0.11169079483900572 * tmp_66;
   real_t tmp_138 = tmp_38 + tmp_40 + tmp_51;
   real_t tmp_139 = tmp_50 * tmp_58;
   real_t tmp_140 = tmp_39 * tmp_60;
   real_t tmp_141 = tmp_37 * tmp_62;
   real_t tmp_142 = tmp_139 + tmp_140 + tmp_141;
   real_t tmp_143 = tmp_52 * tmp_57;
   real_t tmp_144 = tmp_138 * tmp_143;
   real_t tmp_145 = tmp_74 + tmp_75 + tmp_79;
   real_t tmp_146 = tmp_57 * tmp_80;
   real_t tmp_147 = tmp_145 * tmp_146;
   real_t tmp_148 = tmp_88 + tmp_89 + tmp_93;
   real_t tmp_149 = tmp_57 * tmp_94;
   real_t tmp_150 = tmp_148 * tmp_149;
   real_t tmp_151 = tmp_102 + tmp_103 + tmp_107;
   real_t tmp_152 = tmp_108 * tmp_57;
   real_t tmp_153 = tmp_151 * tmp_152;
   real_t tmp_154 = tmp_116 + tmp_117 + tmp_121;
   real_t tmp_155 = tmp_122 * tmp_57;
   real_t tmp_156 = tmp_154 * tmp_155;
   real_t tmp_157 = tmp_130 + tmp_131 + tmp_135;
   real_t tmp_158 = tmp_136 * tmp_57;
   real_t tmp_159 = tmp_157 * tmp_158;
   real_t tmp_160 = tmp_27 + tmp_36 + tmp_49;
   real_t tmp_161 = tmp_48 * tmp_58;
   real_t tmp_162 = tmp_35 * tmp_60;
   real_t tmp_163 = tmp_26 * tmp_62;
   real_t tmp_164 = tmp_161 + tmp_162 + tmp_163;
   real_t tmp_165 = tmp_143 * tmp_160;
   real_t tmp_166 = tmp_70 + tmp_73 + tmp_78;
   real_t tmp_167 = tmp_146 * tmp_166;
   real_t tmp_168 = tmp_84 + tmp_87 + tmp_92;
   real_t tmp_169 = tmp_149 * tmp_168;
   real_t tmp_170 = tmp_101 + tmp_106 + tmp_98;
   real_t tmp_171 = tmp_152 * tmp_170;
   real_t tmp_172 = tmp_112 + tmp_115 + tmp_120;
   real_t tmp_173 = tmp_155 * tmp_172;
   real_t tmp_174 = tmp_126 + tmp_129 + tmp_134;
   real_t tmp_175 = tmp_158 * tmp_174;
   real_t tmp_176 = tmp_25 + tmp_34 + tmp_47;
   real_t tmp_177 = tmp_46 * tmp_58;
   real_t tmp_178 = tmp_33 * tmp_60;
   real_t tmp_179 = tmp_24 * tmp_62;
   real_t tmp_180 = tmp_177 + tmp_178 + tmp_179;
   real_t tmp_181 = tmp_143 * tmp_176;
   real_t tmp_182 = tmp_69 + tmp_72 + tmp_77;
   real_t tmp_183 = tmp_146 * tmp_182;
   real_t tmp_184 = tmp_83 + tmp_86 + tmp_91;
   real_t tmp_185 = tmp_149 * tmp_184;
   real_t tmp_186 = tmp_100 + tmp_105 + tmp_97;
   real_t tmp_187 = tmp_152 * tmp_186;
   real_t tmp_188 = tmp_111 + tmp_114 + tmp_119;
   real_t tmp_189 = tmp_155 * tmp_188;
   real_t tmp_190 = tmp_125 + tmp_128 + tmp_133;
   real_t tmp_191 = tmp_158 * tmp_190;
   real_t tmp_192 = -tmp_139 - tmp_140 - tmp_141;
   real_t tmp_193 = tmp_138 * tmp_57;
   real_t tmp_194 = tmp_160 * tmp_193;
   real_t tmp_195 = tmp_145 * tmp_57;
   real_t tmp_196 = tmp_166 * tmp_195;
   real_t tmp_197 = tmp_148 * tmp_57;
   real_t tmp_198 = tmp_168 * tmp_197;
   real_t tmp_199 = tmp_151 * tmp_57;
   real_t tmp_200 = tmp_170 * tmp_199;
   real_t tmp_201 = tmp_154 * tmp_57;
   real_t tmp_202 = tmp_172 * tmp_201;
   real_t tmp_203 = tmp_157 * tmp_57;
   real_t tmp_204 = tmp_174 * tmp_203;
   real_t tmp_205 = tmp_176 * tmp_193;
   real_t tmp_206 = tmp_182 * tmp_195;
   real_t tmp_207 = tmp_184 * tmp_197;
   real_t tmp_208 = tmp_186 * tmp_199;
   real_t tmp_209 = tmp_188 * tmp_201;
   real_t tmp_210 = tmp_190 * tmp_203;
   real_t tmp_211 = -tmp_161 - tmp_162 - tmp_163;
   real_t tmp_212 = tmp_160 * tmp_176 * tmp_57;
   real_t tmp_213 = tmp_166 * tmp_182 * tmp_57;
   real_t tmp_214 = tmp_168 * tmp_184 * tmp_57;
   real_t tmp_215 = tmp_170 * tmp_186 * tmp_57;
   real_t tmp_216 = tmp_172 * tmp_188 * tmp_57;
   real_t tmp_217 = tmp_174 * tmp_190 * tmp_57;
   real_t tmp_218 = -tmp_177 - tmp_178 - tmp_179;
   real_t a_0_0   = tmp_109 * ( ( tmp_108 * tmp_108 ) * tmp_57 + tmp_108 * tmp_64 - tmp_108 * tmp_65 ) +
                  tmp_123 * ( ( tmp_122 * tmp_122 ) * tmp_57 + tmp_122 * tmp_64 - tmp_122 * tmp_65 ) +
                  tmp_137 * ( ( tmp_136 * tmp_136 ) * tmp_57 + tmp_136 * tmp_64 - tmp_136 * tmp_65 ) +
                  tmp_67 * ( ( tmp_52 * tmp_52 ) * tmp_57 + tmp_52 * tmp_64 - tmp_52 * tmp_65 ) +
                  tmp_81 * ( tmp_57 * ( tmp_80 * tmp_80 ) + tmp_64 * tmp_80 - tmp_65 * tmp_80 ) +
                  tmp_95 * ( tmp_57 * ( tmp_94 * tmp_94 ) + tmp_64 * tmp_94 - tmp_65 * tmp_94 );
   real_t a_0_1 = tmp_109 * ( -tmp_108 * tmp_142 + tmp_151 * tmp_64 + tmp_153 ) +
                  tmp_123 * ( -tmp_122 * tmp_142 + tmp_154 * tmp_64 + tmp_156 ) +
                  tmp_137 * ( -tmp_136 * tmp_142 + tmp_157 * tmp_64 + tmp_159 ) +
                  tmp_67 * ( tmp_138 * tmp_64 - tmp_142 * tmp_52 + tmp_144 ) +
                  tmp_81 * ( -tmp_142 * tmp_80 + tmp_145 * tmp_64 + tmp_147 ) +
                  tmp_95 * ( -tmp_142 * tmp_94 + tmp_148 * tmp_64 + tmp_150 );
   real_t a_0_2 = tmp_109 * ( -tmp_108 * tmp_164 + tmp_170 * tmp_64 + tmp_171 ) +
                  tmp_123 * ( -tmp_122 * tmp_164 + tmp_172 * tmp_64 + tmp_173 ) +
                  tmp_137 * ( -tmp_136 * tmp_164 + tmp_174 * tmp_64 + tmp_175 ) +
                  tmp_67 * ( tmp_160 * tmp_64 - tmp_164 * tmp_52 + tmp_165 ) +
                  tmp_81 * ( -tmp_164 * tmp_80 + tmp_166 * tmp_64 + tmp_167 ) +
                  tmp_95 * ( -tmp_164 * tmp_94 + tmp_168 * tmp_64 + tmp_169 );
   real_t a_0_3 = tmp_109 * ( -tmp_108 * tmp_180 + tmp_186 * tmp_64 + tmp_187 ) +
                  tmp_123 * ( -tmp_122 * tmp_180 + tmp_188 * tmp_64 + tmp_189 ) +
                  tmp_137 * ( -tmp_136 * tmp_180 + tmp_190 * tmp_64 + tmp_191 ) +
                  tmp_67 * ( tmp_176 * tmp_64 - tmp_180 * tmp_52 + tmp_181 ) +
                  tmp_81 * ( -tmp_180 * tmp_80 + tmp_182 * tmp_64 + tmp_183 ) +
                  tmp_95 * ( -tmp_180 * tmp_94 + tmp_184 * tmp_64 + tmp_185 );
   real_t a_1_0 = tmp_109 * ( tmp_108 * tmp_192 - tmp_151 * tmp_65 + tmp_153 ) +
                  tmp_123 * ( tmp_122 * tmp_192 - tmp_154 * tmp_65 + tmp_156 ) +
                  tmp_137 * ( tmp_136 * tmp_192 - tmp_157 * tmp_65 + tmp_159 ) +
                  tmp_67 * ( -tmp_138 * tmp_65 + tmp_144 + tmp_192 * tmp_52 ) +
                  tmp_81 * ( -tmp_145 * tmp_65 + tmp_147 + tmp_192 * tmp_80 ) +
                  tmp_95 * ( -tmp_148 * tmp_65 + tmp_150 + tmp_192 * tmp_94 );
   real_t a_1_1 = tmp_109 * ( -tmp_142 * tmp_151 + ( tmp_151 * tmp_151 ) * tmp_57 + tmp_151 * tmp_192 ) +
                  tmp_123 * ( -tmp_142 * tmp_154 + ( tmp_154 * tmp_154 ) * tmp_57 + tmp_154 * tmp_192 ) +
                  tmp_137 * ( -tmp_142 * tmp_157 + ( tmp_157 * tmp_157 ) * tmp_57 + tmp_157 * tmp_192 ) +
                  tmp_67 * ( ( tmp_138 * tmp_138 ) * tmp_57 - tmp_138 * tmp_142 + tmp_138 * tmp_192 ) +
                  tmp_81 * ( -tmp_142 * tmp_145 + ( tmp_145 * tmp_145 ) * tmp_57 + tmp_145 * tmp_192 ) +
                  tmp_95 * ( -tmp_142 * tmp_148 + ( tmp_148 * tmp_148 ) * tmp_57 + tmp_148 * tmp_192 );
   real_t a_1_2 = tmp_109 * ( -tmp_151 * tmp_164 + tmp_170 * tmp_192 + tmp_200 ) +
                  tmp_123 * ( -tmp_154 * tmp_164 + tmp_172 * tmp_192 + tmp_202 ) +
                  tmp_137 * ( -tmp_157 * tmp_164 + tmp_174 * tmp_192 + tmp_204 ) +
                  tmp_67 * ( -tmp_138 * tmp_164 + tmp_160 * tmp_192 + tmp_194 ) +
                  tmp_81 * ( -tmp_145 * tmp_164 + tmp_166 * tmp_192 + tmp_196 ) +
                  tmp_95 * ( -tmp_148 * tmp_164 + tmp_168 * tmp_192 + tmp_198 );
   real_t a_1_3 = tmp_109 * ( -tmp_151 * tmp_180 + tmp_186 * tmp_192 + tmp_208 ) +
                  tmp_123 * ( -tmp_154 * tmp_180 + tmp_188 * tmp_192 + tmp_209 ) +
                  tmp_137 * ( -tmp_157 * tmp_180 + tmp_190 * tmp_192 + tmp_210 ) +
                  tmp_67 * ( -tmp_138 * tmp_180 + tmp_176 * tmp_192 + tmp_205 ) +
                  tmp_81 * ( -tmp_145 * tmp_180 + tmp_182 * tmp_192 + tmp_206 ) +
                  tmp_95 * ( -tmp_148 * tmp_180 + tmp_184 * tmp_192 + tmp_207 );
   real_t a_2_0 = tmp_109 * ( tmp_108 * tmp_211 - tmp_170 * tmp_65 + tmp_171 ) +
                  tmp_123 * ( tmp_122 * tmp_211 - tmp_172 * tmp_65 + tmp_173 ) +
                  tmp_137 * ( tmp_136 * tmp_211 - tmp_174 * tmp_65 + tmp_175 ) +
                  tmp_67 * ( -tmp_160 * tmp_65 + tmp_165 + tmp_211 * tmp_52 ) +
                  tmp_81 * ( -tmp_166 * tmp_65 + tmp_167 + tmp_211 * tmp_80 ) +
                  tmp_95 * ( -tmp_168 * tmp_65 + tmp_169 + tmp_211 * tmp_94 );
   real_t a_2_1 = tmp_109 * ( -tmp_142 * tmp_170 + tmp_151 * tmp_211 + tmp_200 ) +
                  tmp_123 * ( -tmp_142 * tmp_172 + tmp_154 * tmp_211 + tmp_202 ) +
                  tmp_137 * ( -tmp_142 * tmp_174 + tmp_157 * tmp_211 + tmp_204 ) +
                  tmp_67 * ( tmp_138 * tmp_211 - tmp_142 * tmp_160 + tmp_194 ) +
                  tmp_81 * ( -tmp_142 * tmp_166 + tmp_145 * tmp_211 + tmp_196 ) +
                  tmp_95 * ( -tmp_142 * tmp_168 + tmp_148 * tmp_211 + tmp_198 );
   real_t a_2_2 = tmp_109 * ( -tmp_164 * tmp_170 + ( tmp_170 * tmp_170 ) * tmp_57 + tmp_170 * tmp_211 ) +
                  tmp_123 * ( -tmp_164 * tmp_172 + ( tmp_172 * tmp_172 ) * tmp_57 + tmp_172 * tmp_211 ) +
                  tmp_137 * ( -tmp_164 * tmp_174 + ( tmp_174 * tmp_174 ) * tmp_57 + tmp_174 * tmp_211 ) +
                  tmp_67 * ( ( tmp_160 * tmp_160 ) * tmp_57 - tmp_160 * tmp_164 + tmp_160 * tmp_211 ) +
                  tmp_81 * ( -tmp_164 * tmp_166 + ( tmp_166 * tmp_166 ) * tmp_57 + tmp_166 * tmp_211 ) +
                  tmp_95 * ( -tmp_164 * tmp_168 + ( tmp_168 * tmp_168 ) * tmp_57 + tmp_168 * tmp_211 );
   real_t a_2_3 = tmp_109 * ( -tmp_170 * tmp_180 + tmp_186 * tmp_211 + tmp_215 ) +
                  tmp_123 * ( -tmp_172 * tmp_180 + tmp_188 * tmp_211 + tmp_216 ) +
                  tmp_137 * ( -tmp_174 * tmp_180 + tmp_190 * tmp_211 + tmp_217 ) +
                  tmp_67 * ( -tmp_160 * tmp_180 + tmp_176 * tmp_211 + tmp_212 ) +
                  tmp_81 * ( -tmp_166 * tmp_180 + tmp_182 * tmp_211 + tmp_213 ) +
                  tmp_95 * ( -tmp_168 * tmp_180 + tmp_184 * tmp_211 + tmp_214 );
   real_t a_3_0 = tmp_109 * ( tmp_108 * tmp_218 - tmp_186 * tmp_65 + tmp_187 ) +
                  tmp_123 * ( tmp_122 * tmp_218 - tmp_188 * tmp_65 + tmp_189 ) +
                  tmp_137 * ( tmp_136 * tmp_218 - tmp_190 * tmp_65 + tmp_191 ) +
                  tmp_67 * ( -tmp_176 * tmp_65 + tmp_181 + tmp_218 * tmp_52 ) +
                  tmp_81 * ( -tmp_182 * tmp_65 + tmp_183 + tmp_218 * tmp_80 ) +
                  tmp_95 * ( -tmp_184 * tmp_65 + tmp_185 + tmp_218 * tmp_94 );
   real_t a_3_1 = tmp_109 * ( -tmp_142 * tmp_186 + tmp_151 * tmp_218 + tmp_208 ) +
                  tmp_123 * ( -tmp_142 * tmp_188 + tmp_154 * tmp_218 + tmp_209 ) +
                  tmp_137 * ( -tmp_142 * tmp_190 + tmp_157 * tmp_218 + tmp_210 ) +
                  tmp_67 * ( tmp_138 * tmp_218 - tmp_142 * tmp_176 + tmp_205 ) +
                  tmp_81 * ( -tmp_142 * tmp_182 + tmp_145 * tmp_218 + tmp_206 ) +
                  tmp_95 * ( -tmp_142 * tmp_184 + tmp_148 * tmp_218 + tmp_207 );
   real_t a_3_2 = tmp_109 * ( -tmp_164 * tmp_186 + tmp_170 * tmp_218 + tmp_215 ) +
                  tmp_123 * ( -tmp_164 * tmp_188 + tmp_172 * tmp_218 + tmp_216 ) +
                  tmp_137 * ( -tmp_164 * tmp_190 + tmp_174 * tmp_218 + tmp_217 ) +
                  tmp_67 * ( tmp_160 * tmp_218 - tmp_164 * tmp_176 + tmp_212 ) +
                  tmp_81 * ( -tmp_164 * tmp_182 + tmp_166 * tmp_218 + tmp_213 ) +
                  tmp_95 * ( -tmp_164 * tmp_184 + tmp_168 * tmp_218 + tmp_214 );
   real_t a_3_3 = tmp_109 * ( -tmp_180 * tmp_186 + ( tmp_186 * tmp_186 ) * tmp_57 + tmp_186 * tmp_218 ) +
                  tmp_123 * ( -tmp_180 * tmp_188 + ( tmp_188 * tmp_188 ) * tmp_57 + tmp_188 * tmp_218 ) +
                  tmp_137 * ( -tmp_180 * tmp_190 + ( tmp_190 * tmp_190 ) * tmp_57 + tmp_190 * tmp_218 ) +
                  tmp_67 * ( ( tmp_176 * tmp_176 ) * tmp_57 - tmp_176 * tmp_180 + tmp_176 * tmp_218 ) +
                  tmp_81 * ( -tmp_180 * tmp_182 + ( tmp_182 * tmp_182 ) * tmp_57 + tmp_182 * tmp_218 ) +
                  tmp_95 * ( -tmp_180 * tmp_184 + ( tmp_184 * tmp_184 ) * tmp_57 + tmp_184 * tmp_218 );

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

void DGDiffusionForm_Example::integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&     outwardNormal,
                                                        const DGBasisInfo& trialBasis,
                                                        const DGBasisInfo& testBasis,
                                                        int                trialDegree,
                                                        int                testDegree,
                                                        MatrixXr&          elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 2, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 2, uint_c( trialDegree ) ) ) );

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

   const auto p_affine_10_0 = outwardNormal( 0 );
   const auto p_affine_10_1 = outwardNormal( 1 );

   real_t tmp_0  = -p_affine_0_1;
   real_t tmp_1  = p_affine_2_1 + tmp_0;
   real_t tmp_2  = -p_affine_0_0;
   real_t tmp_3  = p_affine_1_0 + tmp_2;
   real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_0 ) * ( p_affine_2_0 + tmp_2 ) );
   real_t tmp_5  = tmp_1 * tmp_4;
   real_t tmp_6  = tmp_4 * ( p_affine_0_1 - p_affine_1_1 );
   real_t tmp_7  = 0.5 * p_affine_10_0;
   real_t tmp_8  = tmp_3 * tmp_4;
   real_t tmp_9  = tmp_4 * ( p_affine_0_0 - p_affine_2_0 );
   real_t tmp_10 = 0.5 * p_affine_10_1;
   real_t tmp_11 = tmp_10 * ( -tmp_8 - tmp_9 ) + tmp_7 * ( -tmp_5 - tmp_6 );
   real_t tmp_12 = -p_affine_3_1;
   real_t tmp_13 = -p_affine_6_1 + p_affine_7_1;
   real_t tmp_14 = p_affine_6_1 + 0.21132486540518713 * tmp_13;
   real_t tmp_15 = tmp_12 + tmp_14;
   real_t tmp_16 = -p_affine_3_0;
   real_t tmp_17 = p_affine_4_0 + tmp_16;
   real_t tmp_18 = p_affine_5_1 + tmp_12;
   real_t tmp_19 = 1.0 / ( tmp_17 * tmp_18 - ( p_affine_4_1 + tmp_12 ) * ( p_affine_5_0 + tmp_16 ) );
   real_t tmp_20 = tmp_17 * tmp_19;
   real_t tmp_21 = tmp_15 * tmp_20;
   real_t tmp_22 = tmp_19 * ( p_affine_3_0 - p_affine_5_0 );
   real_t tmp_23 = tmp_15 * tmp_22;
   real_t tmp_24 = -p_affine_6_0 + p_affine_7_0;
   real_t tmp_25 = p_affine_6_0 + 0.21132486540518713 * tmp_24;
   real_t tmp_26 = tmp_16 + tmp_25;
   real_t tmp_27 = tmp_18 * tmp_19;
   real_t tmp_28 = tmp_26 * tmp_27;
   real_t tmp_29 = tmp_19 * ( p_affine_3_1 - p_affine_4_1 );
   real_t tmp_30 = tmp_26 * tmp_29;
   real_t tmp_31 = -tmp_21 - tmp_23 - tmp_28 - tmp_30 + 1;
   real_t tmp_32 = tmp_10 * ( -tmp_20 - tmp_22 ) + tmp_7 * ( -tmp_27 - tmp_29 );
   real_t tmp_33 = tmp_0 + tmp_14;
   real_t tmp_34 = tmp_33 * tmp_8;
   real_t tmp_35 = tmp_33 * tmp_9;
   real_t tmp_36 = tmp_2 + tmp_25;
   real_t tmp_37 = tmp_36 * tmp_5;
   real_t tmp_38 = tmp_36 * tmp_6;
   real_t tmp_39 = -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1;
   real_t tmp_40 = std::abs( std::pow( ( tmp_13 * tmp_13 ) + ( tmp_24 * tmp_24 ), 1.0 / 2.0 ) );
   real_t tmp_41 = sigma_0 * std::pow( tmp_40, -beta_0 );
   real_t tmp_42 = tmp_39 * tmp_41;
   real_t tmp_43 = 0.5 * tmp_40;
   real_t tmp_44 = p_affine_6_1 + 0.78867513459481287 * tmp_13;
   real_t tmp_45 = tmp_12 + tmp_44;
   real_t tmp_46 = tmp_20 * tmp_45;
   real_t tmp_47 = tmp_22 * tmp_45;
   real_t tmp_48 = p_affine_6_0 + 0.78867513459481287 * tmp_24;
   real_t tmp_49 = tmp_16 + tmp_48;
   real_t tmp_50 = tmp_27 * tmp_49;
   real_t tmp_51 = tmp_29 * tmp_49;
   real_t tmp_52 = -tmp_46 - tmp_47 - tmp_50 - tmp_51 + 1;
   real_t tmp_53 = tmp_0 + tmp_44;
   real_t tmp_54 = tmp_53 * tmp_8;
   real_t tmp_55 = tmp_53 * tmp_9;
   real_t tmp_56 = tmp_2 + tmp_48;
   real_t tmp_57 = tmp_5 * tmp_56;
   real_t tmp_58 = tmp_56 * tmp_6;
   real_t tmp_59 = -tmp_54 - tmp_55 - tmp_57 - tmp_58 + 1;
   real_t tmp_60 = tmp_41 * tmp_59;
   real_t tmp_61 = 0.5 * tmp_40;
   real_t tmp_62 = tmp_23 + tmp_28;
   real_t tmp_63 = tmp_10 * tmp_22 + tmp_27 * tmp_7;
   real_t tmp_64 = tmp_47 + tmp_50;
   real_t tmp_65 = tmp_21 + tmp_30;
   real_t tmp_66 = tmp_10 * tmp_20 + tmp_29 * tmp_7;
   real_t tmp_67 = tmp_46 + tmp_51;
   real_t tmp_68 = tmp_35 + tmp_37;
   real_t tmp_69 = tmp_10 * tmp_9 + tmp_5 * tmp_7;
   real_t tmp_70 = tmp_41 * tmp_68;
   real_t tmp_71 = tmp_55 + tmp_57;
   real_t tmp_72 = tmp_41 * tmp_71;
   real_t tmp_73 = tmp_34 + tmp_38;
   real_t tmp_74 = tmp_10 * tmp_8 + tmp_6 * tmp_7;
   real_t tmp_75 = tmp_41 * tmp_73;
   real_t tmp_76 = tmp_54 + tmp_58;
   real_t tmp_77 = tmp_41 * tmp_76;
   real_t a_0_0  = tmp_43 * ( tmp_11 * tmp_31 - tmp_31 * tmp_42 - tmp_32 * tmp_39 ) +
                  tmp_61 * ( tmp_11 * tmp_52 - tmp_32 * tmp_59 - tmp_52 * tmp_60 );
   real_t a_0_1 = tmp_43 * ( tmp_11 * tmp_62 - tmp_39 * tmp_63 - tmp_42 * tmp_62 ) +
                  tmp_61 * ( tmp_11 * tmp_64 - tmp_59 * tmp_63 - tmp_60 * tmp_64 );
   real_t a_0_2 = tmp_43 * ( tmp_11 * tmp_65 - tmp_39 * tmp_66 - tmp_42 * tmp_65 ) +
                  tmp_61 * ( tmp_11 * tmp_67 - tmp_59 * tmp_66 - tmp_60 * tmp_67 );
   real_t a_1_0 = tmp_43 * ( tmp_31 * tmp_69 - tmp_31 * tmp_70 - tmp_32 * tmp_68 ) +
                  tmp_61 * ( -tmp_32 * tmp_71 + tmp_52 * tmp_69 - tmp_52 * tmp_72 );
   real_t a_1_1 = tmp_43 * ( tmp_62 * tmp_69 - tmp_62 * tmp_70 - tmp_63 * tmp_68 ) +
                  tmp_61 * ( -tmp_63 * tmp_71 + tmp_64 * tmp_69 - tmp_64 * tmp_72 );
   real_t a_1_2 = tmp_43 * ( tmp_65 * tmp_69 - tmp_65 * tmp_70 - tmp_66 * tmp_68 ) +
                  tmp_61 * ( -tmp_66 * tmp_71 + tmp_67 * tmp_69 - tmp_67 * tmp_72 );
   real_t a_2_0 = tmp_43 * ( tmp_31 * tmp_74 - tmp_31 * tmp_75 - tmp_32 * tmp_73 ) +
                  tmp_61 * ( -tmp_32 * tmp_76 + tmp_52 * tmp_74 - tmp_52 * tmp_77 );
   real_t a_2_1 = tmp_43 * ( tmp_62 * tmp_74 - tmp_62 * tmp_75 - tmp_63 * tmp_73 ) +
                  tmp_61 * ( -tmp_63 * tmp_76 + tmp_64 * tmp_74 - tmp_64 * tmp_77 );
   real_t a_2_2 = tmp_43 * ( tmp_65 * tmp_74 - tmp_65 * tmp_75 - tmp_66 * tmp_73 ) +
                  tmp_61 * ( -tmp_66 * tmp_76 + tmp_67 * tmp_74 - tmp_67 * tmp_77 );

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

void DGDiffusionForm_Example::integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&     outwardNormal,
                                                        const DGBasisInfo& trialBasis,
                                                        const DGBasisInfo& testBasis,
                                                        int                trialDegree,
                                                        int                testDegree,
                                                        MatrixXr&          elMat ) const
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

   real_t tmp_0  = -p_affine_0_1;
   real_t tmp_1  = p_affine_1_1 + tmp_0;
   real_t tmp_2  = -p_affine_0_2;
   real_t tmp_3  = p_affine_2_2 + tmp_2;
   real_t tmp_4  = tmp_1 * tmp_3;
   real_t tmp_5  = p_affine_2_1 + tmp_0;
   real_t tmp_6  = p_affine_1_2 + tmp_2;
   real_t tmp_7  = tmp_5 * tmp_6;
   real_t tmp_8  = -p_affine_0_0;
   real_t tmp_9  = p_affine_1_0 + tmp_8;
   real_t tmp_10 = p_affine_3_2 + tmp_2;
   real_t tmp_11 = tmp_10 * tmp_5;
   real_t tmp_12 = p_affine_2_0 + tmp_8;
   real_t tmp_13 = p_affine_3_1 + tmp_0;
   real_t tmp_14 = tmp_13 * tmp_6;
   real_t tmp_15 = p_affine_3_0 + tmp_8;
   real_t tmp_16 = tmp_13 * tmp_3;
   real_t tmp_17 = tmp_1 * tmp_10;
   real_t tmp_18 =
       1.0 / ( tmp_11 * tmp_9 + tmp_12 * tmp_14 - tmp_12 * tmp_17 + tmp_15 * tmp_4 - tmp_15 * tmp_7 - tmp_16 * tmp_9 );
   real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
   real_t tmp_20 = tmp_18 * ( tmp_14 - tmp_17 );
   real_t tmp_21 = tmp_18 * ( tmp_11 - tmp_16 );
   real_t tmp_22 = 0.5 * p_affine_13_0;
   real_t tmp_23 = tmp_18 * ( tmp_12 * tmp_6 - tmp_3 * tmp_9 );
   real_t tmp_24 = tmp_18 * ( tmp_10 * tmp_9 - tmp_15 * tmp_6 );
   real_t tmp_25 = tmp_18 * ( -tmp_10 * tmp_12 + tmp_15 * tmp_3 );
   real_t tmp_26 = 0.5 * p_affine_13_1;
   real_t tmp_27 = tmp_18 * ( -tmp_1 * tmp_12 + tmp_5 * tmp_9 );
   real_t tmp_28 = tmp_18 * ( tmp_1 * tmp_15 - tmp_13 * tmp_9 );
   real_t tmp_29 = tmp_18 * ( tmp_12 * tmp_13 - tmp_15 * tmp_5 );
   real_t tmp_30 = 0.5 * p_affine_13_2;
   real_t tmp_31 =
       tmp_22 * ( -tmp_19 - tmp_20 - tmp_21 ) + tmp_26 * ( -tmp_23 - tmp_24 - tmp_25 ) + tmp_30 * ( -tmp_27 - tmp_28 - tmp_29 );
   real_t tmp_32 = -p_affine_4_2;
   real_t tmp_33 = p_affine_8_2 + tmp_32;
   real_t tmp_34 = -p_affine_8_2;
   real_t tmp_35 = p_affine_9_2 + tmp_34;
   real_t tmp_36 = p_affine_10_2 + tmp_34;
   real_t tmp_37 = 0.091576213509770743 * tmp_35 + 0.81684757298045851 * tmp_36;
   real_t tmp_38 = tmp_33 + tmp_37;
   real_t tmp_39 = -p_affine_4_0;
   real_t tmp_40 = p_affine_5_0 + tmp_39;
   real_t tmp_41 = -p_affine_4_1;
   real_t tmp_42 = p_affine_6_1 + tmp_41;
   real_t tmp_43 = tmp_40 * tmp_42;
   real_t tmp_44 = p_affine_6_0 + tmp_39;
   real_t tmp_45 = p_affine_5_1 + tmp_41;
   real_t tmp_46 = tmp_44 * tmp_45;
   real_t tmp_47 = p_affine_7_2 + tmp_32;
   real_t tmp_48 = p_affine_7_1 + tmp_41;
   real_t tmp_49 = p_affine_5_2 + tmp_32;
   real_t tmp_50 = tmp_44 * tmp_49;
   real_t tmp_51 = p_affine_7_0 + tmp_39;
   real_t tmp_52 = p_affine_6_2 + tmp_32;
   real_t tmp_53 = tmp_45 * tmp_52;
   real_t tmp_54 = tmp_40 * tmp_48;
   real_t tmp_55 = tmp_49 * tmp_51;
   real_t tmp_56 =
       1.0 / ( -tmp_42 * tmp_55 + tmp_43 * tmp_47 - tmp_46 * tmp_47 + tmp_48 * tmp_50 + tmp_51 * tmp_53 - tmp_52 * tmp_54 );
   real_t tmp_57 = tmp_56 * ( tmp_43 - tmp_46 );
   real_t tmp_58 = tmp_38 * tmp_57;
   real_t tmp_59 = tmp_56 * ( tmp_45 * tmp_51 - tmp_54 );
   real_t tmp_60 = tmp_38 * tmp_59;
   real_t tmp_61 = p_affine_8_1 + tmp_41;
   real_t tmp_62 = -p_affine_8_1;
   real_t tmp_63 = p_affine_9_1 + tmp_62;
   real_t tmp_64 = p_affine_10_1 + tmp_62;
   real_t tmp_65 = 0.091576213509770743 * tmp_63 + 0.81684757298045851 * tmp_64;
   real_t tmp_66 = tmp_61 + tmp_65;
   real_t tmp_67 = tmp_56 * ( -tmp_40 * tmp_52 + tmp_50 );
   real_t tmp_68 = tmp_66 * tmp_67;
   real_t tmp_69 = tmp_56 * ( tmp_40 * tmp_47 - tmp_55 );
   real_t tmp_70 = tmp_66 * tmp_69;
   real_t tmp_71 = tmp_56 * ( -tmp_42 * tmp_51 + tmp_44 * tmp_48 );
   real_t tmp_72 = tmp_38 * tmp_71;
   real_t tmp_73 = tmp_56 * ( -tmp_44 * tmp_47 + tmp_51 * tmp_52 );
   real_t tmp_74 = tmp_66 * tmp_73;
   real_t tmp_75 = p_affine_8_0 + tmp_39;
   real_t tmp_76 = -p_affine_8_0;
   real_t tmp_77 = p_affine_9_0 + tmp_76;
   real_t tmp_78 = p_affine_10_0 + tmp_76;
   real_t tmp_79 = 0.091576213509770743 * tmp_77 + 0.81684757298045851 * tmp_78;
   real_t tmp_80 = tmp_75 + tmp_79;
   real_t tmp_81 = tmp_56 * ( -tmp_42 * tmp_49 + tmp_53 );
   real_t tmp_82 = tmp_80 * tmp_81;
   real_t tmp_83 = tmp_56 * ( -tmp_45 * tmp_47 + tmp_48 * tmp_49 );
   real_t tmp_84 = tmp_80 * tmp_83;
   real_t tmp_85 = tmp_56 * ( tmp_42 * tmp_47 - tmp_48 * tmp_52 );
   real_t tmp_86 = tmp_80 * tmp_85;
   real_t tmp_87 = -tmp_58 - tmp_60 - tmp_68 - tmp_70 - tmp_72 - tmp_74 - tmp_82 - tmp_84 - tmp_86 + 1;
   real_t tmp_88 =
       tmp_22 * ( -tmp_81 - tmp_83 - tmp_85 ) + tmp_26 * ( -tmp_67 - tmp_69 - tmp_73 ) + tmp_30 * ( -tmp_57 - tmp_59 - tmp_71 );
   real_t tmp_89  = p_affine_8_2 + tmp_2;
   real_t tmp_90  = tmp_37 + tmp_89;
   real_t tmp_91  = tmp_27 * tmp_90;
   real_t tmp_92  = tmp_28 * tmp_90;
   real_t tmp_93  = p_affine_8_1 + tmp_0;
   real_t tmp_94  = tmp_65 + tmp_93;
   real_t tmp_95  = tmp_23 * tmp_94;
   real_t tmp_96  = tmp_24 * tmp_94;
   real_t tmp_97  = tmp_29 * tmp_90;
   real_t tmp_98  = tmp_25 * tmp_94;
   real_t tmp_99  = p_affine_8_0 + tmp_8;
   real_t tmp_100 = tmp_79 + tmp_99;
   real_t tmp_101 = tmp_100 * tmp_19;
   real_t tmp_102 = tmp_100 * tmp_20;
   real_t tmp_103 = tmp_100 * tmp_21;
   real_t tmp_104 = -tmp_101 - tmp_102 - tmp_103 - tmp_91 - tmp_92 - tmp_95 - tmp_96 - tmp_97 - tmp_98 + 1;
   real_t tmp_105 = p_affine_8_1 - p_affine_9_1;
   real_t tmp_106 = p_affine_8_0 - p_affine_9_0;
   real_t tmp_107 = p_affine_8_2 - p_affine_9_2;
   real_t tmp_108 =
       std::pow( ( std::abs( tmp_105 * tmp_36 - tmp_107 * tmp_64 ) * std::abs( tmp_105 * tmp_36 - tmp_107 * tmp_64 ) ) +
                     ( std::abs( tmp_105 * tmp_78 - tmp_106 * tmp_64 ) * std::abs( tmp_105 * tmp_78 - tmp_106 * tmp_64 ) ) +
                     ( std::abs( tmp_106 * tmp_36 - tmp_107 * tmp_78 ) * std::abs( tmp_106 * tmp_36 - tmp_107 * tmp_78 ) ),
                 1.0 / 2.0 );
   real_t tmp_109 = sigma_0 * std::pow( 0.5 * tmp_108, -beta_0 );
   real_t tmp_110 = tmp_104 * tmp_109;
   real_t tmp_111 = 1.0 * tmp_108;
   real_t tmp_112 = 0.054975871827660928 * tmp_111;
   real_t tmp_113 = 0.44594849091596489 * tmp_35 + 0.10810301816807022 * tmp_36;
   real_t tmp_114 = tmp_113 + tmp_33;
   real_t tmp_115 = tmp_114 * tmp_57;
   real_t tmp_116 = tmp_114 * tmp_59;
   real_t tmp_117 = 0.44594849091596489 * tmp_63 + 0.10810301816807022 * tmp_64;
   real_t tmp_118 = tmp_117 + tmp_61;
   real_t tmp_119 = tmp_118 * tmp_67;
   real_t tmp_120 = tmp_118 * tmp_69;
   real_t tmp_121 = tmp_114 * tmp_71;
   real_t tmp_122 = tmp_118 * tmp_73;
   real_t tmp_123 = 0.44594849091596489 * tmp_77 + 0.10810301816807022 * tmp_78;
   real_t tmp_124 = tmp_123 + tmp_75;
   real_t tmp_125 = tmp_124 * tmp_81;
   real_t tmp_126 = tmp_124 * tmp_83;
   real_t tmp_127 = tmp_124 * tmp_85;
   real_t tmp_128 = -tmp_115 - tmp_116 - tmp_119 - tmp_120 - tmp_121 - tmp_122 - tmp_125 - tmp_126 - tmp_127 + 1;
   real_t tmp_129 = tmp_113 + tmp_89;
   real_t tmp_130 = tmp_129 * tmp_27;
   real_t tmp_131 = tmp_129 * tmp_28;
   real_t tmp_132 = tmp_117 + tmp_93;
   real_t tmp_133 = tmp_132 * tmp_23;
   real_t tmp_134 = tmp_132 * tmp_24;
   real_t tmp_135 = tmp_129 * tmp_29;
   real_t tmp_136 = tmp_132 * tmp_25;
   real_t tmp_137 = tmp_123 + tmp_99;
   real_t tmp_138 = tmp_137 * tmp_19;
   real_t tmp_139 = tmp_137 * tmp_20;
   real_t tmp_140 = tmp_137 * tmp_21;
   real_t tmp_141 = -tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 - tmp_136 - tmp_138 - tmp_139 - tmp_140 + 1;
   real_t tmp_142 = tmp_109 * tmp_141;
   real_t tmp_143 = 0.11169079483900572 * tmp_111;
   real_t tmp_144 = 0.81684757298045851 * tmp_35 + 0.091576213509770743 * tmp_36;
   real_t tmp_145 = tmp_144 + tmp_33;
   real_t tmp_146 = tmp_145 * tmp_57;
   real_t tmp_147 = tmp_145 * tmp_59;
   real_t tmp_148 = 0.81684757298045851 * tmp_63 + 0.091576213509770743 * tmp_64;
   real_t tmp_149 = tmp_148 + tmp_61;
   real_t tmp_150 = tmp_149 * tmp_67;
   real_t tmp_151 = tmp_149 * tmp_69;
   real_t tmp_152 = tmp_145 * tmp_71;
   real_t tmp_153 = tmp_149 * tmp_73;
   real_t tmp_154 = 0.81684757298045851 * tmp_77 + 0.091576213509770743 * tmp_78;
   real_t tmp_155 = tmp_154 + tmp_75;
   real_t tmp_156 = tmp_155 * tmp_81;
   real_t tmp_157 = tmp_155 * tmp_83;
   real_t tmp_158 = tmp_155 * tmp_85;
   real_t tmp_159 = -tmp_146 - tmp_147 - tmp_150 - tmp_151 - tmp_152 - tmp_153 - tmp_156 - tmp_157 - tmp_158 + 1;
   real_t tmp_160 = tmp_144 + tmp_89;
   real_t tmp_161 = tmp_160 * tmp_27;
   real_t tmp_162 = tmp_160 * tmp_28;
   real_t tmp_163 = tmp_148 + tmp_93;
   real_t tmp_164 = tmp_163 * tmp_23;
   real_t tmp_165 = tmp_163 * tmp_24;
   real_t tmp_166 = tmp_160 * tmp_29;
   real_t tmp_167 = tmp_163 * tmp_25;
   real_t tmp_168 = tmp_154 + tmp_99;
   real_t tmp_169 = tmp_168 * tmp_19;
   real_t tmp_170 = tmp_168 * tmp_20;
   real_t tmp_171 = tmp_168 * tmp_21;
   real_t tmp_172 = -tmp_161 - tmp_162 - tmp_164 - tmp_165 - tmp_166 - tmp_167 - tmp_169 - tmp_170 - tmp_171 + 1;
   real_t tmp_173 = tmp_109 * tmp_172;
   real_t tmp_174 = 0.054975871827660928 * tmp_111;
   real_t tmp_175 = 0.10810301816807022 * tmp_35 + 0.44594849091596489 * tmp_36;
   real_t tmp_176 = tmp_175 + tmp_33;
   real_t tmp_177 = tmp_176 * tmp_57;
   real_t tmp_178 = tmp_176 * tmp_59;
   real_t tmp_179 = 0.10810301816807022 * tmp_63 + 0.44594849091596489 * tmp_64;
   real_t tmp_180 = tmp_179 + tmp_61;
   real_t tmp_181 = tmp_180 * tmp_67;
   real_t tmp_182 = tmp_180 * tmp_69;
   real_t tmp_183 = tmp_176 * tmp_71;
   real_t tmp_184 = tmp_180 * tmp_73;
   real_t tmp_185 = 0.10810301816807022 * tmp_77 + 0.44594849091596489 * tmp_78;
   real_t tmp_186 = tmp_185 + tmp_75;
   real_t tmp_187 = tmp_186 * tmp_81;
   real_t tmp_188 = tmp_186 * tmp_83;
   real_t tmp_189 = tmp_186 * tmp_85;
   real_t tmp_190 = -tmp_177 - tmp_178 - tmp_181 - tmp_182 - tmp_183 - tmp_184 - tmp_187 - tmp_188 - tmp_189 + 1;
   real_t tmp_191 = tmp_175 + tmp_89;
   real_t tmp_192 = tmp_191 * tmp_27;
   real_t tmp_193 = tmp_191 * tmp_28;
   real_t tmp_194 = tmp_179 + tmp_93;
   real_t tmp_195 = tmp_194 * tmp_23;
   real_t tmp_196 = tmp_194 * tmp_24;
   real_t tmp_197 = tmp_191 * tmp_29;
   real_t tmp_198 = tmp_194 * tmp_25;
   real_t tmp_199 = tmp_185 + tmp_99;
   real_t tmp_200 = tmp_19 * tmp_199;
   real_t tmp_201 = tmp_199 * tmp_20;
   real_t tmp_202 = tmp_199 * tmp_21;
   real_t tmp_203 = -tmp_192 - tmp_193 - tmp_195 - tmp_196 - tmp_197 - tmp_198 - tmp_200 - tmp_201 - tmp_202 + 1;
   real_t tmp_204 = tmp_109 * tmp_203;
   real_t tmp_205 = 0.11169079483900572 * tmp_111;
   real_t tmp_206 = 0.091576213509770743 * tmp_35 + 0.091576213509770743 * tmp_36;
   real_t tmp_207 = tmp_206 + tmp_33;
   real_t tmp_208 = tmp_207 * tmp_57;
   real_t tmp_209 = tmp_207 * tmp_59;
   real_t tmp_210 = 0.091576213509770743 * tmp_63 + 0.091576213509770743 * tmp_64;
   real_t tmp_211 = tmp_210 + tmp_61;
   real_t tmp_212 = tmp_211 * tmp_67;
   real_t tmp_213 = tmp_211 * tmp_69;
   real_t tmp_214 = tmp_207 * tmp_71;
   real_t tmp_215 = tmp_211 * tmp_73;
   real_t tmp_216 = 0.091576213509770743 * tmp_77 + 0.091576213509770743 * tmp_78;
   real_t tmp_217 = tmp_216 + tmp_75;
   real_t tmp_218 = tmp_217 * tmp_81;
   real_t tmp_219 = tmp_217 * tmp_83;
   real_t tmp_220 = tmp_217 * tmp_85;
   real_t tmp_221 = -tmp_208 - tmp_209 - tmp_212 - tmp_213 - tmp_214 - tmp_215 - tmp_218 - tmp_219 - tmp_220 + 1;
   real_t tmp_222 = tmp_206 + tmp_89;
   real_t tmp_223 = tmp_222 * tmp_27;
   real_t tmp_224 = tmp_222 * tmp_28;
   real_t tmp_225 = tmp_210 + tmp_93;
   real_t tmp_226 = tmp_225 * tmp_23;
   real_t tmp_227 = tmp_225 * tmp_24;
   real_t tmp_228 = tmp_222 * tmp_29;
   real_t tmp_229 = tmp_225 * tmp_25;
   real_t tmp_230 = tmp_216 + tmp_99;
   real_t tmp_231 = tmp_19 * tmp_230;
   real_t tmp_232 = tmp_20 * tmp_230;
   real_t tmp_233 = tmp_21 * tmp_230;
   real_t tmp_234 = -tmp_223 - tmp_224 - tmp_226 - tmp_227 - tmp_228 - tmp_229 - tmp_231 - tmp_232 - tmp_233 + 1;
   real_t tmp_235 = tmp_109 * tmp_234;
   real_t tmp_236 = 0.054975871827660928 * tmp_111;
   real_t tmp_237 = 0.44594849091596489 * tmp_35 + 0.44594849091596489 * tmp_36;
   real_t tmp_238 = tmp_237 + tmp_33;
   real_t tmp_239 = tmp_238 * tmp_57;
   real_t tmp_240 = tmp_238 * tmp_59;
   real_t tmp_241 = 0.44594849091596489 * tmp_63 + 0.44594849091596489 * tmp_64;
   real_t tmp_242 = tmp_241 + tmp_61;
   real_t tmp_243 = tmp_242 * tmp_67;
   real_t tmp_244 = tmp_242 * tmp_69;
   real_t tmp_245 = tmp_238 * tmp_71;
   real_t tmp_246 = tmp_242 * tmp_73;
   real_t tmp_247 = 0.44594849091596489 * tmp_77 + 0.44594849091596489 * tmp_78;
   real_t tmp_248 = tmp_247 + tmp_75;
   real_t tmp_249 = tmp_248 * tmp_81;
   real_t tmp_250 = tmp_248 * tmp_83;
   real_t tmp_251 = tmp_248 * tmp_85;
   real_t tmp_252 = -tmp_239 - tmp_240 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_249 - tmp_250 - tmp_251 + 1;
   real_t tmp_253 = tmp_237 + tmp_89;
   real_t tmp_254 = tmp_253 * tmp_27;
   real_t tmp_255 = tmp_253 * tmp_28;
   real_t tmp_256 = tmp_241 + tmp_93;
   real_t tmp_257 = tmp_23 * tmp_256;
   real_t tmp_258 = tmp_24 * tmp_256;
   real_t tmp_259 = tmp_253 * tmp_29;
   real_t tmp_260 = tmp_25 * tmp_256;
   real_t tmp_261 = tmp_247 + tmp_99;
   real_t tmp_262 = tmp_19 * tmp_261;
   real_t tmp_263 = tmp_20 * tmp_261;
   real_t tmp_264 = tmp_21 * tmp_261;
   real_t tmp_265 = -tmp_254 - tmp_255 - tmp_257 - tmp_258 - tmp_259 - tmp_260 - tmp_262 - tmp_263 - tmp_264 + 1;
   real_t tmp_266 = tmp_109 * tmp_265;
   real_t tmp_267 = 0.11169079483900572 * tmp_111;
   real_t tmp_268 = tmp_72 + tmp_74 + tmp_86;
   real_t tmp_269 = tmp_22 * tmp_85 + tmp_26 * tmp_73 + tmp_30 * tmp_71;
   real_t tmp_270 = tmp_121 + tmp_122 + tmp_127;
   real_t tmp_271 = tmp_152 + tmp_153 + tmp_158;
   real_t tmp_272 = tmp_183 + tmp_184 + tmp_189;
   real_t tmp_273 = tmp_214 + tmp_215 + tmp_220;
   real_t tmp_274 = tmp_245 + tmp_246 + tmp_251;
   real_t tmp_275 = tmp_60 + tmp_70 + tmp_84;
   real_t tmp_276 = tmp_22 * tmp_83 + tmp_26 * tmp_69 + tmp_30 * tmp_59;
   real_t tmp_277 = tmp_116 + tmp_120 + tmp_126;
   real_t tmp_278 = tmp_147 + tmp_151 + tmp_157;
   real_t tmp_279 = tmp_178 + tmp_182 + tmp_188;
   real_t tmp_280 = tmp_209 + tmp_213 + tmp_219;
   real_t tmp_281 = tmp_240 + tmp_244 + tmp_250;
   real_t tmp_282 = tmp_58 + tmp_68 + tmp_82;
   real_t tmp_283 = tmp_22 * tmp_81 + tmp_26 * tmp_67 + tmp_30 * tmp_57;
   real_t tmp_284 = tmp_115 + tmp_119 + tmp_125;
   real_t tmp_285 = tmp_146 + tmp_150 + tmp_156;
   real_t tmp_286 = tmp_177 + tmp_181 + tmp_187;
   real_t tmp_287 = tmp_208 + tmp_212 + tmp_218;
   real_t tmp_288 = tmp_239 + tmp_243 + tmp_249;
   real_t tmp_289 = tmp_103 + tmp_97 + tmp_98;
   real_t tmp_290 = tmp_21 * tmp_22 + tmp_25 * tmp_26 + tmp_29 * tmp_30;
   real_t tmp_291 = tmp_109 * tmp_289;
   real_t tmp_292 = tmp_135 + tmp_136 + tmp_140;
   real_t tmp_293 = tmp_109 * tmp_292;
   real_t tmp_294 = tmp_166 + tmp_167 + tmp_171;
   real_t tmp_295 = tmp_109 * tmp_294;
   real_t tmp_296 = tmp_197 + tmp_198 + tmp_202;
   real_t tmp_297 = tmp_109 * tmp_296;
   real_t tmp_298 = tmp_228 + tmp_229 + tmp_233;
   real_t tmp_299 = tmp_109 * tmp_298;
   real_t tmp_300 = tmp_259 + tmp_260 + tmp_264;
   real_t tmp_301 = tmp_109 * tmp_300;
   real_t tmp_302 = tmp_102 + tmp_92 + tmp_96;
   real_t tmp_303 = tmp_20 * tmp_22 + tmp_24 * tmp_26 + tmp_28 * tmp_30;
   real_t tmp_304 = tmp_109 * tmp_302;
   real_t tmp_305 = tmp_131 + tmp_134 + tmp_139;
   real_t tmp_306 = tmp_109 * tmp_305;
   real_t tmp_307 = tmp_162 + tmp_165 + tmp_170;
   real_t tmp_308 = tmp_109 * tmp_307;
   real_t tmp_309 = tmp_193 + tmp_196 + tmp_201;
   real_t tmp_310 = tmp_109 * tmp_309;
   real_t tmp_311 = tmp_224 + tmp_227 + tmp_232;
   real_t tmp_312 = tmp_109 * tmp_311;
   real_t tmp_313 = tmp_255 + tmp_258 + tmp_263;
   real_t tmp_314 = tmp_109 * tmp_313;
   real_t tmp_315 = tmp_101 + tmp_91 + tmp_95;
   real_t tmp_316 = tmp_19 * tmp_22 + tmp_23 * tmp_26 + tmp_27 * tmp_30;
   real_t tmp_317 = tmp_109 * tmp_315;
   real_t tmp_318 = tmp_130 + tmp_133 + tmp_138;
   real_t tmp_319 = tmp_109 * tmp_318;
   real_t tmp_320 = tmp_161 + tmp_164 + tmp_169;
   real_t tmp_321 = tmp_109 * tmp_320;
   real_t tmp_322 = tmp_192 + tmp_195 + tmp_200;
   real_t tmp_323 = tmp_109 * tmp_322;
   real_t tmp_324 = tmp_223 + tmp_226 + tmp_231;
   real_t tmp_325 = tmp_109 * tmp_324;
   real_t tmp_326 = tmp_254 + tmp_257 + tmp_262;
   real_t tmp_327 = tmp_109 * tmp_326;
   real_t a_0_0   = tmp_112 * ( -tmp_104 * tmp_88 - tmp_110 * tmp_87 + tmp_31 * tmp_87 ) +
                  tmp_143 * ( -tmp_128 * tmp_142 + tmp_128 * tmp_31 - tmp_141 * tmp_88 ) +
                  tmp_174 * ( -tmp_159 * tmp_173 + tmp_159 * tmp_31 - tmp_172 * tmp_88 ) +
                  tmp_205 * ( -tmp_190 * tmp_204 + tmp_190 * tmp_31 - tmp_203 * tmp_88 ) +
                  tmp_236 * ( -tmp_221 * tmp_235 + tmp_221 * tmp_31 - tmp_234 * tmp_88 ) +
                  tmp_267 * ( -tmp_252 * tmp_266 + tmp_252 * tmp_31 - tmp_265 * tmp_88 );
   real_t a_0_1 = tmp_112 * ( -tmp_104 * tmp_269 - tmp_110 * tmp_268 + tmp_268 * tmp_31 ) +
                  tmp_143 * ( -tmp_141 * tmp_269 - tmp_142 * tmp_270 + tmp_270 * tmp_31 ) +
                  tmp_174 * ( -tmp_172 * tmp_269 - tmp_173 * tmp_271 + tmp_271 * tmp_31 ) +
                  tmp_205 * ( -tmp_203 * tmp_269 - tmp_204 * tmp_272 + tmp_272 * tmp_31 ) +
                  tmp_236 * ( -tmp_234 * tmp_269 - tmp_235 * tmp_273 + tmp_273 * tmp_31 ) +
                  tmp_267 * ( -tmp_265 * tmp_269 - tmp_266 * tmp_274 + tmp_274 * tmp_31 );
   real_t a_0_2 = tmp_112 * ( -tmp_104 * tmp_276 - tmp_110 * tmp_275 + tmp_275 * tmp_31 ) +
                  tmp_143 * ( -tmp_141 * tmp_276 - tmp_142 * tmp_277 + tmp_277 * tmp_31 ) +
                  tmp_174 * ( -tmp_172 * tmp_276 - tmp_173 * tmp_278 + tmp_278 * tmp_31 ) +
                  tmp_205 * ( -tmp_203 * tmp_276 - tmp_204 * tmp_279 + tmp_279 * tmp_31 ) +
                  tmp_236 * ( -tmp_234 * tmp_276 - tmp_235 * tmp_280 + tmp_280 * tmp_31 ) +
                  tmp_267 * ( -tmp_265 * tmp_276 - tmp_266 * tmp_281 + tmp_281 * tmp_31 );
   real_t a_0_3 = tmp_112 * ( -tmp_104 * tmp_283 - tmp_110 * tmp_282 + tmp_282 * tmp_31 ) +
                  tmp_143 * ( -tmp_141 * tmp_283 - tmp_142 * tmp_284 + tmp_284 * tmp_31 ) +
                  tmp_174 * ( -tmp_172 * tmp_283 - tmp_173 * tmp_285 + tmp_285 * tmp_31 ) +
                  tmp_205 * ( -tmp_203 * tmp_283 - tmp_204 * tmp_286 + tmp_286 * tmp_31 ) +
                  tmp_236 * ( -tmp_234 * tmp_283 - tmp_235 * tmp_287 + tmp_287 * tmp_31 ) +
                  tmp_267 * ( -tmp_265 * tmp_283 - tmp_266 * tmp_288 + tmp_288 * tmp_31 );
   real_t a_1_0 = tmp_112 * ( -tmp_289 * tmp_88 + tmp_290 * tmp_87 - tmp_291 * tmp_87 ) +
                  tmp_143 * ( tmp_128 * tmp_290 - tmp_128 * tmp_293 - tmp_292 * tmp_88 ) +
                  tmp_174 * ( tmp_159 * tmp_290 - tmp_159 * tmp_295 - tmp_294 * tmp_88 ) +
                  tmp_205 * ( tmp_190 * tmp_290 - tmp_190 * tmp_297 - tmp_296 * tmp_88 ) +
                  tmp_236 * ( tmp_221 * tmp_290 - tmp_221 * tmp_299 - tmp_298 * tmp_88 ) +
                  tmp_267 * ( tmp_252 * tmp_290 - tmp_252 * tmp_301 - tmp_300 * tmp_88 );
   real_t a_1_1 = tmp_112 * ( tmp_268 * tmp_290 - tmp_268 * tmp_291 - tmp_269 * tmp_289 ) +
                  tmp_143 * ( -tmp_269 * tmp_292 + tmp_270 * tmp_290 - tmp_270 * tmp_293 ) +
                  tmp_174 * ( -tmp_269 * tmp_294 + tmp_271 * tmp_290 - tmp_271 * tmp_295 ) +
                  tmp_205 * ( -tmp_269 * tmp_296 + tmp_272 * tmp_290 - tmp_272 * tmp_297 ) +
                  tmp_236 * ( -tmp_269 * tmp_298 + tmp_273 * tmp_290 - tmp_273 * tmp_299 ) +
                  tmp_267 * ( -tmp_269 * tmp_300 + tmp_274 * tmp_290 - tmp_274 * tmp_301 );
   real_t a_1_2 = tmp_112 * ( tmp_275 * tmp_290 - tmp_275 * tmp_291 - tmp_276 * tmp_289 ) +
                  tmp_143 * ( -tmp_276 * tmp_292 + tmp_277 * tmp_290 - tmp_277 * tmp_293 ) +
                  tmp_174 * ( -tmp_276 * tmp_294 + tmp_278 * tmp_290 - tmp_278 * tmp_295 ) +
                  tmp_205 * ( -tmp_276 * tmp_296 + tmp_279 * tmp_290 - tmp_279 * tmp_297 ) +
                  tmp_236 * ( -tmp_276 * tmp_298 + tmp_280 * tmp_290 - tmp_280 * tmp_299 ) +
                  tmp_267 * ( -tmp_276 * tmp_300 + tmp_281 * tmp_290 - tmp_281 * tmp_301 );
   real_t a_1_3 = tmp_112 * ( tmp_282 * tmp_290 - tmp_282 * tmp_291 - tmp_283 * tmp_289 ) +
                  tmp_143 * ( -tmp_283 * tmp_292 + tmp_284 * tmp_290 - tmp_284 * tmp_293 ) +
                  tmp_174 * ( -tmp_283 * tmp_294 + tmp_285 * tmp_290 - tmp_285 * tmp_295 ) +
                  tmp_205 * ( -tmp_283 * tmp_296 + tmp_286 * tmp_290 - tmp_286 * tmp_297 ) +
                  tmp_236 * ( -tmp_283 * tmp_298 + tmp_287 * tmp_290 - tmp_287 * tmp_299 ) +
                  tmp_267 * ( -tmp_283 * tmp_300 + tmp_288 * tmp_290 - tmp_288 * tmp_301 );
   real_t a_2_0 = tmp_112 * ( -tmp_302 * tmp_88 + tmp_303 * tmp_87 - tmp_304 * tmp_87 ) +
                  tmp_143 * ( tmp_128 * tmp_303 - tmp_128 * tmp_306 - tmp_305 * tmp_88 ) +
                  tmp_174 * ( tmp_159 * tmp_303 - tmp_159 * tmp_308 - tmp_307 * tmp_88 ) +
                  tmp_205 * ( tmp_190 * tmp_303 - tmp_190 * tmp_310 - tmp_309 * tmp_88 ) +
                  tmp_236 * ( tmp_221 * tmp_303 - tmp_221 * tmp_312 - tmp_311 * tmp_88 ) +
                  tmp_267 * ( tmp_252 * tmp_303 - tmp_252 * tmp_314 - tmp_313 * tmp_88 );
   real_t a_2_1 = tmp_112 * ( tmp_268 * tmp_303 - tmp_268 * tmp_304 - tmp_269 * tmp_302 ) +
                  tmp_143 * ( -tmp_269 * tmp_305 + tmp_270 * tmp_303 - tmp_270 * tmp_306 ) +
                  tmp_174 * ( -tmp_269 * tmp_307 + tmp_271 * tmp_303 - tmp_271 * tmp_308 ) +
                  tmp_205 * ( -tmp_269 * tmp_309 + tmp_272 * tmp_303 - tmp_272 * tmp_310 ) +
                  tmp_236 * ( -tmp_269 * tmp_311 + tmp_273 * tmp_303 - tmp_273 * tmp_312 ) +
                  tmp_267 * ( -tmp_269 * tmp_313 + tmp_274 * tmp_303 - tmp_274 * tmp_314 );
   real_t a_2_2 = tmp_112 * ( tmp_275 * tmp_303 - tmp_275 * tmp_304 - tmp_276 * tmp_302 ) +
                  tmp_143 * ( -tmp_276 * tmp_305 + tmp_277 * tmp_303 - tmp_277 * tmp_306 ) +
                  tmp_174 * ( -tmp_276 * tmp_307 + tmp_278 * tmp_303 - tmp_278 * tmp_308 ) +
                  tmp_205 * ( -tmp_276 * tmp_309 + tmp_279 * tmp_303 - tmp_279 * tmp_310 ) +
                  tmp_236 * ( -tmp_276 * tmp_311 + tmp_280 * tmp_303 - tmp_280 * tmp_312 ) +
                  tmp_267 * ( -tmp_276 * tmp_313 + tmp_281 * tmp_303 - tmp_281 * tmp_314 );
   real_t a_2_3 = tmp_112 * ( tmp_282 * tmp_303 - tmp_282 * tmp_304 - tmp_283 * tmp_302 ) +
                  tmp_143 * ( -tmp_283 * tmp_305 + tmp_284 * tmp_303 - tmp_284 * tmp_306 ) +
                  tmp_174 * ( -tmp_283 * tmp_307 + tmp_285 * tmp_303 - tmp_285 * tmp_308 ) +
                  tmp_205 * ( -tmp_283 * tmp_309 + tmp_286 * tmp_303 - tmp_286 * tmp_310 ) +
                  tmp_236 * ( -tmp_283 * tmp_311 + tmp_287 * tmp_303 - tmp_287 * tmp_312 ) +
                  tmp_267 * ( -tmp_283 * tmp_313 + tmp_288 * tmp_303 - tmp_288 * tmp_314 );
   real_t a_3_0 = tmp_112 * ( -tmp_315 * tmp_88 + tmp_316 * tmp_87 - tmp_317 * tmp_87 ) +
                  tmp_143 * ( tmp_128 * tmp_316 - tmp_128 * tmp_319 - tmp_318 * tmp_88 ) +
                  tmp_174 * ( tmp_159 * tmp_316 - tmp_159 * tmp_321 - tmp_320 * tmp_88 ) +
                  tmp_205 * ( tmp_190 * tmp_316 - tmp_190 * tmp_323 - tmp_322 * tmp_88 ) +
                  tmp_236 * ( tmp_221 * tmp_316 - tmp_221 * tmp_325 - tmp_324 * tmp_88 ) +
                  tmp_267 * ( tmp_252 * tmp_316 - tmp_252 * tmp_327 - tmp_326 * tmp_88 );
   real_t a_3_1 = tmp_112 * ( tmp_268 * tmp_316 - tmp_268 * tmp_317 - tmp_269 * tmp_315 ) +
                  tmp_143 * ( -tmp_269 * tmp_318 + tmp_270 * tmp_316 - tmp_270 * tmp_319 ) +
                  tmp_174 * ( -tmp_269 * tmp_320 + tmp_271 * tmp_316 - tmp_271 * tmp_321 ) +
                  tmp_205 * ( -tmp_269 * tmp_322 + tmp_272 * tmp_316 - tmp_272 * tmp_323 ) +
                  tmp_236 * ( -tmp_269 * tmp_324 + tmp_273 * tmp_316 - tmp_273 * tmp_325 ) +
                  tmp_267 * ( -tmp_269 * tmp_326 + tmp_274 * tmp_316 - tmp_274 * tmp_327 );
   real_t a_3_2 = tmp_112 * ( tmp_275 * tmp_316 - tmp_275 * tmp_317 - tmp_276 * tmp_315 ) +
                  tmp_143 * ( -tmp_276 * tmp_318 + tmp_277 * tmp_316 - tmp_277 * tmp_319 ) +
                  tmp_174 * ( -tmp_276 * tmp_320 + tmp_278 * tmp_316 - tmp_278 * tmp_321 ) +
                  tmp_205 * ( -tmp_276 * tmp_322 + tmp_279 * tmp_316 - tmp_279 * tmp_323 ) +
                  tmp_236 * ( -tmp_276 * tmp_324 + tmp_280 * tmp_316 - tmp_280 * tmp_325 ) +
                  tmp_267 * ( -tmp_276 * tmp_326 + tmp_281 * tmp_316 - tmp_281 * tmp_327 );
   real_t a_3_3 = tmp_112 * ( tmp_282 * tmp_316 - tmp_282 * tmp_317 - tmp_283 * tmp_315 ) +
                  tmp_143 * ( -tmp_283 * tmp_318 + tmp_284 * tmp_316 - tmp_284 * tmp_319 ) +
                  tmp_174 * ( -tmp_283 * tmp_320 + tmp_285 * tmp_316 - tmp_285 * tmp_321 ) +
                  tmp_205 * ( -tmp_283 * tmp_322 + tmp_286 * tmp_316 - tmp_286 * tmp_323 ) +
                  tmp_236 * ( -tmp_283 * tmp_324 + tmp_287 * tmp_316 - tmp_287 * tmp_325 ) +
                  tmp_267 * ( -tmp_283 * tmp_326 + tmp_288 * tmp_316 - tmp_288 * tmp_327 );

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

void DGDiffusionForm_Example::integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                                 const std::vector< Point3D >& coordsFacet,
                                                                 const Point3D&,
                                                                 const Point3D&     outwardNormal,
                                                                 const DGBasisInfo& trialBasis,
                                                                 const DGBasisInfo& testBasis,
                                                                 int                trialDegree,
                                                                 int                testDegree,
                                                                 MatrixXr&          elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 2, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 2, uint_c( trialDegree ) ) ) );

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

   real_t tmp_0  = -p_affine_6_1 + p_affine_7_1;
   real_t tmp_1  = -p_affine_0_1;
   real_t tmp_2  = p_affine_6_1 + tmp_1;
   real_t tmp_3  = 0.21132486540518713 * tmp_0 + tmp_2;
   real_t tmp_4  = -p_affine_0_0;
   real_t tmp_5  = p_affine_1_0 + tmp_4;
   real_t tmp_6  = p_affine_2_1 + tmp_1;
   real_t tmp_7  = 1.0 / ( tmp_5 * tmp_6 - ( p_affine_1_1 + tmp_1 ) * ( p_affine_2_0 + tmp_4 ) );
   real_t tmp_8  = tmp_5 * tmp_7;
   real_t tmp_9  = tmp_3 * tmp_8;
   real_t tmp_10 = tmp_7 * ( p_affine_0_0 - p_affine_2_0 );
   real_t tmp_11 = tmp_10 * tmp_3;
   real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
   real_t tmp_13 = p_affine_6_0 + tmp_4;
   real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
   real_t tmp_15 = tmp_6 * tmp_7;
   real_t tmp_16 = tmp_14 * tmp_15;
   real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
   real_t tmp_18 = tmp_14 * tmp_17;
   real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
   real_t tmp_20 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
   real_t tmp_21 = 4 * sigma_0 * std::pow( tmp_20, -beta_0 );
   real_t tmp_22 = p_affine_10_0 * ( -tmp_15 - tmp_17 );
   real_t tmp_23 = p_affine_10_1 * ( -tmp_10 - tmp_8 );
   real_t tmp_24 = tmp_22 + tmp_23;
   real_t tmp_25 = -tmp_22 - tmp_23;
   real_t tmp_26 = 0.5 * tmp_20;
   real_t tmp_27 = 0.78867513459481287 * tmp_0 + tmp_2;
   real_t tmp_28 = tmp_27 * tmp_8;
   real_t tmp_29 = tmp_10 * tmp_27;
   real_t tmp_30 = 0.78867513459481287 * tmp_12 + tmp_13;
   real_t tmp_31 = tmp_15 * tmp_30;
   real_t tmp_32 = tmp_17 * tmp_30;
   real_t tmp_33 = -tmp_28 - tmp_29 - tmp_31 - tmp_32 + 1;
   real_t tmp_34 = 0.5 * tmp_20;
   real_t tmp_35 = tmp_11 + tmp_16;
   real_t tmp_36 = p_affine_10_0 * tmp_15;
   real_t tmp_37 = p_affine_10_1 * tmp_10;
   real_t tmp_38 = tmp_36 + tmp_37;
   real_t tmp_39 = tmp_19 * tmp_21;
   real_t tmp_40 = tmp_35 * tmp_39;
   real_t tmp_41 = tmp_29 + tmp_31;
   real_t tmp_42 = tmp_21 * tmp_33;
   real_t tmp_43 = tmp_41 * tmp_42;
   real_t tmp_44 = tmp_18 + tmp_9;
   real_t tmp_45 = p_affine_10_0 * tmp_17;
   real_t tmp_46 = p_affine_10_1 * tmp_8;
   real_t tmp_47 = tmp_45 + tmp_46;
   real_t tmp_48 = tmp_39 * tmp_44;
   real_t tmp_49 = tmp_28 + tmp_32;
   real_t tmp_50 = tmp_42 * tmp_49;
   real_t tmp_51 = -tmp_36 - tmp_37;
   real_t tmp_52 = tmp_21 * tmp_35 * tmp_44;
   real_t tmp_53 = tmp_21 * tmp_41 * tmp_49;
   real_t tmp_54 = -tmp_45 - tmp_46;
   real_t a_0_0  = tmp_26 * ( ( tmp_19 * tmp_19 ) * tmp_21 - tmp_19 * tmp_24 + tmp_19 * tmp_25 ) +
                  tmp_34 * ( tmp_21 * ( tmp_33 * tmp_33 ) - tmp_24 * tmp_33 + tmp_25 * tmp_33 );
   real_t a_0_1 =
       tmp_26 * ( -tmp_19 * tmp_38 + tmp_25 * tmp_35 + tmp_40 ) + tmp_34 * ( tmp_25 * tmp_41 - tmp_33 * tmp_38 + tmp_43 );
   real_t a_0_2 =
       tmp_26 * ( -tmp_19 * tmp_47 + tmp_25 * tmp_44 + tmp_48 ) + tmp_34 * ( tmp_25 * tmp_49 - tmp_33 * tmp_47 + tmp_50 );
   real_t a_1_0 =
       tmp_26 * ( tmp_19 * tmp_51 - tmp_24 * tmp_35 + tmp_40 ) + tmp_34 * ( -tmp_24 * tmp_41 + tmp_33 * tmp_51 + tmp_43 );
   real_t a_1_1 = tmp_26 * ( tmp_21 * ( tmp_35 * tmp_35 ) - tmp_35 * tmp_38 + tmp_35 * tmp_51 ) +
                  tmp_34 * ( tmp_21 * ( tmp_41 * tmp_41 ) - tmp_38 * tmp_41 + tmp_41 * tmp_51 );
   real_t a_1_2 =
       tmp_26 * ( -tmp_35 * tmp_47 + tmp_44 * tmp_51 + tmp_52 ) + tmp_34 * ( -tmp_41 * tmp_47 + tmp_49 * tmp_51 + tmp_53 );
   real_t a_2_0 =
       tmp_26 * ( tmp_19 * tmp_54 - tmp_24 * tmp_44 + tmp_48 ) + tmp_34 * ( -tmp_24 * tmp_49 + tmp_33 * tmp_54 + tmp_50 );
   real_t a_2_1 =
       tmp_26 * ( tmp_35 * tmp_54 - tmp_38 * tmp_44 + tmp_52 ) + tmp_34 * ( -tmp_38 * tmp_49 + tmp_41 * tmp_54 + tmp_53 );
   real_t a_2_2 = tmp_26 * ( tmp_21 * ( tmp_44 * tmp_44 ) - tmp_44 * tmp_47 + tmp_44 * tmp_54 ) +
                  tmp_34 * ( tmp_21 * ( tmp_49 * tmp_49 ) - tmp_47 * tmp_49 + tmp_49 * tmp_54 );

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

void DGDiffusionForm_Example::integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                                                 const std::vector< Point3D >& coordsFacet,
                                                                 const Point3D&,
                                                                 const Point3D&     outwardNormal,
                                                                 const DGBasisInfo& trialBasis,
                                                                 const DGBasisInfo& testBasis,
                                                                 int                trialDegree,
                                                                 int                testDegree,
                                                                 MatrixXr&          elMat ) const
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

   real_t tmp_0  = -p_affine_8_2;
   real_t tmp_1  = p_affine_9_2 + tmp_0;
   real_t tmp_2  = p_affine_10_2 + tmp_0;
   real_t tmp_3  = -p_affine_0_2;
   real_t tmp_4  = p_affine_8_2 + tmp_3;
   real_t tmp_5  = 0.091576213509770743 * tmp_1 + 0.81684757298045851 * tmp_2 + tmp_4;
   real_t tmp_6  = -p_affine_0_0;
   real_t tmp_7  = p_affine_1_0 + tmp_6;
   real_t tmp_8  = -p_affine_0_1;
   real_t tmp_9  = p_affine_2_1 + tmp_8;
   real_t tmp_10 = p_affine_2_0 + tmp_6;
   real_t tmp_11 = p_affine_1_1 + tmp_8;
   real_t tmp_12 = p_affine_3_2 + tmp_3;
   real_t tmp_13 = tmp_12 * tmp_9;
   real_t tmp_14 = p_affine_3_1 + tmp_8;
   real_t tmp_15 = p_affine_1_2 + tmp_3;
   real_t tmp_16 = tmp_14 * tmp_15;
   real_t tmp_17 = p_affine_3_0 + tmp_6;
   real_t tmp_18 = p_affine_2_2 + tmp_3;
   real_t tmp_19 = tmp_11 * tmp_18;
   real_t tmp_20 = tmp_14 * tmp_18;
   real_t tmp_21 = tmp_11 * tmp_12;
   real_t tmp_22 = tmp_15 * tmp_9;
   real_t tmp_23 =
       1.0 / ( tmp_10 * tmp_16 - tmp_10 * tmp_21 + tmp_13 * tmp_7 + tmp_17 * tmp_19 - tmp_17 * tmp_22 - tmp_20 * tmp_7 );
   real_t tmp_24 = tmp_23 * ( -tmp_10 * tmp_11 + tmp_7 * tmp_9 );
   real_t tmp_25 = tmp_24 * tmp_5;
   real_t tmp_26 = tmp_23 * ( tmp_11 * tmp_17 - tmp_14 * tmp_7 );
   real_t tmp_27 = tmp_26 * tmp_5;
   real_t tmp_28 = -p_affine_8_1;
   real_t tmp_29 = p_affine_9_1 + tmp_28;
   real_t tmp_30 = p_affine_10_1 + tmp_28;
   real_t tmp_31 = p_affine_8_1 + tmp_8;
   real_t tmp_32 = 0.091576213509770743 * tmp_29 + 0.81684757298045851 * tmp_30 + tmp_31;
   real_t tmp_33 = tmp_23 * ( tmp_10 * tmp_15 - tmp_18 * tmp_7 );
   real_t tmp_34 = tmp_32 * tmp_33;
   real_t tmp_35 = tmp_23 * ( tmp_12 * tmp_7 - tmp_15 * tmp_17 );
   real_t tmp_36 = tmp_32 * tmp_35;
   real_t tmp_37 = tmp_23 * ( tmp_10 * tmp_14 - tmp_17 * tmp_9 );
   real_t tmp_38 = tmp_37 * tmp_5;
   real_t tmp_39 = tmp_23 * ( -tmp_10 * tmp_12 + tmp_17 * tmp_18 );
   real_t tmp_40 = tmp_32 * tmp_39;
   real_t tmp_41 = -p_affine_8_0;
   real_t tmp_42 = p_affine_9_0 + tmp_41;
   real_t tmp_43 = p_affine_10_0 + tmp_41;
   real_t tmp_44 = p_affine_8_0 + tmp_6;
   real_t tmp_45 = 0.091576213509770743 * tmp_42 + 0.81684757298045851 * tmp_43 + tmp_44;
   real_t tmp_46 = tmp_23 * ( tmp_19 - tmp_22 );
   real_t tmp_47 = tmp_45 * tmp_46;
   real_t tmp_48 = tmp_23 * ( tmp_16 - tmp_21 );
   real_t tmp_49 = tmp_45 * tmp_48;
   real_t tmp_50 = tmp_23 * ( tmp_13 - tmp_20 );
   real_t tmp_51 = tmp_45 * tmp_50;
   real_t tmp_52 = -tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1;
   real_t tmp_53 = p_affine_8_1 - p_affine_9_1;
   real_t tmp_54 = p_affine_8_0 - p_affine_9_0;
   real_t tmp_55 = p_affine_8_2 - p_affine_9_2;
   real_t tmp_56 =
       std::pow( ( std::abs( tmp_2 * tmp_53 - tmp_30 * tmp_55 ) * std::abs( tmp_2 * tmp_53 - tmp_30 * tmp_55 ) ) +
                     ( std::abs( tmp_2 * tmp_54 - tmp_43 * tmp_55 ) * std::abs( tmp_2 * tmp_54 - tmp_43 * tmp_55 ) ) +
                     ( std::abs( tmp_30 * tmp_54 - tmp_43 * tmp_53 ) * std::abs( tmp_30 * tmp_54 - tmp_43 * tmp_53 ) ),
                 1.0 / 2.0 );
   real_t tmp_57  = 4 * sigma_0 * std::pow( 0.5 * tmp_56, -beta_0 );
   real_t tmp_58  = p_affine_13_0 * ( -tmp_46 - tmp_48 - tmp_50 );
   real_t tmp_59  = p_affine_13_1 * ( -tmp_33 - tmp_35 - tmp_39 );
   real_t tmp_60  = p_affine_13_2 * ( -tmp_24 - tmp_26 - tmp_37 );
   real_t tmp_61  = tmp_58 + tmp_59 + tmp_60;
   real_t tmp_62  = -tmp_58 - tmp_59 - tmp_60;
   real_t tmp_63  = 1.0 * tmp_56;
   real_t tmp_64  = 0.054975871827660928 * tmp_63;
   real_t tmp_65  = 0.44594849091596489 * tmp_1 + 0.10810301816807022 * tmp_2 + tmp_4;
   real_t tmp_66  = tmp_24 * tmp_65;
   real_t tmp_67  = tmp_26 * tmp_65;
   real_t tmp_68  = 0.44594849091596489 * tmp_29 + 0.10810301816807022 * tmp_30 + tmp_31;
   real_t tmp_69  = tmp_33 * tmp_68;
   real_t tmp_70  = tmp_35 * tmp_68;
   real_t tmp_71  = tmp_37 * tmp_65;
   real_t tmp_72  = tmp_39 * tmp_68;
   real_t tmp_73  = 0.44594849091596489 * tmp_42 + 0.10810301816807022 * tmp_43 + tmp_44;
   real_t tmp_74  = tmp_46 * tmp_73;
   real_t tmp_75  = tmp_48 * tmp_73;
   real_t tmp_76  = tmp_50 * tmp_73;
   real_t tmp_77  = -tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 + 1;
   real_t tmp_78  = 0.11169079483900572 * tmp_63;
   real_t tmp_79  = 0.81684757298045851 * tmp_1 + 0.091576213509770743 * tmp_2 + tmp_4;
   real_t tmp_80  = tmp_24 * tmp_79;
   real_t tmp_81  = tmp_26 * tmp_79;
   real_t tmp_82  = 0.81684757298045851 * tmp_29 + 0.091576213509770743 * tmp_30 + tmp_31;
   real_t tmp_83  = tmp_33 * tmp_82;
   real_t tmp_84  = tmp_35 * tmp_82;
   real_t tmp_85  = tmp_37 * tmp_79;
   real_t tmp_86  = tmp_39 * tmp_82;
   real_t tmp_87  = 0.81684757298045851 * tmp_42 + 0.091576213509770743 * tmp_43 + tmp_44;
   real_t tmp_88  = tmp_46 * tmp_87;
   real_t tmp_89  = tmp_48 * tmp_87;
   real_t tmp_90  = tmp_50 * tmp_87;
   real_t tmp_91  = -tmp_80 - tmp_81 - tmp_83 - tmp_84 - tmp_85 - tmp_86 - tmp_88 - tmp_89 - tmp_90 + 1;
   real_t tmp_92  = 0.054975871827660928 * tmp_63;
   real_t tmp_93  = 0.10810301816807022 * tmp_1 + 0.44594849091596489 * tmp_2 + tmp_4;
   real_t tmp_94  = tmp_24 * tmp_93;
   real_t tmp_95  = tmp_26 * tmp_93;
   real_t tmp_96  = 0.10810301816807022 * tmp_29 + 0.44594849091596489 * tmp_30 + tmp_31;
   real_t tmp_97  = tmp_33 * tmp_96;
   real_t tmp_98  = tmp_35 * tmp_96;
   real_t tmp_99  = tmp_37 * tmp_93;
   real_t tmp_100 = tmp_39 * tmp_96;
   real_t tmp_101 = 0.10810301816807022 * tmp_42 + 0.44594849091596489 * tmp_43 + tmp_44;
   real_t tmp_102 = tmp_101 * tmp_46;
   real_t tmp_103 = tmp_101 * tmp_48;
   real_t tmp_104 = tmp_101 * tmp_50;
   real_t tmp_105 = -tmp_100 - tmp_102 - tmp_103 - tmp_104 - tmp_94 - tmp_95 - tmp_97 - tmp_98 - tmp_99 + 1;
   real_t tmp_106 = 0.11169079483900572 * tmp_63;
   real_t tmp_107 = 0.091576213509770743 * tmp_1 + 0.091576213509770743 * tmp_2 + tmp_4;
   real_t tmp_108 = tmp_107 * tmp_24;
   real_t tmp_109 = tmp_107 * tmp_26;
   real_t tmp_110 = 0.091576213509770743 * tmp_29 + 0.091576213509770743 * tmp_30 + tmp_31;
   real_t tmp_111 = tmp_110 * tmp_33;
   real_t tmp_112 = tmp_110 * tmp_35;
   real_t tmp_113 = tmp_107 * tmp_37;
   real_t tmp_114 = tmp_110 * tmp_39;
   real_t tmp_115 = 0.091576213509770743 * tmp_42 + 0.091576213509770743 * tmp_43 + tmp_44;
   real_t tmp_116 = tmp_115 * tmp_46;
   real_t tmp_117 = tmp_115 * tmp_48;
   real_t tmp_118 = tmp_115 * tmp_50;
   real_t tmp_119 = -tmp_108 - tmp_109 - tmp_111 - tmp_112 - tmp_113 - tmp_114 - tmp_116 - tmp_117 - tmp_118 + 1;
   real_t tmp_120 = 0.054975871827660928 * tmp_63;
   real_t tmp_121 = 0.44594849091596489 * tmp_1 + 0.44594849091596489 * tmp_2 + tmp_4;
   real_t tmp_122 = tmp_121 * tmp_24;
   real_t tmp_123 = tmp_121 * tmp_26;
   real_t tmp_124 = 0.44594849091596489 * tmp_29 + 0.44594849091596489 * tmp_30 + tmp_31;
   real_t tmp_125 = tmp_124 * tmp_33;
   real_t tmp_126 = tmp_124 * tmp_35;
   real_t tmp_127 = tmp_121 * tmp_37;
   real_t tmp_128 = tmp_124 * tmp_39;
   real_t tmp_129 = 0.44594849091596489 * tmp_42 + 0.44594849091596489 * tmp_43 + tmp_44;
   real_t tmp_130 = tmp_129 * tmp_46;
   real_t tmp_131 = tmp_129 * tmp_48;
   real_t tmp_132 = tmp_129 * tmp_50;
   real_t tmp_133 = -tmp_122 - tmp_123 - tmp_125 - tmp_126 - tmp_127 - tmp_128 - tmp_130 - tmp_131 - tmp_132 + 1;
   real_t tmp_134 = 0.11169079483900572 * tmp_63;
   real_t tmp_135 = tmp_38 + tmp_40 + tmp_51;
   real_t tmp_136 = p_affine_13_0 * tmp_50;
   real_t tmp_137 = p_affine_13_1 * tmp_39;
   real_t tmp_138 = p_affine_13_2 * tmp_37;
   real_t tmp_139 = tmp_136 + tmp_137 + tmp_138;
   real_t tmp_140 = tmp_52 * tmp_57;
   real_t tmp_141 = tmp_135 * tmp_140;
   real_t tmp_142 = tmp_71 + tmp_72 + tmp_76;
   real_t tmp_143 = tmp_57 * tmp_77;
   real_t tmp_144 = tmp_142 * tmp_143;
   real_t tmp_145 = tmp_85 + tmp_86 + tmp_90;
   real_t tmp_146 = tmp_57 * tmp_91;
   real_t tmp_147 = tmp_145 * tmp_146;
   real_t tmp_148 = tmp_100 + tmp_104 + tmp_99;
   real_t tmp_149 = tmp_105 * tmp_57;
   real_t tmp_150 = tmp_148 * tmp_149;
   real_t tmp_151 = tmp_113 + tmp_114 + tmp_118;
   real_t tmp_152 = tmp_119 * tmp_57;
   real_t tmp_153 = tmp_151 * tmp_152;
   real_t tmp_154 = tmp_127 + tmp_128 + tmp_132;
   real_t tmp_155 = tmp_133 * tmp_57;
   real_t tmp_156 = tmp_154 * tmp_155;
   real_t tmp_157 = tmp_27 + tmp_36 + tmp_49;
   real_t tmp_158 = p_affine_13_0 * tmp_48;
   real_t tmp_159 = p_affine_13_1 * tmp_35;
   real_t tmp_160 = p_affine_13_2 * tmp_26;
   real_t tmp_161 = tmp_158 + tmp_159 + tmp_160;
   real_t tmp_162 = tmp_140 * tmp_157;
   real_t tmp_163 = tmp_67 + tmp_70 + tmp_75;
   real_t tmp_164 = tmp_143 * tmp_163;
   real_t tmp_165 = tmp_81 + tmp_84 + tmp_89;
   real_t tmp_166 = tmp_146 * tmp_165;
   real_t tmp_167 = tmp_103 + tmp_95 + tmp_98;
   real_t tmp_168 = tmp_149 * tmp_167;
   real_t tmp_169 = tmp_109 + tmp_112 + tmp_117;
   real_t tmp_170 = tmp_152 * tmp_169;
   real_t tmp_171 = tmp_123 + tmp_126 + tmp_131;
   real_t tmp_172 = tmp_155 * tmp_171;
   real_t tmp_173 = tmp_25 + tmp_34 + tmp_47;
   real_t tmp_174 = p_affine_13_0 * tmp_46;
   real_t tmp_175 = p_affine_13_1 * tmp_33;
   real_t tmp_176 = p_affine_13_2 * tmp_24;
   real_t tmp_177 = tmp_174 + tmp_175 + tmp_176;
   real_t tmp_178 = tmp_140 * tmp_173;
   real_t tmp_179 = tmp_66 + tmp_69 + tmp_74;
   real_t tmp_180 = tmp_143 * tmp_179;
   real_t tmp_181 = tmp_80 + tmp_83 + tmp_88;
   real_t tmp_182 = tmp_146 * tmp_181;
   real_t tmp_183 = tmp_102 + tmp_94 + tmp_97;
   real_t tmp_184 = tmp_149 * tmp_183;
   real_t tmp_185 = tmp_108 + tmp_111 + tmp_116;
   real_t tmp_186 = tmp_152 * tmp_185;
   real_t tmp_187 = tmp_122 + tmp_125 + tmp_130;
   real_t tmp_188 = tmp_155 * tmp_187;
   real_t tmp_189 = -tmp_136 - tmp_137 - tmp_138;
   real_t tmp_190 = tmp_135 * tmp_57;
   real_t tmp_191 = tmp_157 * tmp_190;
   real_t tmp_192 = tmp_142 * tmp_57;
   real_t tmp_193 = tmp_163 * tmp_192;
   real_t tmp_194 = tmp_145 * tmp_57;
   real_t tmp_195 = tmp_165 * tmp_194;
   real_t tmp_196 = tmp_148 * tmp_57;
   real_t tmp_197 = tmp_167 * tmp_196;
   real_t tmp_198 = tmp_151 * tmp_57;
   real_t tmp_199 = tmp_169 * tmp_198;
   real_t tmp_200 = tmp_154 * tmp_57;
   real_t tmp_201 = tmp_171 * tmp_200;
   real_t tmp_202 = tmp_173 * tmp_190;
   real_t tmp_203 = tmp_179 * tmp_192;
   real_t tmp_204 = tmp_181 * tmp_194;
   real_t tmp_205 = tmp_183 * tmp_196;
   real_t tmp_206 = tmp_185 * tmp_198;
   real_t tmp_207 = tmp_187 * tmp_200;
   real_t tmp_208 = -tmp_158 - tmp_159 - tmp_160;
   real_t tmp_209 = tmp_157 * tmp_173 * tmp_57;
   real_t tmp_210 = tmp_163 * tmp_179 * tmp_57;
   real_t tmp_211 = tmp_165 * tmp_181 * tmp_57;
   real_t tmp_212 = tmp_167 * tmp_183 * tmp_57;
   real_t tmp_213 = tmp_169 * tmp_185 * tmp_57;
   real_t tmp_214 = tmp_171 * tmp_187 * tmp_57;
   real_t tmp_215 = -tmp_174 - tmp_175 - tmp_176;
   real_t a_0_0   = tmp_106 * ( ( tmp_105 * tmp_105 ) * tmp_57 - tmp_105 * tmp_61 + tmp_105 * tmp_62 ) +
                  tmp_120 * ( ( tmp_119 * tmp_119 ) * tmp_57 - tmp_119 * tmp_61 + tmp_119 * tmp_62 ) +
                  tmp_134 * ( ( tmp_133 * tmp_133 ) * tmp_57 - tmp_133 * tmp_61 + tmp_133 * tmp_62 ) +
                  tmp_64 * ( ( tmp_52 * tmp_52 ) * tmp_57 - tmp_52 * tmp_61 + tmp_52 * tmp_62 ) +
                  tmp_78 * ( tmp_57 * ( tmp_77 * tmp_77 ) - tmp_61 * tmp_77 + tmp_62 * tmp_77 ) +
                  tmp_92 * ( tmp_57 * ( tmp_91 * tmp_91 ) - tmp_61 * tmp_91 + tmp_62 * tmp_91 );
   real_t a_0_1 = tmp_106 * ( -tmp_105 * tmp_139 + tmp_148 * tmp_62 + tmp_150 ) +
                  tmp_120 * ( -tmp_119 * tmp_139 + tmp_151 * tmp_62 + tmp_153 ) +
                  tmp_134 * ( -tmp_133 * tmp_139 + tmp_154 * tmp_62 + tmp_156 ) +
                  tmp_64 * ( tmp_135 * tmp_62 - tmp_139 * tmp_52 + tmp_141 ) +
                  tmp_78 * ( -tmp_139 * tmp_77 + tmp_142 * tmp_62 + tmp_144 ) +
                  tmp_92 * ( -tmp_139 * tmp_91 + tmp_145 * tmp_62 + tmp_147 );
   real_t a_0_2 = tmp_106 * ( -tmp_105 * tmp_161 + tmp_167 * tmp_62 + tmp_168 ) +
                  tmp_120 * ( -tmp_119 * tmp_161 + tmp_169 * tmp_62 + tmp_170 ) +
                  tmp_134 * ( -tmp_133 * tmp_161 + tmp_171 * tmp_62 + tmp_172 ) +
                  tmp_64 * ( tmp_157 * tmp_62 - tmp_161 * tmp_52 + tmp_162 ) +
                  tmp_78 * ( -tmp_161 * tmp_77 + tmp_163 * tmp_62 + tmp_164 ) +
                  tmp_92 * ( -tmp_161 * tmp_91 + tmp_165 * tmp_62 + tmp_166 );
   real_t a_0_3 = tmp_106 * ( -tmp_105 * tmp_177 + tmp_183 * tmp_62 + tmp_184 ) +
                  tmp_120 * ( -tmp_119 * tmp_177 + tmp_185 * tmp_62 + tmp_186 ) +
                  tmp_134 * ( -tmp_133 * tmp_177 + tmp_187 * tmp_62 + tmp_188 ) +
                  tmp_64 * ( tmp_173 * tmp_62 - tmp_177 * tmp_52 + tmp_178 ) +
                  tmp_78 * ( -tmp_177 * tmp_77 + tmp_179 * tmp_62 + tmp_180 ) +
                  tmp_92 * ( -tmp_177 * tmp_91 + tmp_181 * tmp_62 + tmp_182 );
   real_t a_1_0 = tmp_106 * ( tmp_105 * tmp_189 - tmp_148 * tmp_61 + tmp_150 ) +
                  tmp_120 * ( tmp_119 * tmp_189 - tmp_151 * tmp_61 + tmp_153 ) +
                  tmp_134 * ( tmp_133 * tmp_189 - tmp_154 * tmp_61 + tmp_156 ) +
                  tmp_64 * ( -tmp_135 * tmp_61 + tmp_141 + tmp_189 * tmp_52 ) +
                  tmp_78 * ( -tmp_142 * tmp_61 + tmp_144 + tmp_189 * tmp_77 ) +
                  tmp_92 * ( -tmp_145 * tmp_61 + tmp_147 + tmp_189 * tmp_91 );
   real_t a_1_1 = tmp_106 * ( -tmp_139 * tmp_148 + ( tmp_148 * tmp_148 ) * tmp_57 + tmp_148 * tmp_189 ) +
                  tmp_120 * ( -tmp_139 * tmp_151 + ( tmp_151 * tmp_151 ) * tmp_57 + tmp_151 * tmp_189 ) +
                  tmp_134 * ( -tmp_139 * tmp_154 + ( tmp_154 * tmp_154 ) * tmp_57 + tmp_154 * tmp_189 ) +
                  tmp_64 * ( ( tmp_135 * tmp_135 ) * tmp_57 - tmp_135 * tmp_139 + tmp_135 * tmp_189 ) +
                  tmp_78 * ( -tmp_139 * tmp_142 + ( tmp_142 * tmp_142 ) * tmp_57 + tmp_142 * tmp_189 ) +
                  tmp_92 * ( -tmp_139 * tmp_145 + ( tmp_145 * tmp_145 ) * tmp_57 + tmp_145 * tmp_189 );
   real_t a_1_2 = tmp_106 * ( -tmp_148 * tmp_161 + tmp_167 * tmp_189 + tmp_197 ) +
                  tmp_120 * ( -tmp_151 * tmp_161 + tmp_169 * tmp_189 + tmp_199 ) +
                  tmp_134 * ( -tmp_154 * tmp_161 + tmp_171 * tmp_189 + tmp_201 ) +
                  tmp_64 * ( -tmp_135 * tmp_161 + tmp_157 * tmp_189 + tmp_191 ) +
                  tmp_78 * ( -tmp_142 * tmp_161 + tmp_163 * tmp_189 + tmp_193 ) +
                  tmp_92 * ( -tmp_145 * tmp_161 + tmp_165 * tmp_189 + tmp_195 );
   real_t a_1_3 = tmp_106 * ( -tmp_148 * tmp_177 + tmp_183 * tmp_189 + tmp_205 ) +
                  tmp_120 * ( -tmp_151 * tmp_177 + tmp_185 * tmp_189 + tmp_206 ) +
                  tmp_134 * ( -tmp_154 * tmp_177 + tmp_187 * tmp_189 + tmp_207 ) +
                  tmp_64 * ( -tmp_135 * tmp_177 + tmp_173 * tmp_189 + tmp_202 ) +
                  tmp_78 * ( -tmp_142 * tmp_177 + tmp_179 * tmp_189 + tmp_203 ) +
                  tmp_92 * ( -tmp_145 * tmp_177 + tmp_181 * tmp_189 + tmp_204 );
   real_t a_2_0 = tmp_106 * ( tmp_105 * tmp_208 - tmp_167 * tmp_61 + tmp_168 ) +
                  tmp_120 * ( tmp_119 * tmp_208 - tmp_169 * tmp_61 + tmp_170 ) +
                  tmp_134 * ( tmp_133 * tmp_208 - tmp_171 * tmp_61 + tmp_172 ) +
                  tmp_64 * ( -tmp_157 * tmp_61 + tmp_162 + tmp_208 * tmp_52 ) +
                  tmp_78 * ( -tmp_163 * tmp_61 + tmp_164 + tmp_208 * tmp_77 ) +
                  tmp_92 * ( -tmp_165 * tmp_61 + tmp_166 + tmp_208 * tmp_91 );
   real_t a_2_1 = tmp_106 * ( -tmp_139 * tmp_167 + tmp_148 * tmp_208 + tmp_197 ) +
                  tmp_120 * ( -tmp_139 * tmp_169 + tmp_151 * tmp_208 + tmp_199 ) +
                  tmp_134 * ( -tmp_139 * tmp_171 + tmp_154 * tmp_208 + tmp_201 ) +
                  tmp_64 * ( tmp_135 * tmp_208 - tmp_139 * tmp_157 + tmp_191 ) +
                  tmp_78 * ( -tmp_139 * tmp_163 + tmp_142 * tmp_208 + tmp_193 ) +
                  tmp_92 * ( -tmp_139 * tmp_165 + tmp_145 * tmp_208 + tmp_195 );
   real_t a_2_2 = tmp_106 * ( -tmp_161 * tmp_167 + ( tmp_167 * tmp_167 ) * tmp_57 + tmp_167 * tmp_208 ) +
                  tmp_120 * ( -tmp_161 * tmp_169 + ( tmp_169 * tmp_169 ) * tmp_57 + tmp_169 * tmp_208 ) +
                  tmp_134 * ( -tmp_161 * tmp_171 + ( tmp_171 * tmp_171 ) * tmp_57 + tmp_171 * tmp_208 ) +
                  tmp_64 * ( ( tmp_157 * tmp_157 ) * tmp_57 - tmp_157 * tmp_161 + tmp_157 * tmp_208 ) +
                  tmp_78 * ( -tmp_161 * tmp_163 + ( tmp_163 * tmp_163 ) * tmp_57 + tmp_163 * tmp_208 ) +
                  tmp_92 * ( -tmp_161 * tmp_165 + ( tmp_165 * tmp_165 ) * tmp_57 + tmp_165 * tmp_208 );
   real_t a_2_3 = tmp_106 * ( -tmp_167 * tmp_177 + tmp_183 * tmp_208 + tmp_212 ) +
                  tmp_120 * ( -tmp_169 * tmp_177 + tmp_185 * tmp_208 + tmp_213 ) +
                  tmp_134 * ( -tmp_171 * tmp_177 + tmp_187 * tmp_208 + tmp_214 ) +
                  tmp_64 * ( -tmp_157 * tmp_177 + tmp_173 * tmp_208 + tmp_209 ) +
                  tmp_78 * ( -tmp_163 * tmp_177 + tmp_179 * tmp_208 + tmp_210 ) +
                  tmp_92 * ( -tmp_165 * tmp_177 + tmp_181 * tmp_208 + tmp_211 );
   real_t a_3_0 = tmp_106 * ( tmp_105 * tmp_215 - tmp_183 * tmp_61 + tmp_184 ) +
                  tmp_120 * ( tmp_119 * tmp_215 - tmp_185 * tmp_61 + tmp_186 ) +
                  tmp_134 * ( tmp_133 * tmp_215 - tmp_187 * tmp_61 + tmp_188 ) +
                  tmp_64 * ( -tmp_173 * tmp_61 + tmp_178 + tmp_215 * tmp_52 ) +
                  tmp_78 * ( -tmp_179 * tmp_61 + tmp_180 + tmp_215 * tmp_77 ) +
                  tmp_92 * ( -tmp_181 * tmp_61 + tmp_182 + tmp_215 * tmp_91 );
   real_t a_3_1 = tmp_106 * ( -tmp_139 * tmp_183 + tmp_148 * tmp_215 + tmp_205 ) +
                  tmp_120 * ( -tmp_139 * tmp_185 + tmp_151 * tmp_215 + tmp_206 ) +
                  tmp_134 * ( -tmp_139 * tmp_187 + tmp_154 * tmp_215 + tmp_207 ) +
                  tmp_64 * ( tmp_135 * tmp_215 - tmp_139 * tmp_173 + tmp_202 ) +
                  tmp_78 * ( -tmp_139 * tmp_179 + tmp_142 * tmp_215 + tmp_203 ) +
                  tmp_92 * ( -tmp_139 * tmp_181 + tmp_145 * tmp_215 + tmp_204 );
   real_t a_3_2 = tmp_106 * ( -tmp_161 * tmp_183 + tmp_167 * tmp_215 + tmp_212 ) +
                  tmp_120 * ( -tmp_161 * tmp_185 + tmp_169 * tmp_215 + tmp_213 ) +
                  tmp_134 * ( -tmp_161 * tmp_187 + tmp_171 * tmp_215 + tmp_214 ) +
                  tmp_64 * ( tmp_157 * tmp_215 - tmp_161 * tmp_173 + tmp_209 ) +
                  tmp_78 * ( -tmp_161 * tmp_179 + tmp_163 * tmp_215 + tmp_210 ) +
                  tmp_92 * ( -tmp_161 * tmp_181 + tmp_165 * tmp_215 + tmp_211 );
   real_t a_3_3 = tmp_106 * ( -tmp_177 * tmp_183 + ( tmp_183 * tmp_183 ) * tmp_57 + tmp_183 * tmp_215 ) +
                  tmp_120 * ( -tmp_177 * tmp_185 + ( tmp_185 * tmp_185 ) * tmp_57 + tmp_185 * tmp_215 ) +
                  tmp_134 * ( -tmp_177 * tmp_187 + ( tmp_187 * tmp_187 ) * tmp_57 + tmp_187 * tmp_215 ) +
                  tmp_64 * ( ( tmp_173 * tmp_173 ) * tmp_57 - tmp_173 * tmp_177 + tmp_173 * tmp_215 ) +
                  tmp_78 * ( -tmp_177 * tmp_179 + ( tmp_179 * tmp_179 ) * tmp_57 + tmp_179 * tmp_215 ) +
                  tmp_92 * ( -tmp_177 * tmp_181 + ( tmp_181 * tmp_181 ) * tmp_57 + tmp_181 * tmp_215 );

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

void DGDiffusionForm_Example::integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                               const std::vector< Point3D >& coordsFacet,
                                                               const Point3D&,
                                                               const Point3D&     outwardNormal,
                                                               const DGBasisInfo& basis,
                                                               int                degree,
                                                               MatrixXr&          elMat ) const
{
   elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, uint_c( degree ) ) ), 1 );

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

   real_t Scalar_Variable_Coefficient_2D_g_out0_id0 = 0;
   real_t Scalar_Variable_Coefficient_2D_g_out0_id1 = 0;
   Scalar_Variable_Coefficient_2D_g( 0.78867513459481287 * p_affine_6_0 + 0.21132486540518713 * p_affine_7_0,
                                     0.78867513459481287 * p_affine_6_1 + 0.21132486540518713 * p_affine_7_1,
                                     &Scalar_Variable_Coefficient_2D_g_out0_id0 );
   Scalar_Variable_Coefficient_2D_g( 0.21132486540518713 * p_affine_6_0 + 0.78867513459481287 * p_affine_7_0,
                                     0.21132486540518713 * p_affine_6_1 + 0.78867513459481287 * p_affine_7_1,
                                     &Scalar_Variable_Coefficient_2D_g_out0_id1 );
   real_t tmp_0  = -p_affine_6_1 + p_affine_7_1;
   real_t tmp_1  = -p_affine_0_1;
   real_t tmp_2  = p_affine_6_1 + tmp_1;
   real_t tmp_3  = 0.21132486540518713 * tmp_0 + tmp_2;
   real_t tmp_4  = -p_affine_0_0;
   real_t tmp_5  = p_affine_1_0 + tmp_4;
   real_t tmp_6  = p_affine_2_1 + tmp_1;
   real_t tmp_7  = 1.0 / ( tmp_5 * tmp_6 - ( p_affine_1_1 + tmp_1 ) * ( p_affine_2_0 + tmp_4 ) );
   real_t tmp_8  = tmp_5 * tmp_7;
   real_t tmp_9  = tmp_3 * tmp_8;
   real_t tmp_10 = tmp_7 * ( p_affine_0_0 - p_affine_2_0 );
   real_t tmp_11 = tmp_10 * tmp_3;
   real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
   real_t tmp_13 = p_affine_6_0 + tmp_4;
   real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
   real_t tmp_15 = tmp_6 * tmp_7;
   real_t tmp_16 = tmp_14 * tmp_15;
   real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
   real_t tmp_18 = tmp_14 * tmp_17;
   real_t tmp_19 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
   real_t tmp_20 = 4 * sigma_0 * std::pow( tmp_19, -beta_0 );
   real_t tmp_21 = -p_affine_10_0 * ( -tmp_15 - tmp_17 ) - p_affine_10_1 * ( -tmp_10 - tmp_8 );
   real_t tmp_22 = 0.5 * Scalar_Variable_Coefficient_2D_g_out0_id0 * tmp_19;
   real_t tmp_23 = 0.78867513459481287 * tmp_0 + tmp_2;
   real_t tmp_24 = tmp_23 * tmp_8;
   real_t tmp_25 = tmp_10 * tmp_23;
   real_t tmp_26 = 0.78867513459481287 * tmp_12 + tmp_13;
   real_t tmp_27 = tmp_15 * tmp_26;
   real_t tmp_28 = tmp_17 * tmp_26;
   real_t tmp_29 = 0.5 * Scalar_Variable_Coefficient_2D_g_out0_id1 * tmp_19;
   real_t tmp_30 = -p_affine_10_0 * tmp_15 - p_affine_10_1 * tmp_10;
   real_t tmp_31 = -p_affine_10_0 * tmp_17 - p_affine_10_1 * tmp_8;
   real_t a_0_0  = tmp_22 * ( tmp_20 * ( -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1 ) + tmp_21 ) +
                  tmp_29 * ( tmp_20 * ( -tmp_24 - tmp_25 - tmp_27 - tmp_28 + 1 ) + tmp_21 );
   real_t a_1_0 = tmp_22 * ( tmp_20 * ( tmp_11 + tmp_16 ) + tmp_30 ) + tmp_29 * ( tmp_20 * ( tmp_25 + tmp_27 ) + tmp_30 );
   real_t a_2_0 = tmp_22 * ( tmp_20 * ( tmp_18 + tmp_9 ) + tmp_31 ) + tmp_29 * ( tmp_20 * ( tmp_24 + tmp_28 ) + tmp_31 );

   elMat( 0, 0 ) = a_0_0;
   elMat( 1, 0 ) = a_1_0;
   elMat( 2, 0 ) = a_2_0;
}

void DGDiffusionForm_Example::integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                                               const std::vector< Point3D >& coordsFacet,
                                                               const Point3D&,
                                                               const Point3D&     outwardNormal,
                                                               const DGBasisInfo& basis,
                                                               int                degree,
                                                               MatrixXr&          elMat ) const
{
   elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, uint_c( degree ) ) ), 1 );

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

   real_t Scalar_Variable_Coefficient_3D_g_out0_id0 = 0;
   real_t Scalar_Variable_Coefficient_3D_g_out0_id1 = 0;
   real_t Scalar_Variable_Coefficient_3D_g_out0_id2 = 0;
   real_t Scalar_Variable_Coefficient_3D_g_out0_id3 = 0;
   real_t Scalar_Variable_Coefficient_3D_g_out0_id4 = 0;
   real_t Scalar_Variable_Coefficient_3D_g_out0_id5 = 0;
   Scalar_Variable_Coefficient_3D_g(
       0.81684757298045851 * p_affine_10_0 + 0.091576213509770743 * p_affine_8_0 + 0.091576213509770743 * p_affine_9_0,
       0.81684757298045851 * p_affine_10_1 + 0.091576213509770743 * p_affine_8_1 + 0.091576213509770743 * p_affine_9_1,
       0.81684757298045851 * p_affine_10_2 + 0.091576213509770743 * p_affine_8_2 + 0.091576213509770743 * p_affine_9_2,
       &Scalar_Variable_Coefficient_3D_g_out0_id0 );
   Scalar_Variable_Coefficient_3D_g(
       0.10810301816807022 * p_affine_10_0 + 0.44594849091596489 * p_affine_8_0 + 0.44594849091596489 * p_affine_9_0,
       0.10810301816807022 * p_affine_10_1 + 0.44594849091596489 * p_affine_8_1 + 0.44594849091596489 * p_affine_9_1,
       0.10810301816807022 * p_affine_10_2 + 0.44594849091596489 * p_affine_8_2 + 0.44594849091596489 * p_affine_9_2,
       &Scalar_Variable_Coefficient_3D_g_out0_id1 );
   Scalar_Variable_Coefficient_3D_g(
       0.091576213509770743 * p_affine_10_0 + 0.091576213509770743 * p_affine_8_0 + 0.81684757298045851 * p_affine_9_0,
       0.091576213509770743 * p_affine_10_1 + 0.091576213509770743 * p_affine_8_1 + 0.81684757298045851 * p_affine_9_1,
       0.091576213509770743 * p_affine_10_2 + 0.091576213509770743 * p_affine_8_2 + 0.81684757298045851 * p_affine_9_2,
       &Scalar_Variable_Coefficient_3D_g_out0_id2 );
   Scalar_Variable_Coefficient_3D_g(
       0.44594849091596489 * p_affine_10_0 + 0.44594849091596489 * p_affine_8_0 + 0.10810301816807022 * p_affine_9_0,
       0.44594849091596489 * p_affine_10_1 + 0.44594849091596489 * p_affine_8_1 + 0.10810301816807022 * p_affine_9_1,
       0.44594849091596489 * p_affine_10_2 + 0.44594849091596489 * p_affine_8_2 + 0.10810301816807022 * p_affine_9_2,
       &Scalar_Variable_Coefficient_3D_g_out0_id3 );
   Scalar_Variable_Coefficient_3D_g(
       0.091576213509770743 * p_affine_10_0 + 0.81684757298045851 * p_affine_8_0 + 0.091576213509770743 * p_affine_9_0,
       0.091576213509770743 * p_affine_10_1 + 0.81684757298045851 * p_affine_8_1 + 0.091576213509770743 * p_affine_9_1,
       0.091576213509770743 * p_affine_10_2 + 0.81684757298045851 * p_affine_8_2 + 0.091576213509770743 * p_affine_9_2,
       &Scalar_Variable_Coefficient_3D_g_out0_id4 );
   Scalar_Variable_Coefficient_3D_g(
       0.44594849091596489 * p_affine_10_0 + 0.10810301816807022 * p_affine_8_0 + 0.44594849091596489 * p_affine_9_0,
       0.44594849091596489 * p_affine_10_1 + 0.10810301816807022 * p_affine_8_1 + 0.44594849091596489 * p_affine_9_1,
       0.44594849091596489 * p_affine_10_2 + 0.10810301816807022 * p_affine_8_2 + 0.44594849091596489 * p_affine_9_2,
       &Scalar_Variable_Coefficient_3D_g_out0_id5 );
   real_t tmp_0  = -p_affine_8_2;
   real_t tmp_1  = p_affine_9_2 + tmp_0;
   real_t tmp_2  = p_affine_10_2 + tmp_0;
   real_t tmp_3  = -p_affine_0_2;
   real_t tmp_4  = p_affine_8_2 + tmp_3;
   real_t tmp_5  = 0.091576213509770743 * tmp_1 + 0.81684757298045851 * tmp_2 + tmp_4;
   real_t tmp_6  = -p_affine_0_0;
   real_t tmp_7  = p_affine_1_0 + tmp_6;
   real_t tmp_8  = -p_affine_0_1;
   real_t tmp_9  = p_affine_2_1 + tmp_8;
   real_t tmp_10 = p_affine_2_0 + tmp_6;
   real_t tmp_11 = p_affine_1_1 + tmp_8;
   real_t tmp_12 = p_affine_3_2 + tmp_3;
   real_t tmp_13 = tmp_12 * tmp_9;
   real_t tmp_14 = p_affine_3_1 + tmp_8;
   real_t tmp_15 = p_affine_1_2 + tmp_3;
   real_t tmp_16 = tmp_14 * tmp_15;
   real_t tmp_17 = p_affine_3_0 + tmp_6;
   real_t tmp_18 = p_affine_2_2 + tmp_3;
   real_t tmp_19 = tmp_11 * tmp_18;
   real_t tmp_20 = tmp_14 * tmp_18;
   real_t tmp_21 = tmp_11 * tmp_12;
   real_t tmp_22 = tmp_15 * tmp_9;
   real_t tmp_23 =
       1.0 / ( tmp_10 * tmp_16 - tmp_10 * tmp_21 + tmp_13 * tmp_7 + tmp_17 * tmp_19 - tmp_17 * tmp_22 - tmp_20 * tmp_7 );
   real_t tmp_24 = tmp_23 * ( -tmp_10 * tmp_11 + tmp_7 * tmp_9 );
   real_t tmp_25 = tmp_24 * tmp_5;
   real_t tmp_26 = tmp_23 * ( tmp_11 * tmp_17 - tmp_14 * tmp_7 );
   real_t tmp_27 = tmp_26 * tmp_5;
   real_t tmp_28 = -p_affine_8_1;
   real_t tmp_29 = p_affine_9_1 + tmp_28;
   real_t tmp_30 = p_affine_10_1 + tmp_28;
   real_t tmp_31 = p_affine_8_1 + tmp_8;
   real_t tmp_32 = 0.091576213509770743 * tmp_29 + 0.81684757298045851 * tmp_30 + tmp_31;
   real_t tmp_33 = tmp_23 * ( tmp_10 * tmp_15 - tmp_18 * tmp_7 );
   real_t tmp_34 = tmp_32 * tmp_33;
   real_t tmp_35 = tmp_23 * ( tmp_12 * tmp_7 - tmp_15 * tmp_17 );
   real_t tmp_36 = tmp_32 * tmp_35;
   real_t tmp_37 = tmp_23 * ( tmp_10 * tmp_14 - tmp_17 * tmp_9 );
   real_t tmp_38 = tmp_37 * tmp_5;
   real_t tmp_39 = tmp_23 * ( -tmp_10 * tmp_12 + tmp_17 * tmp_18 );
   real_t tmp_40 = tmp_32 * tmp_39;
   real_t tmp_41 = -p_affine_8_0;
   real_t tmp_42 = p_affine_9_0 + tmp_41;
   real_t tmp_43 = p_affine_10_0 + tmp_41;
   real_t tmp_44 = p_affine_8_0 + tmp_6;
   real_t tmp_45 = 0.091576213509770743 * tmp_42 + 0.81684757298045851 * tmp_43 + tmp_44;
   real_t tmp_46 = tmp_23 * ( tmp_19 - tmp_22 );
   real_t tmp_47 = tmp_45 * tmp_46;
   real_t tmp_48 = tmp_23 * ( tmp_16 - tmp_21 );
   real_t tmp_49 = tmp_45 * tmp_48;
   real_t tmp_50 = tmp_23 * ( tmp_13 - tmp_20 );
   real_t tmp_51 = tmp_45 * tmp_50;
   real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
   real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
   real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
   real_t tmp_55 =
       std::pow( ( std::abs( tmp_2 * tmp_52 - tmp_30 * tmp_54 ) * std::abs( tmp_2 * tmp_52 - tmp_30 * tmp_54 ) ) +
                     ( std::abs( tmp_2 * tmp_53 - tmp_43 * tmp_54 ) * std::abs( tmp_2 * tmp_53 - tmp_43 * tmp_54 ) ) +
                     ( std::abs( tmp_30 * tmp_53 - tmp_43 * tmp_52 ) * std::abs( tmp_30 * tmp_53 - tmp_43 * tmp_52 ) ),
                 1.0 / 2.0 );
   real_t tmp_56 = 4 * sigma_0 * std::pow( 0.5 * tmp_55, -beta_0 );
   real_t tmp_57 = -p_affine_13_0 * ( -tmp_46 - tmp_48 - tmp_50 ) - p_affine_13_1 * ( -tmp_33 - tmp_35 - tmp_39 ) -
                   p_affine_13_2 * ( -tmp_24 - tmp_26 - tmp_37 );
   real_t tmp_58  = 1.0 * tmp_55;
   real_t tmp_59  = 0.054975871827660928 * Scalar_Variable_Coefficient_3D_g_out0_id0 * tmp_58;
   real_t tmp_60  = 0.44594849091596489 * tmp_1 + 0.10810301816807022 * tmp_2 + tmp_4;
   real_t tmp_61  = tmp_24 * tmp_60;
   real_t tmp_62  = tmp_26 * tmp_60;
   real_t tmp_63  = 0.44594849091596489 * tmp_29 + 0.10810301816807022 * tmp_30 + tmp_31;
   real_t tmp_64  = tmp_33 * tmp_63;
   real_t tmp_65  = tmp_35 * tmp_63;
   real_t tmp_66  = tmp_37 * tmp_60;
   real_t tmp_67  = tmp_39 * tmp_63;
   real_t tmp_68  = 0.44594849091596489 * tmp_42 + 0.10810301816807022 * tmp_43 + tmp_44;
   real_t tmp_69  = tmp_46 * tmp_68;
   real_t tmp_70  = tmp_48 * tmp_68;
   real_t tmp_71  = tmp_50 * tmp_68;
   real_t tmp_72  = 0.11169079483900572 * Scalar_Variable_Coefficient_3D_g_out0_id1 * tmp_58;
   real_t tmp_73  = 0.81684757298045851 * tmp_1 + 0.091576213509770743 * tmp_2 + tmp_4;
   real_t tmp_74  = tmp_24 * tmp_73;
   real_t tmp_75  = tmp_26 * tmp_73;
   real_t tmp_76  = 0.81684757298045851 * tmp_29 + 0.091576213509770743 * tmp_30 + tmp_31;
   real_t tmp_77  = tmp_33 * tmp_76;
   real_t tmp_78  = tmp_35 * tmp_76;
   real_t tmp_79  = tmp_37 * tmp_73;
   real_t tmp_80  = tmp_39 * tmp_76;
   real_t tmp_81  = 0.81684757298045851 * tmp_42 + 0.091576213509770743 * tmp_43 + tmp_44;
   real_t tmp_82  = tmp_46 * tmp_81;
   real_t tmp_83  = tmp_48 * tmp_81;
   real_t tmp_84  = tmp_50 * tmp_81;
   real_t tmp_85  = 0.054975871827660928 * Scalar_Variable_Coefficient_3D_g_out0_id2 * tmp_58;
   real_t tmp_86  = 0.10810301816807022 * tmp_1 + 0.44594849091596489 * tmp_2 + tmp_4;
   real_t tmp_87  = tmp_24 * tmp_86;
   real_t tmp_88  = tmp_26 * tmp_86;
   real_t tmp_89  = 0.10810301816807022 * tmp_29 + 0.44594849091596489 * tmp_30 + tmp_31;
   real_t tmp_90  = tmp_33 * tmp_89;
   real_t tmp_91  = tmp_35 * tmp_89;
   real_t tmp_92  = tmp_37 * tmp_86;
   real_t tmp_93  = tmp_39 * tmp_89;
   real_t tmp_94  = 0.10810301816807022 * tmp_42 + 0.44594849091596489 * tmp_43 + tmp_44;
   real_t tmp_95  = tmp_46 * tmp_94;
   real_t tmp_96  = tmp_48 * tmp_94;
   real_t tmp_97  = tmp_50 * tmp_94;
   real_t tmp_98  = 0.11169079483900572 * Scalar_Variable_Coefficient_3D_g_out0_id3 * tmp_58;
   real_t tmp_99  = 0.091576213509770743 * tmp_1 + 0.091576213509770743 * tmp_2 + tmp_4;
   real_t tmp_100 = tmp_24 * tmp_99;
   real_t tmp_101 = tmp_26 * tmp_99;
   real_t tmp_102 = 0.091576213509770743 * tmp_29 + 0.091576213509770743 * tmp_30 + tmp_31;
   real_t tmp_103 = tmp_102 * tmp_33;
   real_t tmp_104 = tmp_102 * tmp_35;
   real_t tmp_105 = tmp_37 * tmp_99;
   real_t tmp_106 = tmp_102 * tmp_39;
   real_t tmp_107 = 0.091576213509770743 * tmp_42 + 0.091576213509770743 * tmp_43 + tmp_44;
   real_t tmp_108 = tmp_107 * tmp_46;
   real_t tmp_109 = tmp_107 * tmp_48;
   real_t tmp_110 = tmp_107 * tmp_50;
   real_t tmp_111 = 0.054975871827660928 * Scalar_Variable_Coefficient_3D_g_out0_id4 * tmp_58;
   real_t tmp_112 = 0.44594849091596489 * tmp_1 + 0.44594849091596489 * tmp_2 + tmp_4;
   real_t tmp_113 = tmp_112 * tmp_24;
   real_t tmp_114 = tmp_112 * tmp_26;
   real_t tmp_115 = 0.44594849091596489 * tmp_29 + 0.44594849091596489 * tmp_30 + tmp_31;
   real_t tmp_116 = tmp_115 * tmp_33;
   real_t tmp_117 = tmp_115 * tmp_35;
   real_t tmp_118 = tmp_112 * tmp_37;
   real_t tmp_119 = tmp_115 * tmp_39;
   real_t tmp_120 = 0.44594849091596489 * tmp_42 + 0.44594849091596489 * tmp_43 + tmp_44;
   real_t tmp_121 = tmp_120 * tmp_46;
   real_t tmp_122 = tmp_120 * tmp_48;
   real_t tmp_123 = tmp_120 * tmp_50;
   real_t tmp_124 = 0.11169079483900572 * Scalar_Variable_Coefficient_3D_g_out0_id5 * tmp_58;
   real_t tmp_125 = -p_affine_13_0 * tmp_50 - p_affine_13_1 * tmp_39 - p_affine_13_2 * tmp_37;
   real_t tmp_126 = -p_affine_13_0 * tmp_48 - p_affine_13_1 * tmp_35 - p_affine_13_2 * tmp_26;
   real_t tmp_127 = -p_affine_13_0 * tmp_46 - p_affine_13_1 * tmp_33 - p_affine_13_2 * tmp_24;
   real_t a_0_0 =
       tmp_111 * ( tmp_56 * ( -tmp_100 - tmp_101 - tmp_103 - tmp_104 - tmp_105 - tmp_106 - tmp_108 - tmp_109 - tmp_110 + 1 ) +
                   tmp_57 ) +
       tmp_124 * ( tmp_56 * ( -tmp_113 - tmp_114 - tmp_116 - tmp_117 - tmp_118 - tmp_119 - tmp_121 - tmp_122 - tmp_123 + 1 ) +
                   tmp_57 ) +
       tmp_59 * ( tmp_56 * ( -tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1 ) + tmp_57 ) +
       tmp_72 * ( tmp_56 * ( -tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1 ) + tmp_57 ) +
       tmp_85 * ( tmp_56 * ( -tmp_74 - tmp_75 - tmp_77 - tmp_78 - tmp_79 - tmp_80 - tmp_82 - tmp_83 - tmp_84 + 1 ) + tmp_57 ) +
       tmp_98 * ( tmp_56 * ( -tmp_87 - tmp_88 - tmp_90 - tmp_91 - tmp_92 - tmp_93 - tmp_95 - tmp_96 - tmp_97 + 1 ) + tmp_57 );
   real_t a_1_0 = tmp_111 * ( tmp_125 + tmp_56 * ( tmp_105 + tmp_106 + tmp_110 ) ) +
                  tmp_124 * ( tmp_125 + tmp_56 * ( tmp_118 + tmp_119 + tmp_123 ) ) +
                  tmp_59 * ( tmp_125 + tmp_56 * ( tmp_38 + tmp_40 + tmp_51 ) ) +
                  tmp_72 * ( tmp_125 + tmp_56 * ( tmp_66 + tmp_67 + tmp_71 ) ) +
                  tmp_85 * ( tmp_125 + tmp_56 * ( tmp_79 + tmp_80 + tmp_84 ) ) +
                  tmp_98 * ( tmp_125 + tmp_56 * ( tmp_92 + tmp_93 + tmp_97 ) );
   real_t a_2_0 = tmp_111 * ( tmp_126 + tmp_56 * ( tmp_101 + tmp_104 + tmp_109 ) ) +
                  tmp_124 * ( tmp_126 + tmp_56 * ( tmp_114 + tmp_117 + tmp_122 ) ) +
                  tmp_59 * ( tmp_126 + tmp_56 * ( tmp_27 + tmp_36 + tmp_49 ) ) +
                  tmp_72 * ( tmp_126 + tmp_56 * ( tmp_62 + tmp_65 + tmp_70 ) ) +
                  tmp_85 * ( tmp_126 + tmp_56 * ( tmp_75 + tmp_78 + tmp_83 ) ) +
                  tmp_98 * ( tmp_126 + tmp_56 * ( tmp_88 + tmp_91 + tmp_96 ) );
   real_t a_3_0 = tmp_111 * ( tmp_127 + tmp_56 * ( tmp_100 + tmp_103 + tmp_108 ) ) +
                  tmp_124 * ( tmp_127 + tmp_56 * ( tmp_113 + tmp_116 + tmp_121 ) ) +
                  tmp_59 * ( tmp_127 + tmp_56 * ( tmp_25 + tmp_34 + tmp_47 ) ) +
                  tmp_72 * ( tmp_127 + tmp_56 * ( tmp_61 + tmp_64 + tmp_69 ) ) +
                  tmp_85 * ( tmp_127 + tmp_56 * ( tmp_74 + tmp_77 + tmp_82 ) ) +
                  tmp_98 * ( tmp_127 + tmp_56 * ( tmp_87 + tmp_90 + tmp_95 ) );

   elMat( 0, 0 ) = a_0_0;
   elMat( 1, 0 ) = a_1_0;
   elMat( 2, 0 ) = a_2_0;
   elMat( 3, 0 ) = a_3_0;
}

} // namespace dg
} // namespace hyteg
