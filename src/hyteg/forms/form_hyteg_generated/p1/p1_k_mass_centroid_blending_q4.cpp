/*
* Copyright (c) 2023-2024 Andreas Burkhart
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

// This file has been generated with the AHFC. If buggy try fixing the generator itself.

#include "p1_k_mass_centroid_blending_q4.hpp"

namespace hyteg {
namespace forms {

   void p1_k_mass_centroid_blending_q4::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      Point3D PoI_affine_0((1.0/3.0)*coords[0][0] + (1.0/3.0)*coords[1][0] + (1.0/3.0)*coords[2][0], (1.0/3.0)*coords[0][1] + (1.0/3.0)*coords[1][1] + (1.0/3.0)*coords[2][1], 0);
      Point3D PoI_blending_0;
      geometryMap_->evalF(Point3D(PoI_affine_0[0], PoI_affine_0[1], 0), PoI_blending_0);
      Matrix2r J_B_1;
      geometryMap_->evalDF(Point3D(0.44594849091596495*coords[0][0] + 0.44594849091596495*coords[1][0] + 0.10810301816807005*coords[2][0], 0.44594849091596495*coords[0][1] + 0.44594849091596495*coords[1][1] + 0.10810301816807005*coords[2][1], 0), J_B_1);
      Matrix2r J_B_2;
      geometryMap_->evalDF(Point3D(0.10810301816807011*coords[0][0] + 0.44594849091596517*coords[1][0] + 0.44594849091596467*coords[2][0], 0.10810301816807011*coords[0][1] + 0.44594849091596517*coords[1][1] + 0.44594849091596467*coords[2][1], 0), J_B_2);
      Matrix2r J_B_3;
      geometryMap_->evalDF(Point3D(0.44594849091596483*coords[0][0] + 0.10810301816806989*coords[1][0] + 0.44594849091596528*coords[2][0], 0.44594849091596483*coords[0][1] + 0.10810301816806989*coords[1][1] + 0.44594849091596528*coords[2][1], 0), J_B_3);
      Matrix2r J_B_4;
      geometryMap_->evalDF(Point3D(0.81684757298045829*coords[0][0] + 0.091576213509771298*coords[1][0] + 0.09157621350977041*coords[2][0], 0.81684757298045829*coords[0][1] + 0.091576213509771298*coords[1][1] + 0.09157621350977041*coords[2][1], 0), J_B_4);
      Matrix2r J_B_5;
      geometryMap_->evalDF(Point3D(0.091576213509770521*coords[0][0] + 0.81684757298045851*coords[1][0] + 0.091576213509770965*coords[2][0], 0.091576213509770521*coords[0][1] + 0.81684757298045851*coords[1][1] + 0.091576213509770965*coords[2][1], 0), J_B_5);
      Matrix2r J_B_6;
      geometryMap_->evalDF(Point3D(0.091576213509770632*coords[0][0] + 0.091576213509770632*coords[1][0] + 0.81684757298045874*coords[2][0], 0.091576213509770632*coords[0][1] + 0.091576213509770632*coords[1][1] + 0.81684757298045874*coords[2][1], 0), J_B_6);
      Point3D B_1;
      geometryMap_->evalF(Point3D(0.44594849091596495*coords[0][0] + 0.44594849091596495*coords[1][0] + 0.10810301816807005*coords[2][0], 0.44594849091596495*coords[0][1] + 0.44594849091596495*coords[1][1] + 0.10810301816807005*coords[2][1], 0), B_1);
      Point3D B_2;
      geometryMap_->evalF(Point3D(0.10810301816807011*coords[0][0] + 0.44594849091596517*coords[1][0] + 0.44594849091596467*coords[2][0], 0.10810301816807011*coords[0][1] + 0.44594849091596517*coords[1][1] + 0.44594849091596467*coords[2][1], 0), B_2);
      Point3D B_3;
      geometryMap_->evalF(Point3D(0.44594849091596483*coords[0][0] + 0.10810301816806989*coords[1][0] + 0.44594849091596528*coords[2][0], 0.44594849091596483*coords[0][1] + 0.10810301816806989*coords[1][1] + 0.44594849091596528*coords[2][1], 0), B_3);
      Point3D B_4;
      geometryMap_->evalF(Point3D(0.81684757298045829*coords[0][0] + 0.091576213509771298*coords[1][0] + 0.09157621350977041*coords[2][0], 0.81684757298045829*coords[0][1] + 0.091576213509771298*coords[1][1] + 0.09157621350977041*coords[2][1], 0), B_4);
      Point3D B_5;
      geometryMap_->evalF(Point3D(0.091576213509770521*coords[0][0] + 0.81684757298045851*coords[1][0] + 0.091576213509770965*coords[2][0], 0.091576213509770521*coords[0][1] + 0.81684757298045851*coords[1][1] + 0.091576213509770965*coords[2][1], 0), B_5);
      Point3D B_6;
      geometryMap_->evalF(Point3D(0.091576213509770632*coords[0][0] + 0.091576213509770632*coords[1][0] + 0.81684757298045874*coords[2][0], 0.091576213509770632*coords[0][1] + 0.091576213509770632*coords[1][1] + 0.81684757298045874*coords[2][1], 0), B_6);
      real_t ScalarCoeff0_1 = k_(Point3D(B_1[0], B_1[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_2 = k_(Point3D(B_2[0], B_2[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_3 = k_(Point3D(B_3[0], B_3[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_4 = k_(Point3D(B_4[0], B_4[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_5 = k_(Point3D(B_5[0], B_5[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_6 = k_(Point3D(B_6[0], B_6[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      const real_t tmp0 = std::abs(coords[0][0]*coords[1][1] - coords[0][0]*coords[2][1] - coords[0][1]*coords[1][0] + coords[0][1]*coords[2][0] + coords[1][0]*coords[2][1] - coords[1][1]*coords[2][0]);
      const real_t tmp1 = tmp0*std::abs(J_B_2(0, 0)*J_B_2(1, 1) - J_B_2(0, 1)*J_B_2(1, 0))*ScalarCoeff0_2;
      const real_t tmp2 = tmp0*std::abs(J_B_5(0, 0)*J_B_5(1, 1) - J_B_5(0, 1)*J_B_5(1, 0))*ScalarCoeff0_5;
      const real_t tmp3 = tmp0*std::abs(J_B_4(0, 0)*J_B_4(1, 1) - J_B_4(0, 1)*J_B_4(1, 0))*ScalarCoeff0_4;
      const real_t tmp4 = tmp0*std::abs(J_B_3(0, 0)*J_B_3(1, 1) - J_B_3(0, 1)*J_B_3(1, 0))*ScalarCoeff0_3;
      const real_t tmp5 = tmp0*std::abs(J_B_6(0, 0)*J_B_6(1, 1) - J_B_6(0, 1)*J_B_6(1, 0))*ScalarCoeff0_6;
      const real_t tmp6 = tmp0*std::abs(J_B_1(0, 0)*J_B_1(1, 1) - J_B_1(0, 1)*J_B_1(1, 0))*ScalarCoeff0_1;
      const real_t tmp7 = 0.00046103881469491135*tmp5 + 0.022211954685772795*tmp6;
      const real_t tmp8 = 0.0053844320361136249*tmp1 + 0.0041124045469858221*tmp2 + 0.004112404546985856*tmp3 + 0.0053844320361136101*tmp4 + tmp7;
      const real_t tmp9 = 0.0041124045469858282*tmp5 + 0.0053844320361136197*tmp6;
      const real_t tmp10 = 0.0053844320361136188*tmp1 + 0.00046103881469491248*tmp2 + 0.0041124045469858161*tmp3 + 0.022211954685772806*tmp4 + tmp9;
      const real_t tmp11 = 0.022211954685772792*tmp1 + 0.0041124045469858421*tmp2 + 0.00046103881469491362*tmp3 + 0.0053844320361136153*tmp4 + tmp9;
      elMat(0,0) = 0.0013052479514599731*tmp1 + 0.00046103881469491026*tmp2 + 0.036682098380937819*tmp3 + 0.022211954685772785*tmp4 + tmp7;
      elMat(0,1) = tmp8;
      elMat(0,2) = tmp10;
      elMat(1,0) = tmp8;
      elMat(1,1) = 0.022211954685772816*tmp1 + 0.03668209838093784*tmp2 + 0.00046103881469491812*tmp3 + 0.0013052479514599677*tmp4 + tmp7;
      elMat(1,2) = tmp11;
      elMat(2,0) = tmp10;
      elMat(2,1) = tmp11;
      elMat(2,2) = 0.022211954685772764*tmp1 + 0.00046103881469491476*tmp2 + 0.00046103881469490912*tmp3 + 0.022211954685772826*tmp4 + 0.036682098380937854*tmp5 + 0.0013052479514599716*tmp6;
   }

   void p1_k_mass_centroid_blending_q4::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const   {
      Point3D PoI_affine_0((1.0/3.0)*coords[0][0] + (1.0/3.0)*coords[1][0] + (1.0/3.0)*coords[2][0], (1.0/3.0)*coords[0][1] + (1.0/3.0)*coords[1][1] + (1.0/3.0)*coords[2][1], 0);
      Point3D PoI_blending_0;
      geometryMap_->evalF(Point3D(PoI_affine_0[0], PoI_affine_0[1], 0), PoI_blending_0);
      Matrix2r J_B_1;
      geometryMap_->evalDF(Point3D(0.44594849091596495*coords[0][0] + 0.44594849091596495*coords[1][0] + 0.10810301816807005*coords[2][0], 0.44594849091596495*coords[0][1] + 0.44594849091596495*coords[1][1] + 0.10810301816807005*coords[2][1], 0), J_B_1);
      Matrix2r J_B_2;
      geometryMap_->evalDF(Point3D(0.10810301816807011*coords[0][0] + 0.44594849091596517*coords[1][0] + 0.44594849091596467*coords[2][0], 0.10810301816807011*coords[0][1] + 0.44594849091596517*coords[1][1] + 0.44594849091596467*coords[2][1], 0), J_B_2);
      Matrix2r J_B_3;
      geometryMap_->evalDF(Point3D(0.44594849091596483*coords[0][0] + 0.10810301816806989*coords[1][0] + 0.44594849091596528*coords[2][0], 0.44594849091596483*coords[0][1] + 0.10810301816806989*coords[1][1] + 0.44594849091596528*coords[2][1], 0), J_B_3);
      Matrix2r J_B_4;
      geometryMap_->evalDF(Point3D(0.81684757298045829*coords[0][0] + 0.091576213509771298*coords[1][0] + 0.09157621350977041*coords[2][0], 0.81684757298045829*coords[0][1] + 0.091576213509771298*coords[1][1] + 0.09157621350977041*coords[2][1], 0), J_B_4);
      Matrix2r J_B_5;
      geometryMap_->evalDF(Point3D(0.091576213509770521*coords[0][0] + 0.81684757298045851*coords[1][0] + 0.091576213509770965*coords[2][0], 0.091576213509770521*coords[0][1] + 0.81684757298045851*coords[1][1] + 0.091576213509770965*coords[2][1], 0), J_B_5);
      Matrix2r J_B_6;
      geometryMap_->evalDF(Point3D(0.091576213509770632*coords[0][0] + 0.091576213509770632*coords[1][0] + 0.81684757298045874*coords[2][0], 0.091576213509770632*coords[0][1] + 0.091576213509770632*coords[1][1] + 0.81684757298045874*coords[2][1], 0), J_B_6);
      Point3D B_1;
      geometryMap_->evalF(Point3D(0.44594849091596495*coords[0][0] + 0.44594849091596495*coords[1][0] + 0.10810301816807005*coords[2][0], 0.44594849091596495*coords[0][1] + 0.44594849091596495*coords[1][1] + 0.10810301816807005*coords[2][1], 0), B_1);
      Point3D B_2;
      geometryMap_->evalF(Point3D(0.10810301816807011*coords[0][0] + 0.44594849091596517*coords[1][0] + 0.44594849091596467*coords[2][0], 0.10810301816807011*coords[0][1] + 0.44594849091596517*coords[1][1] + 0.44594849091596467*coords[2][1], 0), B_2);
      Point3D B_3;
      geometryMap_->evalF(Point3D(0.44594849091596483*coords[0][0] + 0.10810301816806989*coords[1][0] + 0.44594849091596528*coords[2][0], 0.44594849091596483*coords[0][1] + 0.10810301816806989*coords[1][1] + 0.44594849091596528*coords[2][1], 0), B_3);
      Point3D B_4;
      geometryMap_->evalF(Point3D(0.81684757298045829*coords[0][0] + 0.091576213509771298*coords[1][0] + 0.09157621350977041*coords[2][0], 0.81684757298045829*coords[0][1] + 0.091576213509771298*coords[1][1] + 0.09157621350977041*coords[2][1], 0), B_4);
      Point3D B_5;
      geometryMap_->evalF(Point3D(0.091576213509770521*coords[0][0] + 0.81684757298045851*coords[1][0] + 0.091576213509770965*coords[2][0], 0.091576213509770521*coords[0][1] + 0.81684757298045851*coords[1][1] + 0.091576213509770965*coords[2][1], 0), B_5);
      Point3D B_6;
      geometryMap_->evalF(Point3D(0.091576213509770632*coords[0][0] + 0.091576213509770632*coords[1][0] + 0.81684757298045874*coords[2][0], 0.091576213509770632*coords[0][1] + 0.091576213509770632*coords[1][1] + 0.81684757298045874*coords[2][1], 0), B_6);
      real_t ScalarCoeff0_1 = k_(Point3D(B_1[0], B_1[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_2 = k_(Point3D(B_2[0], B_2[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_3 = k_(Point3D(B_3[0], B_3[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_4 = k_(Point3D(B_4[0], B_4[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_5 = k_(Point3D(B_5[0], B_5[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      real_t ScalarCoeff0_6 = k_(Point3D(B_6[0], B_6[1], 0), PoI_blending_0[0], PoI_blending_0[1], 0);
      const real_t tmp0 = std::abs(coords[0][0]*coords[1][1] - coords[0][0]*coords[2][1] - coords[0][1]*coords[1][0] + coords[0][1]*coords[2][0] + coords[1][0]*coords[2][1] - coords[1][1]*coords[2][0]);
      const real_t tmp1 = tmp0*std::abs(J_B_2(0, 0)*J_B_2(1, 1) - J_B_2(0, 1)*J_B_2(1, 0))*ScalarCoeff0_2;
      const real_t tmp2 = tmp0*std::abs(J_B_5(0, 0)*J_B_5(1, 1) - J_B_5(0, 1)*J_B_5(1, 0))*ScalarCoeff0_5;
      const real_t tmp3 = tmp0*std::abs(J_B_4(0, 0)*J_B_4(1, 1) - J_B_4(0, 1)*J_B_4(1, 0))*ScalarCoeff0_4;
      const real_t tmp4 = tmp0*std::abs(J_B_3(0, 0)*J_B_3(1, 1) - J_B_3(0, 1)*J_B_3(1, 0))*ScalarCoeff0_3;
      const real_t tmp5 = tmp0*std::abs(J_B_6(0, 0)*J_B_6(1, 1) - J_B_6(0, 1)*J_B_6(1, 0))*ScalarCoeff0_6;
      const real_t tmp6 = tmp0*std::abs(J_B_1(0, 0)*J_B_1(1, 1) - J_B_1(0, 1)*J_B_1(1, 0))*ScalarCoeff0_1;
      const real_t tmp7 = 0.00046103881469491135*tmp5 + 0.022211954685772795*tmp6;
      elMat(0,0) = 0.0013052479514599731*tmp1 + 0.00046103881469491026*tmp2 + 0.036682098380937819*tmp3 + 0.022211954685772785*tmp4 + tmp7;
      elMat(0,1) = 0.0053844320361136249*tmp1 + 0.0041124045469858221*tmp2 + 0.004112404546985856*tmp3 + 0.0053844320361136101*tmp4 + tmp7;
      elMat(0,2) = 0.0053844320361136188*tmp1 + 0.00046103881469491248*tmp2 + 0.0041124045469858161*tmp3 + 0.022211954685772806*tmp4 + 0.0041124045469858282*tmp5 + 0.0053844320361136197*tmp6;
   }

   void p1_k_mass_centroid_blending_q4::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
   {
      Point3D PoI_affine_0((1.0/4.0)*coords[0][0] + (1.0/4.0)*coords[1][0] + (1.0/4.0)*coords[2][0] + (1.0/4.0)*coords[3][0], (1.0/4.0)*coords[0][1] + (1.0/4.0)*coords[1][1] + (1.0/4.0)*coords[2][1] + (1.0/4.0)*coords[3][1], (1.0/4.0)*coords[0][2] + (1.0/4.0)*coords[1][2] + (1.0/4.0)*coords[2][2] + (1.0/4.0)*coords[3][2]);
      Point3D PoI_blending_0;
      geometryMap_->evalF(Point3D(PoI_affine_0[0], PoI_affine_0[1], PoI_affine_0[2]), PoI_blending_0);
      Matrix3r J_B_1;
      geometryMap_->evalDF(Point3D(0.11180076739738332*coords[0][0] + 0.097204644587583267*coords[1][0] + 0.10660417256199362*coords[2][0] + 0.6843904154530398*coords[3][0], 0.11180076739738332*coords[0][1] + 0.097204644587583267*coords[1][1] + 0.10660417256199362*coords[2][1] + 0.6843904154530398*coords[3][1], 0.11180076739738332*coords[0][2] + 0.097204644587583267*coords[1][2] + 0.10660417256199362*coords[2][2] + 0.6843904154530398*coords[3][2]), J_B_1);
      Matrix3r J_B_2;
      geometryMap_->evalDF(Point3D(0.32329398483747895*coords[0][0] + 0.02956949520647989*coords[1][0] + 0.32923295974264682*coords[2][0] + 0.31790356021339439*coords[3][0], 0.32329398483747895*coords[0][1] + 0.02956949520647989*coords[1][1] + 0.32923295974264682*coords[2][1] + 0.31790356021339439*coords[3][1], 0.32329398483747895*coords[0][2] + 0.02956949520647989*coords[1][2] + 0.32923295974264682*coords[2][2] + 0.31790356021339439*coords[3][2]), J_B_2);
      Matrix3r J_B_3;
      geometryMap_->evalDF(Point3D(0.10962240533194151*coords[0][0] + 0.43271023904776851*coords[1][0] + 0.10384411641099289*coords[2][0] + 0.35382323920929709*coords[3][0], 0.10962240533194151*coords[0][1] + 0.43271023904776851*coords[1][1] + 0.10384411641099289*coords[2][1] + 0.35382323920929709*coords[3][1], 0.10962240533194151*coords[0][2] + 0.43271023904776851*coords[1][2] + 0.10384411641099289*coords[2][2] + 0.35382323920929709*coords[3][2]), J_B_3);
      Matrix3r J_B_4;
      geometryMap_->evalDF(Point3D(0.3284732067220385*coords[0][0] + 0.24027666492807276*coords[1][0] + 0.30444840243449695*coords[2][0] + 0.12680172591539179*coords[3][0], 0.3284732067220385*coords[0][1] + 0.24027666492807276*coords[1][1] + 0.30444840243449695*coords[2][1] + 0.12680172591539179*coords[3][1], 0.3284732067220385*coords[0][2] + 0.24027666492807276*coords[1][2] + 0.30444840243449695*coords[2][2] + 0.12680172591539179*coords[3][2]), J_B_4);
      Matrix3r J_B_5;
      geometryMap_->evalDF(Point3D(0.0023910074574391427*coords[0][0] + 0.12941137378891016*coords[1][0] + 0.53800720391618584*coords[2][0] + 0.33019041483746481*coords[3][0], 0.0023910074574391427*coords[0][1] + 0.12941137378891016*coords[1][1] + 0.53800720391618584*coords[2][1] + 0.33019041483746481*coords[3][1], 0.0023910074574391427*coords[0][2] + 0.12941137378891016*coords[1][2] + 0.53800720391618584*coords[2][2] + 0.33019041483746481*coords[3][2]), J_B_5);
      Matrix3r J_B_6;
      geometryMap_->evalDF(Point3D(0.56297276014304598*coords[0][0] + 0.12154199133392796*coords[1][0] + 0.008991260093335729*coords[2][0] + 0.30649398842969033*coords[3][0], 0.56297276014304598*coords[0][1] + 0.12154199133392796*coords[1][1] + 0.008991260093335729*coords[2][1] + 0.30649398842969033*coords[3][1], 0.56297276014304598*coords[0][2] + 0.12154199133392796*coords[1][2] + 0.008991260093335729*coords[2][2] + 0.30649398842969033*coords[3][2]), J_B_6);
      Matrix3r J_B_7;
      geometryMap_->evalDF(Point3D(0.056824017127933668*coords[0][0] + 0.45076587609127672*coords[1][0] + 0.43295349048135567*coords[2][0] + 0.05945661629943394*coords[3][0], 0.056824017127933668*coords[0][1] + 0.45076587609127672*coords[1][1] + 0.43295349048135567*coords[2][1] + 0.05945661629943394*coords[3][1], 0.056824017127933668*coords[0][2] + 0.45076587609127672*coords[1][2] + 0.43295349048135567*coords[2][2] + 0.05945661629943394*coords[3][2]), J_B_7);
      Matrix3r J_B_8;
      geometryMap_->evalDF(Point3D(0.47961101102565518*coords[0][0] + 0.41926631387951319*coords[1][0] + 0.053341239535745211*coords[2][0] + 0.047781435559086427*coords[3][0], 0.47961101102565518*coords[0][1] + 0.41926631387951319*coords[1][1] + 0.053341239535745211*coords[2][1] + 0.047781435559086427*coords[3][1], 0.47961101102565518*coords[0][2] + 0.41926631387951319*coords[1][2] + 0.053341239535745211*coords[2][2] + 0.047781435559086427*coords[3][2]), J_B_8);
      Matrix3r J_B_9;
      geometryMap_->evalDF(Point3D(0.15636389323939504*coords[0][0] + 0.067223294893383301*coords[1][0] + 0.74122888209362303*coords[2][0] + 0.035183929773598632*coords[3][0], 0.15636389323939504*coords[0][1] + 0.067223294893383301*coords[1][1] + 0.74122888209362303*coords[2][1] + 0.035183929773598632*coords[3][1], 0.15636389323939504*coords[0][2] + 0.067223294893383301*coords[1][2] + 0.74122888209362303*coords[2][2] + 0.035183929773598632*coords[3][2]), J_B_9);
      Matrix3r J_B_10;
      geometryMap_->evalDF(Point3D(0.097987203649279098*coords[0][0] + 0.7525085070096551*coords[1][0] + 0.08140491840285935*coords[2][0] + 0.068099370938206449*coords[3][0], 0.097987203649279098*coords[0][1] + 0.7525085070096551*coords[1][1] + 0.08140491840285935*coords[2][1] + 0.068099370938206449*coords[3][1], 0.097987203649279098*coords[0][2] + 0.7525085070096551*coords[1][2] + 0.08140491840285935*coords[2][2] + 0.068099370938206449*coords[3][2]), J_B_10);
      Matrix3r J_B_11;
      geometryMap_->evalDF(Point3D(0.77125473269537637*coords[0][0] + 0.040490506727590303*coords[1][0] + 0.17469405869723076*coords[2][0] + 0.013560701879802628*coords[3][0], 0.77125473269537637*coords[0][1] + 0.040490506727590303*coords[1][1] + 0.17469405869723076*coords[2][1] + 0.013560701879802628*coords[3][1], 0.77125473269537637*coords[0][2] + 0.040490506727590303*coords[1][2] + 0.17469405869723076*coords[2][2] + 0.013560701879802628*coords[3][2]), J_B_11);
      Point3D B_1;
      geometryMap_->evalF(Point3D(0.11180076739738332*coords[0][0] + 0.097204644587583267*coords[1][0] + 0.10660417256199362*coords[2][0] + 0.6843904154530398*coords[3][0], 0.11180076739738332*coords[0][1] + 0.097204644587583267*coords[1][1] + 0.10660417256199362*coords[2][1] + 0.6843904154530398*coords[3][1], 0.11180076739738332*coords[0][2] + 0.097204644587583267*coords[1][2] + 0.10660417256199362*coords[2][2] + 0.6843904154530398*coords[3][2]), B_1);
      Point3D B_2;
      geometryMap_->evalF(Point3D(0.32329398483747895*coords[0][0] + 0.02956949520647989*coords[1][0] + 0.32923295974264682*coords[2][0] + 0.31790356021339439*coords[3][0], 0.32329398483747895*coords[0][1] + 0.02956949520647989*coords[1][1] + 0.32923295974264682*coords[2][1] + 0.31790356021339439*coords[3][1], 0.32329398483747895*coords[0][2] + 0.02956949520647989*coords[1][2] + 0.32923295974264682*coords[2][2] + 0.31790356021339439*coords[3][2]), B_2);
      Point3D B_3;
      geometryMap_->evalF(Point3D(0.10962240533194151*coords[0][0] + 0.43271023904776851*coords[1][0] + 0.10384411641099289*coords[2][0] + 0.35382323920929709*coords[3][0], 0.10962240533194151*coords[0][1] + 0.43271023904776851*coords[1][1] + 0.10384411641099289*coords[2][1] + 0.35382323920929709*coords[3][1], 0.10962240533194151*coords[0][2] + 0.43271023904776851*coords[1][2] + 0.10384411641099289*coords[2][2] + 0.35382323920929709*coords[3][2]), B_3);
      Point3D B_4;
      geometryMap_->evalF(Point3D(0.3284732067220385*coords[0][0] + 0.24027666492807276*coords[1][0] + 0.30444840243449695*coords[2][0] + 0.12680172591539179*coords[3][0], 0.3284732067220385*coords[0][1] + 0.24027666492807276*coords[1][1] + 0.30444840243449695*coords[2][1] + 0.12680172591539179*coords[3][1], 0.3284732067220385*coords[0][2] + 0.24027666492807276*coords[1][2] + 0.30444840243449695*coords[2][2] + 0.12680172591539179*coords[3][2]), B_4);
      Point3D B_5;
      geometryMap_->evalF(Point3D(0.0023910074574391427*coords[0][0] + 0.12941137378891016*coords[1][0] + 0.53800720391618584*coords[2][0] + 0.33019041483746481*coords[3][0], 0.0023910074574391427*coords[0][1] + 0.12941137378891016*coords[1][1] + 0.53800720391618584*coords[2][1] + 0.33019041483746481*coords[3][1], 0.0023910074574391427*coords[0][2] + 0.12941137378891016*coords[1][2] + 0.53800720391618584*coords[2][2] + 0.33019041483746481*coords[3][2]), B_5);
      Point3D B_6;
      geometryMap_->evalF(Point3D(0.56297276014304598*coords[0][0] + 0.12154199133392796*coords[1][0] + 0.008991260093335729*coords[2][0] + 0.30649398842969033*coords[3][0], 0.56297276014304598*coords[0][1] + 0.12154199133392796*coords[1][1] + 0.008991260093335729*coords[2][1] + 0.30649398842969033*coords[3][1], 0.56297276014304598*coords[0][2] + 0.12154199133392796*coords[1][2] + 0.008991260093335729*coords[2][2] + 0.30649398842969033*coords[3][2]), B_6);
      Point3D B_7;
      geometryMap_->evalF(Point3D(0.056824017127933668*coords[0][0] + 0.45076587609127672*coords[1][0] + 0.43295349048135567*coords[2][0] + 0.05945661629943394*coords[3][0], 0.056824017127933668*coords[0][1] + 0.45076587609127672*coords[1][1] + 0.43295349048135567*coords[2][1] + 0.05945661629943394*coords[3][1], 0.056824017127933668*coords[0][2] + 0.45076587609127672*coords[1][2] + 0.43295349048135567*coords[2][2] + 0.05945661629943394*coords[3][2]), B_7);
      Point3D B_8;
      geometryMap_->evalF(Point3D(0.47961101102565518*coords[0][0] + 0.41926631387951319*coords[1][0] + 0.053341239535745211*coords[2][0] + 0.047781435559086427*coords[3][0], 0.47961101102565518*coords[0][1] + 0.41926631387951319*coords[1][1] + 0.053341239535745211*coords[2][1] + 0.047781435559086427*coords[3][1], 0.47961101102565518*coords[0][2] + 0.41926631387951319*coords[1][2] + 0.053341239535745211*coords[2][2] + 0.047781435559086427*coords[3][2]), B_8);
      Point3D B_9;
      geometryMap_->evalF(Point3D(0.15636389323939504*coords[0][0] + 0.067223294893383301*coords[1][0] + 0.74122888209362303*coords[2][0] + 0.035183929773598632*coords[3][0], 0.15636389323939504*coords[0][1] + 0.067223294893383301*coords[1][1] + 0.74122888209362303*coords[2][1] + 0.035183929773598632*coords[3][1], 0.15636389323939504*coords[0][2] + 0.067223294893383301*coords[1][2] + 0.74122888209362303*coords[2][2] + 0.035183929773598632*coords[3][2]), B_9);
      Point3D B_10;
      geometryMap_->evalF(Point3D(0.097987203649279098*coords[0][0] + 0.7525085070096551*coords[1][0] + 0.08140491840285935*coords[2][0] + 0.068099370938206449*coords[3][0], 0.097987203649279098*coords[0][1] + 0.7525085070096551*coords[1][1] + 0.08140491840285935*coords[2][1] + 0.068099370938206449*coords[3][1], 0.097987203649279098*coords[0][2] + 0.7525085070096551*coords[1][2] + 0.08140491840285935*coords[2][2] + 0.068099370938206449*coords[3][2]), B_10);
      Point3D B_11;
      geometryMap_->evalF(Point3D(0.77125473269537637*coords[0][0] + 0.040490506727590303*coords[1][0] + 0.17469405869723076*coords[2][0] + 0.013560701879802628*coords[3][0], 0.77125473269537637*coords[0][1] + 0.040490506727590303*coords[1][1] + 0.17469405869723076*coords[2][1] + 0.013560701879802628*coords[3][1], 0.77125473269537637*coords[0][2] + 0.040490506727590303*coords[1][2] + 0.17469405869723076*coords[2][2] + 0.013560701879802628*coords[3][2]), B_11);
      real_t ScalarCoeff0_1 = k_(Point3D(B_1[0], B_1[1], B_1[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_2 = k_(Point3D(B_2[0], B_2[1], B_2[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_3 = k_(Point3D(B_3[0], B_3[1], B_3[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_4 = k_(Point3D(B_4[0], B_4[1], B_4[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_5 = k_(Point3D(B_5[0], B_5[1], B_5[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_6 = k_(Point3D(B_6[0], B_6[1], B_6[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_7 = k_(Point3D(B_7[0], B_7[1], B_7[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_8 = k_(Point3D(B_8[0], B_8[1], B_8[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_9 = k_(Point3D(B_9[0], B_9[1], B_9[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_10 = k_(Point3D(B_10[0], B_10[1], B_10[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_11 = k_(Point3D(B_11[0], B_11[1], B_11[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      const real_t tmp0 = coords[0][0]*coords[1][1];
      const real_t tmp1 = coords[0][0]*coords[1][2];
      const real_t tmp2 = coords[2][1]*coords[3][2];
      const real_t tmp3 = coords[0][1]*coords[1][0];
      const real_t tmp4 = coords[0][1]*coords[1][2];
      const real_t tmp5 = coords[2][2]*coords[3][0];
      const real_t tmp6 = coords[0][2]*coords[1][0];
      const real_t tmp7 = coords[0][2]*coords[1][1];
      const real_t tmp8 = coords[2][0]*coords[3][1];
      const real_t tmp9 = coords[2][2]*coords[3][1];
      const real_t tmp10 = coords[2][0]*coords[3][2];
      const real_t tmp11 = coords[2][1]*coords[3][0];
      const real_t tmp12 = std::abs(coords[0][0]*tmp2 - coords[0][0]*tmp9 - coords[0][1]*tmp10 + coords[0][1]*tmp5 - coords[0][2]*tmp11 + coords[0][2]*tmp8 - coords[1][0]*tmp2 + coords[1][0]*tmp9 + coords[1][1]*tmp10 - coords[1][1]*tmp5 + coords[1][2]*tmp11 - coords[1][2]*tmp8 + coords[2][0]*tmp4 - coords[2][0]*tmp7 - coords[2][1]*tmp1 + coords[2][1]*tmp6 + coords[2][2]*tmp0 - coords[2][2]*tmp3 - coords[3][0]*tmp4 + coords[3][0]*tmp7 + coords[3][1]*tmp1 - coords[3][1]*tmp6 - coords[3][2]*tmp0 + coords[3][2]*tmp3);
      const real_t tmp13 = tmp12*std::abs(J_B_1(0, 0)*J_B_1(1, 1)*J_B_1(2, 2) - J_B_1(0, 0)*J_B_1(1, 2)*J_B_1(2, 1) - J_B_1(0, 1)*J_B_1(1, 0)*J_B_1(2, 2) + J_B_1(0, 1)*J_B_1(1, 2)*J_B_1(2, 0) + J_B_1(0, 2)*J_B_1(1, 0)*J_B_1(2, 1) - J_B_1(0, 2)*J_B_1(1, 1)*J_B_1(2, 0))*ScalarCoeff0_1;
      const real_t tmp14 = tmp12*std::abs(J_B_10(0, 0)*J_B_10(1, 1)*J_B_10(2, 2) - J_B_10(0, 0)*J_B_10(1, 2)*J_B_10(2, 1) - J_B_10(0, 1)*J_B_10(1, 0)*J_B_10(2, 2) + J_B_10(0, 1)*J_B_10(1, 2)*J_B_10(2, 0) + J_B_10(0, 2)*J_B_10(1, 0)*J_B_10(2, 1) - J_B_10(0, 2)*J_B_10(1, 1)*J_B_10(2, 0))*ScalarCoeff0_10;
      const real_t tmp15 = tmp12*std::abs(J_B_11(0, 0)*J_B_11(1, 1)*J_B_11(2, 2) - J_B_11(0, 0)*J_B_11(1, 2)*J_B_11(2, 1) - J_B_11(0, 1)*J_B_11(1, 0)*J_B_11(2, 2) + J_B_11(0, 1)*J_B_11(1, 2)*J_B_11(2, 0) + J_B_11(0, 2)*J_B_11(1, 0)*J_B_11(2, 1) - J_B_11(0, 2)*J_B_11(1, 1)*J_B_11(2, 0))*ScalarCoeff0_11;
      const real_t tmp16 = tmp12*std::abs(J_B_2(0, 0)*J_B_2(1, 1)*J_B_2(2, 2) - J_B_2(0, 0)*J_B_2(1, 2)*J_B_2(2, 1) - J_B_2(0, 1)*J_B_2(1, 0)*J_B_2(2, 2) + J_B_2(0, 1)*J_B_2(1, 2)*J_B_2(2, 0) + J_B_2(0, 2)*J_B_2(1, 0)*J_B_2(2, 1) - J_B_2(0, 2)*J_B_2(1, 1)*J_B_2(2, 0))*ScalarCoeff0_2;
      const real_t tmp17 = tmp12*std::abs(J_B_3(0, 0)*J_B_3(1, 1)*J_B_3(2, 2) - J_B_3(0, 0)*J_B_3(1, 2)*J_B_3(2, 1) - J_B_3(0, 1)*J_B_3(1, 0)*J_B_3(2, 2) + J_B_3(0, 1)*J_B_3(1, 2)*J_B_3(2, 0) + J_B_3(0, 2)*J_B_3(1, 0)*J_B_3(2, 1) - J_B_3(0, 2)*J_B_3(1, 1)*J_B_3(2, 0))*ScalarCoeff0_3;
      const real_t tmp18 = tmp12*std::abs(J_B_4(0, 0)*J_B_4(1, 1)*J_B_4(2, 2) - J_B_4(0, 0)*J_B_4(1, 2)*J_B_4(2, 1) - J_B_4(0, 1)*J_B_4(1, 0)*J_B_4(2, 2) + J_B_4(0, 1)*J_B_4(1, 2)*J_B_4(2, 0) + J_B_4(0, 2)*J_B_4(1, 0)*J_B_4(2, 1) - J_B_4(0, 2)*J_B_4(1, 1)*J_B_4(2, 0))*ScalarCoeff0_4;
      const real_t tmp19 = tmp12*std::abs(J_B_5(0, 0)*J_B_5(1, 1)*J_B_5(2, 2) - J_B_5(0, 0)*J_B_5(1, 2)*J_B_5(2, 1) - J_B_5(0, 1)*J_B_5(1, 0)*J_B_5(2, 2) + J_B_5(0, 1)*J_B_5(1, 2)*J_B_5(2, 0) + J_B_5(0, 2)*J_B_5(1, 0)*J_B_5(2, 1) - J_B_5(0, 2)*J_B_5(1, 1)*J_B_5(2, 0))*ScalarCoeff0_5;
      const real_t tmp20 = tmp12*std::abs(J_B_6(0, 0)*J_B_6(1, 1)*J_B_6(2, 2) - J_B_6(0, 0)*J_B_6(1, 2)*J_B_6(2, 1) - J_B_6(0, 1)*J_B_6(1, 0)*J_B_6(2, 2) + J_B_6(0, 1)*J_B_6(1, 2)*J_B_6(2, 0) + J_B_6(0, 2)*J_B_6(1, 0)*J_B_6(2, 1) - J_B_6(0, 2)*J_B_6(1, 1)*J_B_6(2, 0))*ScalarCoeff0_6;
      const real_t tmp21 = tmp12*std::abs(J_B_7(0, 0)*J_B_7(1, 1)*J_B_7(2, 2) - J_B_7(0, 0)*J_B_7(1, 2)*J_B_7(2, 1) - J_B_7(0, 1)*J_B_7(1, 0)*J_B_7(2, 2) + J_B_7(0, 1)*J_B_7(1, 2)*J_B_7(2, 0) + J_B_7(0, 2)*J_B_7(1, 0)*J_B_7(2, 1) - J_B_7(0, 2)*J_B_7(1, 1)*J_B_7(2, 0))*ScalarCoeff0_7;
      const real_t tmp22 = tmp12*std::abs(J_B_8(0, 0)*J_B_8(1, 1)*J_B_8(2, 2) - J_B_8(0, 0)*J_B_8(1, 2)*J_B_8(2, 1) - J_B_8(0, 1)*J_B_8(1, 0)*J_B_8(2, 2) + J_B_8(0, 1)*J_B_8(1, 2)*J_B_8(2, 0) + J_B_8(0, 2)*J_B_8(1, 0)*J_B_8(2, 1) - J_B_8(0, 2)*J_B_8(1, 1)*J_B_8(2, 0))*ScalarCoeff0_8;
      const real_t tmp23 = tmp12*std::abs(J_B_9(0, 0)*J_B_9(1, 1)*J_B_9(2, 2) - J_B_9(0, 0)*J_B_9(1, 2)*J_B_9(2, 1) - J_B_9(0, 1)*J_B_9(1, 0)*J_B_9(2, 2) + J_B_9(0, 1)*J_B_9(1, 2)*J_B_9(2, 0) + J_B_9(0, 2)*J_B_9(1, 0)*J_B_9(2, 1) - J_B_9(0, 2)*J_B_9(1, 1)*J_B_9(2, 0))*ScalarCoeff0_9;
      const real_t tmp24 = 0.00019284118258298789*tmp13 + 0.00067927474029457304*tmp14 + 0.0002042920824082004*tmp15 + 0.0001756332617534064*tmp16 + 0.0012252085563357041*tmp17 + 0.0025441402188600688*tmp18 + 3.927756412637813e-6*tmp19 + 0.00090579285128178464*tmp20 + 0.00029657108070109601*tmp21 + 0.0020086081656411564*tmp22 + 9.7043437061721164e-5*tmp23;
      const real_t tmp25 = 0.0002114886052241351*tmp13 + 7.3482630816416435e-5*tmp14 + 0.00088140692522574384*tmp15 + 0.0019555375630375097*tmp16 + 0.00029403209924464934*tmp17 + 0.0032236148501277274*tmp18 + 1.6329022584012045e-5*tmp19 + 6.7007451722450668e-5*tmp20 + 0.00028485182968767357*tmp21 + 0.00025554557032146538*tmp22 + 0.001070036785341543*tmp23;
      const real_t tmp26 = 0.0013577402358126134*tmp13 + 6.1471972844661093e-5*tmp14 + 6.8419593871221171e-5*tmp15 + 0.0018882445849485883*tmp16 + 0.0010018419278999056*tmp17 + 0.0013426246398866516*tmp18 + 1.0021588375878364e-5*tmp19 + 0.0022841493761420629*tmp20 + 3.9118118486819423e-5*tmp21 + 0.00022890983237355513*tmp22 + 5.0791462691369921e-5*tmp23;
      const real_t tmp27 = 0.0001838777602667592*tmp13 + 0.00056432169454210179*tmp14 + 4.6273444457025656e-5*tmp15 + 0.00017885967975988623*tmp16 + 0.0011606267857981148*tmp17 + 0.0023580596814301185*tmp18 + 0.00088379533850990657*tmp19 + 1.4466453251644587e-5*tmp20 + 0.0022596340606522779*tmp21 + 0.00022339280549000479*tmp22 + 0.0004600256291754938*tmp23;
      const real_t tmp28 = 0.001180480779665163*tmp13 + 0.00047208391285298519*tmp14 + 3.5919961440753758e-6*tmp15 + 0.00017270485014240849*tmp16 + 0.0039545497911395592*tmp17 + 0.00098212385095753632*tmp18 + 0.00054241048694111049*tmp19 + 0.00049313232066488426*tmp20 + 0.00031031091855146937*tmp21 + 0.00020010837829764484*tmp22 + 2.1836047976503086e-5*tmp23;
      const real_t tmp29 = 0.0012946313139199203*tmp13 + 5.1069126856537646e-5*tmp14 + 1.5497469306941285e-5*tmp15 + 0.0019229320141331003*tmp16 + 0.00094903399967577397*tmp17 + 0.001244423953971335*tmp18 + 0.0022549853301921377*tmp19 + 3.6480239519413834e-5*tmp20 + 0.0002980487265059288*tmp21 + 2.5458827925182763e-5*tmp22 + 0.00024077233132705636*tmp23;
      elMat(0,0) = 0.00022179796335936567*tmp13 + 8.8451136021777227e-5*tmp14 + 0.0038913130050347313*tmp15 + 0.0019202619681454574*tmp16 + 0.0003103931843035691*tmp17 + 0.0034779985659016809*tmp18 + 7.2569315962444466e-8*tmp19 + 0.0041955598720028473*tmp20 + 3.7386060177271278e-5*tmp21 + 0.0022977056853520192*tmp22 + 0.00022572665705198103*tmp23;
      elMat(0,1) = tmp24;
      elMat(0,2) = tmp25;
      elMat(0,3) = tmp26;
      elMat(1,0) = tmp24;
      elMat(1,1) = 0.00016766484748893877*tmp13 + 0.0052165997357982675*tmp14 + 1.0725237183613929e-5*tmp15 + 1.6063976241706176e-5*tmp16 + 0.0048362402347407445*tmp17 + 0.0018610270621383387*tmp18 + 0.00021258668670655885*tmp19 + 0.00019555451822012007*tmp20 + 0.0023525989497467184*tmp21 + 0.0017558849198139254*tmp22 + 4.1720498587739075e-5*tmp23;
      elMat(1,2) = tmp27;
      elMat(1,3) = tmp28;
      elMat(2,0) = tmp25;
      elMat(2,1) = tmp27;
      elMat(2,2) = 0.00020165843483053514*tmp13 + 6.104723211663799e-5*tmp14 + 0.00019964422466934559*tmp15 + 0.0019914611776350134*tmp16 + 0.00027853342070056913*tmp17 + 0.0029878369714824581*tmp18 + 0.003674238554035199*tmp19 + 1.0701786467876464e-6*tmp20 + 0.0021703427558741411*tmp21 + 2.8421193770479881e-5*tmp22 + 0.0050724125229055011*tmp23;
      elMat(2,3) = tmp29;
      elMat(3,0) = tmp26;
      elMat(3,1) = tmp28;
      elMat(3,2) = tmp29;
      elMat(3,3) = 0.0083114313586263554*tmp13 + 4.2721932304909983e-5*tmp14 + 1.2029977592257612e-6*tmp15 + 0.0018567610418442587*tmp16 + 0.0032335995094421228*tmp17 + 0.00051829835161631868*tmp18 + 0.0013839490181706453*tmp19 + 0.0012435380573033172*tmp20 + 4.0930421303904594e-5*tmp21 + 2.2805232058805314e-5*tmp22 + 1.1428746236802438e-5*tmp23;
   }

   void p1_k_mass_centroid_blending_q4::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const   {
      Point3D PoI_affine_0((1.0/4.0)*coords[0][0] + (1.0/4.0)*coords[1][0] + (1.0/4.0)*coords[2][0] + (1.0/4.0)*coords[3][0], (1.0/4.0)*coords[0][1] + (1.0/4.0)*coords[1][1] + (1.0/4.0)*coords[2][1] + (1.0/4.0)*coords[3][1], (1.0/4.0)*coords[0][2] + (1.0/4.0)*coords[1][2] + (1.0/4.0)*coords[2][2] + (1.0/4.0)*coords[3][2]);
      Point3D PoI_blending_0;
      geometryMap_->evalF(Point3D(PoI_affine_0[0], PoI_affine_0[1], PoI_affine_0[2]), PoI_blending_0);
      Matrix3r J_B_1;
      geometryMap_->evalDF(Point3D(0.11180076739738332*coords[0][0] + 0.097204644587583267*coords[1][0] + 0.10660417256199362*coords[2][0] + 0.6843904154530398*coords[3][0], 0.11180076739738332*coords[0][1] + 0.097204644587583267*coords[1][1] + 0.10660417256199362*coords[2][1] + 0.6843904154530398*coords[3][1], 0.11180076739738332*coords[0][2] + 0.097204644587583267*coords[1][2] + 0.10660417256199362*coords[2][2] + 0.6843904154530398*coords[3][2]), J_B_1);
      Matrix3r J_B_2;
      geometryMap_->evalDF(Point3D(0.32329398483747895*coords[0][0] + 0.02956949520647989*coords[1][0] + 0.32923295974264682*coords[2][0] + 0.31790356021339439*coords[3][0], 0.32329398483747895*coords[0][1] + 0.02956949520647989*coords[1][1] + 0.32923295974264682*coords[2][1] + 0.31790356021339439*coords[3][1], 0.32329398483747895*coords[0][2] + 0.02956949520647989*coords[1][2] + 0.32923295974264682*coords[2][2] + 0.31790356021339439*coords[3][2]), J_B_2);
      Matrix3r J_B_3;
      geometryMap_->evalDF(Point3D(0.10962240533194151*coords[0][0] + 0.43271023904776851*coords[1][0] + 0.10384411641099289*coords[2][0] + 0.35382323920929709*coords[3][0], 0.10962240533194151*coords[0][1] + 0.43271023904776851*coords[1][1] + 0.10384411641099289*coords[2][1] + 0.35382323920929709*coords[3][1], 0.10962240533194151*coords[0][2] + 0.43271023904776851*coords[1][2] + 0.10384411641099289*coords[2][2] + 0.35382323920929709*coords[3][2]), J_B_3);
      Matrix3r J_B_4;
      geometryMap_->evalDF(Point3D(0.3284732067220385*coords[0][0] + 0.24027666492807276*coords[1][0] + 0.30444840243449695*coords[2][0] + 0.12680172591539179*coords[3][0], 0.3284732067220385*coords[0][1] + 0.24027666492807276*coords[1][1] + 0.30444840243449695*coords[2][1] + 0.12680172591539179*coords[3][1], 0.3284732067220385*coords[0][2] + 0.24027666492807276*coords[1][2] + 0.30444840243449695*coords[2][2] + 0.12680172591539179*coords[3][2]), J_B_4);
      Matrix3r J_B_5;
      geometryMap_->evalDF(Point3D(0.0023910074574391427*coords[0][0] + 0.12941137378891016*coords[1][0] + 0.53800720391618584*coords[2][0] + 0.33019041483746481*coords[3][0], 0.0023910074574391427*coords[0][1] + 0.12941137378891016*coords[1][1] + 0.53800720391618584*coords[2][1] + 0.33019041483746481*coords[3][1], 0.0023910074574391427*coords[0][2] + 0.12941137378891016*coords[1][2] + 0.53800720391618584*coords[2][2] + 0.33019041483746481*coords[3][2]), J_B_5);
      Matrix3r J_B_6;
      geometryMap_->evalDF(Point3D(0.56297276014304598*coords[0][0] + 0.12154199133392796*coords[1][0] + 0.008991260093335729*coords[2][0] + 0.30649398842969033*coords[3][0], 0.56297276014304598*coords[0][1] + 0.12154199133392796*coords[1][1] + 0.008991260093335729*coords[2][1] + 0.30649398842969033*coords[3][1], 0.56297276014304598*coords[0][2] + 0.12154199133392796*coords[1][2] + 0.008991260093335729*coords[2][2] + 0.30649398842969033*coords[3][2]), J_B_6);
      Matrix3r J_B_7;
      geometryMap_->evalDF(Point3D(0.056824017127933668*coords[0][0] + 0.45076587609127672*coords[1][0] + 0.43295349048135567*coords[2][0] + 0.05945661629943394*coords[3][0], 0.056824017127933668*coords[0][1] + 0.45076587609127672*coords[1][1] + 0.43295349048135567*coords[2][1] + 0.05945661629943394*coords[3][1], 0.056824017127933668*coords[0][2] + 0.45076587609127672*coords[1][2] + 0.43295349048135567*coords[2][2] + 0.05945661629943394*coords[3][2]), J_B_7);
      Matrix3r J_B_8;
      geometryMap_->evalDF(Point3D(0.47961101102565518*coords[0][0] + 0.41926631387951319*coords[1][0] + 0.053341239535745211*coords[2][0] + 0.047781435559086427*coords[3][0], 0.47961101102565518*coords[0][1] + 0.41926631387951319*coords[1][1] + 0.053341239535745211*coords[2][1] + 0.047781435559086427*coords[3][1], 0.47961101102565518*coords[0][2] + 0.41926631387951319*coords[1][2] + 0.053341239535745211*coords[2][2] + 0.047781435559086427*coords[3][2]), J_B_8);
      Matrix3r J_B_9;
      geometryMap_->evalDF(Point3D(0.15636389323939504*coords[0][0] + 0.067223294893383301*coords[1][0] + 0.74122888209362303*coords[2][0] + 0.035183929773598632*coords[3][0], 0.15636389323939504*coords[0][1] + 0.067223294893383301*coords[1][1] + 0.74122888209362303*coords[2][1] + 0.035183929773598632*coords[3][1], 0.15636389323939504*coords[0][2] + 0.067223294893383301*coords[1][2] + 0.74122888209362303*coords[2][2] + 0.035183929773598632*coords[3][2]), J_B_9);
      Matrix3r J_B_10;
      geometryMap_->evalDF(Point3D(0.097987203649279098*coords[0][0] + 0.7525085070096551*coords[1][0] + 0.08140491840285935*coords[2][0] + 0.068099370938206449*coords[3][0], 0.097987203649279098*coords[0][1] + 0.7525085070096551*coords[1][1] + 0.08140491840285935*coords[2][1] + 0.068099370938206449*coords[3][1], 0.097987203649279098*coords[0][2] + 0.7525085070096551*coords[1][2] + 0.08140491840285935*coords[2][2] + 0.068099370938206449*coords[3][2]), J_B_10);
      Matrix3r J_B_11;
      geometryMap_->evalDF(Point3D(0.77125473269537637*coords[0][0] + 0.040490506727590303*coords[1][0] + 0.17469405869723076*coords[2][0] + 0.013560701879802628*coords[3][0], 0.77125473269537637*coords[0][1] + 0.040490506727590303*coords[1][1] + 0.17469405869723076*coords[2][1] + 0.013560701879802628*coords[3][1], 0.77125473269537637*coords[0][2] + 0.040490506727590303*coords[1][2] + 0.17469405869723076*coords[2][2] + 0.013560701879802628*coords[3][2]), J_B_11);
      Point3D B_1;
      geometryMap_->evalF(Point3D(0.11180076739738332*coords[0][0] + 0.097204644587583267*coords[1][0] + 0.10660417256199362*coords[2][0] + 0.6843904154530398*coords[3][0], 0.11180076739738332*coords[0][1] + 0.097204644587583267*coords[1][1] + 0.10660417256199362*coords[2][1] + 0.6843904154530398*coords[3][1], 0.11180076739738332*coords[0][2] + 0.097204644587583267*coords[1][2] + 0.10660417256199362*coords[2][2] + 0.6843904154530398*coords[3][2]), B_1);
      Point3D B_2;
      geometryMap_->evalF(Point3D(0.32329398483747895*coords[0][0] + 0.02956949520647989*coords[1][0] + 0.32923295974264682*coords[2][0] + 0.31790356021339439*coords[3][0], 0.32329398483747895*coords[0][1] + 0.02956949520647989*coords[1][1] + 0.32923295974264682*coords[2][1] + 0.31790356021339439*coords[3][1], 0.32329398483747895*coords[0][2] + 0.02956949520647989*coords[1][2] + 0.32923295974264682*coords[2][2] + 0.31790356021339439*coords[3][2]), B_2);
      Point3D B_3;
      geometryMap_->evalF(Point3D(0.10962240533194151*coords[0][0] + 0.43271023904776851*coords[1][0] + 0.10384411641099289*coords[2][0] + 0.35382323920929709*coords[3][0], 0.10962240533194151*coords[0][1] + 0.43271023904776851*coords[1][1] + 0.10384411641099289*coords[2][1] + 0.35382323920929709*coords[3][1], 0.10962240533194151*coords[0][2] + 0.43271023904776851*coords[1][2] + 0.10384411641099289*coords[2][2] + 0.35382323920929709*coords[3][2]), B_3);
      Point3D B_4;
      geometryMap_->evalF(Point3D(0.3284732067220385*coords[0][0] + 0.24027666492807276*coords[1][0] + 0.30444840243449695*coords[2][0] + 0.12680172591539179*coords[3][0], 0.3284732067220385*coords[0][1] + 0.24027666492807276*coords[1][1] + 0.30444840243449695*coords[2][1] + 0.12680172591539179*coords[3][1], 0.3284732067220385*coords[0][2] + 0.24027666492807276*coords[1][2] + 0.30444840243449695*coords[2][2] + 0.12680172591539179*coords[3][2]), B_4);
      Point3D B_5;
      geometryMap_->evalF(Point3D(0.0023910074574391427*coords[0][0] + 0.12941137378891016*coords[1][0] + 0.53800720391618584*coords[2][0] + 0.33019041483746481*coords[3][0], 0.0023910074574391427*coords[0][1] + 0.12941137378891016*coords[1][1] + 0.53800720391618584*coords[2][1] + 0.33019041483746481*coords[3][1], 0.0023910074574391427*coords[0][2] + 0.12941137378891016*coords[1][2] + 0.53800720391618584*coords[2][2] + 0.33019041483746481*coords[3][2]), B_5);
      Point3D B_6;
      geometryMap_->evalF(Point3D(0.56297276014304598*coords[0][0] + 0.12154199133392796*coords[1][0] + 0.008991260093335729*coords[2][0] + 0.30649398842969033*coords[3][0], 0.56297276014304598*coords[0][1] + 0.12154199133392796*coords[1][1] + 0.008991260093335729*coords[2][1] + 0.30649398842969033*coords[3][1], 0.56297276014304598*coords[0][2] + 0.12154199133392796*coords[1][2] + 0.008991260093335729*coords[2][2] + 0.30649398842969033*coords[3][2]), B_6);
      Point3D B_7;
      geometryMap_->evalF(Point3D(0.056824017127933668*coords[0][0] + 0.45076587609127672*coords[1][0] + 0.43295349048135567*coords[2][0] + 0.05945661629943394*coords[3][0], 0.056824017127933668*coords[0][1] + 0.45076587609127672*coords[1][1] + 0.43295349048135567*coords[2][1] + 0.05945661629943394*coords[3][1], 0.056824017127933668*coords[0][2] + 0.45076587609127672*coords[1][2] + 0.43295349048135567*coords[2][2] + 0.05945661629943394*coords[3][2]), B_7);
      Point3D B_8;
      geometryMap_->evalF(Point3D(0.47961101102565518*coords[0][0] + 0.41926631387951319*coords[1][0] + 0.053341239535745211*coords[2][0] + 0.047781435559086427*coords[3][0], 0.47961101102565518*coords[0][1] + 0.41926631387951319*coords[1][1] + 0.053341239535745211*coords[2][1] + 0.047781435559086427*coords[3][1], 0.47961101102565518*coords[0][2] + 0.41926631387951319*coords[1][2] + 0.053341239535745211*coords[2][2] + 0.047781435559086427*coords[3][2]), B_8);
      Point3D B_9;
      geometryMap_->evalF(Point3D(0.15636389323939504*coords[0][0] + 0.067223294893383301*coords[1][0] + 0.74122888209362303*coords[2][0] + 0.035183929773598632*coords[3][0], 0.15636389323939504*coords[0][1] + 0.067223294893383301*coords[1][1] + 0.74122888209362303*coords[2][1] + 0.035183929773598632*coords[3][1], 0.15636389323939504*coords[0][2] + 0.067223294893383301*coords[1][2] + 0.74122888209362303*coords[2][2] + 0.035183929773598632*coords[3][2]), B_9);
      Point3D B_10;
      geometryMap_->evalF(Point3D(0.097987203649279098*coords[0][0] + 0.7525085070096551*coords[1][0] + 0.08140491840285935*coords[2][0] + 0.068099370938206449*coords[3][0], 0.097987203649279098*coords[0][1] + 0.7525085070096551*coords[1][1] + 0.08140491840285935*coords[2][1] + 0.068099370938206449*coords[3][1], 0.097987203649279098*coords[0][2] + 0.7525085070096551*coords[1][2] + 0.08140491840285935*coords[2][2] + 0.068099370938206449*coords[3][2]), B_10);
      Point3D B_11;
      geometryMap_->evalF(Point3D(0.77125473269537637*coords[0][0] + 0.040490506727590303*coords[1][0] + 0.17469405869723076*coords[2][0] + 0.013560701879802628*coords[3][0], 0.77125473269537637*coords[0][1] + 0.040490506727590303*coords[1][1] + 0.17469405869723076*coords[2][1] + 0.013560701879802628*coords[3][1], 0.77125473269537637*coords[0][2] + 0.040490506727590303*coords[1][2] + 0.17469405869723076*coords[2][2] + 0.013560701879802628*coords[3][2]), B_11);
      real_t ScalarCoeff0_1 = k_(Point3D(B_1[0], B_1[1], B_1[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_2 = k_(Point3D(B_2[0], B_2[1], B_2[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_3 = k_(Point3D(B_3[0], B_3[1], B_3[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_4 = k_(Point3D(B_4[0], B_4[1], B_4[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_5 = k_(Point3D(B_5[0], B_5[1], B_5[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_6 = k_(Point3D(B_6[0], B_6[1], B_6[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_7 = k_(Point3D(B_7[0], B_7[1], B_7[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_8 = k_(Point3D(B_8[0], B_8[1], B_8[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_9 = k_(Point3D(B_9[0], B_9[1], B_9[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_10 = k_(Point3D(B_10[0], B_10[1], B_10[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      real_t ScalarCoeff0_11 = k_(Point3D(B_11[0], B_11[1], B_11[2]), PoI_blending_0[0], PoI_blending_0[1], PoI_blending_0[2]);
      const real_t tmp0 = coords[0][0]*coords[1][1];
      const real_t tmp1 = coords[0][0]*coords[1][2];
      const real_t tmp2 = coords[2][1]*coords[3][2];
      const real_t tmp3 = coords[0][1]*coords[1][0];
      const real_t tmp4 = coords[0][1]*coords[1][2];
      const real_t tmp5 = coords[2][2]*coords[3][0];
      const real_t tmp6 = coords[0][2]*coords[1][0];
      const real_t tmp7 = coords[0][2]*coords[1][1];
      const real_t tmp8 = coords[2][0]*coords[3][1];
      const real_t tmp9 = coords[2][2]*coords[3][1];
      const real_t tmp10 = coords[2][0]*coords[3][2];
      const real_t tmp11 = coords[2][1]*coords[3][0];
      const real_t tmp12 = std::abs(coords[0][0]*tmp2 - coords[0][0]*tmp9 - coords[0][1]*tmp10 + coords[0][1]*tmp5 - coords[0][2]*tmp11 + coords[0][2]*tmp8 - coords[1][0]*tmp2 + coords[1][0]*tmp9 + coords[1][1]*tmp10 - coords[1][1]*tmp5 + coords[1][2]*tmp11 - coords[1][2]*tmp8 + coords[2][0]*tmp4 - coords[2][0]*tmp7 - coords[2][1]*tmp1 + coords[2][1]*tmp6 + coords[2][2]*tmp0 - coords[2][2]*tmp3 - coords[3][0]*tmp4 + coords[3][0]*tmp7 + coords[3][1]*tmp1 - coords[3][1]*tmp6 - coords[3][2]*tmp0 + coords[3][2]*tmp3);
      const real_t tmp13 = tmp12*std::abs(J_B_1(0, 0)*J_B_1(1, 1)*J_B_1(2, 2) - J_B_1(0, 0)*J_B_1(1, 2)*J_B_1(2, 1) - J_B_1(0, 1)*J_B_1(1, 0)*J_B_1(2, 2) + J_B_1(0, 1)*J_B_1(1, 2)*J_B_1(2, 0) + J_B_1(0, 2)*J_B_1(1, 0)*J_B_1(2, 1) - J_B_1(0, 2)*J_B_1(1, 1)*J_B_1(2, 0))*ScalarCoeff0_1;
      const real_t tmp14 = tmp12*std::abs(J_B_10(0, 0)*J_B_10(1, 1)*J_B_10(2, 2) - J_B_10(0, 0)*J_B_10(1, 2)*J_B_10(2, 1) - J_B_10(0, 1)*J_B_10(1, 0)*J_B_10(2, 2) + J_B_10(0, 1)*J_B_10(1, 2)*J_B_10(2, 0) + J_B_10(0, 2)*J_B_10(1, 0)*J_B_10(2, 1) - J_B_10(0, 2)*J_B_10(1, 1)*J_B_10(2, 0))*ScalarCoeff0_10;
      const real_t tmp15 = tmp12*std::abs(J_B_11(0, 0)*J_B_11(1, 1)*J_B_11(2, 2) - J_B_11(0, 0)*J_B_11(1, 2)*J_B_11(2, 1) - J_B_11(0, 1)*J_B_11(1, 0)*J_B_11(2, 2) + J_B_11(0, 1)*J_B_11(1, 2)*J_B_11(2, 0) + J_B_11(0, 2)*J_B_11(1, 0)*J_B_11(2, 1) - J_B_11(0, 2)*J_B_11(1, 1)*J_B_11(2, 0))*ScalarCoeff0_11;
      const real_t tmp16 = tmp12*std::abs(J_B_2(0, 0)*J_B_2(1, 1)*J_B_2(2, 2) - J_B_2(0, 0)*J_B_2(1, 2)*J_B_2(2, 1) - J_B_2(0, 1)*J_B_2(1, 0)*J_B_2(2, 2) + J_B_2(0, 1)*J_B_2(1, 2)*J_B_2(2, 0) + J_B_2(0, 2)*J_B_2(1, 0)*J_B_2(2, 1) - J_B_2(0, 2)*J_B_2(1, 1)*J_B_2(2, 0))*ScalarCoeff0_2;
      const real_t tmp17 = tmp12*std::abs(J_B_3(0, 0)*J_B_3(1, 1)*J_B_3(2, 2) - J_B_3(0, 0)*J_B_3(1, 2)*J_B_3(2, 1) - J_B_3(0, 1)*J_B_3(1, 0)*J_B_3(2, 2) + J_B_3(0, 1)*J_B_3(1, 2)*J_B_3(2, 0) + J_B_3(0, 2)*J_B_3(1, 0)*J_B_3(2, 1) - J_B_3(0, 2)*J_B_3(1, 1)*J_B_3(2, 0))*ScalarCoeff0_3;
      const real_t tmp18 = tmp12*std::abs(J_B_4(0, 0)*J_B_4(1, 1)*J_B_4(2, 2) - J_B_4(0, 0)*J_B_4(1, 2)*J_B_4(2, 1) - J_B_4(0, 1)*J_B_4(1, 0)*J_B_4(2, 2) + J_B_4(0, 1)*J_B_4(1, 2)*J_B_4(2, 0) + J_B_4(0, 2)*J_B_4(1, 0)*J_B_4(2, 1) - J_B_4(0, 2)*J_B_4(1, 1)*J_B_4(2, 0))*ScalarCoeff0_4;
      const real_t tmp19 = tmp12*std::abs(J_B_5(0, 0)*J_B_5(1, 1)*J_B_5(2, 2) - J_B_5(0, 0)*J_B_5(1, 2)*J_B_5(2, 1) - J_B_5(0, 1)*J_B_5(1, 0)*J_B_5(2, 2) + J_B_5(0, 1)*J_B_5(1, 2)*J_B_5(2, 0) + J_B_5(0, 2)*J_B_5(1, 0)*J_B_5(2, 1) - J_B_5(0, 2)*J_B_5(1, 1)*J_B_5(2, 0))*ScalarCoeff0_5;
      const real_t tmp20 = tmp12*std::abs(J_B_6(0, 0)*J_B_6(1, 1)*J_B_6(2, 2) - J_B_6(0, 0)*J_B_6(1, 2)*J_B_6(2, 1) - J_B_6(0, 1)*J_B_6(1, 0)*J_B_6(2, 2) + J_B_6(0, 1)*J_B_6(1, 2)*J_B_6(2, 0) + J_B_6(0, 2)*J_B_6(1, 0)*J_B_6(2, 1) - J_B_6(0, 2)*J_B_6(1, 1)*J_B_6(2, 0))*ScalarCoeff0_6;
      const real_t tmp21 = tmp12*std::abs(J_B_7(0, 0)*J_B_7(1, 1)*J_B_7(2, 2) - J_B_7(0, 0)*J_B_7(1, 2)*J_B_7(2, 1) - J_B_7(0, 1)*J_B_7(1, 0)*J_B_7(2, 2) + J_B_7(0, 1)*J_B_7(1, 2)*J_B_7(2, 0) + J_B_7(0, 2)*J_B_7(1, 0)*J_B_7(2, 1) - J_B_7(0, 2)*J_B_7(1, 1)*J_B_7(2, 0))*ScalarCoeff0_7;
      const real_t tmp22 = tmp12*std::abs(J_B_8(0, 0)*J_B_8(1, 1)*J_B_8(2, 2) - J_B_8(0, 0)*J_B_8(1, 2)*J_B_8(2, 1) - J_B_8(0, 1)*J_B_8(1, 0)*J_B_8(2, 2) + J_B_8(0, 1)*J_B_8(1, 2)*J_B_8(2, 0) + J_B_8(0, 2)*J_B_8(1, 0)*J_B_8(2, 1) - J_B_8(0, 2)*J_B_8(1, 1)*J_B_8(2, 0))*ScalarCoeff0_8;
      const real_t tmp23 = tmp12*std::abs(J_B_9(0, 0)*J_B_9(1, 1)*J_B_9(2, 2) - J_B_9(0, 0)*J_B_9(1, 2)*J_B_9(2, 1) - J_B_9(0, 1)*J_B_9(1, 0)*J_B_9(2, 2) + J_B_9(0, 1)*J_B_9(1, 2)*J_B_9(2, 0) + J_B_9(0, 2)*J_B_9(1, 0)*J_B_9(2, 1) - J_B_9(0, 2)*J_B_9(1, 1)*J_B_9(2, 0))*ScalarCoeff0_9;
      elMat(0,0) = 0.00022179796335936567*tmp13 + 8.8451136021777227e-5*tmp14 + 0.0038913130050347313*tmp15 + 0.0019202619681454574*tmp16 + 0.0003103931843035691*tmp17 + 0.0034779985659016809*tmp18 + 7.2569315962444466e-8*tmp19 + 0.0041955598720028473*tmp20 + 3.7386060177271278e-5*tmp21 + 0.0022977056853520192*tmp22 + 0.00022572665705198103*tmp23;
      elMat(0,1) = 0.00019284118258298789*tmp13 + 0.00067927474029457304*tmp14 + 0.0002042920824082004*tmp15 + 0.0001756332617534064*tmp16 + 0.0012252085563357041*tmp17 + 0.0025441402188600688*tmp18 + 3.927756412637813e-6*tmp19 + 0.00090579285128178464*tmp20 + 0.00029657108070109601*tmp21 + 0.0020086081656411564*tmp22 + 9.7043437061721164e-5*tmp23;
      elMat(0,2) = 0.0002114886052241351*tmp13 + 7.3482630816416435e-5*tmp14 + 0.00088140692522574384*tmp15 + 0.0019555375630375097*tmp16 + 0.00029403209924464934*tmp17 + 0.0032236148501277274*tmp18 + 1.6329022584012045e-5*tmp19 + 6.7007451722450668e-5*tmp20 + 0.00028485182968767357*tmp21 + 0.00025554557032146538*tmp22 + 0.001070036785341543*tmp23;
      elMat(0,3) = 0.0013577402358126134*tmp13 + 6.1471972844661093e-5*tmp14 + 6.8419593871221171e-5*tmp15 + 0.0018882445849485883*tmp16 + 0.0010018419278999056*tmp17 + 0.0013426246398866516*tmp18 + 1.0021588375878364e-5*tmp19 + 0.0022841493761420629*tmp20 + 3.9118118486819423e-5*tmp21 + 0.00022890983237355513*tmp22 + 5.0791462691369921e-5*tmp23;
   }

} // namespace forms
} // namespace hyteg