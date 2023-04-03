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

/*
 * The entire file was generated with the HyTeG form generator.
 *
 * Avoid modifying this file. If buggy, consider fixing the generator itself.
 */

#include "p0_linear_form_affine_q6.hpp"

namespace hyteg {
namespace forms {

   void p0_linear_form_affine_q6::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id9 = 0;
      Scalar_Variable_Coefficient_2D_k( 1.2886799075734023*p_affine_0_0 - 0.0087919971442063111*p_affine_1_0 - 0.27988791042919603*p_affine_2_0, 1.2886799075734023*p_affine_0_1 - 0.0087919971442063111*p_affine_1_1 - 0.27988791042919603*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      Scalar_Variable_Coefficient_2D_k( 0.78822920956686404*p_affine_0_0 + 0.064688806056010906*p_affine_1_0 + 0.14708198437712505*p_affine_2_0, 0.78822920956686404*p_affine_0_1 + 0.064688806056010906*p_affine_1_1 + 0.14708198437712505*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id1 );
      Scalar_Variable_Coefficient_2D_k( 0.40172877323475986*p_affine_0_0 + 0.076825716849082099*p_affine_1_0 + 0.52144550991615801*p_affine_2_0, 0.40172877323475986*p_affine_0_1 + 0.076825716849082099*p_affine_1_1 + 0.52144550991615801*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id2 );
      Scalar_Variable_Coefficient_2D_k( 0.64933214716985033*p_affine_0_0 + 0.29502511936941517*p_affine_1_0 + 0.055642733460734617*p_affine_2_0, 0.64933214716985033*p_affine_0_1 + 0.29502511936941517*p_affine_1_1 + 0.055642733460734617*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id3 );
      Scalar_Variable_Coefficient_2D_k( -0.0087919971442063094*p_affine_0_0 + 1.2886799075734023*p_affine_1_0 - 0.27988791042919603*p_affine_2_0, -0.0087919971442063094*p_affine_0_1 + 1.2886799075734023*p_affine_1_1 - 0.27988791042919603*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id4 );
      Scalar_Variable_Coefficient_2D_k( 0.064688806056010906*p_affine_0_0 + 0.78822920956686404*p_affine_1_0 + 0.14708198437712505*p_affine_2_0, 0.064688806056010906*p_affine_0_1 + 0.78822920956686404*p_affine_1_1 + 0.14708198437712505*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id5 );
      Scalar_Variable_Coefficient_2D_k( 0.076825716849082126*p_affine_0_0 + 0.40172877323475981*p_affine_1_0 + 0.52144550991615801*p_affine_2_0, 0.076825716849082126*p_affine_0_1 + 0.40172877323475981*p_affine_1_1 + 0.52144550991615801*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id6 );
      Scalar_Variable_Coefficient_2D_k( 0.29502511936941517*p_affine_0_0 + 0.64933214716985022*p_affine_1_0 + 0.055642733460734617*p_affine_2_0, 0.29502511936941517*p_affine_0_1 + 0.64933214716985022*p_affine_1_1 + 0.055642733460734617*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id7 );
      Scalar_Variable_Coefficient_2D_k( 0.085333161170310645*p_affine_0_0 + 0.085333161170310701*p_affine_1_0 + 0.8293336776593786*p_affine_2_0, 0.085333161170310645*p_affine_0_1 + 0.085333161170310701*p_affine_1_1 + 0.8293336776593786*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id8 );
      Scalar_Variable_Coefficient_2D_k( 0.35369054666996447*p_affine_0_0 + 0.35369054666996447*p_affine_1_0 + 0.29261890666007107*p_affine_2_0, 0.35369054666996447*p_affine_0_1 + 0.35369054666996447*p_affine_1_1 + 0.29261890666007107*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id9 );
      real_t tmp_0 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.00025771437964227723*Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_0 + 0.052397656566423402*Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_0 + 0.074275554650521658*Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_0 + 0.047488619588783712*Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_0 + 0.00025771437964227723*Scalar_Variable_Coefficient_2D_k_out0_id4*tmp_0 + 0.052397656566423402*Scalar_Variable_Coefficient_2D_k_out0_id5*tmp_0 + 0.074275554650521658*Scalar_Variable_Coefficient_2D_k_out0_id6*tmp_0 + 0.047488619588783712*Scalar_Variable_Coefficient_2D_k_out0_id7*tmp_0 + 0.045496761795224737*Scalar_Variable_Coefficient_2D_k_out0_id8*tmp_0 + 0.10566414783403316*Scalar_Variable_Coefficient_2D_k_out0_id9*tmp_0;
      (elMat(0, 0)) = a_0_0;
   }

   void p0_linear_form_affine_q6::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id9 = 0;
      Scalar_Variable_Coefficient_2D_k( 1.2886799075734023*p_affine_0_0 - 0.0087919971442063111*p_affine_1_0 - 0.27988791042919603*p_affine_2_0, 1.2886799075734023*p_affine_0_1 - 0.0087919971442063111*p_affine_1_1 - 0.27988791042919603*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      Scalar_Variable_Coefficient_2D_k( 0.78822920956686404*p_affine_0_0 + 0.064688806056010906*p_affine_1_0 + 0.14708198437712505*p_affine_2_0, 0.78822920956686404*p_affine_0_1 + 0.064688806056010906*p_affine_1_1 + 0.14708198437712505*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id1 );
      Scalar_Variable_Coefficient_2D_k( 0.40172877323475986*p_affine_0_0 + 0.076825716849082099*p_affine_1_0 + 0.52144550991615801*p_affine_2_0, 0.40172877323475986*p_affine_0_1 + 0.076825716849082099*p_affine_1_1 + 0.52144550991615801*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id2 );
      Scalar_Variable_Coefficient_2D_k( 0.64933214716985033*p_affine_0_0 + 0.29502511936941517*p_affine_1_0 + 0.055642733460734617*p_affine_2_0, 0.64933214716985033*p_affine_0_1 + 0.29502511936941517*p_affine_1_1 + 0.055642733460734617*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id3 );
      Scalar_Variable_Coefficient_2D_k( -0.0087919971442063094*p_affine_0_0 + 1.2886799075734023*p_affine_1_0 - 0.27988791042919603*p_affine_2_0, -0.0087919971442063094*p_affine_0_1 + 1.2886799075734023*p_affine_1_1 - 0.27988791042919603*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id4 );
      Scalar_Variable_Coefficient_2D_k( 0.064688806056010906*p_affine_0_0 + 0.78822920956686404*p_affine_1_0 + 0.14708198437712505*p_affine_2_0, 0.064688806056010906*p_affine_0_1 + 0.78822920956686404*p_affine_1_1 + 0.14708198437712505*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id5 );
      Scalar_Variable_Coefficient_2D_k( 0.076825716849082126*p_affine_0_0 + 0.40172877323475981*p_affine_1_0 + 0.52144550991615801*p_affine_2_0, 0.076825716849082126*p_affine_0_1 + 0.40172877323475981*p_affine_1_1 + 0.52144550991615801*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id6 );
      Scalar_Variable_Coefficient_2D_k( 0.29502511936941517*p_affine_0_0 + 0.64933214716985022*p_affine_1_0 + 0.055642733460734617*p_affine_2_0, 0.29502511936941517*p_affine_0_1 + 0.64933214716985022*p_affine_1_1 + 0.055642733460734617*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id7 );
      Scalar_Variable_Coefficient_2D_k( 0.085333161170310645*p_affine_0_0 + 0.085333161170310701*p_affine_1_0 + 0.8293336776593786*p_affine_2_0, 0.085333161170310645*p_affine_0_1 + 0.085333161170310701*p_affine_1_1 + 0.8293336776593786*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id8 );
      Scalar_Variable_Coefficient_2D_k( 0.35369054666996447*p_affine_0_0 + 0.35369054666996447*p_affine_1_0 + 0.29261890666007107*p_affine_2_0, 0.35369054666996447*p_affine_0_1 + 0.35369054666996447*p_affine_1_1 + 0.29261890666007107*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id9 );
      real_t tmp_0 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.00025771437964227723*Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_0 + 0.052397656566423402*Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_0 + 0.074275554650521658*Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_0 + 0.047488619588783712*Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_0 + 0.00025771437964227723*Scalar_Variable_Coefficient_2D_k_out0_id4*tmp_0 + 0.052397656566423402*Scalar_Variable_Coefficient_2D_k_out0_id5*tmp_0 + 0.074275554650521658*Scalar_Variable_Coefficient_2D_k_out0_id6*tmp_0 + 0.047488619588783712*Scalar_Variable_Coefficient_2D_k_out0_id7*tmp_0 + 0.045496761795224737*Scalar_Variable_Coefficient_2D_k_out0_id8*tmp_0 + 0.10566414783403316*Scalar_Variable_Coefficient_2D_k_out0_id9*tmp_0;
      (elMat(0, 0)) = a_0_0;
   }

   void p0_linear_form_affine_q6::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_0_2 = coords[0][2];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_1_2 = coords[1][2];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t p_affine_2_2 = coords[2][2];
      real_t p_affine_3_0 = coords[3][0];
      real_t p_affine_3_1 = coords[3][1];
      real_t p_affine_3_2 = coords[3][2];
      real_t Scalar_Variable_Coefficient_3D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id13 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id14 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id15 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id16 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id17 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id18 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id19 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id20 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id21 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id22 = 0;
      Scalar_Variable_Coefficient_3D_k( 0.024318974248142888*p_affine_0_0 + 0.038836084344884447*p_affine_1_0 + 0.90292879901361134*p_affine_2_0 + 0.03391614239336127*p_affine_3_0, 0.024318974248142888*p_affine_0_1 + 0.038836084344884447*p_affine_1_1 + 0.90292879901361134*p_affine_2_1 + 0.03391614239336127*p_affine_3_1, 0.024318974248142888*p_affine_0_2 + 0.038836084344884447*p_affine_1_2 + 0.90292879901361134*p_affine_2_2 + 0.03391614239336127*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id0 );
      Scalar_Variable_Coefficient_3D_k( 0.022865823814023127*p_affine_0_0 + 0.90377000133218188*p_affine_1_0 + 0.029335721083178661*p_affine_2_0 + 0.044028453770616323*p_affine_3_0, 0.022865823814023127*p_affine_0_1 + 0.90377000133218188*p_affine_1_1 + 0.029335721083178661*p_affine_2_1 + 0.044028453770616323*p_affine_3_1, 0.022865823814023127*p_affine_0_2 + 0.90377000133218188*p_affine_1_2 + 0.029335721083178661*p_affine_2_2 + 0.044028453770616323*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id1 );
      Scalar_Variable_Coefficient_3D_k( 0.62432136355342938*p_affine_0_0 + 0.065924923160009946*p_affine_1_0 + 0.25359367474320033*p_affine_2_0 + 0.05616003854336038*p_affine_3_0, 0.62432136355342938*p_affine_0_1 + 0.065924923160009946*p_affine_1_1 + 0.25359367474320033*p_affine_2_1 + 0.05616003854336038*p_affine_3_1, 0.62432136355342938*p_affine_0_2 + 0.065924923160009946*p_affine_1_2 + 0.25359367474320033*p_affine_2_2 + 0.05616003854336038*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id2 );
      Scalar_Variable_Coefficient_3D_k( 0.063199980942569456*p_affine_0_0 + 0.61745572014726879*p_affine_1_0 + 0.25844914898392563*p_affine_2_0 + 0.060895149926236132*p_affine_3_0, 0.063199980942569456*p_affine_0_1 + 0.61745572014726879*p_affine_1_1 + 0.25844914898392563*p_affine_2_1 + 0.060895149926236132*p_affine_3_1, 0.063199980942569456*p_affine_0_2 + 0.61745572014726879*p_affine_1_2 + 0.25844914898392563*p_affine_2_2 + 0.060895149926236132*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id3 );
      Scalar_Variable_Coefficient_3D_k( 0.24844954011889497*p_affine_0_0 + 0.62654020170888236*p_affine_1_0 + 0.062115533183598748*p_affine_2_0 + 0.062894724988623918*p_affine_3_0, 0.24844954011889497*p_affine_0_1 + 0.62654020170888236*p_affine_1_1 + 0.062115533183598748*p_affine_2_1 + 0.062894724988623918*p_affine_3_1, 0.24844954011889497*p_affine_0_2 + 0.62654020170888236*p_affine_1_2 + 0.062115533183598748*p_affine_2_2 + 0.062894724988623918*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id4 );
      Scalar_Variable_Coefficient_3D_k( 0.063732895294997705*p_affine_0_0 + 0.27790366933007798*p_affine_1_0 + 0.59490968902179553*p_affine_2_0 + 0.063453746353128726*p_affine_3_0, 0.063732895294997705*p_affine_0_1 + 0.27790366933007798*p_affine_1_1 + 0.59490968902179553*p_affine_2_1 + 0.063453746353128726*p_affine_3_1, 0.063732895294997705*p_affine_0_2 + 0.27790366933007798*p_affine_1_2 + 0.59490968902179553*p_affine_2_2 + 0.063453746353128726*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id5 );
      Scalar_Variable_Coefficient_3D_k( 0.06517799276337044*p_affine_0_0 + 0.59471730187579563*p_affine_1_0 + 0.066603298007603146*p_affine_2_0 + 0.27350140735323081*p_affine_3_0, 0.06517799276337044*p_affine_0_1 + 0.59471730187579563*p_affine_1_1 + 0.066603298007603146*p_affine_2_1 + 0.27350140735323081*p_affine_3_1, 0.06517799276337044*p_affine_0_2 + 0.59471730187579563*p_affine_1_2 + 0.066603298007603146*p_affine_2_2 + 0.27350140735323081*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id6 );
      Scalar_Variable_Coefficient_3D_k( 0.083678814060055035*p_affine_0_0 + 0.066098662414680515*p_affine_1_0 + 0.63005455511098962*p_affine_2_0 + 0.2201679684142748*p_affine_3_0, 0.083678814060055035*p_affine_0_1 + 0.066098662414680515*p_affine_1_1 + 0.63005455511098962*p_affine_2_1 + 0.2201679684142748*p_affine_3_1, 0.083678814060055035*p_affine_0_2 + 0.066098662414680515*p_affine_1_2 + 0.63005455511098962*p_affine_2_2 + 0.2201679684142748*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id7 );
      Scalar_Variable_Coefficient_3D_k( 0.5773457813897267*p_affine_0_0 + 0.28772509482646419*p_affine_1_0 + 0.064620638073368528*p_affine_2_0 + 0.07030848571044053*p_affine_3_0, 0.5773457813897267*p_affine_0_1 + 0.28772509482646419*p_affine_1_1 + 0.064620638073368528*p_affine_2_1 + 0.07030848571044053*p_affine_3_1, 0.5773457813897267*p_affine_0_2 + 0.28772509482646419*p_affine_1_2 + 0.064620638073368528*p_affine_2_2 + 0.07030848571044053*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id8 );
      Scalar_Variable_Coefficient_3D_k( 0.038288670738244879*p_affine_0_0 + 0.32838817123122171*p_affine_1_0 + 0.32028743369769253*p_affine_2_0 + 0.31303572433284083*p_affine_3_0, 0.038288670738244879*p_affine_0_1 + 0.32838817123122171*p_affine_1_1 + 0.32028743369769253*p_affine_2_1 + 0.31303572433284083*p_affine_3_1, 0.038288670738244879*p_affine_0_2 + 0.32838817123122171*p_affine_1_2 + 0.32028743369769253*p_affine_2_2 + 0.31303572433284083*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id9 );
      Scalar_Variable_Coefficient_3D_k( 0.35193919733470436*p_affine_0_0 + 0.05509902249072568*p_affine_1_0 + 0.38108430890631018*p_affine_2_0 + 0.2118774712682597*p_affine_3_0, 0.35193919733470436*p_affine_0_1 + 0.05509902249072568*p_affine_1_1 + 0.38108430890631018*p_affine_2_1 + 0.2118774712682597*p_affine_3_1, 0.35193919733470436*p_affine_0_2 + 0.05509902249072568*p_affine_1_2 + 0.38108430890631018*p_affine_2_2 + 0.2118774712682597*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id10 );
      Scalar_Variable_Coefficient_3D_k( 0.53006327548101662*p_affine_0_0 + 0.066789599781738124*p_affine_1_0 + 0.076992717100967251*p_affine_2_0 + 0.32615440763627801*p_affine_3_0, 0.53006327548101662*p_affine_0_1 + 0.066789599781738124*p_affine_1_1 + 0.076992717100967251*p_affine_2_1 + 0.32615440763627801*p_affine_3_1, 0.53006327548101662*p_affine_0_2 + 0.066789599781738124*p_affine_1_2 + 0.076992717100967251*p_affine_2_2 + 0.32615440763627801*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id11 );
      Scalar_Variable_Coefficient_3D_k( 0.0087819577775188895*p_affine_0_0 + 0.0405760510668179*p_affine_1_0 + 0.088600350468910211*p_affine_2_0 + 0.86204164068675293*p_affine_3_0, 0.0087819577775188895*p_affine_0_1 + 0.0405760510668179*p_affine_1_1 + 0.088600350468910211*p_affine_2_1 + 0.86204164068675293*p_affine_3_1, 0.0087819577775188895*p_affine_0_2 + 0.0405760510668179*p_affine_1_2 + 0.088600350468910211*p_affine_2_2 + 0.86204164068675293*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id12 );
      Scalar_Variable_Coefficient_3D_k( 0.15210381130993078*p_affine_0_0 + 0.1246499636374863*p_affine_1_0 + 0.201234567364421*p_affine_2_0 + 0.52201165768816193*p_affine_3_0, 0.15210381130993078*p_affine_0_1 + 0.1246499636374863*p_affine_1_1 + 0.201234567364421*p_affine_2_1 + 0.52201165768816193*p_affine_3_1, 0.15210381130993078*p_affine_0_2 + 0.1246499636374863*p_affine_1_2 + 0.201234567364421*p_affine_2_2 + 0.52201165768816193*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id13 );
      Scalar_Variable_Coefficient_3D_k( 0.30416926534978178*p_affine_0_0 + 0.31919428034893121*p_affine_1_0 + 0.044383344357208213*p_affine_2_0 + 0.33225310994407881*p_affine_3_0, 0.30416926534978178*p_affine_0_1 + 0.31919428034893121*p_affine_1_1 + 0.044383344357208213*p_affine_2_1 + 0.33225310994407881*p_affine_3_1, 0.30416926534978178*p_affine_0_2 + 0.31919428034893121*p_affine_1_2 + 0.044383344357208213*p_affine_2_2 + 0.33225310994407881*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id14 );
      Scalar_Variable_Coefficient_3D_k( 0.25582078426498617*p_affine_0_0 + 0.27942005294598821*p_affine_1_0 + 0.26956992963327198*p_affine_2_0 + 0.19518923315575359*p_affine_3_0, 0.25582078426498617*p_affine_0_1 + 0.27942005294598821*p_affine_1_1 + 0.26956992963327198*p_affine_2_1 + 0.19518923315575359*p_affine_3_1, 0.25582078426498617*p_affine_0_2 + 0.27942005294598821*p_affine_1_2 + 0.26956992963327198*p_affine_2_2 + 0.19518923315575359*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id15 );
      Scalar_Variable_Coefficient_3D_k( 0.84113895166231845*p_affine_0_0 + 0.051325206165202961*p_affine_1_0 + 0.037264752138355499*p_affine_2_0 + 0.070271090034123102*p_affine_3_0, 0.84113895166231845*p_affine_0_1 + 0.051325206165202961*p_affine_1_1 + 0.037264752138355499*p_affine_2_1 + 0.070271090034123102*p_affine_3_1, 0.84113895166231845*p_affine_0_2 + 0.051325206165202961*p_affine_1_2 + 0.037264752138355499*p_affine_2_2 + 0.070271090034123102*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id16 );
      Scalar_Variable_Coefficient_3D_k( 0.21129765858158644*p_affine_0_0 + 0.0073545238380693517*p_affine_1_0 + 0.25118449527752967*p_affine_2_0 + 0.5301633223028146*p_affine_3_0, 0.21129765858158644*p_affine_0_1 + 0.0073545238380693517*p_affine_1_1 + 0.25118449527752967*p_affine_2_1 + 0.5301633223028146*p_affine_3_1, 0.21129765858158644*p_affine_0_2 + 0.0073545238380693517*p_affine_1_2 + 0.25118449527752967*p_affine_2_2 + 0.5301633223028146*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id17 );
      Scalar_Variable_Coefficient_3D_k( 0.023467795573054673*p_affine_0_0 + 0.064775160447105046*p_affine_1_0 + 0.39086205067101182*p_affine_2_0 + 0.52089499330882849*p_affine_3_0, 0.023467795573054673*p_affine_0_1 + 0.064775160447105046*p_affine_1_1 + 0.39086205067101182*p_affine_2_1 + 0.52089499330882849*p_affine_3_1, 0.023467795573054673*p_affine_0_2 + 0.064775160447105046*p_affine_1_2 + 0.39086205067101182*p_affine_2_2 + 0.52089499330882849*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id18 );
      Scalar_Variable_Coefficient_3D_k( 0.21304118323618559*p_affine_0_0 + 0.060010583020269119*p_affine_1_0 + 0.02584268626070331*p_affine_2_0 + 0.70110554748284204*p_affine_3_0, 0.21304118323618559*p_affine_0_1 + 0.060010583020269119*p_affine_1_1 + 0.02584268626070331*p_affine_2_1 + 0.70110554748284204*p_affine_3_1, 0.21304118323618559*p_affine_0_2 + 0.060010583020269119*p_affine_1_2 + 0.02584268626070331*p_affine_2_2 + 0.70110554748284204*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id19 );
      Scalar_Variable_Coefficient_3D_k( 0.26784419818357552*p_affine_0_0 + 0.064769436930052882*p_affine_1_0 + 0.63676750855851394*p_affine_2_0 + 0.030618856327857661*p_affine_3_0, 0.26784419818357552*p_affine_0_1 + 0.064769436930052882*p_affine_1_1 + 0.63676750855851394*p_affine_2_1 + 0.030618856327857661*p_affine_3_1, 0.26784419818357552*p_affine_0_2 + 0.064769436930052882*p_affine_1_2 + 0.63676750855851394*p_affine_2_2 + 0.030618856327857661*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id20 );
      Scalar_Variable_Coefficient_3D_k( 0.053996140835914419*p_affine_0_0 + 0.27578630046985059*p_affine_1_0 + 0.060016149166168681*p_affine_2_0 + 0.61020140952806623*p_affine_3_0, 0.053996140835914419*p_affine_0_1 + 0.27578630046985059*p_affine_1_1 + 0.060016149166168681*p_affine_2_1 + 0.61020140952806623*p_affine_3_1, 0.053996140835914419*p_affine_0_2 + 0.27578630046985059*p_affine_1_2 + 0.060016149166168681*p_affine_2_2 + 0.61020140952806623*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id21 );
      Scalar_Variable_Coefficient_3D_k( 0.3293797185491984*p_affine_0_0 + 0.32511965857702518*p_affine_1_0 + 0.32683350461904581*p_affine_2_0 + 0.018667118254730612*p_affine_3_0, 0.3293797185491984*p_affine_0_1 + 0.32511965857702518*p_affine_1_1 + 0.32683350461904581*p_affine_2_1 + 0.018667118254730612*p_affine_3_1, 0.3293797185491984*p_affine_0_2 + 0.32511965857702518*p_affine_1_2 + 0.32683350461904581*p_affine_2_2 + 0.018667118254730612*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id22 );
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
      real_t a_0_0 = 0.0011826324752765877*Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_12 + 0.0012068794819778289*Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_12 + 0.010435745880218542*Scalar_Variable_Coefficient_3D_k_out0_id10*tmp_12 + 0.010511060314253419*Scalar_Variable_Coefficient_3D_k_out0_id11*tmp_12 + 0.0017372226206159912*Scalar_Variable_Coefficient_3D_k_out0_id12*tmp_12 + 0.010722336995514586*Scalar_Variable_Coefficient_3D_k_out0_id13*tmp_12 + 0.011189302702092835*Scalar_Variable_Coefficient_3D_k_out0_id14*tmp_12 + 0.018766567415677993*Scalar_Variable_Coefficient_3D_k_out0_id15*tmp_12 + 0.0026542465308339578*Scalar_Variable_Coefficient_3D_k_out0_id16*tmp_12 + 0.0037609445463571375*Scalar_Variable_Coefficient_3D_k_out0_id17*tmp_12 + 0.0040385478129073907*Scalar_Variable_Coefficient_3D_k_out0_id18*tmp_12 + 0.0042507207111737392*Scalar_Variable_Coefficient_3D_k_out0_id19*tmp_12 + 0.0072650663434381943*Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_12 + 0.0052515683137844054*Scalar_Variable_Coefficient_3D_k_out0_id20*tmp_12 + 0.0066190162748470446*Scalar_Variable_Coefficient_3D_k_out0_id21*tmp_12 + 0.0072065494492455648*Scalar_Variable_Coefficient_3D_k_out0_id22*tmp_12 + 0.0077688556877634502*Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_12 + 0.0078580050787102011*Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_12 + 0.0081483459837403596*Scalar_Variable_Coefficient_3D_k_out0_id5*tmp_12 + 0.0082947716819190519*Scalar_Variable_Coefficient_3D_k_out0_id6*tmp_12 + 0.0088388873180282289*Scalar_Variable_Coefficient_3D_k_out0_id7*tmp_12 + 0.0089891684380519962*Scalar_Variable_Coefficient_3D_k_out0_id8*tmp_12 + 0.0099702246102381931*Scalar_Variable_Coefficient_3D_k_out0_id9*tmp_12;
      (elMat(0, 0)) = a_0_0;
   }

   void p0_linear_form_affine_q6::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_0_2 = coords[0][2];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_1_2 = coords[1][2];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t p_affine_2_2 = coords[2][2];
      real_t p_affine_3_0 = coords[3][0];
      real_t p_affine_3_1 = coords[3][1];
      real_t p_affine_3_2 = coords[3][2];
      real_t Scalar_Variable_Coefficient_3D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id13 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id14 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id15 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id16 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id17 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id18 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id19 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id20 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id21 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id22 = 0;
      Scalar_Variable_Coefficient_3D_k( 0.024318974248142888*p_affine_0_0 + 0.038836084344884447*p_affine_1_0 + 0.90292879901361134*p_affine_2_0 + 0.03391614239336127*p_affine_3_0, 0.024318974248142888*p_affine_0_1 + 0.038836084344884447*p_affine_1_1 + 0.90292879901361134*p_affine_2_1 + 0.03391614239336127*p_affine_3_1, 0.024318974248142888*p_affine_0_2 + 0.038836084344884447*p_affine_1_2 + 0.90292879901361134*p_affine_2_2 + 0.03391614239336127*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id0 );
      Scalar_Variable_Coefficient_3D_k( 0.022865823814023127*p_affine_0_0 + 0.90377000133218188*p_affine_1_0 + 0.029335721083178661*p_affine_2_0 + 0.044028453770616323*p_affine_3_0, 0.022865823814023127*p_affine_0_1 + 0.90377000133218188*p_affine_1_1 + 0.029335721083178661*p_affine_2_1 + 0.044028453770616323*p_affine_3_1, 0.022865823814023127*p_affine_0_2 + 0.90377000133218188*p_affine_1_2 + 0.029335721083178661*p_affine_2_2 + 0.044028453770616323*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id1 );
      Scalar_Variable_Coefficient_3D_k( 0.62432136355342938*p_affine_0_0 + 0.065924923160009946*p_affine_1_0 + 0.25359367474320033*p_affine_2_0 + 0.05616003854336038*p_affine_3_0, 0.62432136355342938*p_affine_0_1 + 0.065924923160009946*p_affine_1_1 + 0.25359367474320033*p_affine_2_1 + 0.05616003854336038*p_affine_3_1, 0.62432136355342938*p_affine_0_2 + 0.065924923160009946*p_affine_1_2 + 0.25359367474320033*p_affine_2_2 + 0.05616003854336038*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id2 );
      Scalar_Variable_Coefficient_3D_k( 0.063199980942569456*p_affine_0_0 + 0.61745572014726879*p_affine_1_0 + 0.25844914898392563*p_affine_2_0 + 0.060895149926236132*p_affine_3_0, 0.063199980942569456*p_affine_0_1 + 0.61745572014726879*p_affine_1_1 + 0.25844914898392563*p_affine_2_1 + 0.060895149926236132*p_affine_3_1, 0.063199980942569456*p_affine_0_2 + 0.61745572014726879*p_affine_1_2 + 0.25844914898392563*p_affine_2_2 + 0.060895149926236132*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id3 );
      Scalar_Variable_Coefficient_3D_k( 0.24844954011889497*p_affine_0_0 + 0.62654020170888236*p_affine_1_0 + 0.062115533183598748*p_affine_2_0 + 0.062894724988623918*p_affine_3_0, 0.24844954011889497*p_affine_0_1 + 0.62654020170888236*p_affine_1_1 + 0.062115533183598748*p_affine_2_1 + 0.062894724988623918*p_affine_3_1, 0.24844954011889497*p_affine_0_2 + 0.62654020170888236*p_affine_1_2 + 0.062115533183598748*p_affine_2_2 + 0.062894724988623918*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id4 );
      Scalar_Variable_Coefficient_3D_k( 0.063732895294997705*p_affine_0_0 + 0.27790366933007798*p_affine_1_0 + 0.59490968902179553*p_affine_2_0 + 0.063453746353128726*p_affine_3_0, 0.063732895294997705*p_affine_0_1 + 0.27790366933007798*p_affine_1_1 + 0.59490968902179553*p_affine_2_1 + 0.063453746353128726*p_affine_3_1, 0.063732895294997705*p_affine_0_2 + 0.27790366933007798*p_affine_1_2 + 0.59490968902179553*p_affine_2_2 + 0.063453746353128726*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id5 );
      Scalar_Variable_Coefficient_3D_k( 0.06517799276337044*p_affine_0_0 + 0.59471730187579563*p_affine_1_0 + 0.066603298007603146*p_affine_2_0 + 0.27350140735323081*p_affine_3_0, 0.06517799276337044*p_affine_0_1 + 0.59471730187579563*p_affine_1_1 + 0.066603298007603146*p_affine_2_1 + 0.27350140735323081*p_affine_3_1, 0.06517799276337044*p_affine_0_2 + 0.59471730187579563*p_affine_1_2 + 0.066603298007603146*p_affine_2_2 + 0.27350140735323081*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id6 );
      Scalar_Variable_Coefficient_3D_k( 0.083678814060055035*p_affine_0_0 + 0.066098662414680515*p_affine_1_0 + 0.63005455511098962*p_affine_2_0 + 0.2201679684142748*p_affine_3_0, 0.083678814060055035*p_affine_0_1 + 0.066098662414680515*p_affine_1_1 + 0.63005455511098962*p_affine_2_1 + 0.2201679684142748*p_affine_3_1, 0.083678814060055035*p_affine_0_2 + 0.066098662414680515*p_affine_1_2 + 0.63005455511098962*p_affine_2_2 + 0.2201679684142748*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id7 );
      Scalar_Variable_Coefficient_3D_k( 0.5773457813897267*p_affine_0_0 + 0.28772509482646419*p_affine_1_0 + 0.064620638073368528*p_affine_2_0 + 0.07030848571044053*p_affine_3_0, 0.5773457813897267*p_affine_0_1 + 0.28772509482646419*p_affine_1_1 + 0.064620638073368528*p_affine_2_1 + 0.07030848571044053*p_affine_3_1, 0.5773457813897267*p_affine_0_2 + 0.28772509482646419*p_affine_1_2 + 0.064620638073368528*p_affine_2_2 + 0.07030848571044053*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id8 );
      Scalar_Variable_Coefficient_3D_k( 0.038288670738244879*p_affine_0_0 + 0.32838817123122171*p_affine_1_0 + 0.32028743369769253*p_affine_2_0 + 0.31303572433284083*p_affine_3_0, 0.038288670738244879*p_affine_0_1 + 0.32838817123122171*p_affine_1_1 + 0.32028743369769253*p_affine_2_1 + 0.31303572433284083*p_affine_3_1, 0.038288670738244879*p_affine_0_2 + 0.32838817123122171*p_affine_1_2 + 0.32028743369769253*p_affine_2_2 + 0.31303572433284083*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id9 );
      Scalar_Variable_Coefficient_3D_k( 0.35193919733470436*p_affine_0_0 + 0.05509902249072568*p_affine_1_0 + 0.38108430890631018*p_affine_2_0 + 0.2118774712682597*p_affine_3_0, 0.35193919733470436*p_affine_0_1 + 0.05509902249072568*p_affine_1_1 + 0.38108430890631018*p_affine_2_1 + 0.2118774712682597*p_affine_3_1, 0.35193919733470436*p_affine_0_2 + 0.05509902249072568*p_affine_1_2 + 0.38108430890631018*p_affine_2_2 + 0.2118774712682597*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id10 );
      Scalar_Variable_Coefficient_3D_k( 0.53006327548101662*p_affine_0_0 + 0.066789599781738124*p_affine_1_0 + 0.076992717100967251*p_affine_2_0 + 0.32615440763627801*p_affine_3_0, 0.53006327548101662*p_affine_0_1 + 0.066789599781738124*p_affine_1_1 + 0.076992717100967251*p_affine_2_1 + 0.32615440763627801*p_affine_3_1, 0.53006327548101662*p_affine_0_2 + 0.066789599781738124*p_affine_1_2 + 0.076992717100967251*p_affine_2_2 + 0.32615440763627801*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id11 );
      Scalar_Variable_Coefficient_3D_k( 0.0087819577775188895*p_affine_0_0 + 0.0405760510668179*p_affine_1_0 + 0.088600350468910211*p_affine_2_0 + 0.86204164068675293*p_affine_3_0, 0.0087819577775188895*p_affine_0_1 + 0.0405760510668179*p_affine_1_1 + 0.088600350468910211*p_affine_2_1 + 0.86204164068675293*p_affine_3_1, 0.0087819577775188895*p_affine_0_2 + 0.0405760510668179*p_affine_1_2 + 0.088600350468910211*p_affine_2_2 + 0.86204164068675293*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id12 );
      Scalar_Variable_Coefficient_3D_k( 0.15210381130993078*p_affine_0_0 + 0.1246499636374863*p_affine_1_0 + 0.201234567364421*p_affine_2_0 + 0.52201165768816193*p_affine_3_0, 0.15210381130993078*p_affine_0_1 + 0.1246499636374863*p_affine_1_1 + 0.201234567364421*p_affine_2_1 + 0.52201165768816193*p_affine_3_1, 0.15210381130993078*p_affine_0_2 + 0.1246499636374863*p_affine_1_2 + 0.201234567364421*p_affine_2_2 + 0.52201165768816193*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id13 );
      Scalar_Variable_Coefficient_3D_k( 0.30416926534978178*p_affine_0_0 + 0.31919428034893121*p_affine_1_0 + 0.044383344357208213*p_affine_2_0 + 0.33225310994407881*p_affine_3_0, 0.30416926534978178*p_affine_0_1 + 0.31919428034893121*p_affine_1_1 + 0.044383344357208213*p_affine_2_1 + 0.33225310994407881*p_affine_3_1, 0.30416926534978178*p_affine_0_2 + 0.31919428034893121*p_affine_1_2 + 0.044383344357208213*p_affine_2_2 + 0.33225310994407881*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id14 );
      Scalar_Variable_Coefficient_3D_k( 0.25582078426498617*p_affine_0_0 + 0.27942005294598821*p_affine_1_0 + 0.26956992963327198*p_affine_2_0 + 0.19518923315575359*p_affine_3_0, 0.25582078426498617*p_affine_0_1 + 0.27942005294598821*p_affine_1_1 + 0.26956992963327198*p_affine_2_1 + 0.19518923315575359*p_affine_3_1, 0.25582078426498617*p_affine_0_2 + 0.27942005294598821*p_affine_1_2 + 0.26956992963327198*p_affine_2_2 + 0.19518923315575359*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id15 );
      Scalar_Variable_Coefficient_3D_k( 0.84113895166231845*p_affine_0_0 + 0.051325206165202961*p_affine_1_0 + 0.037264752138355499*p_affine_2_0 + 0.070271090034123102*p_affine_3_0, 0.84113895166231845*p_affine_0_1 + 0.051325206165202961*p_affine_1_1 + 0.037264752138355499*p_affine_2_1 + 0.070271090034123102*p_affine_3_1, 0.84113895166231845*p_affine_0_2 + 0.051325206165202961*p_affine_1_2 + 0.037264752138355499*p_affine_2_2 + 0.070271090034123102*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id16 );
      Scalar_Variable_Coefficient_3D_k( 0.21129765858158644*p_affine_0_0 + 0.0073545238380693517*p_affine_1_0 + 0.25118449527752967*p_affine_2_0 + 0.5301633223028146*p_affine_3_0, 0.21129765858158644*p_affine_0_1 + 0.0073545238380693517*p_affine_1_1 + 0.25118449527752967*p_affine_2_1 + 0.5301633223028146*p_affine_3_1, 0.21129765858158644*p_affine_0_2 + 0.0073545238380693517*p_affine_1_2 + 0.25118449527752967*p_affine_2_2 + 0.5301633223028146*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id17 );
      Scalar_Variable_Coefficient_3D_k( 0.023467795573054673*p_affine_0_0 + 0.064775160447105046*p_affine_1_0 + 0.39086205067101182*p_affine_2_0 + 0.52089499330882849*p_affine_3_0, 0.023467795573054673*p_affine_0_1 + 0.064775160447105046*p_affine_1_1 + 0.39086205067101182*p_affine_2_1 + 0.52089499330882849*p_affine_3_1, 0.023467795573054673*p_affine_0_2 + 0.064775160447105046*p_affine_1_2 + 0.39086205067101182*p_affine_2_2 + 0.52089499330882849*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id18 );
      Scalar_Variable_Coefficient_3D_k( 0.21304118323618559*p_affine_0_0 + 0.060010583020269119*p_affine_1_0 + 0.02584268626070331*p_affine_2_0 + 0.70110554748284204*p_affine_3_0, 0.21304118323618559*p_affine_0_1 + 0.060010583020269119*p_affine_1_1 + 0.02584268626070331*p_affine_2_1 + 0.70110554748284204*p_affine_3_1, 0.21304118323618559*p_affine_0_2 + 0.060010583020269119*p_affine_1_2 + 0.02584268626070331*p_affine_2_2 + 0.70110554748284204*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id19 );
      Scalar_Variable_Coefficient_3D_k( 0.26784419818357552*p_affine_0_0 + 0.064769436930052882*p_affine_1_0 + 0.63676750855851394*p_affine_2_0 + 0.030618856327857661*p_affine_3_0, 0.26784419818357552*p_affine_0_1 + 0.064769436930052882*p_affine_1_1 + 0.63676750855851394*p_affine_2_1 + 0.030618856327857661*p_affine_3_1, 0.26784419818357552*p_affine_0_2 + 0.064769436930052882*p_affine_1_2 + 0.63676750855851394*p_affine_2_2 + 0.030618856327857661*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id20 );
      Scalar_Variable_Coefficient_3D_k( 0.053996140835914419*p_affine_0_0 + 0.27578630046985059*p_affine_1_0 + 0.060016149166168681*p_affine_2_0 + 0.61020140952806623*p_affine_3_0, 0.053996140835914419*p_affine_0_1 + 0.27578630046985059*p_affine_1_1 + 0.060016149166168681*p_affine_2_1 + 0.61020140952806623*p_affine_3_1, 0.053996140835914419*p_affine_0_2 + 0.27578630046985059*p_affine_1_2 + 0.060016149166168681*p_affine_2_2 + 0.61020140952806623*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id21 );
      Scalar_Variable_Coefficient_3D_k( 0.3293797185491984*p_affine_0_0 + 0.32511965857702518*p_affine_1_0 + 0.32683350461904581*p_affine_2_0 + 0.018667118254730612*p_affine_3_0, 0.3293797185491984*p_affine_0_1 + 0.32511965857702518*p_affine_1_1 + 0.32683350461904581*p_affine_2_1 + 0.018667118254730612*p_affine_3_1, 0.3293797185491984*p_affine_0_2 + 0.32511965857702518*p_affine_1_2 + 0.32683350461904581*p_affine_2_2 + 0.018667118254730612*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id22 );
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
      real_t a_0_0 = 0.0011826324752765877*Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_12 + 0.0012068794819778289*Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_12 + 0.010435745880218542*Scalar_Variable_Coefficient_3D_k_out0_id10*tmp_12 + 0.010511060314253419*Scalar_Variable_Coefficient_3D_k_out0_id11*tmp_12 + 0.0017372226206159912*Scalar_Variable_Coefficient_3D_k_out0_id12*tmp_12 + 0.010722336995514586*Scalar_Variable_Coefficient_3D_k_out0_id13*tmp_12 + 0.011189302702092835*Scalar_Variable_Coefficient_3D_k_out0_id14*tmp_12 + 0.018766567415677993*Scalar_Variable_Coefficient_3D_k_out0_id15*tmp_12 + 0.0026542465308339578*Scalar_Variable_Coefficient_3D_k_out0_id16*tmp_12 + 0.0037609445463571375*Scalar_Variable_Coefficient_3D_k_out0_id17*tmp_12 + 0.0040385478129073907*Scalar_Variable_Coefficient_3D_k_out0_id18*tmp_12 + 0.0042507207111737392*Scalar_Variable_Coefficient_3D_k_out0_id19*tmp_12 + 0.0072650663434381943*Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_12 + 0.0052515683137844054*Scalar_Variable_Coefficient_3D_k_out0_id20*tmp_12 + 0.0066190162748470446*Scalar_Variable_Coefficient_3D_k_out0_id21*tmp_12 + 0.0072065494492455648*Scalar_Variable_Coefficient_3D_k_out0_id22*tmp_12 + 0.0077688556877634502*Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_12 + 0.0078580050787102011*Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_12 + 0.0081483459837403596*Scalar_Variable_Coefficient_3D_k_out0_id5*tmp_12 + 0.0082947716819190519*Scalar_Variable_Coefficient_3D_k_out0_id6*tmp_12 + 0.0088388873180282289*Scalar_Variable_Coefficient_3D_k_out0_id7*tmp_12 + 0.0089891684380519962*Scalar_Variable_Coefficient_3D_k_out0_id8*tmp_12 + 0.0099702246102381931*Scalar_Variable_Coefficient_3D_k_out0_id9*tmp_12;
      (elMat(0, 0)) = a_0_0;
   }

   void p0_linear_form_affine_q6::Scalar_Variable_Coefficient_2D_k( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_k( Point3D( {in_0, in_1, 0} ) );
   }

   void p0_linear_form_affine_q6::Scalar_Variable_Coefficient_3D_k( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_k( Point3D( {in_0, in_1, in_2} ) );
   }

} // namespace forms
} // namespace hyteg
