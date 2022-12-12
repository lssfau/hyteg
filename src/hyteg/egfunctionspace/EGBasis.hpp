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

#include "hyteg/types/pointnd.hpp"

#include "Eigen/Eigen"

namespace hyteg {

    class EGBasis {
    public:

        void integrateBasisFunction(uint_t degree,
                                    const std::array<Eigen::Matrix<real_t, 2, 1>, 3> &coords,
                                    const std::function<real_t(const Point3D &)> &f0,
                                    const std::function<real_t(const Point3D &)> &f1,
                                    std::vector<real_t> &values) {
            WALBERLA_CHECK(degree == 0);
            WALBERLA_CHECK(values.size() == 1);

            const auto p_affine_0_0 = coords[0](0);
            const auto p_affine_0_1 = coords[0](1);

            const auto p_affine_1_0 = coords[1](0);
            const auto p_affine_1_1 = coords[1](1);

            const auto p_affine_2_0 = coords[2](0);
            const auto p_affine_2_1 = coords[2](1);

            callback_Scalar_Variable_Coefficient_2D_f0 = f0;
            callback_Scalar_Variable_Coefficient_2D_f1 = f1;

            real_t Scalar_Variable_Coefficient_2D_f0_out0_id0 = 0;
            real_t Scalar_Variable_Coefficient_2D_f1_out0_id1 = 0;
            real_t Scalar_Variable_Coefficient_2D_f0_out0_id2 = 0;
            real_t Scalar_Variable_Coefficient_2D_f1_out0_id3 = 0;
            real_t Scalar_Variable_Coefficient_2D_f0_out0_id4 = 0;
            real_t Scalar_Variable_Coefficient_2D_f1_out0_id5 = 0;
            real_t Scalar_Variable_Coefficient_2D_f0_out0_id6 = 0;
            real_t Scalar_Variable_Coefficient_2D_f1_out0_id7 = 0;
            real_t Scalar_Variable_Coefficient_2D_f0_out0_id8 = 0;
            real_t Scalar_Variable_Coefficient_2D_f1_out0_id9 = 0;
            real_t Scalar_Variable_Coefficient_2D_f0_out0_id10 = 0;
            real_t Scalar_Variable_Coefficient_2D_f1_out0_id11 = 0;

            Scalar_Variable_Coefficient_2D_f0(
                    0.091576213509770743 * p_affine_0_0 + 0.091576213509770743 * p_affine_1_0 +
                    0.81684757298045851 * p_affine_2_0,
                    0.091576213509770743 * p_affine_0_1 + 0.091576213509770743 * p_affine_1_1 +
                    0.81684757298045851 * p_affine_2_1, &Scalar_Variable_Coefficient_2D_f0_out0_id0);
            Scalar_Variable_Coefficient_2D_f1(
                    0.091576213509770743 * p_affine_0_0 + 0.091576213509770743 * p_affine_1_0 +
                    0.81684757298045851 * p_affine_2_0,
                    0.091576213509770743 * p_affine_0_1 + 0.091576213509770743 * p_affine_1_1 +
                    0.81684757298045851 * p_affine_2_1, &Scalar_Variable_Coefficient_2D_f1_out0_id1);
            Scalar_Variable_Coefficient_2D_f0(0.44594849091596489 * p_affine_0_0 + 0.44594849091596489 * p_affine_1_0 +
                                              0.10810301816807022 * p_affine_2_0,
                                              0.44594849091596489 * p_affine_0_1 + 0.44594849091596489 * p_affine_1_1 +
                                              0.10810301816807022 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f0_out0_id2);
            Scalar_Variable_Coefficient_2D_f1(0.44594849091596489 * p_affine_0_0 + 0.44594849091596489 * p_affine_1_0 +
                                              0.10810301816807022 * p_affine_2_0,
                                              0.44594849091596489 * p_affine_0_1 + 0.44594849091596489 * p_affine_1_1 +
                                              0.10810301816807022 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f1_out0_id3);
            Scalar_Variable_Coefficient_2D_f0(0.091576213509770743 * p_affine_0_0 + 0.81684757298045851 * p_affine_1_0 +
                                              0.091576213509770743 * p_affine_2_0,
                                              0.091576213509770743 * p_affine_0_1 + 0.81684757298045851 * p_affine_1_1 +
                                              0.091576213509770743 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f0_out0_id4);
            Scalar_Variable_Coefficient_2D_f1(0.091576213509770743 * p_affine_0_0 + 0.81684757298045851 * p_affine_1_0 +
                                              0.091576213509770743 * p_affine_2_0,
                                              0.091576213509770743 * p_affine_0_1 + 0.81684757298045851 * p_affine_1_1 +
                                              0.091576213509770743 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f1_out0_id5);
            Scalar_Variable_Coefficient_2D_f0(0.44594849091596489 * p_affine_0_0 + 0.10810301816807022 * p_affine_1_0 +
                                              0.44594849091596489 * p_affine_2_0,
                                              0.44594849091596489 * p_affine_0_1 + 0.10810301816807022 * p_affine_1_1 +
                                              0.44594849091596489 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f0_out0_id6);
            Scalar_Variable_Coefficient_2D_f1(0.44594849091596489 * p_affine_0_0 + 0.10810301816807022 * p_affine_1_0 +
                                              0.44594849091596489 * p_affine_2_0,
                                              0.44594849091596489 * p_affine_0_1 + 0.10810301816807022 * p_affine_1_1 +
                                              0.44594849091596489 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f1_out0_id7);
            Scalar_Variable_Coefficient_2D_f0(0.81684757298045851 * p_affine_0_0 + 0.091576213509770743 * p_affine_1_0 +
                                              0.091576213509770743 * p_affine_2_0,
                                              0.81684757298045851 * p_affine_0_1 + 0.091576213509770743 * p_affine_1_1 +
                                              0.091576213509770743 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f0_out0_id8);
            Scalar_Variable_Coefficient_2D_f1(0.81684757298045851 * p_affine_0_0 + 0.091576213509770743 * p_affine_1_0 +
                                              0.091576213509770743 * p_affine_2_0,
                                              0.81684757298045851 * p_affine_0_1 + 0.091576213509770743 * p_affine_1_1 +
                                              0.091576213509770743 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f1_out0_id9);
            Scalar_Variable_Coefficient_2D_f0(0.10810301816807022 * p_affine_0_0 + 0.44594849091596489 * p_affine_1_0 +
                                              0.44594849091596489 * p_affine_2_0,
                                              0.10810301816807022 * p_affine_0_1 + 0.44594849091596489 * p_affine_1_1 +
                                              0.44594849091596489 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f0_out0_id10);
            Scalar_Variable_Coefficient_2D_f1(0.10810301816807022 * p_affine_0_0 + 0.44594849091596489 * p_affine_1_0 +
                                              0.44594849091596489 * p_affine_2_0,
                                              0.10810301816807022 * p_affine_0_1 + 0.44594849091596489 * p_affine_1_1 +
                                              0.44594849091596489 * p_affine_2_1,
                                              &Scalar_Variable_Coefficient_2D_f1_out0_id11);

            real_t tmp_0 = std::abs(
                    p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                    p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0);
            real_t a_0_0 = 0.054975871827660928 * tmp_0 *
                           (-0.24175711982356257 * Scalar_Variable_Coefficient_2D_f0_out0_id0 +
                            0.4835142396471252 * Scalar_Variable_Coefficient_2D_f1_out0_id1) +
                           0.11169079483900572 * tmp_0 *
                           (0.11261515758263158 * Scalar_Variable_Coefficient_2D_f0_out0_id10 +
                            0.11261515758263158 * Scalar_Variable_Coefficient_2D_f1_out0_id11) +
                           0.11169079483900572 * tmp_0 *
                           (0.11261515758263158 * Scalar_Variable_Coefficient_2D_f0_out0_id2 -
                            0.2252303151652631 * Scalar_Variable_Coefficient_2D_f1_out0_id3) +
                           0.054975871827660928 * tmp_0 *
                           (0.4835142396471252 * Scalar_Variable_Coefficient_2D_f0_out0_id4 -
                            0.24175711982356257 * Scalar_Variable_Coefficient_2D_f1_out0_id5) +
                           0.11169079483900572 * tmp_0 *
                           (-0.2252303151652631 * Scalar_Variable_Coefficient_2D_f0_out0_id6 +
                            0.11261515758263158 * Scalar_Variable_Coefficient_2D_f1_out0_id7) +
                           0.054975871827660928 * tmp_0 *
                           (-0.24175711982356257 * Scalar_Variable_Coefficient_2D_f0_out0_id8 -
                            0.24175711982356257 * Scalar_Variable_Coefficient_2D_f1_out0_id9);

            values[0] = a_0_0;
        }

        void integrateBasisFunction(uint_t degree,
                                    const std::array<Eigen::Matrix<real_t, 3, 1>, 4> &coords,
                                    const std::function<real_t(const Point3D &)> &f0,
                                    const std::function<real_t(const Point3D &)> &f1,
                                    const std::function<real_t(const Point3D &)> &f2,
                                    std::vector<real_t> &values) {
            WALBERLA_CHECK(degree == 0);
            WALBERLA_CHECK(values.size() == 1);
            const auto p_affine_0_0 = coords[0](0);
            const auto p_affine_0_1 = coords[0](1);
            const auto p_affine_0_2 = coords[0](2);

            const auto p_affine_1_0 = coords[1](0);
            const auto p_affine_1_1 = coords[1](1);
            const auto p_affine_1_2 = coords[1](2);

            const auto p_affine_2_0 = coords[2](0);
            const auto p_affine_2_1 = coords[2](1);
            const auto p_affine_2_2 = coords[2](2);

            const auto p_affine_3_0 = coords[3](0);
            const auto p_affine_3_1 = coords[3](1);
            const auto p_affine_3_2 = coords[3](2);


            callback_Scalar_Variable_Coefficient_3D_f1 = f0;
            callback_Scalar_Variable_Coefficient_3D_f2 = f1;
            callback_Scalar_Variable_Coefficient_3D_f3 = f2;

            real_t Scalar_Variable_Coefficient_3D_f1_out0_id0 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id1 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id2 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id3 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id4 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id5 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id6 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id7 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id8 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id9 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id10 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id11 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id12 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id13 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id14 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id15 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id16 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id17 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id18 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id19 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id20 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id21 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id22 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id23 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id24 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id25 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id26 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id27 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id28 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id29 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id30 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id31 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id32 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id33 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id34 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id35 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id36 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id37 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id38 = 0;
            real_t Scalar_Variable_Coefficient_3D_f1_out0_id39 = 0;
            real_t Scalar_Variable_Coefficient_3D_f2_out0_id40 = 0;
            real_t Scalar_Variable_Coefficient_3D_f3_out0_id41 = 0;
            Scalar_Variable_Coefficient_3D_f1( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id0 );
            Scalar_Variable_Coefficient_3D_f2( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id1 );
            Scalar_Variable_Coefficient_3D_f3( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id2 );
            Scalar_Variable_Coefficient_3D_f1( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id3 );
            Scalar_Variable_Coefficient_3D_f2( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id4 );
            Scalar_Variable_Coefficient_3D_f3( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id5 );
            Scalar_Variable_Coefficient_3D_f1( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id6 );
            Scalar_Variable_Coefficient_3D_f2( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id7 );
            Scalar_Variable_Coefficient_3D_f3( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id8 );
            Scalar_Variable_Coefficient_3D_f1( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id9 );
            Scalar_Variable_Coefficient_3D_f2( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id10 );
            Scalar_Variable_Coefficient_3D_f3( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id11 );
            Scalar_Variable_Coefficient_3D_f1( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id12 );
            Scalar_Variable_Coefficient_3D_f2( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id13 );
            Scalar_Variable_Coefficient_3D_f3( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id14 );
            Scalar_Variable_Coefficient_3D_f1( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id15 );
            Scalar_Variable_Coefficient_3D_f2( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id16 );
            Scalar_Variable_Coefficient_3D_f3( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id17 );
            Scalar_Variable_Coefficient_3D_f1( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id18 );
            Scalar_Variable_Coefficient_3D_f2( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id19 );
            Scalar_Variable_Coefficient_3D_f3( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id20 );
            Scalar_Variable_Coefficient_3D_f1( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id21 );
            Scalar_Variable_Coefficient_3D_f2( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id22 );
            Scalar_Variable_Coefficient_3D_f3( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id23 );
            Scalar_Variable_Coefficient_3D_f1( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id24 );
            Scalar_Variable_Coefficient_3D_f2( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id25 );
            Scalar_Variable_Coefficient_3D_f3( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id26 );
            Scalar_Variable_Coefficient_3D_f1( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id27 );
            Scalar_Variable_Coefficient_3D_f2( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id28 );
            Scalar_Variable_Coefficient_3D_f3( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id29 );
            Scalar_Variable_Coefficient_3D_f1( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id30 );
            Scalar_Variable_Coefficient_3D_f2( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id31 );
            Scalar_Variable_Coefficient_3D_f3( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id32 );
            Scalar_Variable_Coefficient_3D_f1( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id33 );
            Scalar_Variable_Coefficient_3D_f2( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id34 );
            Scalar_Variable_Coefficient_3D_f3( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id35 );
            Scalar_Variable_Coefficient_3D_f1( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id36 );
            Scalar_Variable_Coefficient_3D_f2( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id37 );
            Scalar_Variable_Coefficient_3D_f3( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id38 );
            Scalar_Variable_Coefficient_3D_f1( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f1_out0_id39 );
            Scalar_Variable_Coefficient_3D_f2( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f2_out0_id40 );
            Scalar_Variable_Coefficient_3D_f3( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_f3_out0_id41 );
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
            real_t a_0_0 = 0.018781320953002646*tmp_12*(0.060885919263300614*Scalar_Variable_Coefficient_3D_f1_out0_id0 + 0.060885919263300614*Scalar_Variable_Coefficient_3D_f2_out0_id1 - 0.18265775778990179*Scalar_Variable_Coefficient_3D_f3_out0_id2) + 0.0070910034628469103*tmp_12*(-0.20449629587435036*Scalar_Variable_Coefficient_3D_f1_out0_id12 + 0.20449629587435036*Scalar_Variable_Coefficient_3D_f2_out0_id13 - 0.20449629587435036*Scalar_Variable_Coefficient_3D_f3_out0_id14) + 0.0070910034628469103*tmp_12*(0.20449629587435036*Scalar_Variable_Coefficient_3D_f1_out0_id15 - 0.20449629587435036*Scalar_Variable_Coefficient_3D_f2_out0_id16 - 0.20449629587435036*Scalar_Variable_Coefficient_3D_f3_out0_id17) + 0.018781320953002646*tmp_12*(0.060885919263300614*Scalar_Variable_Coefficient_3D_f1_out0_id18 - 0.18265775778990179*Scalar_Variable_Coefficient_3D_f2_out0_id19 + 0.060885919263300614*Scalar_Variable_Coefficient_3D_f3_out0_id20) + 0.012248840519393657*tmp_12*(-0.15726474968910875*Scalar_Variable_Coefficient_3D_f1_out0_id21 + 0.47179424906732625*Scalar_Variable_Coefficient_3D_f2_out0_id22 - 0.15726474968910875*Scalar_Variable_Coefficient_3D_f3_out0_id23) + 0.018781320953002646*tmp_12*(-0.18265775778990179*Scalar_Variable_Coefficient_3D_f1_out0_id24 + 0.060885919263300614*Scalar_Variable_Coefficient_3D_f2_out0_id25 + 0.060885919263300614*Scalar_Variable_Coefficient_3D_f3_out0_id26) + 0.012248840519393657*tmp_12*(0.47179424906732625*Scalar_Variable_Coefficient_3D_f1_out0_id27 - 0.15726474968910875*Scalar_Variable_Coefficient_3D_f2_out0_id28 - 0.15726474968910875*Scalar_Variable_Coefficient_3D_f3_out0_id29) + 0.012248840519393657*tmp_12*(-0.15726474968910875*Scalar_Variable_Coefficient_3D_f1_out0_id3 - 0.15726474968910875*Scalar_Variable_Coefficient_3D_f2_out0_id4 + 0.47179424906732625*Scalar_Variable_Coefficient_3D_f3_out0_id5) + 0.018781320953002646*tmp_12*(0.060885919263300614*Scalar_Variable_Coefficient_3D_f1_out0_id30 + 0.060885919263300614*Scalar_Variable_Coefficient_3D_f2_out0_id31 + 0.060885919263300614*Scalar_Variable_Coefficient_3D_f3_out0_id32) + 0.012248840519393657*tmp_12*(-0.15726474968910875*Scalar_Variable_Coefficient_3D_f1_out0_id33 - 0.15726474968910875*Scalar_Variable_Coefficient_3D_f2_out0_id34 - 0.15726474968910875*Scalar_Variable_Coefficient_3D_f3_out0_id35) + 0.0070910034628469103*tmp_12*(-0.20449629587435036*Scalar_Variable_Coefficient_3D_f1_out0_id36 + 0.20449629587435036*Scalar_Variable_Coefficient_3D_f2_out0_id37 + 0.20449629587435036*Scalar_Variable_Coefficient_3D_f3_out0_id38) + 0.0070910034628469103*tmp_12*(0.20449629587435036*Scalar_Variable_Coefficient_3D_f1_out0_id39 - 0.20449629587435036*Scalar_Variable_Coefficient_3D_f2_out0_id40 + 0.20449629587435036*Scalar_Variable_Coefficient_3D_f3_out0_id41) + 0.0070910034628469103*tmp_12*(-0.20449629587435036*Scalar_Variable_Coefficient_3D_f1_out0_id6 - 0.20449629587435036*Scalar_Variable_Coefficient_3D_f2_out0_id7 + 0.20449629587435036*Scalar_Variable_Coefficient_3D_f3_out0_id8) + 0.0070910034628469103*tmp_12*(0.20449629587435036*Scalar_Variable_Coefficient_3D_f1_out0_id9 + 0.20449629587435036*Scalar_Variable_Coefficient_3D_f2_out0_id10 - 0.20449629587435036*Scalar_Variable_Coefficient_3D_f3_out0_id11);

            values[0] = a_0_0;
        }

    private:
        void Scalar_Variable_Coefficient_2D_f0(real_t in_0, real_t in_1, real_t *out_0) const {
            *out_0 = callback_Scalar_Variable_Coefficient_2D_f0(Point3D({in_0, in_1, 0}));
        }

        void Scalar_Variable_Coefficient_2D_f1(real_t in_0, real_t in_1, real_t *out_0) const {
            *out_0 = callback_Scalar_Variable_Coefficient_2D_f1(Point3D({in_0, in_1, 0}));
        }
        std::function<real_t(const Point3D &)> callback_Scalar_Variable_Coefficient_2D_f0;
        std::function<real_t(const Point3D &)> callback_Scalar_Variable_Coefficient_2D_f1;


        void Scalar_Variable_Coefficient_3D_f1(real_t in_0, real_t in_1, real_t in_2, real_t *out_0) const {
            *out_0 = callback_Scalar_Variable_Coefficient_3D_f1(Point3D({in_0, in_1, in_2}));
        }
        void Scalar_Variable_Coefficient_3D_f2(real_t in_0, real_t in_1, real_t in_2, real_t *out_0) const {
            *out_0 = callback_Scalar_Variable_Coefficient_3D_f2(Point3D({in_0, in_1, in_2}));
        }
        void Scalar_Variable_Coefficient_3D_f3(real_t in_0, real_t in_1, real_t in_2, real_t *out_0) const {
            *out_0 = callback_Scalar_Variable_Coefficient_3D_f3(Point3D({in_0, in_1, in_2}));
        }
        std::function<real_t(const Point3D &)> callback_Scalar_Variable_Coefficient_3D_f1;
        std::function<real_t(const Point3D &)> callback_Scalar_Variable_Coefficient_3D_f2;
        std::function<real_t(const Point3D &)> callback_Scalar_Variable_Coefficient_3D_f3;
    };

}