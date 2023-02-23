/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

// Some selected quadrature/cubature formulas from the FEM book by Lang & Junger, pp 235

namespace hyteg {

namespace quadrature {

/// 2DD-5: exact for polynomial integrands up to order 3
static const std::array< Point2D, 7 > D5_points = { Point2D( 0.0, 0.0 ),
                                                    Point2D( 1.0, 0.0 ),
                                                    Point2D( 0.0, 1.0 ),
                                                    Point2D( 0.5, 0.0 ),
                                                    Point2D( 0.5, 0.5 ),
                                                    Point2D( 0.0, 0.5 ),
                                                    Point2D( 1.0 / 3.0, 1.0 / 3.0 ) };

static const std::array< real_t, 7 > D5_weights =
    { 3.0 / 120.0, 3.0 / 120.0, 3.0 / 120.0, 8.0 / 120.0, 8.0 / 120.0, 8.0 / 120.0, 27.0 / 120.0 };

/// 2DD-6: exact for polynomial integrands up to order 5
static const std::array< real_t, 7 > D6_weights = { 0.06296959027241,
                                                    0.06296959027241,
                                                    0.06296959027241,
                                                    0.06619707639425,
                                                    0.06619707639425,
                                                    0.06619707639425,
                                                    0.11250000000000 };

static const std::array< Point2D, 7 > D6_points = { Point2D( 0.10128650732346, 0.10128650732346 ),
                                                    Point2D( 0.79742698535309, 0.10128650732346 ),
                                                    Point2D( 0.10128650732346, 0.79742698535309 ),
                                                    Point2D( 0.47014206410512, 0.47014206410512 ),
                                                    Point2D( 0.05971587178977, 0.47014206410512 ),
                                                    Point2D( 0.47014206410512, 0.05971587178977 ),
                                                    Point2D( 0.33333333333333, 0.33333333333333 ) };

} // namespace quadrature

namespace cubature {

/// 3DT-1: exact for polynomial integrands up to order 1
static const std::array< Point3D, 1 > T1_points  = { Point3D( 0.25, 0.25, 0.25 ) };
static const std::array< real_t, 1 >  T1_weights = { 1.0 / 6.0 };

/// 3DT-2: exact for polynomial integrands up to order 1
static const std::array< Point3D, 4 > T2_points  = { Point3D( 0.0, 0.0, 0.0 ),
                                                     Point3D( 1.0, 0.0, 0.0 ),
                                                     Point3D( 0.0, 1.0, 0.0 ),
                                                     Point3D( 0.0, 0.0, 1.0 ) };
static const std::array< real_t, 4 >  T2_weights = { 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0 };

/// 3DT-3: exact for polynomial integrands up to order 2
static const std::array< Point3D, 4 > T3_points = { Point3D( { ( 5.0 - 1.0 * std::sqrt( 5.0 ) ) / 20.0,
                                                               ( 5.0 - 1.0 * std::sqrt( 5.0 ) ) / 20.0,
                                                               ( 5.0 - 1.0 * std::sqrt( 5.0 ) ) / 20.0 } ),
                                                    Point3D( { ( 5.0 + 3.0 * std::sqrt( 5.0 ) ) / 20.0,
                                                               ( 5.0 - 1.0 * std::sqrt( 5.0 ) ) / 20.0,
                                                               ( 5.0 - 1.0 * std::sqrt( 5.0 ) ) / 20.0 } ),
                                                    Point3D( { ( 5.0 - 1.0 * std::sqrt( 5.0 ) ) / 20.0,
                                                               ( 5.0 + 3.0 * std::sqrt( 5.0 ) ) / 20.0,
                                                               ( 5.0 - 1.0 * std::sqrt( 5.0 ) ) / 20.0 } ),
                                                    Point3D( { ( 5.0 - 1.0 * std::sqrt( 5.0 ) ) / 20.0,
                                                               ( 5.0 - 1.0 * std::sqrt( 5.0 ) ) / 20.0,
                                                               ( 5.0 + 3.0 * std::sqrt( 5.0 ) ) / 20.0 } ) };

static const std::array< real_t, 4 > T3_weights = { 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0 };

/// 3DT-4: exact for polynomial integrands up to order 5
/// Generated using quadpy 0.14.11 and the tetrahedron.stroud_t3_5_1() scheme
static const std::array< real_t, 15 > T4_weights = { 0.019753086419753287,
                                                     0.011989513963169754,
                                                     0.011989513963169754,
                                                     0.011989513963169754,
                                                     0.011989513963169754,
                                                     0.011511367871045284,
                                                     0.011511367871045284,
                                                     0.011511367871045284,
                                                     0.011511367871045284,
                                                     0.008818342151675533,
                                                     0.008818342151675533,
                                                     0.008818342151675533,
                                                     0.008818342151675533,
                                                     0.008818342151675533,
                                                     0.008818342151675533 };

static const std::array< Point3D, 15 > T4_points = { Point3D( 0.25, 0.25, 0.25 ),
                                                     Point3D( 0.09197107805272303, 0.09197107805272303, 0.7240867658418308 ),
                                                     Point3D( 0.09197107805272303, 0.7240867658418308, 0.09197107805272303 ),
                                                     Point3D( 0.7240867658418308, 0.09197107805272303, 0.09197107805272303 ),
                                                     Point3D( 0.09197107805272303, 0.09197107805272303, 0.09197107805272303 ),
                                                     Point3D( 0.3197936278296299, 0.3197936278296299, 0.040619116511110276 ),
                                                     Point3D( 0.3197936278296299, 0.040619116511110276, 0.3197936278296299 ),
                                                     Point3D( 0.040619116511110276, 0.3197936278296299, 0.3197936278296299 ),
                                                     Point3D( 0.3197936278296299, 0.3197936278296299, 0.3197936278296299 ),
                                                     Point3D( 0.44364916731037085, 0.056350832689629156, 0.056350832689629156 ),
                                                     Point3D( 0.056350832689629156, 0.44364916731037085, 0.056350832689629156 ),
                                                     Point3D( 0.056350832689629156, 0.056350832689629156, 0.44364916731037085 ),
                                                     Point3D( 0.44364916731037085, 0.44364916731037085, 0.056350832689629156 ),
                                                     Point3D( 0.44364916731037085, 0.056350832689629156, 0.44364916731037085 ),
                                                     Point3D( 0.056350832689629156, 0.44364916731037085, 0.44364916731037085 ) };

} // namespace cubature
} // namespace hyteg
