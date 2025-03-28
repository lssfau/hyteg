/*
* Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl,
 * Marcus Mohr.
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

namespace hyteg {
template< typename ValueType, int N >
class PointND;


using Point2D = PointND<walberla::real_t, 2>;
using Point3D = PointND<walberla::real_t, 3>;
using Point4D = PointND<walberla::real_t, 4>;
using Point6D = PointND<walberla::real_t, 6>;
using Point10D = PointND<walberla::real_t, 10>;

template< typename ValueType, int M, int N >
class Matrix;

template< int M, int N >
using Matrixr = Matrix<walberla::real_t, M, N>;

using Matrix2r = Matrixr<2, 2>;
using Matrix3r = Matrixr<3, 3>;
using Matrix4r = Matrixr<4, 4>;
using Matrix6r = Matrixr<6, 6>;
using Matrix10r = Matrixr<10, 10>;
}
