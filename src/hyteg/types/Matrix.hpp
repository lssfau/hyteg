/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include "hyteg/eigen/EigenWrapper.hpp"

namespace hyteg {

using walberla::int_c;
using walberla::real_t;
using walberla::uint_t;

template < typename ValueType, int M, int N >
using Matrix = Eigen::Matrix< ValueType, M, N, N == 1 ? Eigen::ColMajor : Eigen::RowMajor >;

template < int M, int N >
using Matrixr = Matrix< real_t, M, N >;

using Matrix2r  = Matrixr< 2, 2 >;
using Matrix3r  = Matrixr< 3, 3 >;
using Matrix4r  = Matrixr< 4, 4 >;
using Matrix6r  = Matrixr< 6, 6 >;
using Matrix10r = Matrixr< 10, 10 >;
using MatrixXr  = Matrixr< Eigen::Dynamic, Eigen::Dynamic >;
using VectorXr  = Matrixr< Eigen::Dynamic, 1 >;

} // namespace hyteg
