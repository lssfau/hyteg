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

template < uint_t M, uint_t N >
using Matrixr = Eigen::Matrix< real_t, M, N, N == 1 ? Eigen::ColMajor : Eigen::RowMajor >;
template < typename T, uint_t M, uint_t N >
using Matrix = Eigen::Matrix< T, M, N, N == 1 ? Eigen::ColMajor : Eigen::RowMajor >;
typedef Eigen::Matrix< real_t, 2, 2, Eigen::RowMajor >   Matrix2r;
typedef Eigen::Matrix< real_t, 3, 3, Eigen::RowMajor >   Matrix3r;
typedef Eigen::Matrix< real_t, 4, 4, Eigen::RowMajor >   Matrix4r;
typedef Eigen::Matrix< real_t, 6, 6, Eigen::RowMajor >   Matrix6r;
typedef Eigen::Matrix< real_t, 10, 10, Eigen::RowMajor > Matrix10r;

} // namespace hyteg
