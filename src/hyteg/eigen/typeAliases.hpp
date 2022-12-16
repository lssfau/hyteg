/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "hyteg/eigen/EigenWrapper.hpp"

namespace Eigen {

// all "standard" Eigen typedefs for real_t
using Matrix2r    = Matrix< walberla::real_t, 2, 2 >;
using Vector2r    = Matrix< walberla::real_t, 2, 1 >;
using RowVector2r = Matrix< walberla::real_t, 1, 2 >;
using Matrix3r    = Matrix< walberla::real_t, 3, 3 >;
using Vector3r    = Matrix< walberla::real_t, 3, 1 >;
using RowVector3r = Matrix< walberla::real_t, 1, 3 >;
using Matrix4r    = Matrix< walberla::real_t, 4, 4 >;
using Vector4r    = Matrix< walberla::real_t, 4, 1 >;
using RowVector4r = Matrix< walberla::real_t, 1, 4 >;
using MatrixXr    = Matrix< walberla::real_t, Dynamic, Dynamic >;
using VectorXr    = Matrix< walberla::real_t, Dynamic, 1 >;
using RowVectorXr = Matrix< walberla::real_t, 1, Dynamic >;
using Matrix2Xr   = Matrix< walberla::real_t, 2, Dynamic >;
using MatrixX2r   = Matrix< walberla::real_t, Dynamic, 2 >;
using Matrix3Xr   = Matrix< walberla::real_t, 3, Dynamic >;
using MatrixX3r   = Matrix< walberla::real_t, Dynamic, 3 >;
using Matrix4Xr   = Matrix< walberla::real_t, 4, Dynamic >;
using MatrixX4r   = Matrix< walberla::real_t, Dynamic, 4 >;

// some more which might come in handy (the same as above for 6 and 10)
using Matrix6r     = Matrix< walberla::real_t, 6, 6 >;
using Vector6r     = Matrix< walberla::real_t, 6, 1 >;
using RowVector6r  = Matrix< walberla::real_t, 1, 6 >;
using Matrix10r    = Matrix< walberla::real_t, 10, 10 >;
using Vector10r    = Matrix< walberla::real_t, 10, 1 >;
using RowVector10r = Matrix< walberla::real_t, 1, 10 >;
using Matrix6Xr    = Matrix< walberla::real_t, 6, Dynamic >;
using MatrixX6r    = Matrix< walberla::real_t, Dynamic, 6 >;
using Matrix10Xr   = Matrix< walberla::real_t, 10, Dynamic >;
using MatrixX10r   = Matrix< walberla::real_t, Dynamic, 10 >;

// partial template spacialization for real_t
template < int Rows, int Cols >
using Matrixr = Matrix< walberla::real_t, Rows, Cols >;

} // namespace Eigen
