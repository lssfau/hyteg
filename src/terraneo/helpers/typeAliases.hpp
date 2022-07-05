/*
 * Copyright (c) 2022 Berta Vilacis, Marcus Mohr.
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

namespace terraneo {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

#ifdef WALBERLA_DOUBLE_ACCURACY
using vec3D = Eigen::Vector3d;
using mat3D = Eigen::Matrix3d;
#else
using vec3D = Eigen::Vector3f;
using mat3D = Eigen::Matrix3f;
#endif

} // namespace terraneo
