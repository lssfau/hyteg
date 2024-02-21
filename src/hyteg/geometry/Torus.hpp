/*
 * Copyright (c) 2017-2021 Nils Kohl.
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

#include <cmath>
#include <vector>

#include "core/math/Constants.h"

#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::math::pi;

inline Point3D torusCoordinates( real_t torusRadiusToTubeCenter, real_t tubeRadius, real_t toroidalAngle, real_t poloidalAngle )
{
   Point3D x( { ( torusRadiusToTubeCenter + tubeRadius * std::cos( poloidalAngle ) ) * std::cos( toroidalAngle ),
                ( torusRadiusToTubeCenter + tubeRadius * std::cos( poloidalAngle ) ) * std::sin( toroidalAngle ),
                tubeRadius * std::sin( poloidalAngle ) } );
   return x;
}

} // namespace hyteg