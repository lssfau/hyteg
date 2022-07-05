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

#include "terraneo/helpers/typeAliases.hpp"

namespace terraneo {
namespace plates {

/// Struct to bundle information on the rotation of a plate
struct RotationInfo
{
   real_t time{ real_c( 0 ) };      //!< time for which the finite rotation is given
   real_t longitude{ real_c( 0 ) }; //!< longitude of finite rotation
   real_t latitude{ real_c( 0 ) };  //!< latitude of finite rotation
   real_t angle{ real_c( 0 ) };     //!< finite rotation angle
   uint_t plateID{ 0 };             //!< ID of the plate itself
   uint_t conjugateID{ 0 };         //!< ID of the conjugate plate (needed to "hop" to reference plate)
};

/// Struct to bundle information on a finite rotation
struct FiniteRotation
{
   real_t time;
   vec3D  lonLatAng; // entries give: longitude, latitude, angle
};

/// Alias for iterating over a vector of RotationInfo elements
using rotIter_t = std::vector< RotationInfo >::const_iterator;

/// A "polygon" is describing a plate boundary and given by a vector of points in 3D
using Polygon = std::vector< vec3D >;

/// Hardcoded radius of our planet
///
/// \todo Find a way to set this and give it a more generic name,
/// as we might want to simulate other planets at some point.
struct constants
{
   // static constexpr real_t earthRadiusInKm{real_c( 6371 )};
   // fails with:
   // call to non-‘constexpr’ function ‘walberla::real_t walberla::real_c(T)
   // [with T = int; walberla::real_t = double]’
   static constexpr real_t earthRadiusInKm{ 6371.0 };
};

} // namespace plates
} // namespace terraneo
