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

namespace terraneo::plates {

/// \defgroup SmoothingStrategies SmoothingStrategies
///
/// In order to avoid problems with discontinuous Dirichlet boundary conditions
/// we can "smooth" plate velocities such that the become zero at the interfaces
/// between plates.
/// @{

/// Linearly decrease smoothing factor from 1.0 to 0.0
///
/// The velocity of a point is scaled by a smoothing factor that linearly
/// increases with distance of the point from the boundary from 0.0 with
/// the given slope and is capped at 1.0.
class LinearDistanceSmoother
{
 public:
   LinearDistanceSmoother( real_t slope )
   : slope_( slope ){};

   real_t operator()( const real_t distanceFromBoundary ) const
   {
      real_t smoothing = distanceFromBoundary * slope_;
      return smoothing > real_c( 1 ) ? real_c( 1 ) : smoothing;
   }

 private:
   real_t slope_{ real_c( 0 ) };
};

/// @}
} // namespace terraneo::plates
